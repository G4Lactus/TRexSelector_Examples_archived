// ============================================================================
// demo_ts_perf_03_memory_vs_np.cpp
// ============================================================================
/**
 * @file demo_ts_perf_03_memory_vs_np.cpp
 *
 * @brief Solver memory vs. the ratio n/p, for all demo T-solvers.
 *
 * @details p = 1,000, L = p, n in {50, 100, 200, 500, 1,000, 2,000}, 10 true
 *          coefficients equal to 3 on a fixed random support, noise sigma
 *          about 1, T = 10. Per run and n the parent draws and scales the data
 *          once; every solver then runs in a forked child as a precentered
 *          reader. The child first reads every page of X, D and y and solves
 *          a small warm-up problem, so its baseline counts the data and the
 *          solver code, and reports its peak resident set minus that baseline
 *          as the solver memory.
 *
 *          Usage: demo_ts_perf_03_memory_vs_np [num_MC] [parallel children]   (default 500, 6)
 */
// ============================================================================

// std includes
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// POSIX includes
#include <sys/wait.h>
#include <unistd.h>

// Eigen includes
#include <Eigen/Dense>

// ml_methods includes
#include <ml_methods/scaler_methods/lp_norm_scaler.hpp>

// utils includes
#include <utils/datageneration/utils_datagen.hpp>
#include <utils/eval_metrics/utils_eval_cdiagnostics.hpp>
#include <utils/logging/logger.hpp>
#include <utils/openmp/utils_openmp.hpp>

// demo suite includes
#include "ts_perf_solvers.hpp"
#include "ts_perf_utils.hpp"

// ============================================================================
// Namespace aliases
// ============================================================================

namespace datagen = trex::utils::datageneration::datagen;
namespace scaler = trex::ml_methods::scaler_methods;
namespace omp_utils = trex::utils::openmp;
namespace cdiagnost = trex::utils::eval::cdiagnostics;


// ============================================================================
// Setup
// ============================================================================

namespace {

constexpr std::size_t p = 1'000;
constexpr std::size_t num_dummies = p;
constexpr std::size_t num_active = 10;
constexpr double beta_active = 3.0;
constexpr double snr = 90.0;             // signal variance 10 * 3^2 over noise variance 1
constexpr std::size_t T_stop = 10;
constexpr unsigned support_seed = 12345;
constexpr int data_seed_offset = 1000;   // run mc uses seed data_seed_offset + mc
const std::vector<std::size_t> n_grid = {50, 100, 200, 500, 1'000, 2'000};

/**
 * @brief Fixed random support: num_active of the p indices, sorted.
 * @return 0-based indices of the true coefficients.
 */
std::vector<std::size_t> draw_support() {
    std::vector<std::size_t> all(p);
    std::iota(all.begin(), all.end(), 0);
    std::vector<std::size_t> support;
    std::mt19937 rng(support_seed);
    std::sample(all.begin(), all.end(), std::back_inserter(support), num_active, rng);
    std::sort(support.begin(), support.end());
    return support;
}

/** @brief What a child reports through its pipe. */
struct ChildReport {
    double solver_mib = 0.0;  ///< Solver memory: peak resident set minus the baseline with the data.
    double steps = 0.0;      ///< Steps to T.
    int ok = 0;              ///< 1 when the solver finished.
};

/** @brief A running child: its pid, the read end of its pipe, and its solver. */
struct Pending {
    pid_t pid;
    int fd;
    std::size_t solver;
};

/**
 * @brief Solve a small problem with the solver of @p spec.
 *
 * @details Loads the solver's code pages and warms the allocator, so the
 *          baseline taken afterwards holds them and they do not count as
 *          solver memory.
 */
void warm_up(const ts_perf::SolverSpec& spec) {
    constexpr Eigen::Index n_w = 20, p_w = 30;
    Eigen::MatrixXd X_w = Eigen::MatrixXd::Random(n_w, p_w);
    Eigen::MatrixXd D_w = Eigen::MatrixXd::Random(n_w, p_w);
    Eigen::VectorXd y_w = X_w.leftCols(2).rowwise().sum() + 0.1 * Eigen::VectorXd::Random(n_w);
    X_w.rowwise() -= X_w.colwise().mean();
    X_w.colwise().normalize();
    D_w.rowwise() -= D_w.colwise().mean();
    D_w.colwise().normalize();
    y_w.array() -= y_w.mean();
    Eigen::Map<Eigen::MatrixXd> X_map(X_w.data(), n_w, p_w);
    Eigen::Map<Eigen::MatrixXd> D_map(D_w.data(), n_w, p_w);
    Eigen::Map<Eigen::VectorXd> y_map(y_w.data(), n_w);
    auto solver = spec.make(X_map, D_map, y_map);
    solver->setTieSeed(1);
    solver->executeStep(2, /*early_stop=*/true);
}

/**
 * @brief Run one solver in a forked child and return the child.
 *
 * @details The child reads all data pages and solves a small warm-up problem
 *          before its baseline, so the peak minus the baseline is the memory
 *          the solver allocates for this problem.
 */
Pending launch(const ts_perf::SolverSpec& spec, std::size_t s, Eigen::Map<Eigen::MatrixXd>& X,
               Eigen::Map<Eigen::MatrixXd>& D, Eigen::Map<Eigen::VectorXd>& y, std::uint32_t seed) {
    int fds[2];
    if (pipe(fds) != 0) throw std::runtime_error("pipe() failed");
    const pid_t pid = fork();
    if (pid < 0) throw std::runtime_error("fork() failed");
    if (pid == 0) {
        close(fds[0]);
        ChildReport report;
        volatile double sink = X.sum() + D.sum() + y.sum();  // touch every data page
        (void)sink;
        try {
            warm_up(spec);
        } catch (...) {
            // the measured run below reports a failing solver
        }
        const double baseline = ts_perf::peak_rss_mib();
        try {
            auto solver = spec.make(X, D, y);
            solver->setTieSeed(seed);
            solver->executeStep(T_stop, /*early_stop=*/true);
            report.steps = static_cast<double>(solver->getNumSteps());
            report.solver_mib = ts_perf::peak_rss_mib() - baseline;
            report.ok = 1;
        } catch (...) {
            report.ok = 0;
        }
        const ssize_t written = write(fds[1], &report, sizeof(report));
        (void)written;
        _exit(0);
    }
    close(fds[1]);
    return {pid, fds[0], s};
}

/** @brief Wait for one of the children and read its report. */
std::pair<std::size_t, ChildReport> collect(std::vector<Pending>& running) {
    int status = 0;
    const pid_t pid = waitpid(-1, &status, 0);
    const auto it = std::find_if(running.begin(), running.end(), [pid](const Pending& c) { return c.pid == pid; });
    if (it == running.end()) throw std::runtime_error("waitpid() returned an unknown child");
    ChildReport report;
    if (read(it->fd, &report, sizeof(report)) != static_cast<ssize_t>(sizeof(report))) report.ok = 0;
    close(it->fd);
    const std::size_t s = it->solver;
    running.erase(it);
    return {s, report};
}

} // namespace


// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {

    std::cout << std::unitbuf << "\n";
    cdiagnost::print_section_header("T-Solver Performance: Memory vs. n/p");

    try {
        const auto num_MC = static_cast<std::size_t>(ts_perf::parse_positive(argc > 1 ? argv[1] : nullptr, 500));
        const auto parallel = static_cast<std::size_t>(ts_perf::parse_positive(argc > 2 ? argv[2] : nullptr, 6));

        // Single-threaded parent and solvers, so fork() stays safe with OpenMP
        omp_utils::set_num_threads(1);
        trex::utils::logging::setWarningsEnabled(false);

        const auto support = draw_support();
        const std::vector<double> coefs(num_active, beta_active);
        const auto solvers = ts_perf::demo_solvers();
        const std::size_t S = solvers.size(), K = n_grid.size();

        std::cout << "p = " << p << ", L = p, n in {50, 100, 200, 500, 1000, 2000}, T = " << T_stop
                  << ", true coefficients: " << num_active << " x " << beta_active << ", SNR = " << snr << "\n";
        std::cout << S << " solvers, " << num_MC << " runs, " << parallel
                  << " children at a time (one forked child per solver)\n\n";

        // memory[s][k] and steps[s][k] collect one value per run
        std::vector<std::vector<std::vector<double>>> memory(S, std::vector<std::vector<double>>(K));
        std::vector<std::vector<std::vector<double>>> steps(S, std::vector<std::vector<double>>(K));
        std::size_t failures = 0;
        const auto t_start = std::chrono::steady_clock::now();

        for (std::size_t mc = 0; mc < num_MC; ++mc) {
            for (std::size_t k = 0; k < K; ++k) {
                const std::size_t n = n_grid[k];
                datagen::SyntheticData data(
                    static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(p), static_cast<Eigen::Index>(num_dummies),
                    support, coefs, snr, /*seed=*/data_seed_offset + static_cast<int>(mc), /*X_seed=*/-1,
                    /*dummy_seed=*/-1, datagen::predictor_policy::Normal(), datagen::dummygen::Distribution::Normal(),
                    datagen::noisegen::noise_policy::Normal());
                Eigen::Map<Eigen::MatrixXd> X_map(data.getX().data(), data.rows(), data.cols());
                Eigen::Map<Eigen::MatrixXd> D_map(data.getD().data(), data.rows(), data.num_dummies());
                Eigen::Map<Eigen::VectorXd> y_map(data.getY().data(), data.rows());

                // Center and scale once in front; the solvers are precentered readers
                scaler::LpNormScaler scaler_X(scaler::LpNormScaler::NormType::L2, true);
                scaler_X.fit_transform_inplace(X_map);
                scaler::LpNormScaler scaler_D(scaler::LpNormScaler::NormType::L2, true);
                scaler_D.fit_transform_inplace(D_map);
                y_map.array() -= y_map.mean();

                std::vector<Pending> running;
                const auto seed = static_cast<std::uint32_t>(data_seed_offset + mc);
                auto store = [&](const std::pair<std::size_t, ChildReport>& result) {
                    if (result.second.ok) {
                        memory[result.first][k].push_back(result.second.solver_mib);
                        steps[result.first][k].push_back(result.second.steps);
                    } else {
                        ++failures;
                    }
                };
                for (std::size_t s = 0; s < S; ++s) {
                    if (running.size() >= parallel) store(collect(running));
                    running.push_back(launch(solvers[s], s, X_map, D_map, y_map, seed));
                }
                while (!running.empty()) store(collect(running));
            }
            if ((mc + 1) % 50 == 0 || mc + 1 == num_MC)
                std::cout << "  runs done: " << (mc + 1) << " of " << num_MC << "\n";
        }
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();

        // Data size per n: X, D and y as doubles
        std::vector<double> data_mib(K);
        for (std::size_t k = 0; k < K; ++k)
            data_mib[k] = 8.0 * static_cast<double>(n_grid[k] * (p + num_dummies + 1)) / (1024.0 * 1024.0);

        std::ostringstream out;
        out << "\n";
        auto header = [&](const std::string& title) {
            out << title << "\n" << std::left << std::setw(22) << "solver" << std::right;
            for (const auto n : n_grid) out << std::setw(9) << ("n=" + std::to_string(n));
            out << "\n" << std::string(22 + 9 * K, '-') << "\n";
        };
        header("Data X, D and y (MiB)");
        out << std::left << std::setw(22) << "all solvers" << std::right << std::fixed << std::setprecision(2);
        for (std::size_t k = 0; k < K; ++k) out << std::setw(9) << data_mib[k];
        out << "\n\n";
        auto table = [&](const std::string& title, auto value, int precision) {
            header(title);
            out << std::fixed << std::setprecision(precision);
            for (std::size_t s = 0; s < S; ++s) {
                out << std::left << std::setw(22) << solvers[s].name << std::right;
                for (std::size_t k = 0; k < K; ++k) out << std::setw(9) << value(s, k);
                out << "\n";
            }
            out << "\n";
        };
        table("Median solver memory (MiB)", [&](std::size_t s, std::size_t k) { return ts_perf::median(memory[s][k]); }, 2);
        table("Largest solver memory over the runs (MiB)", [&](std::size_t s, std::size_t k) {
            return memory[s][k].empty() ? 0.0 : *std::max_element(memory[s][k].begin(), memory[s][k].end()); }, 2);
        table("Median solver memory relative to the data (%)",
              [&](std::size_t s, std::size_t k) { return 100.0 * ts_perf::median(memory[s][k]) / data_mib[k]; }, 1);
        table("Median steps to T = 10", [&](std::size_t s, std::size_t k) { return ts_perf::median(steps[s][k]); }, 1);

        out << "Failed solver runs: " << failures << "\n";
        out << "Wall time: " << std::setprecision(1) << wall << " s\n";
        std::cout << out.str();

        // Result files
        const std::string stem = "memory_vs_np_mc" + std::to_string(num_MC);
        const auto csv_path = ts_perf::output_dir() / (stem + ".csv");
        std::ofstream csv(csv_path);
        csv << "solver,family,n,n_over_p,data_mib,solver_mib,solver_max_mib,steps\n" << std::setprecision(6);
        for (std::size_t s = 0; s < S; ++s) {
            for (std::size_t k = 0; k < K; ++k) {
                const double max_memory = memory[s][k].empty() ? 0.0 : *std::max_element(memory[s][k].begin(), memory[s][k].end());
                csv << solvers[s].name << ',' << solvers[s].family << ',' << n_grid[k] << ','
                    << static_cast<double>(n_grid[k]) / static_cast<double>(p) << ',' << data_mib[k] << ','
                    << ts_perf::median(memory[s][k]) << ',' << max_memory << ',' << ts_perf::median(steps[s][k]) << '\n';
            }
        }
        const auto txt_path = ts_perf::output_dir() / (stem + ".txt");
        std::ofstream(txt_path) << out.str();

        std::cout << "\nResults written to:\n  " << ts_perf::display_path(csv_path)
                  << "\n  " << ts_perf::display_path(txt_path) << "\n\n";

    } catch (const std::exception& e) {
        std::cerr << "\n[ERROR] " << e.what() << "\n"
                  << "Usage: demo_ts_perf_03_memory_vs_np [num_MC] [parallel children]\n\n";
        return 1;
    }

    return 0;
}
