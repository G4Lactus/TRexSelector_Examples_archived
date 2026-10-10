// ==============================================================================
// demo_trex_08_mc_sim_scalability.cpp
// ==============================================================================
/**
 * @file demo_trex_08_mc_sim_scalability.cpp
 *
 * @brief Scalability benchmark for the T-Rex Selector: runtime and memory usage
 *        as the problem scales, across three base solvers (TLARS, TOMP, TAFS).
 *
 * @details Two aims:
 *   (a) high-dimensionality — sweep p up to 200,000;
 *   (b) volumes beyond RAM  — the largest p columns push the design [X|D] past
 *       physical RAM, feasible only because the pipeline is fully out-of-core.
 *
 *  Design of the study:
 *   - 2-D grid  n in {1000,1500,2000,5000} x p in {1k,5k,10k,20k,50k,100k,150k,200k}
 *   - the top-right corner is out-of-core: at n=5000, p=200k the logical
 *     [X|D] = 5000*11*200k*8 ~ 82 GiB > 64 GiB RAM (the real on-disk dummy
 *     footprint is far larger still), yet resident RSS stays a few GB.
 *   - SNR = 1 fixed (SNR is not a scaling axis; sweeping it only burns time),
 *     demo-03 DGP (variable support s=10 redrawn per trial, unit coefficients,
 *     Normal predictors/noise).
 *   - num_MC = 10 SEQUENTIAL trials per (n, p, solver). Trials are NOT run in
 *     parallel: this is a timing/memory benchmark and concurrent trials would
 *     corrupt per-trial wall-clock, confound RSS, and multiply the on-disk
 *     dummy backing (one trial already peaks at >100 GB of D, since the K
 *     random experiments run concurrently via inner OpenMP — which is the
 *     parallelism this benchmark relies on).
 *
 *  Fully memory-mapped pipeline in every point: X/y are generated into mmap
 *  backing files (SyntheticDataMapped, removed per trial by an RAII
 *  MmapFileGuard), and use_memory_mapping = true memory-maps the dummy matrices
 *  D and serializes solver warm-start state between T-loop iterations.
 *
 *  Per grid point we record: timing (DGP, select() median/min/max, per-T-iter),
 *  memory (post-select RSS, physical footprint, process peak RSS, X+y bytes,
 *  theoretical design volume), and FDR / TPR / Avg L / Avg T.
 *
 *  Self-cleanup: the process sweeps its own memory-mapped orphans
 *  (demo08_scalability_* backing files and trex_memmap_* dummy dirs in the temp
 *  dir) at startup and on exit / SIGINT / SIGTERM, so a cancelled run does not
 *  strand tens of GB. (SIGKILL cannot be trapped; the next startup sweep
 *  reclaims it.)
 *
 *  Expandable by design. Every trial is appended to a raw per-trial cache
 *  (<stem>_trials.csv, one row per (n, p, solver, mc) keyed by its seed). A run
 *  computes only the trials it is missing and derives the summary table/CSV from
 *  the cache, so: num_MC top-up (run at 2 now, later 10 computes only the 8 new
 *  trials per cell), grid growth (append new n/p values later; existing cells are
 *  reused), and trial-level resume (a crash/OOM mid-point keeps finished trials).
 *  Per-trial seeds are a stable function of (n, p, mc) (grid-position independent
 *  and identical to the Python demo), so reuse and reproducibility survive
 *  reordering or extending the grid.
 *
 *  Results go to simulation_results/data/ as an aligned table (.txt) and a tidy
 *  summary CSV, both re-derived from the cache and re-saved after each (n, p) point.
 *
 *  Run:  ./demo_trex_08_mc_sim_scalability            (full benchmark, num_MC=10)
 *        ./demo_trex_08_mc_sim_scalability 2          (full benchmark, num_MC=2)
 *        ./demo_trex_08_mc_sim_scalability smoke      (tiny end-to-end check)
 *        ./demo_trex_08_mc_sim_scalability smoke 3    (smoke grid, num_MC=3)
 *
 * @note This simulation takes a substantial amount of time (the largest points
 *       are minutes per trial); the top-right corner needs ~hundreds of GB of
 *       free disk for the memory-mapped dummy backing.
 */
// ==============================================================================

// std includes
#include <algorithm>
#include <cctype>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// System includes (memory profiling hooks + physical RAM query)
#include <sys/resource.h>
#if defined(__APPLE__)
#include <mach/mach.h>
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

// Eigen includes
#include <Eigen/Dense>

// OpenMP compatibility layer
#include <utils/openmp/utils_openmp.hpp>

// T-Rex Selector includes
#include <trex_selector_methods/ad/trex_core/trex.hpp>
#include <utils/datageneration/utils_datagen.hpp>
#include <utils/eval_metrics/utils_eval_cdiagnostics.hpp>
#include <utils/eval_metrics/utils_eval_rates.hpp>
#include <utils/memmap/memory_mapped_matrix.hpp>

// Demo utilities
#include "trex_sim_utils.hpp"
using namespace trex_sim;


// ==============================================================================
// Namespace aliases
// ==============================================================================

namespace cdianostics = trex::utils::eval::cdiagnostics;
namespace datagen     = trex::utils::datageneration::datagen;
namespace dummygen    = trex::utils::datageneration::dummygen;
namespace rates       = trex::utils::eval::rates;

// T-Rex types
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;
using trex::trex_selector_methods::ad::trex_core::TRexSelector;
using trex::trex_selector_methods::ad::trex_core::LLoopStrategy;
using trex::trex_selector_methods::ad::utils::solver_dispatch::SolverTypeForTRex;


// ==============================================================================
// Study configuration
// ==============================================================================

static const std::vector<Eigen::Index> N_VALUES =
    {1'000, 1'500, 2'000, 5'000};
static const std::vector<Eigen::Index> P_VALUES =
    {1'000, 5'000, 10'000, 20'000, 50'000, 100'000, 150'000, 200'000};
static constexpr std::size_t NUM_MC     = 10;      // sequential trials per point
static constexpr double      SNR        = 1.0;     // fixed; not a scaling axis
static constexpr double      TFDR       = 0.1;
static constexpr std::size_t S_TRUE     = 10;      // true support cardinality
static constexpr int         DUMMY_MULT = 10;      // max_dummy_multiplier (10 x p)

/** @brief Theoretical design volume [X|D] in GiB (X plus 10x-p dummies, f64). */
static double design_gib(Eigen::Index n, Eigen::Index p) {
    return static_cast<double>(n) * (1.0 + DUMMY_MULT) *
           static_cast<double>(p) * 8.0 / (1024.0 * 1024.0 * 1024.0);
}

/** @brief Deterministic, grid-position-independent, 32-bit-safe per-trial seed.
 *
 *  Distinct (n, p) already yield different data (the X draw has a different
 *  shape), so we only need determinism in (n, p, mc), not cross-cell seed
 *  uniqueness. Matches the Python demo's _trial_seed so both languages seed
 *  identically. This is the key to expandability: reuse and reproducibility
 *  survive reordering or extending the grid because the seed depends on the
 *  (n, p) values, never on the grid position. */
static unsigned trial_seed(Eigen::Index n, Eigen::Index p, std::size_t mc) {
    const long long v = 24LL + static_cast<long long>(n) * 100003LL +
                        static_cast<long long>(p) * 13LL +
                        static_cast<long long>(mc);
    return static_cast<unsigned>(v % 2147483647LL);   // mod 2^31 - 1
}

/** @brief Physical RAM in GiB (NaN if unavailable) — used to flag out-of-core points. */
static double physical_ram_gib() {
#if defined(__APPLE__)
    int64_t mem = 0;
    std::size_t len = sizeof(mem);
    if (sysctlbyname("hw.memsize", &mem, &len, nullptr, 0) == 0 && mem > 0)
        return static_cast<double>(mem) / (1024.0 * 1024.0 * 1024.0);
#elif defined(__linux__)
    std::ifstream mi("/proc/meminfo");
    std::string line;
    while (std::getline(mi, line)) {
        if (line.rfind("MemTotal:", 0) == 0)
            return std::stod(line.substr(line.find(':') + 1)) /
                   (1024.0 * 1024.0);   // kB -> GiB
    }
#endif
    return std::numeric_limits<double>::quiet_NaN();
}


// ==============================================================================
// Memory profiling hooks
// ==============================================================================

/** @brief Current resident set size (RSS) of this process in MiB. */
static double current_rss_mib() {
#if defined(__APPLE__)
    mach_task_basic_info   info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS)
        return static_cast<double>(info.resident_size) / (1024.0 * 1024.0);
    return 0.0;
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm");
    long total_pages = 0, resident_pages = 0;
    if (statm >> total_pages >> resident_pages)
        return static_cast<double>(resident_pages) *
               static_cast<double>(sysconf(_SC_PAGESIZE)) / (1024.0 * 1024.0);
    return 0.0;
#else
    return 0.0;
#endif
}

/** @brief Physical footprint of this process in MiB (heap incl. compressed/swap). */
static double current_footprint_mib() {
#if defined(__APPLE__)
    task_vm_info_data_t    info;
    mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_VM_INFO,
                  reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS)
        return static_cast<double>(info.phys_footprint) / (1024.0 * 1024.0);
    return 0.0;
#elif defined(__linux__)
    std::ifstream status("/proc/self/status");
    std::string line;
    double kib = 0.0;
    while (std::getline(status, line)) {
        if (line.rfind("VmRSS:", 0) == 0 || line.rfind("VmSwap:", 0) == 0)
            kib += std::stod(line.substr(line.find(':') + 1));
    }
    return kib / 1024.0;
#else
    return 0.0;
#endif
}

/** @brief Peak resident set size of this process in MiB (monotone, process-wide). */
static double peak_rss_mib() {
    struct rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
#if defined(__APPLE__)
    return static_cast<double>(ru.ru_maxrss) / (1024.0 * 1024.0);  // bytes
#else
    return static_cast<double>(ru.ru_maxrss) / 1024.0;             // KiB
#endif
}


// ==============================================================================
// Self-cleanup of memory-mapped orphans
// ==============================================================================

/** @brief Recursive on-disk size of a file or directory in bytes (best-effort). */
static std::uintmax_t path_size_bytes(const std::filesystem::path& p) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (fs::is_directory(p, ec)) {
        std::uintmax_t total = 0;
        for (fs::recursive_directory_iterator it(p, ec), end; it != end;
             it.increment(ec)) {
            std::error_code e2;
            if (it->is_regular_file(e2)) total += it->file_size(e2);
        }
        return total;
    }
    return fs::file_size(p, ec);
}

/** @brief Remove this demo's stray mmap backing (demo08_scalability_* and
 *         trex_memmap_*) from the temp dir and CWD. Assumes SOLO execution. */
static void sweep_mmap_orphans(const char* reason) {
    namespace fs = std::filesystem;
    std::error_code ec;
    std::uintmax_t freed = 0;
    for (const fs::path& dir : {fs::temp_directory_path(ec), fs::current_path(ec)}) {
        if (dir.empty()) continue;
        for (fs::directory_iterator it(dir, ec), end; it != end && !ec;
             it.increment(ec)) {
            const std::string name = it->path().filename().string();
            if (name.rfind("demo08_scalability_", 0) == 0 ||
                name.rfind("trex_memmap_", 0) == 0) {
                freed += path_size_bytes(it->path());
                std::error_code e2;
                fs::remove_all(it->path(), e2);
            }
        }
    }
    if (freed > 0)
        std::cout << "[cleanup] swept "
                  << std::fixed << std::setprecision(1)
                  << (static_cast<double>(freed) / (1024.0 * 1024.0 * 1024.0))
                  << " GiB of mmap orphans (" << reason << ")\n" << std::flush;
}

extern "C" void on_terminating_signal(int sig) {
    sweep_mmap_orphans("signal");
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}

static void install_cleanup_handlers() {
    std::atexit([]() { sweep_mmap_orphans("atexit"); });
    std::signal(SIGINT,  on_terminating_signal);
    std::signal(SIGTERM, on_terminating_signal);
}


// ==============================================================================
// Grid point descriptor and result types
// ==============================================================================

/** @brief One (n, p) grid point on the fully memory-mapped pipeline. */
struct ScalabilityPoint {
    Eigen::Index n;             ///< Number of observations
    Eigen::Index p;            ///< Number of predictors
    std::size_t  num_MC;       ///< Number of sequential MC trials
    bool         out_of_core;  ///< design [X|D] exceeds physical RAM
    std::string label() const {
        return "n" + std::to_string(n) + "_p" + std::to_string(p);
    }
};

/** @brief Measurements from a single MC trial. */
struct ScalabilityTrialResult {
    double dgp_s    = 0.0;   ///< Data generation wall-clock [s]
    double select_s = 0.0;   ///< Selector construction + select() wall-clock [s]
    double fdp      = 0.0;   ///< False discovery proportion
    double tpp      = 0.0;   ///< True positive proportion
    double L        = 0.0;   ///< Calibrated dummy multiplier
    double T        = 0.0;   ///< Stopping time T_stop
    double data_mib = 0.0;   ///< X + y mmap backing-file sizes on disk [MiB]
    double rss_mib  = 0.0;   ///< Current RSS right after select() [MiB]
    double fp_mib   = 0.0;   ///< Physical footprint right after select() [MiB]
    double peak_rss_mib = 0.0;  ///< Process peak RSS up to and incl. this trial [MiB]
};

/** @brief Aggregate result for one grid point (n, p, solver). */
struct ScalabilityGridPointResult {
    double design_gib   = 0.0;   ///< Theoretical [X|D] volume [GiB]
    double avg_fdr      = 0.0;   ///< Mean false discovery rate
    double avg_tpr      = 0.0;   ///< Mean true positive rate
    double avg_L        = 0.0;   ///< Mean calibrated dummy multiplier
    double avg_T        = 0.0;   ///< Mean stopping time T_stop
    double avg_dgp_s    = 0.0;   ///< Mean data-generation time [s]
    double avg_titer_s  = 0.0;   ///< Mean select time per T-loop iteration [s]
    double sel_med_s    = 0.0;   ///< Median select() time across trials [s]
    double sel_min_s    = 0.0;   ///< Min select() time across trials [s]
    double sel_max_s    = 0.0;   ///< Max select() time across trials [s]
    double max_rss_mib  = 0.0;   ///< Max post-select() RSS across trials [MiB]
    double max_fp_mib   = 0.0;   ///< Max post-select() footprint across trials [MiB]
    double peak_rss_mib = 0.0;   ///< Process peak RSS at grid-point end [MiB]
    double data_mib     = 0.0;   ///< X + y mmap backing-file sizes [MiB]
    std::size_t num_mc  = 0;     ///< Number of cached trials this aggregate covers
};


// ==============================================================================
// Grid construction
// ==============================================================================

/** @brief Full n x p benchmark grid; out_of_core flags points whose logical
 *         design [X|D] exceeds this machine's physical RAM. */
static std::vector<ScalabilityPoint> make_benchmark_grid(double ram_gib,
                                                         std::size_t num_MC) {
    std::vector<ScalabilityPoint> grid;
    grid.reserve(N_VALUES.size() * P_VALUES.size());
    for (const auto n : N_VALUES)
        for (const auto p : P_VALUES)
            grid.push_back({n, p, num_MC, design_gib(n, p) > ram_gib});
    return grid;
}

/** @brief Tiny end-to-end smoke grid (pipeline verification only). */
static std::vector<ScalabilityPoint> make_smoke_grid(std::size_t num_MC) {
    return {{150, 500, num_MC, false}, {200, 1'000, num_MC, false}};
}


// ==============================================================================
// Configuration helpers
// ==============================================================================

/** @brief T-Rex control parameters (base config of demo 02/03) for one point. */
static TRexControlParameter make_scalability_trex_ctrl(const DemoSolverInfo& solver)
{
    // Every control option is set explicitly and identically to the Python
    // demo's _make_scalability_ctrl_dict(), so the two demos are provably
    // comparable — nothing is left to a struct/helper default (opt_threshold and
    // parallel_rnd_experiments in particular differ between the C++ struct and
    // the Python helper, which is what made the demos diverge).
    TRexControlParameter ctrl;
    ctrl.K                        = 20;
    ctrl.max_dummy_multiplier     = DUMMY_MULT;
    ctrl.use_max_T_stop           = true;
    ctrl.opt_threshold            = 0.75;
    ctrl.dummy_distribution       = dummygen::Distribution::Normal();
    ctrl.lloop_strategy           = LLoopStrategy::STANDARD;
    ctrl.tloop_stagnation_stop    = true;
    ctrl.tloop_max_stagnant_steps = 5;
    // Sequential random experiments (mirrors R use_openmp=FALSE, as in demos
    // 02/03/07): one experiment working set in RAM at a time (Eigen uses the
    // inner threads), not K concurrent ones — essential for the out-of-core
    // corner. This is the C++ struct default, but pinned here so it is explicit.
    ctrl.parallel_rnd_experiments = false;
    // Always on: memory-mapped dummy matrices D + solver warm-start state
    // serialized to disk between T-loop iterations. Together with the
    // mmap-generated X/y this makes the whole pipeline memory-mapped.
    ctrl.use_memory_mapping       = true;
    ctrl.solver_type              = solver.solver_type;
    ctrl.solver_params.lambda2    = solver.lambda2;
    ctrl.solver_params.rho_afs    = solver.rho_afs;
    return ctrl;
}

/** @brief Draw a variable active support (demo 03: uniform w/o replacement). */
static std::vector<std::size_t> draw_support(unsigned seed, Eigen::Index p,
                                             std::size_t cardinality) {
    std::mt19937 rng_sup(seed + 500000u);
    std::vector<std::size_t> pool(static_cast<std::size_t>(p));
    std::iota(pool.begin(), pool.end(), std::size_t{0});
    std::shuffle(pool.begin(), pool.end(), rng_sup);
    pool.resize(cardinality);
    std::sort(pool.begin(), pool.end());
    return pool;
}

/** @brief Unit or N(0,1) coefficients for the active support (demo 03 pattern). */
static std::vector<double> make_coefs(unsigned seed, std::size_t cardinality,
                                      bool rnd_coef) {
    std::vector<double> coefs;
    coefs.reserve(cardinality);
    if (rnd_coef) {
        std::mt19937 rng_coef(seed + 600000u);
        std::normal_distribution<double> nd(0.0, 1.0);
        for (std::size_t i = 0; i < cardinality; ++i) coefs.push_back(nd(rng_coef));
    } else {
        for (std::size_t i = 0; i < cardinality; ++i) coefs.push_back(1.0);
    }
    return coefs;
}


// ==============================================================================
// Single-trial runner (sequential; timing + memory instrumented)
// ==============================================================================

using bench_clock = std::chrono::steady_clock;

static double seconds_since(const bench_clock::time_point& t0) {
    return std::chrono::duration<double>(bench_clock::now() - t0).count();
}

/** @brief One MC trial on the fully memory-mapped pipeline. MmapFileGuard is
 *         declared BEFORE SyntheticDataMapped so reverse-LIFO destruction closes
 *         the mmap handles before removing the backing files (exception-safe). */
static ScalabilityTrialResult run_scalability_trial(
    const ScalabilityPoint&         point,
    const TRexControlParameter&     trex_ctrl,
    double                          tFDR,
    double                          snr,
    unsigned                        seed,
    const std::vector<std::size_t>& support,
    const std::vector<double>&      coefs,
    const std::string&              mmap_filestem)
{
    ScalabilityTrialResult m;

    const std::string X_path = mmap_filestem + "_X.dat";
    const std::string y_path = mmap_filestem + "_y.dat";

    MmapFileGuard guard{{X_path, y_path}};

    const auto t0 = bench_clock::now();
    datagen::SyntheticDataMapped data(
        X_path, y_path,
        point.n, point.p,
        support, coefs, snr,
        static_cast<int>(seed));
    m.dgp_s = seconds_since(t0);

    std::error_code ec;
    const auto X_bytes = std::filesystem::file_size(X_path, ec);
    const auto y_bytes = std::filesystem::file_size(y_path, ec);
    if (!ec)
        m.data_mib = static_cast<double>(X_bytes + y_bytes) / (1024.0 * 1024.0);

    // Map views over the mmap backing files — X never fully resides in RAM
    Eigen::Map<Eigen::MatrixXd> X_map(data.getX().data(), data.rows(), data.cols());
    Eigen::Map<Eigen::VectorXd> y_map(data.getY().data(), data.rows());

    const auto t1 = bench_clock::now();
    TRexSelector selector(X_map, y_map, tFDR, trex_ctrl, -1, false);
    selector.select();
    m.select_s = seconds_since(t1);

    m.rss_mib      = current_rss_mib();
    m.fp_mib       = current_footprint_mib();
    m.peak_rss_mib = peak_rss_mib();   // process-lifetime max up to this trial

    const auto sel = selector.getSelectedIndices();
    m.fdp = rates::compute_fdp(sel, support);
    m.tpp = rates::compute_tpp(sel, support);
    m.L   = static_cast<double>(selector.getDummyMultiplierL());
    m.T   = static_cast<double>(selector.getTStop());
    return m;
}


// ==============================================================================
// Grid-point loop (n, p, solver): sequential MC trials
// ==============================================================================

static double vmedian(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const std::size_t mid = v.size() / 2;
    return (v.size() % 2) ? v[mid] : 0.5 * (v[mid - 1] + v[mid]);
}


// ==============================================================================
// Expandable per-trial cache (top-up num_MC / grid; trial-level resume)
// ==============================================================================
//
// Every trial is appended to a raw per-trial store (<stem>_trials.csv, one row
// per (n, p, solver, mc)); the summary .txt/.csv are derived from it, so a later
// run computes only the trials it is missing and re-aggregates. This makes both
// the num_MC axis and the grid expandable, and gives trial-level crash resume.

/** @brief One cached trial (a row of the raw per-trial store). */
struct TrialRow {
    std::size_t mc       = 0;
    unsigned    seed     = 0;
    double dgp_s = 0, select_s = 0, rss_mib = 0, fp_mib = 0, peak_rss_mib = 0;
    double data_mib = 0, fdp = 0, tpp = 0, L = 0, T = 0;
};

/** @brief Trial cache: (n, p, solver) -> {mc -> TrialRow}. */
using TrialCache =
    std::map<std::tuple<Eigen::Index, Eigen::Index, std::string>,
             std::map<std::size_t, TrialRow>>;

static const char* const RAW_HEADER =
    "n,p,solver,mc,trial_seed,dgp_s,select_s,rss_mib,fp_mib,peak_rss_mib,"
    "data_mib,fdp,tpp,L,T";

/** @brief Load the raw per-trial store into a TrialCache (empty if absent). */
static TrialCache load_trial_cache(const std::string& raw_path) {
    TrialCache cache;
    std::ifstream in(raw_path);
    if (!in.is_open()) return cache;

    std::string line;
    std::getline(in, line);   // header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) f.push_back(cell);
        if (f.size() < 15) continue;
        const Eigen::Index n = std::stoll(f[0]);
        const Eigen::Index p = std::stoll(f[1]);
        const std::string  solver = f[2];
        TrialRow r;
        r.mc           = static_cast<std::size_t>(std::stoull(f[3]));
        r.seed         = static_cast<unsigned>(std::stoul(f[4]));
        r.dgp_s        = std::stod(f[5]);
        r.select_s     = std::stod(f[6]);
        r.rss_mib      = std::stod(f[7]);
        r.fp_mib       = std::stod(f[8]);
        r.peak_rss_mib = std::stod(f[9]);
        r.data_mib     = std::stod(f[10]);
        r.fdp          = std::stod(f[11]);
        r.tpp          = std::stod(f[12]);
        r.L            = std::stod(f[13]);
        r.T            = std::stod(f[14]);
        cache[{n, p, solver}][r.mc] = r;
    }
    return cache;
}

/** @brief Append one trial to the raw store (write header if the file is new). */
static void append_trial_row(const std::string& raw_path, Eigen::Index n,
                             Eigen::Index p, const std::string& solver,
                             const TrialRow& r) {
    const bool exists = std::filesystem::exists(raw_path);
    std::ofstream out(raw_path, std::ios::app);
    if (!out.is_open()) return;
    if (!exists) out << RAW_HEADER << "\n";
    out << std::fixed << std::setprecision(6)
        << n << "," << p << "," << solver << "," << r.mc << "," << r.seed << ","
        << r.dgp_s << "," << r.select_s << "," << r.rss_mib << "," << r.fp_mib
        << "," << r.peak_rss_mib << "," << r.data_mib << "," << r.fdp << ","
        << r.tpp << "," << r.L << "," << r.T << "\n";
}

/** @brief Aggregate every cached trial of one (n, p, solver) cell into a summary.
 *
 *  Aggregates over all cached trials (monotone: more trials -> a better estimate),
 *  so re-running at a lower num_MC never discards work. */
static ScalabilityGridPointResult aggregate_cell(
    const std::map<std::size_t, TrialRow>& trials, Eigen::Index n, Eigen::Index p)
{
    ScalabilityGridPointResult agg;
    agg.design_gib = design_gib(n, p);
    if (trials.empty()) return agg;

    std::vector<double> selects;
    selects.reserve(trials.size());
    for (const auto& [mc, r] : trials) {
        (void)mc;
        selects.push_back(r.select_s);
        agg.avg_fdr     += r.fdp;
        agg.avg_tpr     += r.tpp;
        agg.avg_L       += r.L;
        agg.avg_T       += r.T;
        agg.avg_dgp_s   += r.dgp_s;
        agg.avg_titer_s += r.select_s / std::max(r.T, 1.0);
        agg.max_rss_mib   = std::max(agg.max_rss_mib, r.rss_mib);
        agg.max_fp_mib    = std::max(agg.max_fp_mib, r.fp_mib);
        agg.peak_rss_mib  = std::max(agg.peak_rss_mib, r.peak_rss_mib);
        agg.data_mib      = r.data_mib;   // identical across a cell's trials
    }
    const double k = static_cast<double>(trials.size());
    agg.avg_fdr     /= k;
    agg.avg_tpr     /= k;
    agg.avg_L       /= k;
    agg.avg_T       /= k;
    agg.avg_dgp_s   /= k;
    agg.avg_titer_s /= k;
    agg.sel_med_s = vmedian(selects);
    agg.sel_min_s = *std::min_element(selects.begin(), selects.end());
    agg.sel_max_s = *std::max_element(selects.begin(), selects.end());
    agg.num_mc    = trials.size();
    return agg;
}

/** @brief Ensure >= num_MC trials exist for this (n, p, solver) cell, then
 *         aggregate. Only the missing trial indices are computed; each is
 *         appended to the raw store immediately (trial-level durability). */
static ScalabilityGridPointResult run_scalability_grid_point(
    const ScalabilityPoint&          point,
    const DemoSolverInfo&            solver,
    double                           tFDR,
    double                           snr,
    bool                             rnd_coef,
    std::size_t                      cardinality_true_support,
    std::map<std::size_t, TrialRow>& cell_cache,
    const std::string&               raw_path)
{
    const TRexControlParameter trex_ctrl = make_scalability_trex_ctrl(solver);
    const double design = design_gib(point.n, point.p);

    const std::string mmap_filestem =
        (std::filesystem::temp_directory_path() /
         ("demo08_scalability_" + point.label())).string();

    const std::size_t have = cell_cache.size();
    std::vector<std::size_t> todo;
    for (std::size_t mc = 0; mc < point.num_MC; ++mc)
        if (cell_cache.find(mc) == cell_cache.end()) todo.push_back(mc);

    std::cout << "  [" << solver.solver_name << "] "
              << "n=" << point.n << ", p=" << point.p
              << "  (design [X|D] = " << std::fixed << std::setprecision(1)
              << design << " GiB"
              << (point.out_of_core ? "  (> RAM — out-of-core)" : "")
              << ") — " << todo.size() << " new / " << point.num_MC
              << " MC trials (mmap)"
              << (have ? "  (reusing " + std::to_string(have) + " cached)" : "")
              << " ...\n" << std::flush;

    for (const std::size_t mc : todo) {

        const unsigned seed = trial_seed(point.n, point.p, mc);
        const auto support = draw_support(seed, point.p,
                                          cardinality_true_support);
        const auto coefs   = make_coefs(seed, cardinality_true_support,
                                        rnd_coef);

        const ScalabilityTrialResult m = run_scalability_trial(
            point, trex_ctrl, tFDR, snr, seed, support, coefs, mmap_filestem);

        TrialRow r;
        r.mc = mc; r.seed = seed;
        r.dgp_s = m.dgp_s; r.select_s = m.select_s; r.rss_mib = m.rss_mib;
        r.fp_mib = m.fp_mib; r.peak_rss_mib = m.peak_rss_mib;
        r.data_mib = m.data_mib; r.fdp = m.fdp; r.tpp = m.tpp; r.L = m.L; r.T = m.T;
        append_trial_row(raw_path, point.n, point.p, solver.solver_name, r);
        cell_cache[mc] = r;   // durable before we continue

        std::cout << "    trial " << (mc + 1) << "/" << point.num_MC
                  << ": dgp=" << std::fixed << std::setprecision(2) << m.dgp_s
                  << "s  select=" << m.select_s
                  << "s  FDP=" << std::setprecision(3) << m.fdp
                  << "  TPP=" << m.tpp
                  << "  L=" << std::setprecision(0) << m.L
                  << "  T=" << m.T
                  << "  RSS=" << std::setprecision(1) << m.rss_mib << " MiB"
                  << "  FP=" << m.fp_mib << " MiB\n"
                  << std::flush;
    }

    const ScalabilityGridPointResult agg =
        aggregate_cell(cell_cache, point.n, point.p);

    std::cout << "  [" << solver.solver_name << "] "
              << "n=" << point.n << ", p=" << point.p << " done ("
              << agg.num_mc << " trials). "
              << "select median=" << std::fixed << std::setprecision(2)
              << agg.sel_med_s << "s [" << agg.sel_min_s << ", "
              << agg.sel_max_s << "]  peakRSS=" << std::setprecision(0)
              << agg.peak_rss_mib << " MiB  FDR=" << std::setprecision(3)
              << agg.avg_fdr << "  TPR=" << agg.avg_tpr << "\n\n" << std::flush;

    return agg;
}


// ==============================================================================
// Output: aligned table (console + .txt) and tidy CSV
// ==============================================================================

/** @brief Metric column descriptor for the results table / CSV. */
struct ScalabilityMetricSpec {
    std::string name;                                  ///< Column / CSV label
    int         precision;                             ///< Table decimals
    double (*get)(const ScalabilityGridPointResult&);  ///< Value extractor
};

/** @brief Per-grid-point metric columns (table order and CSV metric set). */
static const std::vector<ScalabilityMetricSpec>& scalability_metrics() {
    static const std::vector<ScalabilityMetricSpec> specs = {
        {"DesignGiB",   1, [](const ScalabilityGridPointResult& r) { return r.design_gib;   }},
        {"SelMed_s",    3, [](const ScalabilityGridPointResult& r) { return r.sel_med_s;     }},
        {"SelMin_s",    3, [](const ScalabilityGridPointResult& r) { return r.sel_min_s;     }},
        {"SelMax_s",    3, [](const ScalabilityGridPointResult& r) { return r.sel_max_s;     }},
        {"DGP_s",       3, [](const ScalabilityGridPointResult& r) { return r.avg_dgp_s;     }},
        {"s/Titer",     4, [](const ScalabilityGridPointResult& r) { return r.avg_titer_s;   }},
        {"PeakRSS_MiB", 1, [](const ScalabilityGridPointResult& r) { return r.peak_rss_mib;  }},
        {"RSS_MiB",     1, [](const ScalabilityGridPointResult& r) { return r.max_rss_mib;   }},
        {"Foot_MiB",    1, [](const ScalabilityGridPointResult& r) { return r.max_fp_mib;    }},
        {"X+y_MiB",     1, [](const ScalabilityGridPointResult& r) { return r.data_mib;      }},
        {"FDR",         4, [](const ScalabilityGridPointResult& r) { return r.avg_fdr;       }},
        {"TPR",         4, [](const ScalabilityGridPointResult& r) { return r.avg_tpr;       }},
        {"Avg_L",       2, [](const ScalabilityGridPointResult& r) { return r.avg_L;         }},
        {"Avg_T",       2, [](const ScalabilityGridPointResult& r) { return r.avg_T;         }},
    };
    return specs;
}

/** @brief Print the scalability results as an aligned per-solver table (rows =
 *         (n, p) grid points, columns = metrics) to console + .txt, and write a
 *         tidy long-format CSV: n,p,is_barrier,solver,num_mc,snr,metric,value. */
static void save_and_print_scalability_results(
    const std::string&                    file_stem,
    double                                tFDR,
    double                                snr,
    double                                ram_gib,
    std::size_t                           num_MC_target,
    const std::vector<ScalabilityPoint>&  grid,
    const std::vector<DemoSolverInfo>&    solvers,
    const std::map<std::string, std::map<std::pair<Eigen::Index, Eigen::Index>,
                                         ScalabilityGridPointResult>>& results)
{
    const std::string folder = DEMO_OUTPUT_DIR;
    std::filesystem::create_directories(folder);

    std::ofstream out_file(folder + file_stem + ".txt");
    auto print_dual = [&](const std::string& text) {
        std::cout << text;
        if (out_file.is_open()) out_file << text;
    };

    {
        std::ostringstream ss;
        ss << "\n"
           << "======================================================================\n"
           << "=== T-Rex Scalability Results (tFDR = " << std::fixed
           << std::setprecision(2) << tFDR << ", SNR = " << std::setprecision(1)
           << snr << ", num_MC = " << num_MC_target << ") ===\n"
           << "======================================================================\n\n"
           << "Fully memory-mapped pipeline (X/y mmap + D mmap + solver "
              "serialization).\n"
           << "Physical RAM: " << std::setprecision(0) << ram_gib
           << " GiB. Grid: n x p; the barrier row (p*) is the out-of-core "
              "point.\n\n";
        print_dual(ss.str());
    }

    const int n_w = 8, p_w = 9, met_w = 13;
    for (const auto& solver : solvers) {
        const auto& grid_res = results.at(solver.solver_name);
        print_dual(solver.solver_name + "\n");

        std::ostringstream hdr;
        hdr << std::left << std::setw(n_w) << "n" << std::setw(p_w) << "p";
        for (const auto& spec : scalability_metrics())
            hdr << std::right << std::setw(met_w) << spec.name;
        const std::string hdr_str = hdr.str();
        print_dual(hdr_str + "\n" + std::string(hdr_str.size(), '-') + "\n");

        for (const auto& pt : grid) {
            const auto it = grid_res.find({pt.n, pt.p});
            if (it == grid_res.end()) continue;
            const auto& agg = it->second;
            std::ostringstream row;
            const std::string lbl_p =
                std::to_string(pt.p) + (pt.out_of_core ? "*" : "");
            row << std::left << std::setw(n_w) << pt.n << std::setw(p_w) << lbl_p;
            for (const auto& spec : scalability_metrics())
                row << std::right << std::fixed << std::setprecision(spec.precision)
                    << std::setw(met_w) << spec.get(agg);
            print_dual(row.str() + "\n");
        }
        print_dual("\n");
    }
    print_dual("(* = out-of-core point, logical design [X|D] > physical RAM)\n\n");

    // Tidy long-format CSV: n,p,is_barrier,solver,num_mc,snr,metric,value
    std::ofstream csv(folder + file_stem + ".csv");
    if (csv.is_open()) {
        csv << "n,p,is_barrier,solver,num_mc,snr,metric,value\n"
            << std::fixed << std::setprecision(6);
        for (const auto& solver : solvers) {
            const auto& grid_res = results.at(solver.solver_name);
            for (const auto& pt : grid) {
                const auto it = grid_res.find({pt.n, pt.p});
                if (it == grid_res.end()) continue;
                // num_mc = actual cached trial count for this cell.
                for (const auto& spec : scalability_metrics()) {
                    csv << pt.n << "," << pt.p << ","
                        << (pt.out_of_core ? 1 : 0) << ","
                        << solver.solver_name << "," << it->second.num_mc << ","
                        << snr << "," << spec.name << ","
                        << spec.get(it->second) << "\n";
                }
            }
        }
        std::cout << "[Info] CSV results saved to:          "
                  << folder + file_stem + ".csv\n";
    } else {
        std::cout << "[Warning] Could not open CSV file: "
                  << folder + file_stem + ".csv\n";
    }

    if (out_file.is_open()) {
        std::cout << "[Info] Table results saved to:        "
                  << folder + file_stem + ".txt\n\n";
        out_file.close();
    } else {
        std::cout << "[Warning] Could not open output file: "
                  << folder + file_stem + ".txt\n\n";
    }
}


// ==============================================================================
// Scalability benchmark driver
// ==============================================================================

void demo_TRexSelector_scalability(const std::vector<ScalabilityPoint>& grid,
                                   const std::string& file_stem,
                                   double ram_gib,
                                   std::size_t num_MC_target,
                                   bool rnd_coef) {

    std::cout.setf(std::ios::unitbuf);

    // Clear any orphans stranded by a previous cancelled run before we start.
    sweep_mmap_orphans("startup");

    // Load the raw per-trial cache; only missing trials will be computed.
    std::filesystem::create_directories(std::string(DEMO_OUTPUT_DIR));
    const std::string raw_path =
        std::string(DEMO_OUTPUT_DIR) + file_stem + "_trials.csv";
    TrialCache cache = load_trial_cache(raw_path);
    std::size_t cached_trials = 0;
    for (const auto& [key, cell] : cache) { (void)key; cached_trials += cell.size(); }

    cdianostics::print_section_header(
        "Demo: T-Rex Selector Scalability Benchmark  |  TLARS / TOMP / TAFS");

    // Base solvers to compare (TOMP default, TAFS rho_afs = 0.3).
    const std::vector<DemoSolverInfo> solvers = {
        {SolverTypeForTRex::TLARS, "TLARS"},
        {SolverTypeForTRex::TOMP,  "TOMP"},
        {SolverTypeForTRex::TAFS,  "TAFS", 0.0, 0.3}
    };

    const double max_design = design_gib(N_VALUES.back(), P_VALUES.back());
    std::size_t n_over = 0;
    for (const auto n : N_VALUES)
        for (const auto p : P_VALUES)
            if (design_gib(n, p) > ram_gib) ++n_over;

    std::cout << "Fully memory-mapped pipeline (X/y mmap + D mmap + solver "
                 "serialization)\n"
              << "Physical RAM: " << std::fixed << std::setprecision(0)
              << ram_gib << " GiB\n"
              << "Grid: " << N_VALUES.size() << " n x " << P_VALUES.size()
              << " p = " << grid.size() << " points; largest design [X|D] = "
              << std::setprecision(1) << max_design << " GiB; " << n_over
              << " point(s) exceed RAM (out-of-core)\n"
              << "SNR = " << std::setprecision(1) << SNR << " (fixed),  s = "
              << S_TRUE << ",  tFDR = " << std::setprecision(2) << TFDR
              << ",  num_MC = " << num_MC_target
              << " sequential per (n, p, solver)\n"
              << "Trial cache: " << raw_path << "  (" << cached_trials
              << " trial(s) already present — only missing trials are computed)"
              << "\n\n";

    // results: solver -> {(n, p) -> aggregate}
    std::map<std::string, std::map<std::pair<Eigen::Index, Eigen::Index>,
                                   ScalabilityGridPointResult>> results;
    for (const auto& solver : solvers) results[solver.solver_name];

    // (n, p)-outer loop: every solver finishes the current point before the
    // next (larger) one starts, and results are re-saved after each point so a
    // crash on a bigger point preserves the finished ones.
    for (std::size_t gi = 0; gi < grid.size(); ++gi) {
        const auto& pt = grid[gi];

        cdianostics::print_section_header(
            "Grid point: n=" + std::to_string(pt.n) + ", p=" +
            std::to_string(pt.p) + (pt.out_of_core ? "   [out-of-core]" : ""));

        const double X_gib = static_cast<double>(pt.n) *
                             static_cast<double>(pt.p) * 8.0 /
                             (1024.0 * 1024.0 * 1024.0);
        std::cout << "  design [X|D] = " << std::fixed << std::setprecision(1)
                  << design_gib(pt.n, pt.p) << " GiB  (X backing "
                  << std::setprecision(2) << X_gib << " GiB + 10x-p dummies)\n\n"
                  << std::flush;

        // Per-(n, p) seeds are a stable function of (n, p, mc), shared across
        // solvers -> identical data, and independent of grid position.
        for (const auto& solver : solvers) {
            auto& cell_cache = cache[{pt.n, pt.p, solver.solver_name}];
            results[solver.solver_name][{pt.n, pt.p}] =
                run_scalability_grid_point(pt, solver, TFDR, SNR, rnd_coef,
                                           S_TRUE, cell_cache, raw_path);
        }

        // Incremental save over completed points (derived from the trial cache).
        const std::vector<ScalabilityPoint> completed(
            grid.begin(),
            grid.begin() + static_cast<std::ptrdiff_t>(gi) + 1);
        save_and_print_scalability_results(file_stem, TFDR, SNR, ram_gib,
                                           num_MC_target, completed, solvers,
                                           results);
    }
    std::cout << "\n\n";
}


// ==============================================================================
// Main
// ==============================================================================

int main(int argc, char** argv) {

    std::cout.setf(std::ios::unitbuf);
    install_cleanup_handlers();

    // Inner OpenMP (the K random experiments per select()) runs with a fixed
    // thread count (TREX_NUM_WORKERS, default 8); the outer MC loop stays
    // sequential for undistorted timing.
    omp_set_num_threads(trex_mc_workers::count());
    std::cout << "Running with " << omp_get_max_threads() << " threads\n\n";

    const double ram_gib = physical_ram_gib();

    // CLI:  [smoke] [num_MC]  (order-independent; num_MC is any integer token).
    bool smoke = false;
    std::size_t num_MC = 0;   // 0 => use the default for the chosen mode
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "smoke") {
            smoke = true;
        } else if (!a.empty() &&
                   std::all_of(a.begin(), a.end(),
                               [](unsigned char c) { return std::isdigit(c); })) {
            num_MC = static_cast<std::size_t>(std::stoul(a));
        }
        // other tokens are ignored
    }

    if (smoke) {
        const std::size_t mc = num_MC ? num_MC : std::size_t{2};
        demo_TRexSelector_scalability(make_smoke_grid(mc),
                                      "demo_trex_08_scalability_smoke",
                                      ram_gib, mc, /*rnd_coef=*/false);
    } else {
        const std::size_t mc = num_MC ? num_MC : NUM_MC;
        demo_TRexSelector_scalability(make_benchmark_grid(ram_gib, mc),
                                      "demo_trex_08_scalability_results",
                                      ram_gib, mc, /*rnd_coef=*/false);
    }

    return 0;
}
