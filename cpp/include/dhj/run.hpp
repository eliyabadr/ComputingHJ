// run.hpp - the shared main() body; each of the three binaries picks a Mode.
#pragma once

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

#include "args.hpp"
#include "cell.hpp"
#include "checkpoint.hpp"
#include "env.hpp"
#include "pool.hpp"
#include "reach.hpp"
#include "refine.hpp"
#include "timing.hpp"
#include "vi.hpp"

namespace dhj {

inline std::unique_ptr<Environment> make_env(const Args& args) {
#if DHJ_DIM == 3
    if (args.dynamics == Dynamics::DubinsSpeed)
        throw std::invalid_argument("dubins_speed is 4D; rebuild with DHJ_DIM=4 (make dim4)");
    if (args.dynamics == Dynamics::Dubins)
        return std::unique_ptr<Environment>(
            new DubinsCarEnvironment(args.velocity, args.dt, args.tau, /*obstacle_radius=*/1.3));
    return std::unique_ptr<Environment>(
        new EvasionEnvironment(args.velocity, args.dt, args.tau, /*obstacle_radius=*/1.0));
#else
    if (args.dynamics != Dynamics::DubinsSpeed)
        throw std::invalid_argument("the DHJ_DIM=4 build only provides --dynamics dubins_speed");
    return std::unique_ptr<Environment>(
        new DubinsSpeedEnvironment(args.dt, args.tau, /*obstacle_radius=*/1.3));
#endif
}

inline void run_algorithm_1(Mode mode, const Args& args, const Environment& env, ThreadPool& pool,
                            Timers& timers) {
    std::printf("%s\n", std::string(70, '=').c_str());
    std::printf("ALGORITHM 1: Discretization Routine\n");
    std::printf("%s\n", std::string(70, '=').c_str());
    std::printf("Environment: %s\n", env.name());
    std::printf("Initializing grid with resolution %d^%d...\n", args.resolution, static_cast<int>(kDim));

    const double grid0 = now_seconds();
    CellTree tree(env.state_bounds(), args.resolution);
    timers.grid_build += now_seconds() - grid0;
    timers.spatial_index += tree.index_build_seconds();

    GronwallReachabilityAnalyzer reach(env);
    const std::string out_dir = args.out_dir.empty() ? algorithm1_output_dir(args, env) : args.out_dir;
    SafetyValueIterator vi(mode, args, env, tree, reach, out_dir, pool, timers);

    const double t0 = now_seconds();
    vi.value_iteration(args.iterations, args.tolerance, args.plot_freq, args.conservative,
                       args.delta_max);
    std::printf("\nALGORITHM 1 COMPLETE\n");
    std::printf("Total time: %.2f seconds\n", now_seconds() - t0);
    std::printf("Results saved to: %s/\n", out_dir.c_str());
    timers.final_cells = tree.num_leaves();
}

inline void run_algorithm_2(Mode mode, const Args& args, const Environment& env, ThreadPool& pool,
                            Timers& timers) {
    std::printf("%s\n", std::string(70, '=').c_str());
    std::printf("ALGORITHM 2/3: Adaptive Refinement\n");
    std::printf("%s\n", std::string(70, '=').c_str());
    std::printf("Environment: %s\n", env.name());

    // Resolved here as well as in AdaptiveRefinement, because `--resume auto`
    // has to find a checkpoint before the tree it restores into exists.
    const std::string out_dir =
        args.out_dir.empty() ? refine_output_dir(mode, args, env) : args.out_dir;
    const std::string ckpt_dir =
        args.checkpoint_dir.empty() ? out_dir + "/checkpoints" : args.checkpoint_dir;

    std::unique_ptr<CellTree> tree;
    CheckpointMeta resume_meta;
    bool resuming = false;

    if (!args.resume.empty()) {
        std::string path = args.resume;
        if (path == "auto") {
            path = latest_checkpoint(ckpt_dir);
            if (path.empty())
                throw std::runtime_error("--resume auto: no checkpoint found in " + ckpt_dir);
            std::printf("--resume auto selected %s\n", path.c_str());
        }
        const double load0 = now_seconds();
        resume_meta = load_checkpoint(path, mode, tree);
        timers.grid_build += now_seconds() - load0;
        timers.spatial_index += tree->index_build_seconds();
        resuming = true;

        if (resume_meta.initial_resolution != static_cast<std::uint32_t>(args.initial_resolution))
            std::fprintf(stderr,
                         "  warning: checkpoint was made with --initial-resolution %u but this run "
                         "passes %d; the checkpoint's grid wins\n",
                         resume_meta.initial_resolution, args.initial_resolution);
    } else {
        const double grid0 = now_seconds();
        tree.reset(new CellTree(env.state_bounds(), args.initial_resolution));
        timers.grid_build += now_seconds() - grid0;
        timers.spatial_index += tree->index_build_seconds();
    }

    GronwallReachabilityAnalyzer reach(env);
    AdaptiveRefinement adaptive(mode, args, env, *tree, reach, pool, timers);

    const double t0 = now_seconds();
    adaptive.refine(args.epsilon, args.refinements, args.vi_iterations,
                    resuming ? &resume_meta : nullptr);
    std::printf("ALGORITHM 2 COMPLETE\n");
    std::printf("Total time: %.2f seconds\n", now_seconds() - t0);
    timers.final_cells = tree->num_leaves();
}

inline int run_main(Mode mode, int argc, char** argv) {
    const Args args = parse_args(argc, argv, mode);

    std::printf("OPTIMIZED SAFETY VALUE FUNCTION - C++ PORT OF %s.py\n", mode_script(mode));
    std::printf("Algorithm: %d, Workers: %d\n", args.algorithm, args.workers);
    std::printf("Dynamics: %s\n", args.dynamics_name.c_str());
    std::printf("Conservative mode: %s\n", py_bool(args.conservative).c_str());

    std::unique_ptr<Environment> env;
    try {
        env = make_env(args);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    std::printf("  Time discretization: %d steps of dt=%g s over tau=%g s\n", env->n_steps(),
                env->dt(), env->tau());

    ThreadPool pool(static_cast<std::size_t>(args.workers));
    Timers timers;
    const double t0 = now_seconds();

    try {
        if (args.algorithm == 1) run_algorithm_1(mode, args, *env, pool, timers);
        else run_algorithm_2(mode, args, *env, pool, timers);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }

    timers.total = now_seconds() - t0;
    timers.report(mode_script(mode));
    return 0;
}

}  // namespace dhj
