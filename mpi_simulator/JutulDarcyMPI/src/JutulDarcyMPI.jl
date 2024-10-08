module JutulDarcyMPI
    using MPI, PartitionedArrays, HYPRE
    using Jutul
    using JutulDarcy
    using LinearAlgebra
    using ArgParse
    using ThreadPinning

    function to_nothing_or_positive(x)
        if x < 0
            x = nothing
        end
        return x
    end

    function parse_commandline()
        s = ArgParseSettings(
            description = "JutulDarcy.jl MPI compiled reservoir simulator demonstrator",
            version = "0.0.2",
            add_version = true
        )
        s.autofix_names = true

        @add_arg_table s begin
            "filename"
            help = "Either a .mat file (exported from MRST) or a .data file (industry standard reservoir simulator input) to be simulated"
            required = true
        end

        add_arg_group(s, "tolerances and model setup");
        @add_arg_table s begin
            "--tol-cnv"
                help = "CNV tolerance (Inf norm of mass balance)"
                arg_type = Float64
                default = 0.001
            "--tol-mb"
                help = "mass balance tolerance (Scaled L1 norm of mass balance)"
                arg_type = Float64
                default = 1e-7
            "--tol-dp-well"
                help = "tolerance for well connection pressure drop (in bar)"
                arg_type = Float64
                default = 1e-3
            "--multisegment-wells"
                help = "always use multisegment wells (MRST input only)"
                arg_type = Bool
                default = false
            "--number-of-steps"
                help = "number of steps in input file to solve include in solution"
                arg_type = Float64
                default = Inf
            "--can-shut-wells"
                help = "wells can be shut by solver"
                arg_type = Bool
                default = true
        end

        add_arg_group(s, "nonlinear solver");
        @add_arg_table s begin
            "--relaxation"
                help = "use relaxation/dampening for nonlinear solver"
                arg_type = Bool
                default = true
            "--tol-factor-final-iteration"
                help = "factor to relax convergence check for final iteration before cutting time-steps. Set to > 1 to loosen tolerances at last iteration."
                arg_type = Float64
                default = 1.0
            "--max-nonlinear-iterations"
                help = "maximum number of nonlinear iterations before time-step is cut"
                arg_type = Int
                default = -1
            "--max-timestep-cuts"
                help = "maximum number of time-step cuts in a single report step solve"
                arg_type = Int
                default = 20
            "--ds-max"
                help = "maximum change in saturations during iteration"
                arg_type = Float64
                default = 0.2
            "--dz-max"
                help = "maximum change in compositions during iteration"
                arg_type = Float64
                default = 0.2
            "--dp-max-abs"
                help = "maximum absolute change in pressure during iteration (units of bar)"
                arg_type = Float64
                default = Inf
            "--dp-max-rel"
                help = "maximum change in relative pressure during iteration"
                arg_type = Float64
                default = 0.2
            "--presolve-wells"
                help = "solve well system before each Newton iteration"
                arg_type = Bool
                default = false
        end

        add_arg_group(s, "linear solver");
        @add_arg_table s begin
            "--rtol"
                help = "relative tolerance for linear solver"
                arg_type = Float64
                default = 5e-3
            "--backend"
                help = "matrix format to use (either csr or csc)"
                arg_type = Symbol
                default = :csr
                range_tester = x -> x in (:csr, :csc)
            "--linear-solver"
                help = "linear solver must be (bicgstab or gmres)"
                arg_type = Symbol
                default = :gmres
                range_tester = x -> x in (:bicgstab, :gmres)
            "--precond"
                help = "preconditioner to use for linear solver"
                arg_type = Symbol
                default = :cpr
                range_tester = x -> x in (:cpr, :ilu0, :jacobi, :spai0)
        end

        add_arg_group(s, "timestepping");
        @add_arg_table s begin
            "--timesteps"
                help = "type of timestepping (auto, iteration, none)"
                arg_type = Symbol
                default = :auto
                range_tester = x -> x in (:auto, :iteration, :none)
            "--timesteps-target-ds"
                help = "target saturation change for timestepping"
                arg_type = Float64
                default = Inf
            "--timesteps-target-iterations"
                help = "target iterations for timestepping (default: 8 for method=newton, 3 for method=nldd)"
                arg_type = Int
                default = -1
            "--max-timestep"
                help = "maximum time-step in days"
                arg_type = Float64
                default = 90.0
            "--initial-timestep"
                help = "initial time-step in days"
                arg_type = Float64
                default = 1.0
        end

        add_arg_group(s, "output and printing");
        @add_arg_table s begin
            "--info-level"
                help = "level out output. Set to -1 for no output."
                arg_type = Int
                default = 1
            "--report-level"
                help = "amount of convergence/numerical statistics to store. Set to 1 for addditional output."
                arg_type = Int
                default = 0
            "--verbose"
                help = "extra output from the app itself and the parser. For simulation convergence reporting, see --info-level"
                action = :store_true
            "--no-header"
                help = "do not print ASCII header with version numbers and URL"
                action = :store_true
            "--output-path"
                help = "path where output results are to be written. A random temporary folder will be created if not provided."
                arg_type = String
                default = ""
            "--output-substates"
                help = "include all substeps in output states (can drastically increase file sizes!)"
                action = :store_true
            "--presolve-one-step"
                help = "solve one step before running the simulation to ensure everything is compiled. This is generally slower, but makes timings more accurate if a model contains uncompiled features."
                action = :store_true
            "--print-wells"
                help = "print tables of well results after simulation is complete."
                action = :store_true
        end

        add_arg_group(s, "numerical scheme");
        @add_arg_table s begin
            "--method"
                help = "type of method to use: Can be newton or nldd"
                arg_type = String
                default = "newton"
            "--nldd-number-of-subdomains"
                help = "number of domains (per process) to use in NLDD (not compatible with --nldd-cells-per-subdomain)"
                arg_type = Int
                default = -1
            "--nldd-cells-per-subdomain"
                help = "number of cells per subdomain to use in NLDD (not compatible with --nldd-number-of-subdomains)"
                arg_type = Int
                default = -1
            "--nldd-solve-tol-pressure"
                help = "tolerance for maximum value of change in pressure before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-pressure-mean"
                help = "tolerance for average value of change in pressure before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-saturations"
                help = "tolerance for maximum value of change in saturations before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-saturations-mean"
                help = "tolerance for average value of change in saturations before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-densities"
                help = "tolerance for maximum value of change in densities before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-densities-mean"
                help = "tolerance for average value of change in densities before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-mobility"
                help = "tolerance for maximum value of change in mobility before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-mobility-mean"
                help = "tolerance for average value of change in mobility before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-composition"
                help = "tolerance for maximum value of change in composition before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-composition-mean"
                help = "tolerance for average value of change in composition before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-phase-mass-fractions"
                help = "tolerance for maximum value of change in phase-mass-fractions before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-solve-tol-phase-mass-fractions-mean"
                help = "tolerance for average value of change in phase-mass-fractions before a local subdomain is solved"
                arg_type = Float64
                default = -1.0
            "--nldd-subdomain-failure-cuts"
                help = "cut a timestep if any subdomain solver fails"
                arg_type = Bool
                default = false
            "--nldd-always-solve-wells"
                help = "always solve domains containing wells"
                arg_type = Bool
                default = false
            "--nldd-solve-tol-first-newton"
                help = "use newton for first iteration when solve tols in changes are enabled"
                arg_type = Bool
                default = true
            "--nldd-mpi-sync-after-solve"
                action = :store_true
            "--nldd-subdomain-precond"
                help = "preconditioner to use for linear solver in NLDD subdomains"
                arg_type = Symbol
                default = :cpr
                range_tester = x -> x in (:cpr, :ilu0, :jacobi, :spai0)
        end

        add_arg_group(s, "HPC configuration");
        @add_arg_table s begin
            "--partitioner-conn-type"
                help = "Connection type for partitioned (trans, unit)"
                arg_type = Symbol
                default = :unit
            "--pin-threads"
                help = "Pin threads (cputhreads, cores, sockets, compact, numa, random, current, firstn, affinitymask)"
                arg_type = Symbol
                default = :none
            "--blas-threads"
                help = "Number of threads for BLAS"
                arg_type = Int
                default = 1
        end

        return parse_args(s)
    end

    function julia_main()::Cint
        MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        is_main = rank == 0

        args = parse_commandline()
        if isnothing(args)
            return 0
        end
        process_is_verbose = rank == 0 && args["verbose"]
        function verbose_print(arg...; kwarg...)
            if process_is_verbose
                Jutul.jutul_message(arg...; kwarg...)
            end
        end
        if rank == 0 && !args["no_header"]
            print_big_logo()
        end
        if process_is_verbose
            verbose_print("JutulDarcyMPI", "Case is $(args["filename"]).\nSuccessfully parsed $(length(args)) keyword arguments:")
            print_arg = Jutul.OrderedCollections.OrderedDict()
            argkeys = sort(collect(keys(args)))
            for k in argkeys
                if k == "filename"
                    continue
                end
                print_arg[k] = args[k]
            end
            Jutul.pretty_table(print_arg, crop = :none, show_header = false, alignment = :l)
        end
        pt = args["pin_threads"]
        if pt != :none
            verbose_print("JutulDarcyMPI", "Pinning to threads: $pt")
            mpi_pinthreads(pt)
        end

        if args["relaxation"]
            r = SimpleRelaxation()
        else
            r = NoRelaxation()
        end

        if args["multisegment_wells"]
            w = :ms
        else
            w = :simple
        end

        extra_arg = Dict{Symbol, Any}()

        method = Symbol(lowercase(args["method"]))
        if method == :nldd
            verbose_print("Numerical scheme", "NLDD solver selected.")
            cells_per_block = args["nldd_cells_per_subdomain"]
            if cells_per_block < 1
                cells_per_block = missing
            end
            no_blocks = args["nldd_number_of_subdomains"]
            if no_blocks < 1
                no_blocks = missing
            end
            extra_arg[:nldd_arg] = Dict(
                :cells_per_block => cells_per_block,
                :no_blocks => no_blocks,
                :mpi_sync_after_solve => args["nldd_mpi_sync_after_solve"]
            )
            extra_arg[:solve_tol_pressure] = to_nothing_or_positive(args["nldd_solve_tol_pressure"])
            extra_arg[:solve_tol_pressure_mean] = to_nothing_or_positive(args["nldd_solve_tol_pressure_mean"])
            extra_arg[:solve_tol_saturations] = to_nothing_or_positive(args["nldd_solve_tol_saturations"])
            extra_arg[:solve_tol_saturations_mean] = to_nothing_or_positive(args["nldd_solve_tol_saturations_mean"])
            extra_arg[:solve_tol_densities] = to_nothing_or_positive(args["nldd_solve_tol_densities"])
            extra_arg[:solve_tol_densities_mean] = to_nothing_or_positive(args["nldd_solve_tol_densities_mean"])
            extra_arg[:solve_tol_mobility] = to_nothing_or_positive(args["nldd_solve_tol_mobility"])
            extra_arg[:solve_tol_mobility_mean] = to_nothing_or_positive(args["nldd_solve_tol_mobility_mean"])
            extra_arg[:solve_tol_composition] = to_nothing_or_positive(args["nldd_solve_tol_composition"])
            extra_arg[:solve_tol_composition_mean] = to_nothing_or_positive(args["nldd_solve_tol_composition_mean"])
            extra_arg[:solve_tol_phase_mass_fractions] = to_nothing_or_positive(args["nldd_solve_tol_phase_mass_fractions"])
            extra_arg[:solve_tol_phase_mass_fractions_mean] = to_nothing_or_positive(args["nldd_solve_tol_phase_mass_fractions_mean"])
            extra_arg[:solve_tol_first_newton] = args["nldd_solve_tol_first_newton"]
            extra_arg[:subdomain_failure_cuts] = args["nldd_subdomain_failure_cuts"]
            extra_arg[:always_solve_wells] = args["nldd_always_solve_wells"]
            extra_arg[:subdomain_precond] = args["nldd_subdomain_precond"]

            target_its = 3
            max_its = 10
        else
            verbose_print("Numerical scheme", "Newton solver selected.")
            target_its = 8
            max_its = 15
        end
        if args["timesteps_target_iterations"] < 1
            args["timesteps_target_iterations"] = target_its
        end
        if args["max_nonlinear_iterations"] < 1
            args["max_nonlinear_iterations"] = max_its
        end

        num_blas_threads = args["blas_threads"]
        BLAS.set_num_threads(num_blas_threads)
        ENV["OPENBLAS_NUM_THREADS"] = num_blas_threads

        pth = args["filename"]
        @assert isfile(pth) "Input file $pth must exist."

        basepath, ext = splitext(pth)
        folder_pth, name = splitdir(basepath)
        ext = lowercase(ext)
        verbose_print("IO", "Reading case from $pth...")
        t_setup = @elapsed if ext == ".mat"
            case, = setup_case_from_mrst(pth,
                backend = args["backend"],
                wells = w,
                ds_max = args["ds_max"],
                dz_max = args["dz_max"],
                dp_max_abs = args["dp_max_abs"],
                dp_max_rel = args["dp_max_rel"],
                split_wells = true
            )
        else
            @assert ext == ".data" "File must have either extension .mat (for MRST export) or .data (for industry standard input format). Was: $ext"
            case = JutulDarcy.setup_case_from_data_file(pth,
                backend = args["backend"],
                split_wells = true,
                ds_max = args["ds_max"],
                dz_max = args["dz_max"],
                dp_max_abs = args["dp_max_abs"],
                dp_max_rel = args["dp_max_rel"],
                can_shut_wells = args["can_shut_wells"],
                parse_arg = (verbose = args["verbose"] && is_main, silent = !is_main)
            )
        end
        verbose_print("IO", "Case $name set up in $(Jutul.autoformat_time(t_setup)).")
        nstep = args["number_of_steps"]
        nstep_in_case = length(case.dt)
        if nstep < nstep_in_case && nstep > 0
            nstep = Int(nstep)
            verbose_print("IO", "Limiting case to $nstep out of $nstep_in_case steps (--number-of-steps=$nstep).")
            nstep = min(nstep_in_case, nstep)
            case = case[nstep]
        end
        outpth = args["output_path"]
        if outpth == ""
            outpth = joinpath(folder_pth, "$(name)_jutul_output")
            verbose_print("Output", "--output-path not provided, writing output to $outpth")
        end
        if is_main
            mkpath(outpth)
        end
        MPI.Barrier(comm)
        if args["verbose"]
            print("Starting simulation at rank $(rank+1) of $(MPI.Comm_size(comm))\n")
        end
        kwarg = Dict(
            :mode => :mpi,
            :parray_arg => (
                conn = args["partitioner_conn_type"],
            ),
            :method => method,
            :max_timestep_cuts => args["max_timestep_cuts"],
            :linear_solver_arg => (
                rtol = args["rtol"],
                solver = args["linear_solver"]
            ),
            :precond => args["precond"],
            :max_nonlinear_iterations => args["max_nonlinear_iterations"],
            :target_its => args["timesteps_target_iterations"],
            :timesteps => args["timesteps"],
            :target_ds => args["timesteps_target_ds"],
            :max_dt => args["max_timestep"]*si_unit(:day),
            :initial_dt => args["initial_timestep"]*si_unit(:day),
            :tol_cnv => args["tol_cnv"],
            :tol_mb => args["tol_mb"],
            :tol_dp_well => args["tol_dp_well"],
            :tol_factor_final_iteration => args["tol_factor_final_iteration"],
            :relaxation => r,
            :report_level => args["report_level"],
            :output_substates => args["output_substates"],
            :failure_cuts_timestep => true,
        )
        if args["presolve_one_step"]
            subcase = case[1]
            subcase.dt[1] *= 0.01
            try
                result = simulate_reservoir(subcase;
                output_path = mktempdir(),
                info_level = -1,
                kwarg...,
                extra_arg...
            )
            catch
                # No need to do anything, this is just for accurate timings and could fail.
            end
        end
        wells, states = simulate_reservoir(case;
            output_path = outpth,
            info_level = args["info_level"],
            kwarg...,
            extra_arg...
        )
        if rank == 0 && args["print_wells"]
            wells()
        end
        MPI.Barrier(comm)
        return 0 # if things finished successfully
    end
    function print_big_logo()
        println("")
        txt =
"       █████             █████               ████  ██████████
      ░░███             ░░███               ░░███ ░░███░░░░███
       ░███  █████ ████ ███████   █████ ████ ░███  ░███   ░░███  ██████   ████████   ██████  █████ ████
       ░███ ░░███ ░███ ░░░███░   ░░███ ░███  ░███  ░███    ░███ ░░░░░███ ░░███░░███ ███░░███░░███ ░███
       ░███  ░███ ░███   ░███     ░███ ░███  ░███  ░███    ░███  ███████  ░███ ░░░ ░███ ░░░  ░███ ░███
 ███   ░███  ░███ ░███   ░███ ███ ░███ ░███  ░███  ░███    ███  ███░░███  ░███     ░███  ███ ░███ ░███
░░████████   ░░████████  ░░█████  ░░████████ █████ ██████████  ░░████████ █████    ░░██████  ░░███████
 ░░░░░░░░     ░░░░░░░░    ░░░░░    ░░░░░░░░ ░░░░░ ░░░░░░░░░░    ░░░░░░░░ ░░░░░      ░░░░░░    ░░░░░███
                                                                                              ███ ░███
    https://github.com/sintefmath/JutulDarcy.jl                                              ░░██████
                                                                                              ░░░░░░
"
        # This doesn't get called properly when compiled? Always print.
        # wdth = displaysize(stdout)[2]
        # if wdth >= 104
        if Sys.iswindows()
            println("JutulDarcy.jl compiled simulator")
            println("https://github.com/sintefmath/JutulDarcy.jl")
        else
            print(txt)
        end
        jver = pkgversion(Jutul)
        jdver = pkgversion(JutulDarcy)
        gver = pkgversion(JutulDarcy.GeoEnergyIO)
        mcfver = pkgversion(JutulDarcy.MultiComponentFlash)
        jutul_message("Packages", "JutulDarcy@$jdver, Jutul@$jver, GeoEnergyIO@$gver, MultiComponentFlash@$mcfver")
    end
end # module
# using MPI
# MPI.install_mpiexecjl()
# mpiexec()
