module JutulDarcyMPI
    using MPI, PartitionedArrays, HYPRE
    using Jutul
    using JutulDarcy
    using LinearAlgebra
    using ArgParse

    function parse_commandline()
        s = ArgParseSettings(
            description = "JutulDarcy.jl MPI compiled reservoir simulator demonstrator",
            version = "0.0.1",
            add_version = true
        )
        s.autofix_names = true

        @add_arg_table s begin
            "filename"
            help = ".mat file exported from MRST to be simulated"
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
                help = "always use multisegment wells"
                arg_type = Bool
                default = false
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
                default = 15
            "--max-timestep-cuts"
                help = "maximum number of time-step cuts in a single report step solve"
                arg_type = Int
                default = 20
        end

        add_arg_group(s, "linear solver");
        @add_arg_table s begin
            "--rtol"
                help = "relative tolerance for linear solver"
                arg_type = Float64
                default = 1e-3
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
                help = "target iterations for timestepping"
                arg_type = Int
                default = 8
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
            "--verbose"
                help = "extra output from the app itself and the parser. For simulation convergence reporting, see --info-level"
                arg_type = Bool
                default = false
            "--output-path"
                help = "path where output results are to be written. A random temporary folder will be created if not provided."
                default = ""
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
        function verbose_print(arg...; kwarg...)
            if rank == 0 && args["verbose"]
                Jutul.jutul_message(arg...; kwarg...)
            end
        end
        verbose_print("Jutul", "Successfully parsed input arguments:")
        i = 1
        for (arg,val) in args
            verbose_print("$i", "$arg =>  $val [$(typeof(val))]", color = :green)
            i += 1
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
        BLAS.set_num_threads(1)

        pth = args["filename"]
        @assert isfile(pth) "Path $pth must be a valid .mat or .data file"

        basepath, ext = splitext(pth)
        folder_pth, name = splitdir(basepath)
        ext = lowercase(ext)
        if ext == ".mat"
            verbose_print("IO", "Reading case from $pth...")
            t_setup = @elapsed case, = setup_case_from_mrst(pth, backend = args["backend"], wells = w, split_wells = true)
            verbose_print("IO", "Case $name set up in $t_setup s.")
        else
            @assert ext == ".data"
            case = JutulDarcy.setup_case_from_data_file(pth, 
                backend = args["backend"],
                split_wells = true,
                parse_arg = (verbose = args["verbose"],)
            )
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
        result = simulate_reservoir(case,
            mode = :mpi,
            output_path = outpth,
            info_level = args["info_level"],
            max_timestep_cuts = args["max_timestep_cuts"],
            linear_solver_arg = (rtol = args["rtol"], solver = args["linear_solver"]),
            precond = args["precond"],
            max_nonlinear_iterations = args["max_nonlinear_iterations"],
            target_its = args["timesteps_target_iterations"],
            timesteps = args["timesteps"],
            target_ds = args["timesteps_target_ds"],
            max_dt = args["max_timestep"]*si_unit(:day),
            initial_dt = args["initial_timestep"]*si_unit(:day),
            tol_cnv = args["tol_cnv"],
            tol_mb = args["tol_mb"],
            tol_dp_well = args["tol_dp_well"],
            tol_factor_final_iteration = args["tol_factor_final_iteration"],
            relaxation = r,
            failure_cuts_timestep = true
        )
        if MPI.Comm_size(comm) > 1
            # Skip this for single process to avoid breaking interactive
            # development by calling finalize.
            MPI.Finalize()
        end
        return 0 # if things finished successfully
    end
end # module
# using MPI
# MPI.install_mpiexecjl()
# mpiexec()
