"""
Defines the command line interface to MCMRSimulator.jl (`mcmr`).

Functions:
- [`MCMRSimulator.run_main`](@ref MCMRSimulator.CLI.run_main)
"""
module CLI
include("run.jl")
include("geometry.jl")

import ArgParse
import Markdown
import Pkg
import .Run
import .Geometry


"""
    install_cli()

Installs the command line interface for 
"""
function install_cli(;
    command::AbstractString="mcmr",
    destdir::AbstractString=joinpath(DEPOT_PATH[1], "bin"),
    force::Bool=false
    )
    # The code in this function is for a major part taken from the jpkg.jl package (/Applications/Julia-1.9.app/Contents/Resources/julia/bin/).
    julia_executable = joinpath(Sys.BINDIR, "julia")
    if Sys.iswindows()
        @warn "Installing the CLI is untested on Windows. The installer will continue, but don't be surprised if it does not work."
    end
    full_destdir = abspath(expanduser(destdir))
    exec = joinpath(full_destdir, command)

    if !isdir(full_destdir)
        if destdir === joinpath(DEPOT_PATH[1], "bin")
            mkdir(destdir)
        else
            error("Destination directory $destdir does not exist. Please create the directory or set `MCMRSimulator.install_cli(destdir=)` to a valid directory.")
        end
    end

    previous_install = Sys.which(command)
    if !isnothing(previous_install)
        if samefile(previous_install, exec)
            if !force
                error("$command is already installed at $(destdir). Use `MCMRSimulator.install_cli(force=true)` to overwrite.")
            end
        else
            error("Command $command is already available at $previous_install. Continuing with this install would conflict with this. Please remove previous installations or change the name of the executable (`MCMRSimulator.install_cli(command=...)`)")
        end
    elseif ispath(exec)
        error("File/directory $exec already exists and does not appear to be an MCMR executable. I'm not going to overwrite it. Please delete it if you want to install MCMR here.")
    end


    open(exec, "w") do f
        write(f, 
        """
        #!/usr/bin/env bash

        $(julia_executable) --startup-file=no --project=$(Pkg.project().path) -e 'import MCMRSimulator.CLI: run_main; run_main()' -- "\$@"
        """
        )
        
    end
    chmod(exec, 0o0100775) # equivalent to -rwxrwxr-x (chmod +x exec)
    
    @info "Installed MCMRSimulator command line tool to '$exec'"

    if isnothing(Sys.which(command))
        @warn "The directory that the MCMRSimulator CLI was installed too does not appear to be on the PATH. Either use the complete executable ('$(Base.contractuser(exec))') or add '$(full_destdir)' to your PATH variable."
    end
end

"""
    run_main(args=ARGS; kwargs...)

Main run command for the command line interface.

Any keyword arguments are passed on to `ArgParse.ArgParseSettings`.
"""
function run_main(args=ARGS; kwargs...)
    if length(args) == 0
        println(stderr, "No mcmr command given.\n")
    else
        if args[1] == "run"
            return Run.run_main(args[2:end]; kwargs...)
        elseif args[1] == "geometry"
            return Geometry.run_main(args[2:end]; kwargs...)
        else
            println(stderr, "Invalid mcmr command $(args[1]) given.\n")
        end
    end
    println(stderr, "usage: mcmr {run/geometry}")
    return Cint(1)
end

function cmdline_debug_handler(settings::ArgParse.ArgParseSettings, err, err_code::Int = 1)
    println(stderr, err.text)
    println(stderr, ArgParse.usage_string(settings))
    error()
end

"""
    run_main_test(cmd::AbstractString)

Returns the stdout and stderr produced by running [`run_main`](@ref) on the given command.
This is used for testing of the command line interface
"""
function run_main_test(cmd::AbstractString)
    original_stderr = stderr
    err_rd, _ = redirect_stderr()
    err_text = @async read(err_rd, String)
    original_stdout = stdout
    out_rd, _ = redirect_stdout()
    out_text = @async read(out_rd, String)
    prev_interactive = Base.is_interactive
    Base.is_interactive = false

    try
        try
            run_main(split(cmd), exc_handler=cmdline_debug_handler, exit_after_help=false)
        catch err
            if hasproperty(err, :msg) && length(err.msg) > 0
                rethrow(err)
            end
        end
    finally
        close(err_rd)
        redirect_stderr(original_stderr)
        close(out_rd)
        redirect_stdout(original_stdout)
        Base.is_interactive = prev_interactive
    end

    return (fetch(out_text), fetch(err_text))
end

"""
    run_main_test(cmd::AbstractString)

Returns markdown with the stdout and stderr produced by running [`run_main`](@ref) on the given command.
This is used for testing of the command line interface
"""
function run_main_docs(cmd::AbstractString)
    (out, err) = run_main_test(cmd)
    text = ""
    if length(out) > 0
        text = text * "```stdout\n$(out)\n```\n"
    end
    if length(err) > 0
        text = text * "```stderr\n$(err)\n```\n"
    end
    return Markdown.parse(text)
end

end