# [Installation](@id installation)
MCMRSimulator is an application written in the [Julia](https://julialang.org) language.
You can run simulations either directly from the Julia REPL, in a [Jupyter notebook](@ref jupyter_install), or using the command line interface.
## Installing the simulator for a specific project
1. First install julia from the [official website](https://julialang.org/downloads/).
2. Create a directory for the project for which you are going to use the simulator. We will install MCMRSimulator in isolation just for this project. This ensures that if we install a newer version of the simulator for another project in the future, it will not interfere with the reproducibility of the results of this project. We will refer to this newly created project directory below as "<project_dir>".
2. Start the julia REPL in a terminal (`$ julia --project=<project_dir>`). The `--project` flag ensures that we create and activate a Julia environemnt just for this project. This will create a "Project.toml" and "Manifest.toml" inside the directory specifying the installed packages (such as the MCMR simulator).
3. Enter the package manager by pressing "]"
   - First install the required MRIBuilder.jl using `pkg> add https://git.fmrib.ox.ac.uk/ndcn0236/mribuilder.jl.git`.
   - Then install MCMRSimulator.jl using `pkg> add https://git.fmrib.ox.ac.uk/ndcn0236/mcmrsimulator.jl.git{install_version}`.
   - (Optional) Install one of the [Makie backends](https://makie.juliaplots.org/stable/documentation/backends/) for plotting (e.g., `pkg> add CairoMakie`).
   - (Optional) If you want to use a Jupyter notebook for this project, you will also have to install an `IJulia` kernel. You can find instructions to do so [below](@ref jupyter_install).
   - Press "\[backspace\]" to leave the package manager.
4. (Optional) To install the MCMRSimulator command line interface (CLI) run the following in the main julia REPL:
   `using MCMRSimulator; MCMRSimulator.install_cli(destdir="...")`, where `destdir` is the target directory for the executable (called `mcmr` by default).
   Ensure that the `destdir` is a directory that is in your $PATH variable.

## Running MCMRSimulator
After this installation process, you can run MCMRSimulator in one of the following ways:
- *Julia REPL*: Start the REPL in a terminal by typing `$ julia --project=<project_dir>`. Afterwards type `using MCMRSimulator` to import the simulator. You can now follow the steps in the [MCMRSimulator tutorial using Julia](@ref tutorial_julia).
- *Jupyter notebook*: Make sure that you install `IJulia` using the instructions [below](@ref jupyter_install). This will allow you to start a notebook in jupyter running in Julia. Within this notebook, you can follow the steps in the [MCMRSimulator tutorial using Julia](@ref tutorial_julia).
- *Command line interface*: If you followed the instructions in step 4 above, you can run the MCMRSimulator command line interface simply by typing `mcmr` in the terminal. If it is not working, you might want to redo step 4 above and pay attention to any warning messages. With this alias set up, you can now follow [the command line tutorial](@ref tutorial_cli).

## Updating MCMRSimulator
First check the [CHANGELOG](https://git.fmrib.ox.ac.uk/ndcn0236/mcmrsimulator.jl/-/blob/main/CHANGELOG.md) to find a list of changes since your current version.
If you decide to update:
1. Start the julia REPL again in a terminal (`$ julia --project=<project_dir>`)
2. Enter the package manager by pressing "]"
3. Update all installed packages using by typing `update` and pressing enter (`pkg> update`).

## Sharing your MCMRSimulator installation
To share the exact environment used by your installation of MCMRSimulator, simply go to the `<project_dir>` directory and locate the files named "Project.toml" and "Manifest.toml". Transfer these files to any other computer, to ensure that they install the exact same version of all Julia packages used (see https://pkgdocs.julialang.org/v1/environments/ for more details).

## [Running MCMRSimulator in a Jupyter notebook](@id jupyter_install)

### Installing the Julia kernel
You only have to run the following once:
```bash
julia -e "import Pkg; Pkg.add(\"IJulia\"); Pkg.build(\"IJulia\")" 
```
This line installs and builds `IJulia` in the main, global environment.
This will create a new jupyter kernel just for Julia.
This global environment is available within any local environments,
so we only have to do this once across any number of Julia projects.

To test whether this has worked, start a `Jupyter notebook` (or `Jupyter lab`).
If you have FSL installed, you can do so using `fslpython -m notebook`.
Within the `Jupyter notebook` interface press "New".
The Julia kernel should show up here with its version number.
![](jupyter_julia_kernel.png)

### Using the simulator with this Julia kernel
If you followed the steps in the previous section, the Julia kernel will be installed in its own environment separate from the simulator.
To use the simulator using this kernel, we simply start a notebook and then switch to an environment that has the MCMR simulator installed (i.e., created as described [above](@ref installation)).
You can do this by including a block at the top of the notebook with the following code:
```julia
import Pkg
Pkg.activate("<project dir>")
```
where `"<project_dir>"` is the directory where we installed the simulator.

