FROM julia:latest
WORKDIR /env

COPY . .

RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
ENTRYPOINT [ "julia", "--project=/env", "-e", "import MCMRSimulator.CLI: run_main; run_main()", "--" ]