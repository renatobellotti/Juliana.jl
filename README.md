![](juliana_logo.svg)

Implementation of a spot (and later probably beam angle) optimiser in Julia.

## Installation

1. Get access to the repo.
2. Clone the repo:
   ```bash
   git clone git@gitlab.psi.ch:bellotti_r/differentiable-planning.git
   ```
3. Install Julia using [the official instructions](https://julialang.org/downloads/).
4. Open a Julia command line from within the JulianA repo by typing ``julia`` in a terminal. Then type the following and hit enter to setup a development version of the ``Juliana`` package (i. e. local edits to the files will be visible when you import the package in a fresh session). The following command will download all dependencies and setup the package for local development. The bracket ``]`` command will switch to the package manager mode of the Julia REPL, which makes it easier to deal with packages and environments. 
   ```julia
   ] dev .
   ```
5. Test your installation in a Julia REPL. The first import will precompile ``Juliana``, so please be a bit patient.
   ```julia
   using Juliana
   ```
6. (Optional) Install a Jupyter kernel for Julia. Open a Julia REPL in the repo directory, then type the following:
   ```julia
   using Pkg
   Pkg.add("IJulia")
   using IJulia
   notebook()
   ```
   The ``notebook()`` command will start a local Jupyter server and open a browser so you can start editing. The first two lines are needed only once for the installation of Jupyter itself.


## TPS integration

This package is built such that it can support any TPS,
in principle. To add a TPS, you have to create a type,
e. g. ``TpsType``, that is a subtype to ``Juliana.AbstractTps``
and implement the following methods for this type:

```julia
place_spots(tps::TpsType,
            ct::ScalarGrid,
            target_structure::Structure,
            preabsorber::String,
            gantry_angles,
            couch_angles,
            nozzle_extractions)

calculate_dose(tps::TpsType,
               resolution::Real,
               ct::ScalarGrid,
               plan::TreatmentPlan)

calculate_Dij_matrix(tps::TpsType,
                     ct::ScalarGrid,
                     plan::TreatmentPlan,
                     optimisation_points)
```

Currently implemented TPS:

- Fiona standalone<br />
  Needs proprietary files.
  A path to these files is given to the constructor.
  Two paths are needed:
  To a directory that contains all the machine parameters
  and to the JAR file that contains the actual code of
  Fiona standalone.


## Test coverage

You can get an overview over test coverage by running the following command from within the project root directory:

```bash
julia --project -E 'using LocalCoverage; coverage = generate_coverage("Juliana"; run_test = true)'
```
