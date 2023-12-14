![](juliana_logo.svg)

This repository contains the code for the paper [JulianA: An automatic treatment planning platform for intensity-modulated proton therapy and its application to intra- and extracerebral neoplasms](https://doi.org/10.48550/arXiv.2305.10211) by R. Bellotti et al. Below are instructions how to setup the development environment and reproduce the results.

> ⚠️ WARNING:<br />
> The code in this repository is only for transparency purposes. It cannot be run as-is because it relies on a binary of the in-house treatment planning system (TPS) FIonA at the Center for Proton Therapy at the Paul Scherrer Institut. These binaries were not approved for publication. However, they are available upon request.<br />
> Please contact one of the authors of the paper if you are interested.


## General

### Disclaimer

> ⚠️ WARNING:<br />
> This code is used for research purposes only. No patient has been treated using any of the software in this repository. The code is still at an early experimental phase.

### Software architecture

The `JulianA` Julia code emerged from a Julia port of our in-house Python library `pyftpp`. `pyftpp` is a wrapper around a command line version of FIonA, whose development name was `ftpp`.

The code for the `JulianA` paper still calls some dependencies from `pyftpp`. Therefore the relevant parts of the `pyftpp` code is included in this repository.

Future versions of `JulianA` will remove the Python dependencies.

> ⚠️ WARNING:<br />
> JulianA relies on a command line version of `FIonA` for dose calculations called `FIonA standalone`. That version is for research purposes only and not used in clinical applications.

### Hardware requirements

The spot weight optimisation uses  GPUs, so please ensure you have an NVIDIA GPU available in order to reproduce the results.

> ⚠️ WARNING:<br />
> JulianA is, in principle, vendor-agnostic when it comes to GPUs. However, we have tested the code exclusively on NVIDIA GPUs. The experiments were run on a DGX A100 machine.<br />
> Support for AMD GPUs can be added by changing all the device allocations to the corresponding library.

The `analyse.ipynb` notebook might need up to 64GB of RAM.

### Operating system

The `Julia` language runs on all major operating systems. However, the `JulianA` code assumes a UNIX system (i. e. Linux or Mac OS) for the path specifications with slashes `/` rather than backslashes `\`. Feel free to run the code on Windows, but be warned that it is not recommended for scientific computing.

All results were obtained on the Merlin/Gwendolen cluster at the Paul Scherrer Institut.

## Installation

### Download the code

You have already obtained the source code of Juliana if you are reading this text.

### Install Julia

First install Julia 1.8.1:

```bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.1-linux-x86_64.tar.gz
tar -xzvf julia-1.8.1-linux-x86_64.tar.gz
```

Setup the environment and install all dependencies by calling `./julia-1.8.1/bin/julia` and enterthe package mode on the REPL by pressing `]`:

```julia
activate .
instantiate
```

### Install the Python dependencies

```bash
conda env create -f pyftpp/environment.yml --prefix <path>
conda activate <path>
pip install pyftpp
python -m ipykernel install --user --name "pyftpp_reproduction"
```
where `<path>` is the path at which the conda environment will be created.


### Setting up `PyCall.jl`

Setup `PyCall.jl` to use the correct Python binary. This can be done on the Julia REPL by calling `./julia-1.8.1/bin/julia` and entering the following commands.

```julia
using Pkg
ENV["PYTHON"] = "<path>/bin/python"
Pkg.build("PyCall")
```


### Setting up the Jupyter kernel

On the Julia REPL:

```julia
using IJulia
IJulia.installkernel("JulianA reproduction", "--project=@.")
```

### Setting up the structure cache

`JulianA` calculates the distance from each voxel to a structure for each structure in order to build the binary mask and the objective function. This calculation is very time consuming. For this reason, `JulianA` uses a cache system and writes the distances to disk after the patient is loaded for the first time. The second loading will be much faster.

The cache directory can be set using the environment variable `JULIANA_CACHE_DIR`. The default value is `~/juliana_cache` if the environment variable is not set.

Some of the SLURM jobs will time out if the patients are not in the cache yet!


## Reproducing the results

### Generating the plans and dose distributions

The code to run an optimisation is contained in `optimise_julia.ipynb`. We decided to run the notebook for each patient on a high-performance computing cluster using the SLURM job manager. The corresponding job script can be run as follows:

```bash
sbatch --array=0-18 optimise_julia_notebook.slurm 
```

This will start a job array of 19 jobs, each running the notebook for a different patient. The SLURM script assumes that there is an environment variable `DATA_DIR`, which points to the directory that contains the dataset.


### Generating the plots for the paper

- Running the notebook `analyse.ipynb` will generate the figure `dose_metrics_boxplot.png` (Fig. 6 in the paper) and `metrics.csv` (Tab. 4 in the paper).
- The notebook `plot_dvhs.ipynb` generates the plots for Fig. 2, 4, 5.
- The notebook `plot_runtime.ipynb` generates the plot for Fig. 3.



## Erratum

We reran the entire evaluation pipeline and compared the resulting DICOM RT dose files. The results agree voxel by voxel, except for patient `test_00`. The irreproducible dose distribution has been used during the analysis that is published in the paper. However, the maximum absolute voxel-wise deviation was 0.75Gy RBE (prescribed dose: 73.8Gy RBE). The DVH curves are overlapping to such an extent that the difference is not visible. Therefore, all results published in the paper are still valid and reproducible.
