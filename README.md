# j_presto

`j_presto` is a program suite for preparing biomolecular systems, performing molecular dynamics (MD) simulations, and analyzing the resulting structural ensembles.

The suite includes the original Adaptive Lambda Square Dynamics (ALSD) method for enhancing conformational sampling in selected molecular regions. It also supports the zero-dipole summation method (ZDM) and related electrostatic treatments for efficient MD calculations.

## Installation

### 1. Obtain the AMBER parameter database

The `prep` workflow uses parameter files distributed with AmberTools. AmberTools is an external package developed by the AMBER project and is not distributed with j_presto.

1. Download and extract AmberTools from the official AMBER website.
2. Locate the directory containing the `leaprc`, `dat`, `lib`, `frcmod`, and `prep` files. For example, AmberTools24 stores these files under `amber24_src/dat/leap`.
3. Set the directory path in your shell configuration:

```bash
export J_PRESTO_AMBER_DATABASE_PATH="<path-to-the-AMBER-database>"
```

### 2. Obtain j_presto

Clone the repository:

```bash
git clone https://github.com/jikebe-AIST/j_presto.git
cd j_presto
```

Alternatively, download and extract the ZIP archive from the GitHub repository.

### 3. Run the installation script

j_presto contains Fortran and Python 3 programs. Before installation, prepare:

- Intel Fortran (`ifort`) or GNU Fortran (`gfortran`)
- Python 3
- `f2py`
- the Python modules checked by `setup.sh`

Run the installer from the repository root:

```bash
./setup.sh
```

The installer checks the required Python modules, compiles the Fortran programs, builds the `calc_Efield` extension, and creates an installed `J_PRESTO` directory. It asks for:

1. the destination directory, defaulting to `~/local/bin`;
2. the Fortran compiler, defaulting to `ifort`.

If a required Python module is missing, install it and run `setup.sh` again. Some subprograms also use optional packages that may not be required by the core installer. For example, external trajectory conversion with `traj_conv` requires MDAnalysis.

After installation, add the following settings to `.bashrc` or `.bash_profile`:

```bash
export J_PRESTO_PATH="<installation-directory>/J_PRESTO"
export J_PRESTO_AMBER_DATABASE_PATH="<path-to-the-AMBER-database>"
export PATH="$PATH:$J_PRESTO_PATH"
```

To use the Japanese interactive manual by default, also add:

```bash
export J_PRESTO_LANG=ja
```

Reload the shell configuration, for example:

```bash
source ~/.bashrc
```

## Usage

Run a subprogram through the top-level launcher:

```bash
j_presto <subprogram> [arguments]
```

Display the current public subprogram list with:

```bash
j_presto -h
```

### Simulation preparation and execution

| Subprogram | Description |
|---|---|
| `initpep` | Generate extended peptide structures for folding simulations. |
| `pdb_alignfit` | Superimpose a query structure onto a target using sequence alignment. |
| `nt_gen` | Generate an atom-name correspondence table from PREP and PDB files. |
| `nt_conv` | Convert residue and atom names in PREP/PDB files using a name table. |
| `gen_db` | Generate a topology database (`*.tpldb`) for creating j_presto topology files (`*.tpl`). |
| `prep` | Prepare PDB, TPL, and SHK input files for MD simulations. |
| `pepchk` | Check peptide-bond geometry and classify cis, trans, and twisted bonds. |
| `md_run` | Run MD simulations or energy-minimization calculations. |
| `batch_set` | Set up files and directories for multiple MD simulation runs. |
| `GEprep` | Prepare generalized-ensemble parameter files (`*.nf`) for McMD, ALSD, and related methods. |

### Analysis, conversion, and visualization

| Subprogram | Description |
|---|---|
| `Ens_Ana` | Analyze conformational ensembles obtained from simulations. |
| `PCAaxis` | Calculate principal-component axes by covariance-matrix diagonalization. |
| `PCAproj` | Project structures onto a PCA subspace. |
| `distrib` | Generate weighted one- or two-dimensional statistical distributions. |
| `pick_conf` | Extract and rank structures using user-defined conditions and event timing. |
| `ttp_se` | Calculate weighted means and standard errors across independent runs. |
| `pdb_movie` | Combine sequential PDB files into a multi-model PDB trajectory. |
| `traj_conv` | Convert structures and trajectories between j_presto and external formats. |
| `DISTIL` | Detect and interpret structural differences between two molecular ensembles. |
| `struct_formula` | Draw two-dimensional molecular structures with atom names mapped from a PDB file. |

### Enzyme engineering

| Subprogram | Description |
|---|---|
| `MSPER` | Predict mutation sites for improving enzyme regioselectivity. |
| `conmut` | Propose stability-enhancing mutations using sequence consensus. |

### Utilities and documentation

| Subprogram | Description |
|---|---|
| `get_tmpl` | Copy templates for j_presto input files and scripts. |
| `genlist` | Generate file-path lists and data ranges for MD and analysis inputs. |
| `manual` | Open the interactive j_presto manual. |

Development-stage or internal utilities may exist in the `sp` directory without appearing in this public list.

## Parallel execution

The MD calculation itself supports OpenMP parallelization. Direct GPU execution is not currently supported. Multiple independent OpenMP-parallelized runs can be grouped into one MPI job using the following launchers:

| Program | Description |
|---|---|
| `j_presto_mpi` | Execute multiple MD runs as an embarrassingly parallel MPI job. |
| `j_presto_master` | Execute multiple MD runs using a master-worker MPI model. |

## Interactive manual

The interactive manual contains command descriptions, input-file keywords, tutorials, and frequently asked questions.

```bash
j_presto manual
```

Select a language explicitly:

```bash
j_presto manual -l en
j_presto manual -l ja
```

Use locale-based automatic selection:

```bash
j_presto manual -l auto
```

The language may also be changed while the manual is running:

```text
manual[en]> lang ja
manual[ja]> lang en
```

New users should open `tutorial` from the manual home page. The tutorials are organized by task, including installation, system preparation, MD simulation, trajectory analysis, PCA, DISTIL, trajectory conversion, peptide-bond inspection, MSPER, and molecular-structure drawing.

## License

j_presto is distributed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0).

The software may be used and modified for personal, research, and academic purposes. Commercial redistribution, sale, or sublicensing of the software itself is prohibited. Commercial use of results generated with the software, including products developed from simulation or analysis results, is permitted provided that the software itself is not sold, redistributed, or sublicensed.

See `LICENSE.md` for the complete license and citation information.

## Contact

j_presto is provided by Jinzen Ikebe, Artificial Intelligence Research Center, National Institute of Advanced Industrial Science and Technology (AIST), AIST Tokyo Waterfront Bio-IT Research Building, 2-4-7 Aomi, Koto-ku, Tokyo 135-0064, Japan.

Questions: `ikebe.jinzen@aist.go.jp`
