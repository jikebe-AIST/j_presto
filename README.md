# j_presto
"j_presto" is a program suite that compiles subprograms for preparing input files and executing molecular dynamics (MD) simulations, specifically for biomolecules such as proteins.
This suite features the original computational method, Adaptive Lambda Square Dynamics (ALSD) simulation, which enables comprehensive conformational sampling by enhancing conformational changes in specific regions of the molecule.
Additionally, it incorporates the Zero Dipole Expansion Method (ZDM), a technique that allows for the efficient and highly accurate calculation of electrostatic interactions.

## Installation

### Before Installation

#### Download AMBERtools
To create simulation parameter files for j_presto, a database provided by AMBERtools is required.
AMBERtools is an external program provided by the AMBER project at the University of California, San Diego (UCSD) in the United States, and the developers of j_presto do not hold its copyright.
Users are responsible for downloading it at their own risk.
	1 : Access the AMBERtools webpage (https://ambermd.org/AmberTools.php).
	2 : Click the "Download AmberTools" link
	3 : Read the description on the linked page carefully and download AmberTools using your preferred method.
	4 : Extract the downloaded file (AmberTools***.tar.bz2, where *** is the version number) and copy the directory containing the parameter files (leaprc, dat, lib, frcmod, prep files) to a location of your choice.
	    For example, in the case of AmberTools24, these files are stored inside the amber24_src/dat/leap directory.
	5 : If you are using bash, add the line
		export J_PRESTO_AMBER_DATABASE_PATH="[the path where you placed the db directory]"
	    to your .bashrc or .bash_profile file, and then run
		source .bashrc (if you added it to the .bashrc file)
	    to set the path to the AMBER database directory.

#### Download j_presto repository
	1 : Access the GitHub Repository:
		Open your web browser and navigate to the j_presto GitHub repository page.
	2 : Clone or Download the Repository:
		To clone the repository, click on the green "Code" button located on the right side of the page.
		Copy the URL provided in the dropdown menu.
		Open your terminal and enter the following command:
			git clone [copied URL]
		Alternatively, you can choose to download the repository as a ZIP file.
		Click the "Code" button, then select "Download ZIP".
		This will download the repository contents as a compressed file to your local machine.
		Extract the Files.

### Install j_presto
	1 : The j_presto program is written in Fortran and Python3.
	    Before installing the program, please ensure that you have installed the following three components:
		1 : Intel Fortran (ifort) or gfortran
		2 : python3
		3 : f2py, which is required to compile Fortran subroutines into Python modules.
	2 : Navigate into the downloaded j_presto directory and execute the installation script named "setup.sh" as follows:
		./setup.sh
	3 : The setup.sh script first checks whether the necessary Python modules for installation are present.
	    If any modules are missing, their names will be displayed, and the installation of j_presto will be aborted.
	    In this case, please install the missing Python modules with "pip" command or similar, and then run setup.sh again.
	    !! CAUTION !!
	    Even when the Python modules are installed, if their versions are outdated, it can cause Python errors and lead to abnormal termination of the j_presto program.
	    In such cases, it is recommended to check the error messages and update the relevant Python modules to their latest versions.
	4 : When you run setup.sh, you will be prompted to specify two things: the location to store the j_presto executable files, and the Fortran compiler to use.
	    For the first, choose an appropriate directory (the default is ~/local/bin/), and for the latter, select either ifort or gfortran (the default is ifort).
	5 : Finally, if you are using bash, add the path where you saved the j_presto directry including executable files to your .bashrc or .bash_profile file as follows:
		export J_PRESTO_PATH="[the path where you placed the J_PRESTO directory]"
	   and please add J_PRESTO_PATH to your PATH as shown below:
		export PATH="$PATH:$J_PRESTO_PATH"
	   After that, run
		source .bashrc (if you added it to the .bashrc file)
	   to set the directory containing the j_presto executables to your system's path.

## Usage
"j_presto" is a program for performing molecular dynamics (MD) simulations.
With j_presto, you can execute the following sub-programs using the command:
	j_presto [sub_program_name]

Sub-program name list:

    For MD simulation,
        pdb_alignfit    : Superimposition of query structure to target one based on sequence alignment.
        nt_gen          : Generation of atomic name correspondence table from a PREP and a PDB file.
        nt_conv         : Convert residue and atom names in a PREP and a PDB file based on a name table.
        gen_db          : Generation of a database file (*.tpldb) for creating j_presto topology files (*.tpl).
        prep            : Preparation of input files (*.pdb, *.tpl, and *.shk files) for MD simulations.
        md_run          : Execution of MD simulation or energy minimization calculation.
        batch_set       : Setting up all the necessary files and directories required to run multiple MD simulations efficiently.
        GEprep          : Preparation of a parameter file (*.nf) for Generalized Ensemble MD simulations such as McMD and ALSD.

    For Analysis,
        Ens_Ana         : Performing analysis of conformational ensembles obtained from simulation results.
        PCAaxis         : Calculating axes for principal component analysis (PCA) through diagonalization.
        PCAproj         : Projecting each structure onto PCA subspace.
        distrib         : Generation of statistical distribution data from weighted input values for one- or two-dimensional variables.
        pick_conf       : Extraction and ranking structures based on user-defined conditions and event timings.
        ttp_se          : Calculation of weighted averages and standard errors of the data

    For Enzyme development,
	MSPER           : Predict mutation sites to improve enzyme selectivity.
        conmut          : Propose mutations to improve enzyme stability from input amino acid sequences using the consensus method

    For user guidance,
        get_tmpl        : Obtaining templates of input files and scripts for excuting j_presto.
        manual          : Viewing the manual.

The simulations performed by j_presto support parallelization using OpenMP, but currently do not support parallelization via MPI or GPU.
However, you can use the programs j_presto_mpi and j_presto_master, which allow you to bundle multiple OpenMP-parallelized MD simulation runs and execute them as a single MPI job.
	j_presto_mpi    : Execute multiple MD simulation runs as a single MPI job using an embarrassingly parallel approach.
	j_presto_master : Execute multiple MD simulation runs as a single MPI job using a master-slave approach.

The most important command here is "j_presto manual".
Running this command allows you to access the tutorial, along with explanations of the options and input files required to run j_presto.
For those who wish to learn how to use j_presto, first execute "j_presto manual," then proceed with "tutorial" to learn the complete workflow.

## License
This program, j_presto, is available under the Creative Commons Attribution-NonCommercial 4.0 International Public License (CC BY-NC 4.0).
The software is freely accessible for personal, research, and academic use, with permission to modify it for personal research needs.
However, commercial redistribution, sale, or sublicensing of the software itself is prohibited.

The license permits commercial use of results generated using the software, including the development and commercialization of products derived from simulation outputs (e.g., protein design or other biomolecular products), as long as the software itself is not sold, redistributed, or sublicensed.

For full details, please refer to the LICENSE.md file.

## Contact
The source codes of "j_presto" are provided by Jinzen Ikebe at Artificial Intelligence Research Center (AIRC), AIST Tokyo Waterfront BIO-IT Research Building, 2-4-7 Aomi, Kyoto-ku, Tokyo 135-0064, JAPAN.
For question, please contact ikebe.jinzen@aist.go.jp.
