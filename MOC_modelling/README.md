# MOC-automation
This repository contains the scripts used for the automated atomistic modelling and analysis for the: Development of an Automated Workflow for Screening the Assembly and Host-Guest Behaviour of Metal-Organic Cages towards Accelerated Discovery.

## Software used
- [Scigress](https://www.fqs.pl/en/chemistry/products/scigress) v2.6
- [ORCA](https://www.faccts.de/orca/) v5.0.4
- [xtb](https://github.com/grimme-lab/xtb) v6.6.1
- [stk](https://github.com/lukasturcani/stk) v2024.3.28.0 (*stk.small.NCore* doesn't work for earlier versions of *stk*)
- [stko](https://github.com/lukasturcani/stko) v2023.11.13.0 
- [gulp](https://gulp.curtin.edu.au/) v5.1
- [OPTIM](https://www-wales.ch.cam.ac.uk/OPTIM/) 
- SLURM scheduler (for job arrays)
- Python 3.11.9 (the used *stk* version requires python 3.11)
- [pywindow](https://github.com/marcinmiklitz/pywindow)

  For installing these tools, refer to their respective documentation or package managers suitable for your system.

## Repository Structure
- Scripts/: contains all bash and python scripts
- Tetrahedrons_stko/: contains all the MD and GFN2-xTB ("vtight") optimized cage structures (*Tetrahedron_l[i]_ald[j]_MD.mol* and *Tetrahedron_l[i]_ald[j]_XTB.mol*) and their vibrational spectra (*Tetrahedron_l[i]_ald[j]_XTB_vibspectrum*) calculated at the GFN2-xTB ("vtight") level.
- Tetrahedrons_ORCA_DFT/: contains all the input files (*Tetrahedron_l[i]_ald[j]_ORCA_DFT.inp* and *Tetrahedron_l[i]_ald[j].xyz*) and output files (*Tetrahedron_l[i]_ald[j]_ORCA_DFT.out*) for the r2SCAN-3c single-point energy calculations using _ORCA_.

Here, "*i*" refers to the type of linker (l0 = A or l1 = B) and "*j*" refers to the type of aldehyde (ald1 - ald4).

## Code
### Tetrahedrons
1) Construct a series of tetrahedral MOCs with _stk_, optimise them to the GFN2-xTB level with stko and create ORCA input files:
   
   >>sbatch ConstructionTet.sh

   With this script, first the building blocks used in the construction process are defined. Four different metal complex building blocks are created by building four different imine-groups (ald1, ald2, ald3 and ald4) around a Zn^{2+} metal ion using the *stk.metal_complex.OctahedralLambda* function. Two different tritopic organic linker building blocks are created from their SMILES strings using the *stk.BuildingBlock* function. Both the Zn metal complexes and organic linkers are placed on a tetrahedral topology graph using the *stk.cage.M4L4Tetrahedron* function. Only one diastereomer is made, the one with all Lambda metal complexes (similar to CCDC [929026](https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid=929026&DatabaseToSearch=Published)). After construction of the eight tetrahedrons, a series of optimization methods is performed on each structure. First, two geometry optimizations are performed with UFF4MOF in GULP. Secondly, a conformer search is performed on the UFF4MOF optimized structures using high-temperature molecular dynamics (MD). The lowest energy cage conformer is saved as *Tetrahedron_l[i]_ald[j]_MD.mol*, with *i* noting the relevant organic linker and *j* noting the relevant aldehyde.
   
   Lastly, GFN2-xTB geometry optimizations are performed at the GFN2-xTB level (opt-level *vtight*). The *xtb* program is called through *stko* by running the following line in the background:

   ```
   /home/pcpt3/.conda/envs/env3/bin/xtb input_structure_1.xyz --gfn 2 --ohess vtight --parallel 32 --etemp 300 --chrg 8 --uhf 0 -I det_control.in
   ```

   The results after the geometry optimization are saved as *Tetrahedron_l[i]_ald[j]_XTB.mol*. The vibrational spectrum that are created by the vibrational analysis is saved as *vibspectrum* in their respective output folder. In this repository these spectra are renamed to match the corresponding structure (*Tetrahedron_l[i]_ald[j]_XTB_vibspectrum*).

    Besides the construction and optimization of the tetrahedrons with *stk* and *stko*, running this script also leads to the creation of the input files for the ORCA calculations run in step 2: *Tetrahedron_l[i]_ald[j]_XTB.xyz* and *Tetrahedron_l[i]_ald[j]_ORCA_DFT.inp*, with the latter looking like this:

   ```
   ! r2SCAN-3c ENERGY TightSCF def2/J 
   %PAL NPROCS 32 END 
   %maxcore 5316 

   *xyzfile 8 1 /rds/user/pcpt3/hpc-work/CageConstruction/GreenawayCollab/Tetrahedrons_ORCA_DFT/Tetrahedron_l[i]_ald[j]_ORCA_DFT/Tetrahedron_l[i]_ald[j]_XTB.xyz 
   ```
   
2) Run DFT single-point energy calculations on GFN2-xTB optimized structures (r<sup>2</sup>SCAN-3c//GFN2-xTB level) using ORCA:
   
   >>sbatch ORCAARRAY_Tet.sh

    Eight output files (*Tetrahedron_l[i]_ald[j]_ORCA_DFT.out*) are created after running this script, each containing the DFT single point energy in the gasphase of the corresponding structure.
   
### Icosahedrons
Five different icosahedrons were constructed from linker l0 and l1 using MM3 modelling in Scigress. The crystal structure of a previously published Fe-based icosahedron was used as a starting point (see CCDC [929025](https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid=929025&DatabaseToSearch=Published)). Four icosahedrons were constructed with linker A and four different aldehydes, and one icosahedron was constructed with linker B and aldehyde ald1. Each structure was optimized at the GFN-FF level of theory. Optimization at the GFN2-xTB level did not lead to convergence for the linker A structures, but did for the linker B structure. The construction and optmization process of the icosahedrons was not automated due to the complexity of these geometries (e.g. the presence of non-C3-symmetric meridional octahedral complexes at the vertices instead of C3-symmetric facial octahedral complexes).

Jobs were run through XTBOPTIM:

>>sbatch sbatch.optim

In this script, the XTBOPTIM program is called:

```
XTBOPTIM > ${dir_path}/OPTIM.out
```

This program opens a file called 'odata' which contains the keywords:

```
BFGSMIN 1.0D-6 
CONVERGE 0.1 1.0D-6 
UPDATES 500 100 5 5 
MAXERISE 1.0D-5 0.02D0 
MAXBFGS 0.2 0.1 
BFGSSTEPS 10000 

RADIUS 2000.0 
DUMPDATA 

XTBAPI <input>.xyz <method> 24 
XTBARG accuracy 0.1 
```

Here, /<method/> is either gfnff or gfn2.

## Pywindow 
Scripts/pywindow.py shows the automated analysis for reading in the ORCA_DFT XYZ file, conversion to a pywindow molecule and subsequent diameter and volume calcualtions. This is then written to a JSON file, indexed via variable name. 
