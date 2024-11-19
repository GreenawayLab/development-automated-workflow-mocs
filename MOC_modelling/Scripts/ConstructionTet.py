import os
import sys
import stk
import stko
import pandas as pd
import re
import shutil

modules_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '/home/pcpt3/Python_Calculations/Scripts'))
sys.path.append(modules_path)

from env_set_Paula_CSD3 import xtb_path, gulp_path

os.chdir(os.path.dirname(__file__))
cd = os.getcwd()
print(cd)

xyzwriter = stk.XyzWriter()

#List of tritopic linkers
Linkers_Br = [
    'BrC(C=C1)=CC=C1C2=CC(C3=CC=C(Br)C=C3)=CC(C4=CC=C(Br)C=C4)=C2',
    'BrC1=CC(Br)=CC(Br)=C1'
]

#List of imine moieties
Aldehydes_Br_dict = {'ald1': 'C1=NC(C=NBr)=CC=C1', 'ald2':'Br/N=C/C1=CC=CC(C)=N1', 'ald3':'Br/N=C/C1=CC=C(C)C=N1',  'ald4':'CC1=CC=NC(/C=N/Br)=C1'}

#Define the metal center buildingblock
metal_atom = stk.BuildingBlock(
    smiles= '[Zn+2]',
    functional_groups=(
        stk.SingleAtom(stk.Zn(0, charge=2))
        for j in range(6)
    ),
    position_matrix=[[0, 0, 0]],
)

#Define octahedral lambda complex building block
def make_complex_bb(a, stereo):
    fp_Br_factories = [
            #CNBr_metal
            stk.SmartsFunctionalGroupFactory(
                smarts='[#6]~[#7X2]~[#35]',
                bonders=(1, ),
                deleters=(),
            ),
            #CNC_metal
	    stk.SmartsFunctionalGroupFactory(
                smarts='[#6]~[#7X2]~[#6]',
                bonders=(1, ),
                deleters=(),
            )]

    fp_Br = stk.BuildingBlock(
             smiles=(Aldehydes_Br_dict[a]),
             functional_groups=fp_Br_factories
          )

    if stereo == "Lambda":
        metal_complex = stk.ConstructedMolecule(
                topology_graph=stk.metal_complex.OctahedralLambda(
                    metals={metal_atom: 0},
                    ligands={
                        fp_Br: (0, 1, 2),
                        },
                    #optimizer = stk.Collapser(),
                    optimizer = stk.MCHammer(),
                )
        )

    if stereo == "Delta":
        metal_complex = stk.ConstructedMolecule(
                topology_graph=stk.metal_complex.OctahedralDelta(
                    metals={metal_atom: 0},
                    ligands={
                        fp_Br: (0, 1, 2),
                        },
                    #optimizer = stk.Collapser(),
                    optimizer = stk.MCHammer(),
                )
        )

  # More than one placer id has to be defined for complexes to ensure alignment
  # Use the NCCN plane (attached to the Br) to define the orientation of the complex.
    fgfactory = stk.SmartsFunctionalGroupFactory(
      smarts='[#7X3]~[#6]~[#6]~[#7X3]~[#35]',
      bonders=(3, ),
      deleters=(4, ),
      placers=(0, 1, 2, 3),
    )

    bb_metal_complex = stk.BuildingBlock.init_from_molecule(
        metal_complex,
        functional_groups=[fgfactory]
      )

    return bb_metal_complex

#Define optimization functions
def MOC_GULP_twice(MOC, name, gulp_file):
  # Unrestrained UFF opt.
  print(f'..doing unrestrained UFF4MOF optimisation of {name}')
  gulp_opt = stko.GulpUFFOptimizer(
          gulp_path=gulp_path(),
          output_dir= f'{gulp_file}_1',
          maxcyc=1000,
          metal_FF= {30: 'Zn4+2'},  # No alternative available for 90 degrees.
          metal_ligand_bond_order='',
          conjugate_gradient=True,
        )
  gulp_opt.assign_FF(MOC)
  gulp_cage = gulp_opt.optimize(mol=MOC)
  #gulp_cage.write(f'{gulp_file}_1.mol')
  shutil.rmtree(f"{gulp_file}_1")

  print(f'..doing second unrestrained UFF4MOF optimisation of {name}')
  gulp_opt2 = stko.GulpUFFOptimizer(
          gulp_path=gulp_path(),
          output_dir= f'{gulp_file}_2',
          maxcyc=1000,
          metal_FF={30: 'Zn4+2'},  # No alternative available for 90 degrees.
          metal_ligand_bond_order='',
          conjugate_gradient=False,
        )
  gulp_opt2.assign_FF(gulp_cage)
  gulp_cage2 = gulp_opt2.optimize(mol=gulp_cage)
  #gulp_cage2.write(f'{gulp_file}_2.mol')
  shutil.rmtree(f"{gulp_file}_2")
  return gulp_cage2

def MOC_GULP_MD(MOC, name, MD_file):
  print(f'..doing MD UFF4MOF optimisation of {name}')
  gulp_MD = stko.GulpUFFMDOptimizer(
                gulp_path= gulp_path(),
                output_dir= f'{MD_file}',
                metal_FF={30: 'Zn4+2'},  # No alternative available for 90 degrees.
                metal_ligand_bond_order='half',
                integrator='leapfrog verlet',
                ensemble='nvt',
                temperature=400,
                equilbration=0.1,
                production=2,
                timestep=0.5,
                N_conformers=10,
                opt_conformers=False,
                save_conformers=False,
            )
  gulp_MD.assign_FF(MOC)
  gulp_MD_cage = gulp_MD.optimize(MOC)
  gulp_MD_cage.write(f'{MD_file}.mol')
  shutil.rmtree(MD_file)
  return gulp_MD_cage

def MOC_XTB(MOC, name, XTB_file,charge):
  print(f'..doing XTB optimisation of {name}')
  xtb_opt = stko.XTB(
          xtb_path=xtb_path(),
          output_dir= f'{XTB_file}',
          gfn_version=2,
          num_cores=32,
          opt_level = 'vtight',
          charge=charge,
          num_unpaired_electrons=0,
          max_runs=1,
          electronic_temperature=300,
          calculate_hessian=True,
          unlimited_memory=True,
        )
  xtb_cage = xtb_opt.optimize(mol=MOC)
  xtb_cage.write(f'{XTB_file}.mol')
  #shutil.rmtree(XTB_file)
  return xtb_cage

def MOC_OptSequence(MOC, cage_name, charge):
    gulp_file = f"{cd}/{cage_name}_gulp"
    gulp_cage = MOC_GULP_twice(MOC, cage_name, gulp_file)
    MD_file = f"{cd}/{cage_name}_MD"
    MD_cage = MOC_GULP_MD(gulp_cage, cage_name, MD_file)
    XTB_file = f"{cd}/{cage_name}_XTB"
    
    try:
        XTB_cage = MOC_XTB(MD_cage, cage_name, XTB_file, charge)
        return XTB_cage
    except Exception as e:
        # Handle the error and continue
        print(f"MOC_XTB failed for {cage_name} with error: {e}")
        print(f"MD result for {cage_name} printed instead.")
        MD_cage.write(f'{MD_file}.mol')
        return MD_cage
    #return MD_cage

#Make linker building blocks, with three bromines (for tetrahedrons) and with two bromines and one amine (for helicates)
bb_planar = [0] * (len(Linkers_Br)*2)
for n in range(0, len(Linkers_Br)):
    bb_planar[n] = stk.BuildingBlock(
             smiles=(Linkers_Br[n]),
             functional_groups=[stk.BromoFactory()],
          )
    bb_planar[n+2] = stk.BuildingBlock(
             #Change one Br into an I group
             smiles=(Linkers_Br[n].replace('Br', 'I', 1)),
             functional_groups=[stk.BromoFactory()],
          )


#Make complex building blocks
metal_facL = [0] * len(Aldehydes_Br_dict)
for i in range(0, len(Aldehydes_Br_dict)):
    metal_facL[i] = make_complex_bb(list(Aldehydes_Br_dict.keys())[i], "Lambda")

#Make and optimise tetrahedrons
for i in [0,1]:
  for j in range(0,4):
    cage_name = f"Tetrahedron_l{i}_{list(Aldehydes_Br_dict.keys())[j]}"
    cage_file = f"{cage_name}_XTB.mol"

    if not os.path.exists(cage_file):
        collaps_cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L4Tetrahedron(
                building_blocks = {bb_planar[i]: (4, 5, 6, 7),
                      metal_facL[j]: (0, 1, 2, 3)},
                optimizer = stk.Collapser(),
        ))

        XTB_cage = MOC_OptSequence(collaps_cage, cage_name, 8)

    mol_stk = stk.BuildingBlock.init_from_file(path = cage_file) 
    DFT_folder = f"/rds/user/pcpt3/hpc-work/CageConstruction/GreenawayCollab/Tetrahedrons_ORCA_DFT/{cage_name}_ORCA_DFT"
    os.makedirs(DFT_folder, exist_ok=True)
    xyzwriter.write(molecule=mol_stk, path=f"{DFT_folder}/{cage_name}_XTB.xyz")
    inp_file_path = f"{DFT_folder}/{cage_name}_ORCA_DFT.inp"

    with open(inp_file_path, 'w', encoding='utf-8') as file:
        file.write("! r2SCAN-3c ENERGY TightSCF def2/J \n") 
        file.write("%PAL NPROCS 32 END \n")
        file.write("%maxcore 5316 \n") #0.75*3544MB on icelake per CPU or 0.75*7088MB on icelake-himem per cpu
        file.write("\n")
        file.write(f"*xyzfile 8 1 {DFT_folder}/{cage_name}_XTB.xyz \n")
    print(f"{cage_name}_ORCA_DFT.inp file has been created.")
