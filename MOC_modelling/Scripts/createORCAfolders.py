from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import argparse
import stk
import shutil

os.chdir(os.path.dirname(__file__))
cd = os.getcwd()
print(cd)

xyzwriter = stk.XyzWriter()

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description="Create ORCA folders")
parser.add_argument('--arraynr', type=int, help='Loop number for the array job')

args = parser.parse_args()

array_nr = args.arraynr

#List of tritopic linkers
Linkers_Br = [
    'BrC(C=C1)=CC=C1C2=CC(C3=CC=C(Br)C=C3)=CC(C4=CC=C(Br)C=C4)=C2',
    'BrC1=CC(Br)=CC(Br)=C1'
]

#List of imine moieties
Aldehydes_Br_dict = {'ald1': 'C1=NC(C=NBr)=CC=C1', 'ald2':'Br/N=C/C1=CC=CC(C)=N1', 'ald3':'Br/N=C/C1=CC=C(C)C=N1', 'ald4':'CC1=CC=NC(/C=N/Br)=C1'}

#List of helicate types
ru_dict = {"M2A3": "AAA", "M2A2B1":"AAB", "M2A1B2":"ABB", "M2B3":"BBB"}

#List of all the Iodo helicate cage names 
cage_names = [0] * (len(Aldehydes_Br_dict)*len(ru_dict)*2)
t = 0

for a in list(Aldehydes_Br_dict.keys()):
  for type in list(ru_dict.keys()):
    for stereo in ["LL", "LD"]:
      for l in [0,1]:
        cage_names[t] = f"Helicate_{type}_{stereo}_l{l}_{a}_collaps4_collaps_MM3"
        t+=1

cage_name = cage_names[array_nr]
cage_file = f"{cd}/{cage_name}.xyz"

folder = f"{cd}/{cage_name}_ORCA_GFN2"
os.makedirs(folder, exist_ok=True)
shutil.copy(cage_file, f"{folder}/{cage_name}.xyz")
inp_file_path = f"{folder}/{cage_name}_ORCA_GFN2.inp"

#Create ORCA.inp file for constrained GFN2-xTB optimisation of helicates with ORCA
with open(inp_file_path, 'w') as file:
    file.write("! XTB2 OPT Freq \n") 
    file.write("%PAL NPROCS 10 END \n")
    file.write("%maxcore 3320 \n")
    file.write("\n")
    file.write(f"*xyzfile 4 1 {cage_name}.xyz \n")
print(f"{cage_name}_ORCA_GFN2.inp file has been created.")
