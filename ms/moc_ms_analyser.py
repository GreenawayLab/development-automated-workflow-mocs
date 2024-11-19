import numpy as np
import pandas as pd
import rdkit 
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from pyopenms import EmpiricalFormula, FineIsotopePatternGenerator
from pyopenms import *
import csv
from math import isclose
from tqdm import tqdm

def _get_precursor_formula(smiles: str) -> EmpiricalFormula:
    precursor = CalcMolFormula(rdkit.Chem.MolFromSmiles(smiles))
    for i in precursor:
        if i == '-':
            precursor = precursor.replace('-', '') # to avoid error of -ve charge in EmpricialFormula
    precursor_formula = EmpiricalFormula(precursor)
    precursor_formula.setCharge(0) # to avoid error of -ve charge in EmpricialFormula
    return EmpiricalFormula(precursor_formula)

def _get_ligand_formula(triamine_smiles, aldehyde_smiles) -> EmpiricalFormula:
    triamine = CalcMolFormula(rdkit.Chem.MolFromSmiles(triamine_smiles))
    triamine_formula = EmpiricalFormula(triamine)
    aldehyde = CalcMolFormula(rdkit.Chem.MolFromSmiles(aldehyde_smiles))
    aldehyde_formula = EmpiricalFormula(aldehyde)
    H2O = EmpiricalFormula('H2O')
    ligand_formula = triamine_formula + aldehyde_formula + aldehyde_formula + aldehyde_formula - H2O - H2O - H2O
    ligand_formula.setCharge(0) # to avoid error of -ve charge in EmpricialFormula
    return ligand_formula

def _get_intermediate_formula(triamine_smiles, aldehyde_smiles) -> EmpiricalFormula:
    triamine = CalcMolFormula(rdkit.Chem.MolFromSmiles(triamine_smiles))
    triamine_formula = EmpiricalFormula(triamine)
    aldehyde = CalcMolFormula(rdkit.Chem.MolFromSmiles(aldehyde_smiles))
    aldehyde_formula = EmpiricalFormula(aldehyde)
    H2O = EmpiricalFormula('H2O')
    intermediate_formula = triamine_formula + aldehyde_formula + aldehyde_formula - H2O - H2O
    intermediate_formula.setCharge(0) # to avoid error of -ve charge in EmpricialFormula
    return intermediate_formula

def find_solutions(max_L):
    solutions = []
    # set number of metals and ligands with inputs
    for z in range(0, max_L+1):  # set number of tri-ligands
        for y in range(0, max_L+1):  # set number of intermediate ligands
            if y+z>0 and y+z<=max_L: # must have at least one ligand, but only up to the maximum number of ligands
                for m in range((y+z)-1,(y+2*z+1)+1): #set number of metals
                    for a in range(0,2*m+1): #set number of counterions, up to double the number of metals
                        if m < 12:
                            c = 2*m-a #set charge as 2 x metals minus number of counterions
                            line = m,y,z,a,c #create tuple of values
                            solutions.append(line)
    solutions.sort(key=lambda x: x[2])  # Sort solutions by zincs
    return solutions

def calc_formulas(solutions, triamine_smiles, aldehyde_smiles, metal_smiles, counter_ion_smiles):
    aldehyde_formula = _get_precursor_formula(aldehyde_smiles)
    triamine_formula = _get_precursor_formula(triamine_smiles)
    metal_counterion_formula = _get_precursor_formula(metal_smiles) + _get_precursor_formula(counter_ion_smiles)
    metal_formula = _get_precursor_formula(metal_smiles)
    intermediate_formula = _get_intermediate_formula(triamine_smiles, aldehyde_smiles)
    ligand_formula = _get_ligand_formula(triamine_smiles, aldehyde_smiles) 
    counter_ion_formula = _get_precursor_formula(counter_ion_smiles)
    H_adduct = EmpiricalFormula('H')

    LH = ligand_formula + H_adduct
    LH2 = LH + H_adduct 
    LH3 = LH2 + H_adduct
    IH = intermediate_formula + H_adduct
    IH2 = IH + H_adduct
    IH3 = IH2 + H_adduct
    aldH = aldehyde_formula + H_adduct
    triamineH = triamine_formula + H_adduct


    formulas = []
    names = []
    for solution in solutions:
        zincs = solution[0]
        y = solution[1] # intermediate
        z = solution[2] # ligand
        counter_ions = solution[3]
        charge = solution[4]
        if charge == 0:
            charge = 1
        compound_formula=metal_formula #start generating Empirical formula with openpyms
        for i in range(0,zincs-1): #for each component start adding on fragments
                compound_formula = compound_formula + metal_formula
        for i in range(0,y):
                compound_formula = compound_formula + intermediate_formula
        for i in range(0,z):
                compound_formula = compound_formula + ligand_formula
        for i in range(0,counter_ions):
                compound_formula = compound_formula + counter_ion_formula
        gen_name = 'M'+str(zincs)+'I'+str(y)+'L'+str(z)+'_Cnt'+str(counter_ions) #generate name for dictionary key
        names.append(gen_name)
        comp_properties = {'formula': compound_formula,'charge': charge}
        formulas.append(comp_properties) 
    data1 = dict(zip(names,formulas))
    data2 = {
        'L+H' : {'formula': LH, 'charge': 1},
        'L+2H' : {'formula': LH2, 'charge': 2},
        'L+3H' : {'formula': LH3, 'charge': 3},
        'I+H' : {'formula': IH, 'charge': 1},
        'I+2H' : {'formula': IH2, 'charge': 2},
        'I+3H' : {'formula': IH3, 'charge': 3},
        'Ald+H' : {'formula': aldH, 'charge': 1},
        'Triamine+H' : {'formula': triamineH, 'charge': 1},
        'M' : {'formula': metal_formula, 'charge': 2},
        'M+Cnt1' : {'formula': metal_counterion_formula, 'charge': 1}
    }
    # concatenate data1 and data2
    data = {**data1, **data2}
    
    return data

def get_top_10_isotopes(formula_dict, error=1e-1): # error previously 1e-3
    top_10_isotopes_dict = {}

    with tqdm(total=len(formula_dict), desc="Calculating isotope peaks") as pbar:
        for key, value in formula_dict.items():
            formula = value['formula']
            charge = value['charge']

            # Generate isotope distribution with the threshold
            isotopes = formula.getIsotopeDistribution(FineIsotopePatternGenerator(error, False))

            # Normalize m/z values by charge
            normalized_isotopes = [(iso.getMZ() / charge, iso.getIntensity()) for iso in isotopes.getContainer()]

            # Extract and sort isotopes by intensity
            normalized_isotopes.sort(key=lambda x: x[1], reverse=True)

            # Get top 10 highest abundance isotopes
            top_10_isotopes = normalized_isotopes[:10]

            # Store the top 10 isotopes in the dictionary
            top_10_isotopes_dict[key] = {
                'predicted_charge': charge,
                'isotopes': top_10_isotopes
            }
            pbar.update(1)

        return top_10_isotopes_dict

def read_csv(file_path):
    mz_intensity_data = []
    with open(file_path, 'r') as csv_file:
        reader = csv.reader(csv_file)
        next(reader)  # Skip header
        for row in reader:
            mz, intensity = map(float, row)
            mz_intensity_data.append((mz, intensity))
    return mz_intensity_data

def get_var_name(var):
    """
    Retrieve the variable name as a string. 
    """
    return [k for k, v in globals().items() if v is var][0]


def find_matching_isotopes_df(mz_intensity_data, top_10_isotopes_dict, mass_error_ppm=10):
    matching_results = []
    
    # Determine the maximum intensity in the data
    max_intensity = max(intensity for mz, intensity in mz_intensity_data)
    intensity_threshold = 0.005* max_intensity  # alter this value to change the intensity threshold
    
    with tqdm(total=len(mz_intensity_data), desc="Processing peaks") as pbar:
        for mz, intensity in mz_intensity_data:
            if intensity >= intensity_threshold:  # Filtering peaks with intensity above 5% of max intensity
                for formula, isotopes_info in top_10_isotopes_dict.items():
                    predicted_charge = isotopes_info['predicted_charge']
                    isotopes = isotopes_info['isotopes']
                    for iso_mz, iso_intensity in isotopes:
                        theoretical_mass = iso_mz
                        mass_error = theoretical_mass * mass_error_ppm * 1e-6 # mass error from mz to ppm
                        if isclose(mz, theoretical_mass, abs_tol=mass_error):
                            matching_results.append({
                                'Formula': formula,
                                'found_mz': mz,
                                'theoretical_mz': iso_mz,
                                'found_intensity': intensity,
                                'theoretical_intensity': iso_intensity,
                                'predicted_charge': predicted_charge
                            })
            pbar.update(1)
    
    # Convert list of dictionaries to DataFrame
    isotopes_df = pd.DataFrame(matching_results)
    return isotopes_df


def calculate_mz_differences(isotopes_df):
    # Sort by 'Formula' and 'found_mz' to ensure correct difference calculation
    df = isotopes_df.sort_values(by=['Formula', 'found_mz'],ascending=False)

    # Create new columns to store the mz differences and charges
    df['mz_difference'] = None
    df['charge'] = None
    df['splitting_error'] = None

    tolerance = 0.03

    # Iterate over the DataFrame grouped by 'Formula'
    for formula, group in df.groupby('Formula'):
        previous_mz = None
        for index, row in group.iterrows():
            if previous_mz is not None:
                mz_diff = row['found_mz'] - previous_mz
                df.at[index, 'mz_difference'] = mz_diff
                mz_diff_charge_working = row['predicted_charge']*mz_diff
                if round(abs(mz_diff_charge_working)) == 1: # product of charge and splitting must be equal to 1
                    if abs(1-abs(mz_diff_charge_working)) < tolerance: 
                        df.at[index, 'charge'] = row['predicted_charge']
                        df.at[index, 'splitting_error'] = abs(1-abs(mz_diff_charge_working))*100 # percentage error away from expected mz diff 
            previous_mz = row['found_mz']

    return df


# The User should provide the smiles of the triamine precursor, aldehyde precursor, counter ion and metal.

triamineA_smiles = 'NC(C=C1)=CC=C1C2=CC(C3=CC=C(N)C=C3)=CC(C4=CC=C(N)C=C4)=C2'
triamineB_smiles = 'NC1=CC(N)=CC(N)=C1'
aldehyde1_smiles = 'O=CC1=CC=CC=N1'
aldehyde2_smiles = 'O=CC1=CC=CC(C)=N1'
aldehyde3_smiles = 'O=CC1=CC=C(C)C=N1'
aldehyde4_smiles = 'CC1=CC=NC(C=O)=C1'
counter_ion1_smiles = 'O=S([N-]S(=O)(C(F)(F)F)=O)(C(F)(F)F)=O'
counter_ion2_smiles = 'O=S(C(F)(F)F)([O-])=O'
counter_ion3_smiles = 'F[B-](F)(F)F'
metal_smiles = '[Zn]' # The metal is uncharged here, as charges will be explicitly added in the assembly of the MOCs


AB_03_002_001 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde1_smiles, 'counter_ion_smiles': counter_ion1_smiles, 'metal_smiles': metal_smiles}
AB_03_002_002 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde2_smiles, 'counter_ion_smiles': counter_ion1_smiles, 'metal_smiles': metal_smiles}
AB_03_002_003 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde3_smiles, 'counter_ion_smiles': counter_ion1_smiles, 'metal_smiles': metal_smiles}
AB_03_002_004 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde4_smiles, 'counter_ion_smiles': counter_ion1_smiles, 'metal_smiles': metal_smiles}
AB_03_002_005 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde1_smiles, 'counter_ion_smiles': counter_ion1_smiles, 'metal_smiles': metal_smiles}
AB_03_002_006 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde2_smiles, 'counter_ion_smiles': counter_ion1_smiles, 'metal_smiles': metal_smiles}
AB_03_002_007 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde3_smiles, 'counter_ion_smiles': counter_ion1_smiles, 'metal_smiles': metal_smiles}
AB_03_002_008 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde4_smiles, 'counter_ion_smiles': counter_ion1_smiles, 'metal_smiles': metal_smiles}
AB_03_002_009 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde1_smiles, 'counter_ion_smiles': counter_ion2_smiles, 'metal_smiles': metal_smiles}
AB_03_002_010 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde2_smiles, 'counter_ion_smiles': counter_ion2_smiles, 'metal_smiles': metal_smiles}
AB_03_002_011 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde3_smiles, 'counter_ion_smiles': counter_ion2_smiles, 'metal_smiles': metal_smiles}
AB_03_002_012 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde4_smiles, 'counter_ion_smiles': counter_ion2_smiles, 'metal_smiles': metal_smiles}
AB_03_002_013 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde1_smiles, 'counter_ion_smiles': counter_ion2_smiles, 'metal_smiles': metal_smiles}
AB_03_002_014 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde2_smiles, 'counter_ion_smiles': counter_ion2_smiles, 'metal_smiles': metal_smiles}
AB_03_002_015 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde3_smiles, 'counter_ion_smiles': counter_ion2_smiles, 'metal_smiles': metal_smiles}
AB_03_002_016 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde4_smiles, 'counter_ion_smiles': counter_ion2_smiles, 'metal_smiles': metal_smiles}
AB_03_002_017 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde1_smiles, 'counter_ion_smiles': counter_ion3_smiles, 'metal_smiles': metal_smiles}
AB_03_002_018 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde2_smiles, 'counter_ion_smiles': counter_ion3_smiles, 'metal_smiles': metal_smiles}
AB_03_002_019 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde3_smiles, 'counter_ion_smiles': counter_ion3_smiles, 'metal_smiles': metal_smiles}
AB_03_002_020 = {'triamine_smiles': triamineA_smiles, 'aldehyde_smiles': aldehyde4_smiles, 'counter_ion_smiles': counter_ion3_smiles, 'metal_smiles': metal_smiles}
AB_03_002_021 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde1_smiles, 'counter_ion_smiles': counter_ion3_smiles, 'metal_smiles': metal_smiles}
AB_03_002_022 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde2_smiles, 'counter_ion_smiles': counter_ion3_smiles, 'metal_smiles': metal_smiles}
AB_03_002_023 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde3_smiles, 'counter_ion_smiles': counter_ion3_smiles, 'metal_smiles': metal_smiles}
AB_03_002_024 = {'triamine_smiles': triamineB_smiles, 'aldehyde_smiles': aldehyde4_smiles, 'counter_ion_smiles': counter_ion3_smiles, 'metal_smiles': metal_smiles}


# create list of variables to iterate through and run the analysis
data = [AB_03_002_001, AB_03_002_002, AB_03_002_003, AB_03_002_004, AB_03_002_005, AB_03_002_006, AB_03_002_007, AB_03_002_008, AB_03_002_009, AB_03_002_010, AB_03_002_011, AB_03_002_012, AB_03_002_013, AB_03_002_014, AB_03_002_015, AB_03_002_016, AB_03_002_017, AB_03_002_018, AB_03_002_019, AB_03_002_020, AB_03_002_021, AB_03_002_022, AB_03_002_023, AB_03_002_024]
#data = [AB_03_002_001] # for testing purposes, only run one MOC ms analysis
max_ligands = 12 # set the maximum number of ligands to consider in the MOC assembly

for i in data:
    name = get_var_name(i)
    precursor_combinations = find_solutions(max_ligands)
    formula_dict = calc_formulas(precursor_combinations, i['triamine_smiles'], i['aldehyde_smiles'], i['metal_smiles'], i['counter_ion_smiles'])
    # provide the path to the csv file containing the mz and intensity data
    csv_file_path = f'/home/abasford/data/MOC_synth/ms/{name}.csv'
    mz_intensity_data = read_csv(csv_file_path)
    isotopes_dict = get_top_10_isotopes(formula_dict)
    matching_peaks_df = find_matching_isotopes_df(mz_intensity_data, isotopes_dict)
    result_df = calculate_mz_differences(matching_peaks_df)
    result_df.to_csv(f'/home/abasford/projects/MOC_paper/mass_spec/{name}_results_unfiltered.csv', index=False) # uncomment if you want all reuslts without charge filtering
    result_filtered_df = result_df[result_df['predicted_charge'] == result_df['charge']]
    result_filtered_df.to_csv(f'/home/abasford/projects/MOC_paper/mass_spec/{name}_results.csv', index=False)






