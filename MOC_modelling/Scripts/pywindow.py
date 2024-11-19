import pywindow as pw
import json

def get_var_name(var):
    """
    Retrieve the variable name as a string.
    """
    return [k for k, v in globals().items() if v is var][0]

cage1 = pw.MolecularSystem.load_file('/path/to/xyz/file')
cage3 = pw.MolecularSystem.load_file('/path/to/xyz/file')

# convert the systems to Molecule objects
cages = [cage1, cage3]
pw_dict = {}
for cage in cages:
    name = get_var_name(cage)
    #print(name)
    mol = cage.system_to_molecule()
    # append each molecule to the dictionary
    pw_dict[name] = {'name' : name, 'cage' : mol}

#print(pw_dict)
# iterate over pw_cages and add new properties to the dictionary
for cage in pw_dict:
    cage_mol = pw_dict[cage]['cage']
    cage_name = pw_dict[cage]['name']
    cage_mass = cage_mol.molecular_weight()
    cage_diam = cage_mol.calculate_pore_diameter_opt()
    cage_vol = cage_mol.calculate_pore_volume_opt()
    pw_dict[cage].update({'mass' : cage_mass,
                          'diameter' : cage_diam,
                          'volume' : cage_vol
                          })
    
    # drop 'cage' and 'name' key from the dictionary
    pw_dict[cage].pop('cage')
    pw_dict[cage].pop('name')
    

# save the dictionary to a JSON file
with open('/path/to/output/cage_data_pwindow_.json', 'w') as f:
    json.dump(pw_dict, f, indent=4)
    