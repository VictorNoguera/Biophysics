import sys
import warnings
from Bio.PDB import PDBParser
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)


# GENERAL PART (Imports and PDB parsing)

# Check if minimum arguments are provided
if len(sys.argv) < 3:
    print("Usage: python practical.py <exercise_number> <file.pdb> [parameters...]")
    sys.exit()

exercise = sys.argv[1]
pdb_file = sys.argv[2]

# Load structure only once to avoid repetition across exercises
structure = PDBParser().get_structure('PDB', pdb_file)


# #1 CODE (CA atom pairs closer than distance)
# Run: python practical.py 1 1UBQ.pdb <distance>

if exercise == '1':
    distance = float(sys.argv[3])
    
    # Store all CA atoms
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
    
    # Nested loop to compare all vs all
    pairs = []
    for i in range(len(ca_atoms)):
        for j in range(i + 1, len(ca_atoms)):
            dist = ca_atoms[i] - ca_atoms[j]
            if dist < distance:
                pairs.append((ca_atoms[i].get_parent(), ca_atoms[j].get_parent(), dist))
    
    # Sort by residue number
    pairs.sort(key=lambda x: x[0].id[1])
    for r1, r2, d in pairs:
        print(f"{r1.get_resname()} {r1.id[1]} - {r2.get_resname()} {r2.id[1]} : {d:.2f} A")


# #2 CODE (All atoms for a given residue number)
# Run: python practical.py 2 1UBQ.pdb <chain> <residue_number>

elif exercise == '2':
    chain_id = sys.argv[3]
    res_num = int(sys.argv[4])
    
    # Direct access to the residue
    residue = structure[0][chain_id][res_num]
    atoms = list(residue.get_atoms())
    atoms.sort(key=lambda a: a.get_serial_number())
    
    for atom in atoms:
        c = atom.get_coord()
        print(f"Atom: {atom.get_name():<4} | Coordinates: X={c[0]:.3f}, Y={c[1]:.3f}, Z={c[2]:.3f}")


# #3 CODE (Possible hydrogen bonds between polar atoms)
# Run: python practical.py 3 1UBQ.pdb <distance_cutoff>

elif exercise == '3':
    cutoff = float(sys.argv[3])
    
    # Filter only polar atoms
    polar_atoms = []
    for atom in structure.get_atoms():
        if atom.element in ['O', 'N', 'S']:
            polar_atoms.append(atom)
            
    h_bonds = []
    for i in range(len(polar_atoms)):
        for j in range(i + 1, len(polar_atoms)):
            r1 = polar_atoms[i].get_parent()
            r2 = polar_atoms[j].get_parent()
            
            # Ensure they belong to different residues
            if r1 != r2:
                dist = polar_atoms[i] - polar_atoms[j]
                if dist < cutoff:
                    h_bonds.append((r1, polar_atoms[i], r2, polar_atoms[j], dist))
                    
    h_bonds.sort(key=lambda x: x[0].id[1])
    for r1, a1, r2, a2, d in h_bonds:
        print(f"{r1.get_resname()}{r1.id[1]} ({a1.get_name()}) -- {r2.get_resname()}{r2.id[1]} ({a2.get_name()}) : {d:.2f} A")


# #4 CODE (CA atoms of a given residue type)
# Run: python practical.py 4 1UBQ.pdb <residue_type>

elif exercise == '4':
    res_type = sys.argv[3].upper()
    aa_dict = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 
               'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 
               'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER', 'T':'THR', 'V':'VAL', 
               'W':'TRP', 'Y':'TYR'}
               
    # Translate to 3-letter code if necessary
    if len(res_type) == 1:
        res_type = aa_dict.get(res_type, res_type)
        
    valid_residues = []
    for res in structure.get_residues():
        if res.get_resname() == res_type and 'CA' in res:
            valid_residues.append(res)
            
    valid_residues.sort(key=lambda r: r.id[1])
    for res in valid_residues:
        c = res['CA'].get_coord()
        print(f"Chain {res.get_parent().id}, Residue {res.id[1]} | Coordinates: X={c[0]:.2f}, Y={c[1]:.2f}, Z={c[2]:.2f}")