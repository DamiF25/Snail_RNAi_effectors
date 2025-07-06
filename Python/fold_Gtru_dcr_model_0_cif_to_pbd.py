from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO

def cif_to_pdb_with_plddt(cif_file, output_pdb):
    """
    Convert AlphaFold .cif to .pdb and map plDDT to B-factor.

    Parameters:
        cif_file (str): Input .cif file path.
        output_pdb (str): Output .pdb file path.
    """
    # Parse the .cif file
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", cif_file)

    # Access the plDDT scores from the .cif file
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # Look for the pLDDT score in the atom's associated information
                    try:
                        # Extract the plDDT value directly from the atom's attributes
                        # In a .cif file, it may be part of the atom-site data (e.g., '_atom_site.pLDDT')
                        plddt = atom.get_bfactor()  # Attempt to access the B-factor if it's stored there
                        if plddt is not None:
                            atom.bfactor = float(plddt)  # Assign the B-factor from the plDDT score
                        else:
                            atom.bfactor = 0.0  # Default to 0.0 if no pLDDT score found
                    except AttributeError:
                        atom.bfactor = 0.0  # If error occurs, default to 0.0

    # Save the updated structure to a PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    print(f"Converted {cif_file} to {output_pdb} with pLDDT mapped to B-factor.")

# Input and output file paths
input_cif = "fold_gtru_dcr_model_0.cif"  # Replace with your .cif file
output_pdb = "fold_gtru_dcr__model_0.pdb"  # Replace with desired .pdb file name

# Convert .cif to .pdb with pLDDT mapping
cif_to_pdb_with_plddt(input_cif, output_pdb)