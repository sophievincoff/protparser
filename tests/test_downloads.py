from protparser.rcsb import download_rcsb
from Bio.PDB.MMCIFParser import MMCIFParser

# test downloading a structure and make sure it can be parsed
pdb_id = "7UGW"
output_dir = "downloads"
cifparser = MMCIFParser()
download_rcsb(pdb_id, struct_format="cif", output_dir=output_dir)

# try parsing
structure = cifparser.get_structure(pdb_id, f'{output_dir}/{pdb_id}.cif')
print("parsed!")