import requests
import json
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser
import requests
import re
import os
import Bio.PDB as PDB
from bs4 import BeautifulSoup

class AlphaFoldStructure:
    '''
    This class processes an mmCIF file, either uploaded or downloaded from the AlphaFold2 database, to provide comprehensive information.
    '''
    def __init__(self, fold_path=None, uniprot_to_download=None, uniprot_output_dir= None, secondary_structure_types=None, af3=False):
        # If the user provided a PDB path, convert their file to mmcif. Isolate the suffix
        if fold_path is not None:
          fold_fname = fold_path.split('/')[-1]
          prefix, suffix = fold_fname.split('.')

          if suffix == 'pdb': # convert to cif
            # make a directory for converted cif files
            conversion_path = 'mmcif_converted_files'
            if not(os.path.exists(conversion_path)):
              os.makedirs(conversion_path)

            fold_path = self._convert_pdb_to_mmcif__(fold_path, f'{conversion_path}/{prefix}.cif')

        self.file_path = fold_path

        # If user provided a uniprot ID to download, download it and save it as the file path so it can be processed
        if uniprot_to_download is not None:
            if fold_path is not None:
              print("WARNING: both a fold_path and a uniprot_to_download were provided. Running default: downloading the CIF file for provided UniProt ID.")
            self.file_path = self._download_mmCIF(uniprot_to_download, output_path=uniprot_output_dir)

        # Either they provide acceptable secondary structure types, or query the internet for them
        if secondary_structure_types is None:
          self.secondary_structure_types = self._pull_secondary_structure_types()
        else:
          self.secondary_structure_types = secondary_structure_types

        # If there's a CIF file, initialize the object
        if self.file_path:
            self.af3 = af3
            self.cif_lines = self._parse_cif()
            self.secondary_structures = self._extract_secondary_structures()
            self.structure_dict = self._calc_pLDDTs(af3=af3)
            self.sequence = self.structure_dict['seq']
            self.plddts = self.structure_dict['res_pLDDTs']
            self.avg_pLDDT = self.structure_dict['avg_pLDDT']
            self.residues_df = self._create_residues_summary_dataframe()
            self.secondary_structures_df = self._create_secondary_structures_summary_dataframe()
        # Otherwise, print an error.
        else:
            print("ERROR: structure could not be created. No CIF file found.")

    def _convert_pdb_to_mmcif__(self, pdb_filename, mmcif_filename):
      parser = PDB.PDBParser()
      structure = parser.get_structure('structure', pdb_filename)

      io = PDB.MMCIFIO()
      io.set_structure(structure)
      io.save(mmcif_filename)
      return mmcif_filename

    def _download_mmCIF(self, uniprot_id, output_path=None):
        '''
        Download mmCIF file with provided uniprot_id and optional output_path for the downloaded file.

        Return: path to downloaded file if successful, None otherwise
        '''
        full_file_name = f"AF-{uniprot_id}-F1-model_v4.cif"     # define file name that will be found on the AlphaFold2 database.
        # if output path not provided, just save locally under full_file_name
        if output_path is None:
            output_path = full_file_name

        # request the URL for the file 
        url = f"https://alphafold.ebi.ac.uk/files/{full_file_name}"
        response = requests.get(url)

        if response.status_code == 200:
            with open(output_path, 'wb') as file:
                file.write(response.content)
            #print(f"File downloaded successfully and saved as {output_path}")
        else:
            print(f"Failed to download file. Status code: {response.status_code}")
            return None
        
        return output_path

    def _pull_secondary_structure_types(self):
        '''
        Pull a dictionary of secondary structure types and their descriptions from the PDB mmCIF website (necessary for annotating the CIF file)
        Only called if the user does not provide such a dictionary themselves.
        '''

        # request the .html tree from the website with all secondary structure terms
        url = "https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_struct_conf_type.id.html"
        response = requests.get(url)
        
        if response.status_code != 200:
            raise Exception("Failed to retrieve mmCIF dictionary")
        
        # Parse the response content
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Debug: Print the soup to understand the structure
        # print(soup.prettify())
        # write the prettified soup to a txt file
        with open('mmcif_dictionary.txt', 'w') as f:
            f.write(soup.prettify())

        # Find the h4 header with the class "panel-title" and text "Controlled Vocabulary"
        header = soup.find('h4', class_='panel-title')
        if header is None or 'Controlled Vocabulary' not in header.text:
            raise Exception("Could not find the 'Controlled Vocabulary' header")
        
        # Debug: Print the found header
        #print(f"Found header: {header}")

        # The table should be the next sibling of the header
        table = header.find_next('table')
        if table is None:
            raise Exception("Could not find the table following the 'Controlled Vocabulary' header")
        
        # Debug: Print the found table (only the opening <table> tag)
        #print(f"Found table (showing header line): {str(table).split('<thead')[0]}")

        # Iterate through rows in the table and process each entry
        secondary_structure_types = {}
        rows = table.find_all('tr')
        for row in rows[1:]:  # Skip the header row
            cols = row.find_all('td')
            if len(cols) > 1:
                type_id = cols[0].text.strip()
                description = cols[1].text.replace('\t', ' ').strip()
                
                # Replace multiple spaces with a single space
                description = re.sub(' +', ' ', description)
                
                # If this is a protein secondary structure (the table also contains nucleic acid structures), add it to teh dictionary
                if '(protein)' in description:
                  secondary_structure_types[type_id] = description
        
        return secondary_structure_types

    def get_secondary_structure_types(self):
      '''
      Display secondary structure types
      '''
      print("Secondary Structure Types in mmCIF files:")
      for ss_type, description in self.secondary_structure_types.items():
          print(f"{ss_type}: {description}")

      return self.secondary_structure_types

    def _parse_cif(self):
        '''
        Read cif file lines from self.file_path
        '''
        with open(self.file_path, 'r') as file:
            lines = file.readlines()
        return lines

    def _extract_secondary_structures(self):
        '''
        Iterate through the lines of the cif files to find each secondary structure.
        Returns a tuple for each amino acid that has a secondary structure annotation. Tuple contains:
          1. Structure Type (e.g. STRN)
          2. Structure ID (e.g. STRN1)
          3. Description (e.g. beta strand)
          4. Position (e.g. 3)
        '''
        secondary_structures = []
        parsing_secondary_structure = False

        # iterate throhugh cif lines 
        for line in self.cif_lines:
            # hone in on the right section of the cif file
            if line.startswith("_struct_conf.conf_type_id"):
                parsing_secondary_structure = True
                continue
            # if we're in the right section...
            if parsing_secondary_structure:
                if line.startswith("#"):
                    parsing_secondary_structure = False   # no longer in the right section
                    continue
                # still in the right section
                columns = line.split()
                # iterate through columns to find each piece of info we need
                if len(columns) >= 7:
                    sec_struc_type = columns[6]
                    sec_struc_id = columns[13]
                    start_res = int(columns[2])
                    end_res = int(columns[9])
                    sec_struc_name = self.secondary_structure_types.get(sec_struc_type, 'Unknown')
                    # make tuple for this position in the sequence
                    for pos in range(start_res, end_res + 1):
                        secondary_structures.append((sec_struc_type, sec_struc_id, sec_struc_name, pos))
        
        return secondary_structures

    def _calc_pLDDTs(self, af3=False):
        '''
        This method iterates through the cif file to return a dictionary with a few key pieces of info:
          1. Sequence
          2. pLDDTs for each residue
          3. Average pLDDT
        If the structure is from AlphaFold3, pLDDT calculation is a bit different.
        '''

        # define dictionary needed to translate into single-letter AA code
        aa_dict = {
            "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
            "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
            "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
            "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
        }
        
        parser = MMCIFParser(QUIET=True)    # create a parser
        data = parser.get_structure("structure", self.file_path)    # parse structure

        # count models and chains (should be 1 model and 1 chain; don't use this class to parse a complex)
        model = data.get_models()
        models = list(model)
        chains = list(models[0].get_chains())

        # iterate through the chains and get amino acid letters and pLDDTs
        all_pLDDTs = []
        for n in range(len(chains)):
            chainname = chr(n + 65)     # turn chain number into letter (e.g. 1 --> "A" so we have Chain A instead of Chain 1)
            residues = list(chains[n].get_residues())   # extract all residues
            seq = ''
            pLDDTs = [0] * len(residues)    # initialize empty pLDDT array for this chain

            # iterate through all residues in this chain
            for i in range(0, len(residues)):
                r = residues[i]
                try:
                    seq += aa_dict[r.get_resname()]
                except KeyError:
                    print('residue name invalid')
                    badchain = True
                    break

                # Find the alpha carbon (CA) atom in the residue
                if af3:
                    alpha_carbon_bfactor = None
                    for atom in r.get_atoms():
                        if atom.get_name() == 'CA':
                            alpha_carbon_bfactor = atom.get_bfactor()
                            break

                    if alpha_carbon_bfactor is None:
                        #print(f'Alpha carbon not found in residue {r}')
                        badchain = True
                        continue

                    # Assign the B-factor of the alpha carbon as the pLDDT value for this residue
                    pLDDTs[i] = alpha_carbon_bfactor
                else:
                    atoms = list(r.get_atoms())

                    bfactor = atoms[0].get_bfactor()
                    badchain = False

                    for a in range(0, len(atoms)):
                        if atoms[a].get_bfactor() != bfactor:
                            badchain = True

                    pLDDTs[i] = bfactor

            all_pLDDTs.extend(pLDDTs) # add pLDDTs for this chain to list of all pLDDTs

        avg_pLDDT = np.mean(all_pLDDTs) # average pLDDTs across all chains
        return_dict = {
            'avg_pLDDT': round(avg_pLDDT, 2),
            'res_pLDDTs': all_pLDDTs,
            'seq': seq
        }
        return return_dict

    def _create_residues_summary_dataframe(self):
        '''
        Create a dataframe that summarizes the secondary structure information for each residue. 
        Columns:
          1. Position: amino acid position (e.g. 3)
          2. Residue: amino acid 1-letter code (e.g. A)
          3. pLDDT: alphafold2's pLDDT score for this residue to 2 decimal places (e.g. 77.54)
          4. Structure Type: type of secondary structure (e.g. STRN)
          5. Structure ID: ID of this secondary structure (e.g. STRN1)
          5. Description: description of this secondary structure (e.g. beta strand)
          6. Disordered: is this residue disordered or not? A residue is not disordered if it's in a HELX or STRN. (True/False)
        
        '''
        # Convert the secondary structures to a dataframe
        df_secondary_structures = pd.DataFrame(self.secondary_structures, columns=['Structure Type', 'Structure ID', 'Description', 'Position'])

        # Add Residue and pLDDT columns to the dataframe
        df_temp = pd.DataFrame(
            data={
                'Position': list(range(1, len(self.sequence) + 1)),
                'Residue': list(self.sequence),
                'pLDDT': self.plddts
            })

        df_secondary_structures = pd.merge(df_secondary_structures, df_temp, on='Position', how='right')
        # Determine if each residue is disordered or not based on what Structure Type it's in. If helix or strand, it's ordered. If anything else or NaN, it's disordered.
        df_secondary_structures['Disordered'] = df_secondary_structures['Structure Type'].apply(
            lambda x: False if (type(x)==str and (('HELX' in x) or ('STRN' in x))) else True
        )

        return df_secondary_structures

    def _create_secondary_structures_summary_dataframe(self):
        '''
        Create a dataframe grouped by each Structure ID, providing a summary of each secondary structure in the chain. 
        Columns:
          1. Structure ID: ID of this secondary structure (e.g. STRN1)
          2. Start: start position of this secondary structure (e.g. 3)
          3. End: end position of this secondary structure (e.g. 12)
          4. Start Residue: amino acid 1-letter code of the start position (e.g. A)
          5. End Residue: amino acid 1-letter code of the end position (e.g. L)
          6. Disordered: is this residue disordered or not? A residue is not disordered if it's in a HELX or STRN. (True/False)
          7. Description: description of this secondary structure (e.g. beta strand)
          8. Structure Type: type of secondary structure (e.g. STRN)
          9. avg_pLDDT: average pLDDT for this secondary structure (e.g. 77.54)
        '''

        # Apply groupby on self.residues_df to reorganize it by Structure ID
        secondary_structures_df = self.residues_df.groupby('Structure ID').agg({
            'Position': ['first', 'last'],
            'Residue': ['first','last'],
            'Disordered': 'first',
            'Description': 'first',
            'Structure Type': 'first',
            'pLDDT': 'mean'
        }).reset_index()

        # Flatten the multi-level columns
        secondary_structures_df.columns = ['Structure ID', 'Start', 'End', 'Start Residue', 'End Residue', 'Disordered', 'Description', 'Structure Type', 'avg_pLDDT']
        secondary_structures_df['avg_pLDDT'] = secondary_structures_df['avg_pLDDT'].round(2)

        # Display the summarized DataFrame
        return secondary_structures_df

    def get_residues_df(self):
        return self.residues_df

    def get_secondary_structures_df(self):
        return self.secondary_structures_df

    def get_full_sequence(self):
        return ''.join([res for res in self.residues_df['Residue'] if res != 'X'])

    def get_average_plddt(self):
        plddt_values = [plddt for plddt in self.residues_df['pLDDT'] if plddt is not None]
        return sum(plddt_values) / len(plddt_values) if plddt_values else None
    
def pull_secondary_structure_types():
    url = "https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_struct_conf_type.id.html"
    response = requests.get(url)
    
    if response.status_code != 200:
        raise Exception("Failed to retrieve mmCIF dictionary")
    
    soup = BeautifulSoup(response.content, 'html.parser')
    
    # Debug: Print the soup to understand the structure
    # print(soup.prettify())
    # write the prettified soup to a txt file
    with open('mmcif_dictionary.txt', 'w') as f:
        f.write(soup.prettify())

    # Find the h4 header with the class "panel-title" and text "Controlled Vocabulary"
    header = soup.find('h4', class_='panel-title')
    if header is None or 'Controlled Vocabulary' not in header.text:
        raise Exception("Could not find the 'Controlled Vocabulary' header")
    
    # Debug: Print the found header
    #print(f"Found header: {header}")

    # The table should be the next sibling of the header
    table = header.find_next('table')
    if table is None:
        raise Exception("Could not find the table following the 'Controlled Vocabulary' header")
    
    # Debug: Print the found table (only the opening <table> tag)
    #print(f"Found table (showing header line): {str(table).split('<thead')[0]}")

    secondary_structure_types = {}
    rows = table.find_all('tr')
    for row in rows[1:]:  # Skip the header row
        cols = row.find_all('td')
        if len(cols) > 1:
            type_id = cols[0].text.strip()
            description = cols[1].text.replace('\t', ' ').strip()
            
            # Replace multiple spaces with a single space
            description = re.sub(' +', ' ', description)
            
            if '(protein)' in description:
              secondary_structure_types[type_id] = description
    
    return secondary_structure_types