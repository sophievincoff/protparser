import requests
import json
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO
import requests
import re
import os
from bs4 import BeautifulSoup

class RCSBStructure:
    '''
    This class processes an ID from the RCSB PDB website to provide comprehensive information. 
    '''
    def __init__(self, pdb_id, download_struct=False, struct_format='cif', struct_save_dir=None, print_progress=False):
        self.pdb_id = pdb_id.upper()    # turn PDB ID to uppercase so user doesn't have to
        self.fasta_content = None       # will hold text from .fasta file
        self.chain_info = None          # will hold info extracted from this ID's webpage
        self._pull_fasta_and_chain_info(print_progress=print_progress)    # call method to initialize self.fasta_content, self.chain_info
        
        if download_struct:
            download_rcsb(pdb_id, struct_format=struct_format, output_dir=struct_save_dir)

    def _pull_fasta_and_chain_info(self,print_progress=False):
        '''
        Retrieves:
          1. FASTA file for the PDB ID
          2. Chain information including: 
            a. chain identifiers
            b. author-assigned chain identifiers
            c. UniProt IDs and alignments for each
            d. sequence (canonical and PDB)
            e. Pfam domains
            f. Disordered regions
            g. Disordered binding sites
        '''
        # Get the FASTA file
        fasta_url = f'https://www.rcsb.org/fasta/entry/{self.pdb_id}'
        fasta_response = requests.get(fasta_url)
        if fasta_response.status_code != 200:
            raise Exception(f"Failed to retrieve FASTA file for PDB ID {self.pdb_id}")
        self.fasta_content = fasta_response.text
        
        ## Download #1: full RCSB PDB page for this ID
        # Get the chain information including author-assigned chain identifiers and UniProt IDs
        summary_url = f'https://data.rcsb.org/rest/v1/core/entry/{self.pdb_id}'
        summary_response = requests.get(summary_url)
        if summary_response.status_code != 200:
            raise Exception(f"Failed to retrieve chain information for PDB ID {self.pdb_id}")
        summary_data = summary_response.json()

        # Write the JSON response to a text file
        #with open(f'{self.pdb_id}_summary.json', 'w') as json_file:
            #json.dump(summary_data, json_file, indent=4)
        
        ## Downloads #2-(# entities + 1): download info on each entity ID and associate it to a chain_auth_id
        self.chain_info = {}
        for polymer_entity in summary_data['rcsb_entry_container_identifiers']['polymer_entity_ids']:
            # Download info for this polymer ID, write it to a JSON
            polymer_url = f'https://data.rcsb.org/rest/v1/core/polymer_entity/{self.pdb_id}/{polymer_entity}'
            polymer_response = requests.get(polymer_url)
            if polymer_response.status_code != 200:
                raise Exception(f"Failed to retrieve polymer entity information for PDB ID {self.pdb_id} entity {polymer_entity}")
            polymer_data = polymer_response.json()

            # Write the JSON response to a text file
            #with open(f'{self.pdb_id}_{polymer_entity}_summary.json', 'w') as json_file:
                #json.dump(polymer_data, json_file, indent=4)
            
            # Collect all desired information
            chains = polymer_data['rcsb_polymer_entity_container_identifiers']['asym_ids']
            chains_str = ','.join(chains)       # turn chains into a string (if there are multiple letters for this chain, they'll be stored as a string, not a list)
            auth_asym_ids = polymer_data['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']
            entity_id = polymer_data['rcsb_polymer_entity_container_identifiers']['entity_id']
            pdb_sequence = polymer_data['entity_poly']['pdbx_seq_one_letter_code']
            canonical_sequence = polymer_data['entity_poly']['pdbx_seq_one_letter_code_can']
            
            # Add newly acquired information to self.chain_info dict for this chain
            self.chain_info[chains_str] = {
                'author_chain': auth_asym_ids,
                'entity_id': entity_id,
                'pdb_sequence': pdb_sequence,
                'canonical_sequence': canonical_sequence
            }

            # print statement for sanity check (optional)
            if print_progress:
              print(f'Analyzing Entity ID {entity_id}, Chains {chains_str} (file: {self.pdb_id}_{polymer_entity}_summary.json)')

            # sequence alignment info: conflicts, deletions, insertions, mutations
            self.chain_info[chains_str]['rcsb_conflicts'] = {
                'rcsb_conflict_count': polymer_data['entity_poly']['rcsb_conflict_count'],
                'rcsb_deletion_count': polymer_data['entity_poly']['rcsb_deletion_count'],
                'rcsb_insertion_count': polymer_data['entity_poly']['rcsb_insertion_count'],
                'rcsb_mutation_count': polymer_data['entity_poly']['rcsb_mutation_count'],
                'rcsb_non_std_monomer_count': polymer_data['entity_poly']['rcsb_non_std_monomer_count'],
                'rcsb_sample_sequence_length': polymer_data['entity_poly']['rcsb_sample_sequence_length']
            }

            # Extract sequence alignment info: UniProt IDs and their aligned regions with the seq from the crystal structure
            alignment_info = {}
            # Iterate through alignments 
            for alignment in polymer_data.get('rcsb_polymer_entity_align', []):
                # Continue iff the alignment is UniProt
                if alignment.get('reference_database_name', None) == 'UniProt':
                    uniprot_id = alignment.get('reference_database_accession', None)

                    # Continue iff there's a UniProt ID
                    if uniprot_id:
                        alignment_info[uniprot_id] = {}

                        # Iterate through each aligned region and save the alignment
                        for i, aligned_region in enumerate(alignment.get('aligned_regions', [])):
                            alignment_info[uniprot_id][f'alignment_{i+1}'] = {
                                'entity_beg_seq_id': aligned_region.get('entity_beg_seq_id', None),
                                'length': aligned_region.get('length', None),
                                'ref_beg_seq_id': aligned_region.get('ref_beg_seq_id', None)
                            }
                    
            # There is no uniprot id if alignment_info is empty 
            if len(alignment_info) == 0:
                alignment_info = None

            # Add UniProt IDs and all their alignments to self.chain_info for this chain
            self.chain_info[chains_str]['uniprot_ids'] = alignment_info

            # get Pfam domains, Disordered binding sites, and disordered regions
            pfam_info = {}
            pfam_counter = 1                        # initialize counter for Pfam regions
            disordered_binding_site_info = {}
            disordered_binding_site_counter = 1     # initialize counter for disordered binding sites
            disordered_region_info = {}
            disordered_region_counter = 1           # initialize counter for disordered regions 
            
            # Iterate through features 
            if 'rcsb_polymer_entity_feature' in polymer_data:
                for feature in polymer_data['rcsb_polymer_entity_feature']:
                    # Pfam
                    if feature['type'] == 'Pfam':
                        pfam_info[f'pfam_{pfam_counter}'] = {
                            'feature_id': feature['feature_id'],
                            'name': feature['name']
                        }
                        # Iterate through feature_positions to save each alignment withiin this Pfam region
                        for j, feature_positions in enumerate(feature['feature_positions']):
                            pfam_info[f'pfam_{pfam_counter}'][f'feature_positions_{j+1}'] = {
                                'beg_seq_id': feature_positions['beg_seq_id'],
                                'end_seq_id': feature_positions['end_seq_id']
                            }
                        pfam_counter += 1   # increment the Pfam region counter

                    # Disorder (IUPred2)
                    if feature['type'] == 'disorder':
                        disordered_region_info[f'disorder_{disordered_region_counter}'] = {
                            'source': feature['provenance_source'],
                            'name': feature['name']
                        }
                        # Iterate through feature_positions to save each disordered stretch
                        for j, feature_positions in enumerate(feature['feature_positions']):
                            disordered_region_info[f'disorder_{disordered_region_counter}'][f'feature_positions_{j+1}'] = {
                                'beg_seq_id': feature_positions['beg_seq_id'],
                                'values': feature_positions['values']
                            }
                        disordered_region_counter += 1      # increment disordered region counter

                    # Disorder binding (Anchor2)
                    if feature['type'] == 'disorder_binding':
                        disordered_binding_site_info[f'disorder_binding_{disordered_binding_site_counter}'] = {
                            'source': feature['provenance_source'],
                            'name': feature['name']
                        }
                        # Iterate through feature_positions to save each stretch of disordered binding site info
                        for j, feature_positions in enumerate(feature['feature_positions']):
                            disordered_binding_site_info[f'disorder_binding_{disordered_binding_site_counter}'][f'feature_positions_{j+1}'] = {
                                'beg_seq_id': feature_positions['beg_seq_id'],
                                'values': feature_positions['values']
                            }
                        disordered_binding_site_counter += 1      # increment disordered binding site counter
                        
                    if feature['type'] == 'mutation':   
                        if 'mutations' not in self.chain_info[chains_str]:
                            self.chain_info[chains_str]['mutations'] = set()
                        # get beginning and end, if available 
                        for subfeature in feature['feature_positions']:
                            end_pos = None
                            start_pos = subfeature['beg_seq_id']
                            if 'end_seq_id' in subfeature:
                                end_pos = subfeature['end_seq_id']
                            else:
                                end_pos = start_pos
                                
                            # now iterate through and add each number
                            for i in range(start_pos, end_pos+1):
                                self.chain_info[chains_str]['mutations'].add(i)
                
            # Add new info to self.chain_info for this chain
            self.chain_info[chains_str]['pfam_info'] = pfam_info
            self.chain_info[chains_str]['disordered_binding_site_info'] = disordered_binding_site_info
            self.chain_info[chains_str]['disordered_region_info'] = disordered_region_info

    def summarize_all_chains(self):
        '''
        Helper method: print a summary for the user of the key information extracted for each chain
        '''
        print("Chain information:")
        for chain, info in self.chain_info.items():
            print("Chain {} (Auth: {})\tUniProt IDs: {}\tPDB Seq ({} AAs): {}\tCanonical Seq ({} AAs): {}\tPfam domains: {}".format(
                chain,
                ','.join(info['author_chain']),
                ','.join(list(info['uniprot_ids'].keys())) if info['uniprot_ids'] != None else "None",
                str(len(info['pdb_sequence'])) + ' ' * (4 - len(str(len(info['pdb_sequence'])))),
                info['pdb_sequence'][0:7] + '...',
                str(len(info['canonical_sequence'])) + ' ' * (4 - len(str(len(info['canonical_sequence'])))),
                info['canonical_sequence'][0:7] + '...',
                len(info['pfam_info'])
            ))

    def _process_item(self, cur_key, item, indent=0):
        """
        Helper method for detailed_chain_summary.
        Recursively process and print items in a dictionary.

        Parameters:
        cur_key (str): The current key being processed.
        item (dict, str, list): The item to process.
        indent (int): The current level of indentation for printing.
        """
        indent_str = '\t' * indent

        if isinstance(item, dict):
            print(f"{indent_str}{cur_key}:")
            for key, value in item.items():
                self._process_item(key, value, indent + 1)
        elif isinstance(item, list):
            print(f"{indent_str}{cur_key}: {','.join([str(x) for x in item])}")
        else:
            print(f"{indent_str}{cur_key}: {item}")

    def detailed_chain_summary(self, chains=None):
        """
        Print out a detailed, organized summary of everything in the chain_info dict.

        Parameters:
        chain_info (dict): The chain_info dictionary containing the data.
        chains (list): The chains whose disordered information you would like to see. 
        """
        chains_to_summarize = self.chain_info
        if chains:
            if type(chains)==list:
                chains_to_summarize = {chain: self.chain_info[chain] for chain in chains}
            else:
                print('Must provide a list [] of chains. Running default: print summary for all chains')

        for key in chains_to_summarize:
            self._process_item(key, chains_to_summarize[key], 0)

    # Getter methods
    def get_chain_list(self):
        '''
        Getter method: prints all the chains and returns a list of them
        '''
        print('|'.join(list(self.chain_info.keys())))
        return list(self.chain_info.keys())

    def get_fasta_content(self):
        '''
        Getter method: returns the fasta text
        '''
        return self.fasta_content

    def get_chain_info(self):
        '''
        Getter method: returns self.chain_info
        '''
        return self.chain_info

    def get_chain_info_for_chain(self, chain):
        ''''
        Getter method: returns self.chain_info for a chain of your choice
        '''
        return self.chain_info.get(chain, None)
    
def download_rcsb_apicall(pdb_id, struct_format='cif',output_dir=None):
    full_file_name = f"{pdb_id}.{struct_format}"     # define file name that will be found on the AlphaFold2 database.
    full_file_name_uppercase = f"{pdb_id.upper()}.{struct_format}"
    # if output path not provided, just save locally under full_file_name
    output_path=full_file_name_uppercase
    if output_dir is not None:
        os.makedirs(output_dir,exist_ok=True)
        output_path = f"{output_dir}/{full_file_name_uppercase}"

    # request the URL for the file 
    url = f"https://files.rcsb.org/download/{full_file_name}"
    response = requests.get(url)
    return response, output_path

def download_rcsb(pdb_id, struct_format='cif', convert_if_fail=False, output_dir=None):
    '''
    Download mmCIF file with provided uniprot_id and optional output_path for the downloaded file.
    Try both the uppercase and lowercase version of the ID, but convert to uppercase if the lowercase version is found.

    Return: path to downloaded file if successful, None otherwise
    '''
    # try both the uppercase and lowercase ID
    uppercase_id = pdb_id.upper()
    lowercase_id = pdb_id.lower()
    
    for formatted_pdb_id in [uppercase_id, lowercase_id]:
        response, output_path = download_rcsb_apicall(pdb_id, struct_format=struct_format,output_dir=output_dir)

        if response.status_code == 200:
            with open(output_path, 'wb') as file:
                file.write(response.content)
            return output_path
            #print(f"File downloaded successfully and saved as {output_path}")
        else:
            # try downloading cif and converting to pdb, if pdb is the format. PDBs aren't available for everything on RCSB
            if struct_format=='pdb' and convert_if_fail:
                response, output_path = download_rcsb_apicall(pdb_id, struct_format='cif',output_dir=output_dir)
                if response.status_code == 200:
                    with open(output_path, 'wb') as file:
                        file.write(response.content)
                    convert_cif_to_pdb(output_path, output_path.replace('.cif','.pdb'))
                    os.remove(output_path)
                    print(f"Deleted original .cif file download. See {output_path.replace('.cif','.pdb')} for the PDB")
                else:
                    print(f"Failed to download {pdb_id} file. Status code: {response.status_code}")
            else:
                print(f"Failed to download {pdb_id} file. Status code: {response.status_code}")
    
    return None   
    
def convert_cif_to_pdb(cif_file: str, pdb_file: str):
    """
    Converts an MMCIF (.cif) file to a PDB (.pdb) file.
    
    Args:
        cif_file (str): Path to the input CIF file.
        pdb_file (str): Path to the output PDB file.
    """
    try:
        # Parse the CIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('structure_id', cif_file)
        
        # Write to PDB format
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)
        
        print(f"Successfully converted {cif_file} to {pdb_file}")
    except Exception as e:
        print(f"Error during conversion: {e}")