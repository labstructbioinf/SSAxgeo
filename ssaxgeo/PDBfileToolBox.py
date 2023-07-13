import os

def get_all_files_in_dir(directory):
    file_paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            file_paths.append(file_path)
    return file_paths


# Sanitize pdb functions
def __add_cryst1_record(input_content):


    # Check if the file already contains a CRYST1 record
    has_cryst1 = any(line.startswith('CRYST1') for line in input_content)

    if not has_cryst1:
        # Create a new CRYST1 record
        cryst1_record = 'CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           \n'.format(0.0, 0.0, 0.0, 90.0, 90.0, 90.0)

        # Add the CRYST1 record to the content
        input_content.insert(0, cryst1_record)
        # Write the modified content to the output file
    return input_content

def __exclude_hetatoms(input_content):
    output_content = []
    for line in input_content:
        if line.startswith("HETATM"):
            continue
        else:
            output_content.append(line)
    return output_content



def __find_residues_with_missing_ca(pdb_file_path):
    with open(pdb_file_path, 'r') as pdb_file:
        atom_lines = [line for line in pdb_file if line.startswith('ATOM')]

    res_CA_dct = {}
    # necessary to compute dihedral angle
    res_C_dct = {}
    res_N_dct = {}

    missing_ca_residues = []

    for line in atom_lines:
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        res_number = line[22:26].strip()
        
        # add keys to dict if absent 
        if res_number not in res_CA_dct.keys():
            res_CA_dct[res_number] = False
        
        if res_number not in res_C_dct.keys():
            res_C_dct[res_number] = False

        if res_number not in res_N_dct.keys():
            res_N_dct[res_number] = False

        # set true if atoms are present
        if atom_name == 'CA':
            res_CA_dct[res_number] = True
        if atom_name == 'C':
            res_C_dct[res_number] = True
        if atom_name == 'N':
            res_N_dct[res_number] = True
        
    for res in res_CA_dct.keys():
        if res_CA_dct[res] == False:
            missing_ca_residues.append(res)
            continue
        if res_C_dct[res] == False:
            missing_ca_residues.append(res)
            continue
        if res_N_dct[res] == False:
            missing_ca_residues.append(res)
            continue
    return missing_ca_residues
def __delete_residues(content, residues_to_exclude):
    output_content = []
    for line in content:
        if line.startswith('ATOM'):
            residue_number = line[22:26].strip()
            if residue_number not in residues_to_exclude:
                output_content.append(line)
        else:
            output_content.append(line)
    return output_content

def __exclude_letter_post_resnumber(input_content, remove=True, fix=False):
    """
    melodia raises index error when a letter is adjacent to
    residue number. so it must be removed
    ex:
    The line bellow crashes Melodia
    "ATOM    377  N   GLY A  80A      24.355 -13.413  10.541  1.00 14.63           N  "
    The line bellow don't
    "ATOM    377  N   GLY A  80       24.355 -13.413  10.541  1.00 14.63           N  "

    """
    output_content = []
    for line in input_content:
        if line.startswith("ATOM"):
            code_for_insert = line[22:29].strip()
            
            if code_for_insert.isdigit() == False:
                if fix == True:
                    new_line = line[:26]+" "+line[27:]
                    output_content.append(new_line)
                    continue
                if remove == True:
                    continue
            else:
                output_content.append(line)
        else:
            output_content.append(line)
    return output_content

def validate_pdb_file(pdb_file):

    errors_found = []
    # get list os res with missing CA atoms
    res_miss_ca_lst = __find_residues_with_missing_ca(pdb_file)

    with open(pdb_file, 'r') as f:
        content = f.readlines()

    # Check if the file starts with the "HEADER" record
    #if not content[0].startswith('HEADER'):
    #    return False

    # Check if the file contains at least one "ATOM" or "HETATM" record
    has_atom_record = any(line.startswith('ATOM') or line.startswith('HETATM') for line in content)
    if not has_atom_record:
        errors_found.append("No ATOM or HETATM records present")

    # Check if the file ends with the "END" record
    if not content[-1].startswith('END'):
        errors_found.append("No END at the end of file")

    # Check if the "MODEL" and "ENDMDL" records, if present, occur in pairs
    model_count = sum(1 for line in content if line.startswith('MODEL'))
    endmdl_count = sum(1 for line in content if line.startswith('ENDMDL'))
    if model_count != endmdl_count:
        errors_found.append("At least one MODEL missing ENDMDL")

    # Check if the "CRYST1" record is present
    has_cryst1_record = any(line.startswith('CRYST1') for line in content)
    if not has_cryst1_record:
        errors_found.append("No CRYST1 record is present")
    return errors_found, res_miss_ca_lst

def write_sanitized_pdb(errors_found_lst, input_pdb, output_pdb, res_miss_ca_lst):
    # check possible errors at input pdb
    if len(errors_found_lst) == 0:
        return False
    
    with open(input_pdb, 'r') as f:
        content = f.readlines()
    # exclude hetatoms
    content = __exclude_hetatoms(content)
    # remove extra letter post resnumber
    content = __exclude_letter_post_resnumber(content) 
    
    if "No CRYST1 record is present" in errors_found_lst:
        # add cryst1 record if needed
        content = __add_cryst1_record(content)
    # remove residues with missing CA
    if len(res_miss_ca_lst) > 0:
        print(f"WARN: residues with missing CA, C or N {res_miss_ca_lst} from {input_pdb} will not be considered")
        content = __delete_residues(content, res_miss_ca_lst)

    with open(output_pdb, 'w') as f:
        f.writelines(content)

def run_dssp(pdb_flpth, dssp_flpth):
    cmd = f"mkdssp --output-forma dssp {pdb_flpth} > {dssp_flpth}"
    os.system(cmd)

#input_pdb = "/home/antonio/Projects/HlxCnt/mypdb/pdb_chain/gp/6gpz_E.pdb"
#output_pdb = "/home/antonio/Projects/HlxCnt/mypdb/pdb_chain/gp/6gpz_E_sntzd.pdb"

#errors_found_lst, res_missing= validate_pdb_file(input_pdb)
#write_sanitized_pdb(errors_found_lst, input_pdb, output_pdb, res_missing)