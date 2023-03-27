import pandas as pd
import multiprocessing as mp
import os

def decompressGZFiles(mylocalpdb_path):
    
    def decompressGZ(gzfile):
        os.system(f"gunzip {gzfile}")

    # get list of files to unzip
    def getGZFilesList(base_dir):
        # get files
        res = []
        for (dir_path, dir_names, file_names) in os.walk(base_dir):
            res.extend([f for f in file_names if f.endswith(".gz")])
        # add paths it assumes pdb_chain subdir structure (ex: 1F09 would be at subdir F0)
        gz_flpaths = [base_dir+f"{fl[1:3]}/"+fl for fl in res]
        return gz_flpaths

    pool = mp.Pool(ncpus)
    gz_list = getGZFilesList(mylocalpdb_path+'pdb_chain/')
    pool.map(decompressGZ, gz_list)

def addPDBChainsFlPaths(sampled_clstrd_pdb_df, mylocalpdb):
    '''
    add a new column with pdb chains path and excludes collumns
    with missing files
    '''

    def getPathFromIdx(i, mylocalpdb):
        '''get pdb chain paths based on indexes'''
        ch_dir = i[1:3]
        ch_pdb = i+".pdb"
        return f"{mylocalpdb}pdb_chain/{ch_dir}/{ch_pdb}"
    # get list of indexes with missing files
    missing_i = []
    chain_flpath_list = []
    for i_fl in [[i, getPathFromIdx(i, mylocalpdb_path)] for i in sampled_clstrd_pdb_df.index]:
        chain_flpath_list.append(i_fl[1])
        if os.path.isfile(i_fl[1]) == False:
            missing_i.append(i_fl[0])
    # create new column
    sampled_clstrd_pdb_df["chain_flpath"] = chain_flpath_list 
    # exclude columns with missing files
    missing_fls_i = sampled_clstrd_pdb_df.loc[missing_i].index
    if len(missing_fls_i) >0:
        print(f"{len(missing_fls_i)} files which chains were not found at {mylocalpdb}")
    return sampled_clstrd_pdb_df.drop(missing_fls_i, axis=0)

mylocalpdb_path = "/home/antonio/Projects/HlxCnt/mypdb/"
sampled_clstrd_path = "/home/antonio/Projects/HlxCnt/sampled_clustered_pdb.csv"
ncpus = 4

# decompress chain gz files if needed
decompressGZFiles(mylocalpdb_path)

# add paths per entry
# load samples clustered pdb
sampled_clstrd_pdb_df = pd.read_csv(sampled_clstrd_path, index_col=0)
addPDBChainsFlPaths(sampled_clstrd_path, mylocalpdb_path)

# TODO: run melodia on each entry, store xgeo.csv files at mylocalpdb dir
#       check if file already exist to avoid recomputing xgeo data!
# TODO: add xgeo file location colum at sampled clustered

# TODO: export write new columns on csv