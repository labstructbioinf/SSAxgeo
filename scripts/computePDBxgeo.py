#!/usr/bin/env python3
__author__ = "Antonio Marinho da Silva Neto"
__license__ = "MIT License"
__version__ = "1.0_dev"
__maintainer__ = "Antonio Marinho da Silva Neto"
__email__ = "amarinhosn@pm.me"
__status__ = "Alpha"

import pandas as pd
import multiprocessing as mp
import os
import ssaxgeo.PDBx as PDBx

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

def compute_xgeo(kwargs):
    #TODO check if xgeo file exists
    # if exists, then just return the file
    pdb_chain_flpth = kwargs["pdb_chain_flpth"]
    output_dir = kwargs["output_dir"]
    try:
        # if not, compute xgeo and write xgeo csv
        prefix = pdb_chain_flpth.split("/")[-1].split(".")[0]
        chain_output_dir = f"{output_dir}/{prefix[1:3]}"
        os.makedirs(chain_output_dir, exist_ok=True)
        entry = PDBx.entry(pdb_chain_flpth)
        entry.xdata_df.to_csv(chain_output_dir+"/"+prefix+".csv")
    # TODO investigate those erorrs
    except:
        print(f"Error: {pdb_chain_flpth}")
    #except(AttributeError):
    #    print(f"AttributeError: {pdb_chain_flpth}")

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
    index = sampled_clstrd_pdb_df.index
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
output_dir = mylocalpdb_path+"/xgeo_chain/"
ncpus = 4

# decompress chain gz files if needed
decompressGZFiles(mylocalpdb_path)

# add paths per entry
# load samples clustered pdb
sampled_clstrd_pdb_df = pd.read_csv(sampled_clstrd_path, index_col=0)
sampled_clstrd_pdb_df = addPDBChainsFlPaths(sampled_clstrd_pdb_df, mylocalpdb_path)

# TODO: run melodia on each entry, store xgeo.csv files at mylocalpdb dir
#       check if file already exist to avoid recomputing xgeo data!
chain_fls = sampled_clstrd_pdb_df["chain_flpath"].values
kwargs_to_proc = []
for i in chain_fls:
    kwargs_to_proc.append({"pdb_chain_flpth":i, "output_dir":output_dir})

pool = mp.Pool(ncpus)
print(len(kwargs_to_proc))
pool.map(compute_xgeo, kwargs_to_proc)

# TODO: add xgeo file location colum at sampled clustered

# TODO: export write new columns on csv