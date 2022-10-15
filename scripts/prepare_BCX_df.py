import ssaxgeo.PDBx as PDBx
import pandas as pd
import numpy as np
import random
import pickle as pck
import seaborn as sns
import matplotlib.pyplot as plt
##
import ssaxgeo.MathToolBox as MathToolBox
from scipy.cluster.hierarchy import fcluster
import hdbscan
##
import multiprocessing as mp
##
import argparse
import os
##
desc = '''
This script load precompute data from a lPDB copy, get a sample from a given
clustered PDB, compute normalized and smoothed values for each entry,
standardize coordinates and run the SSEx assignment on all residues.

'''
script_dir = os.path.dirname(os.path.realpath(__file__))


parser = argparse.ArgumentParser(description=desc)
parser.add_argument('lpdb_csv',type=str,help='path to a localPDB metadata csv')
parser.add_argument('bc_group_col', type=str,
   help="column name with cluster number for each entry (ex. 'bc-90' for BC90)")
parser.add_argument('-ncpus', type=int, default=1,
    help='number of cpus to use (default=1)')
parser.add_argument('-save_dir', type=str, default=os.getcwd(),
    help='set directory for output files (default=Working dir)')
parser.add_argument('-canonical_dir',type=str,default=script_dir+'/canonical/',
    help='''
set directory of canonical regions to be loaded
(default=canonical folder at this script directory)
''')

args = parser.parse_args()
parser.print_help()

NCPUS = args.ncpus
lpdb_csv_pth = args.lpdb_csv
save_dir = args.save_dir
bc_group = args.bc_group_col
CANONICAL_DIR = args.canonical_dir

### FUNCTIONS ##################################################################

def get_bc_sample(lPDB_df, bc_colnm):
    '''
    filter lPDB dataframe for a representative set of BCX0 group
    Parameters
    ---
    lPDB_df: <pandas dataframe>; metadata of a local PDB metadata
    bc_column: <str>; columns of the dataframe to be used

    Return
    ---
    <pandas dataframe> with a sample obtained
    '''
    df = lPDB_df
    # hard filters [res < 2.0; only proteins and crystallographic structures]
    df_filter = df.loc[(df['res']<=2.0)&
                       (df['res']>0)&
                       (df['method']=='diffraction')&
                       (df['mol_type']=='prot')]

    # get one random entrie per group
    groups = df_filter.groupby(bc_colnm)
    exemplars_lst = []

    for name, group in groups:
        try:
            idxs = random.sample(list(group.index), 1)
        except(ValueError):
            idxs = random.sample(list(group.index), 1)
        exemplars_lst.append(idxs)

    flatten = lambda l: [item for sublist in l for item in sublist]
    return df_filter.loc[flatten(exemplars_lst)]

# ~~! paralelization !~~
def norm_entries(entry):
    # ---- custom values that worked before --------------
    # These values were obtained from a smaller dataset
    c_max = 2.0464588
    c_min = 0.0
    t_max = 0.1688680999999992
    t_min = -0.3012673000000001
    # ----------------------------------------------------
    #print(entry.xdata_df.columns)
    # Normalized KT
    entry.get_normalize_xdf(c_norm=True, k_min=c_min, k_max=c_max, t_min=t_min,
                            t_max=t_max)
    return entry

def smooth_entries(entry):
    # Smooth KTW
    entry.add_smoothed('curv_norm')
    entry.add_smoothed('tor_norm')
    entry.add_smoothed('wri')
    return entry

def standardize_coords(entry):
    entry.standardize_coords()
    return entry

def get_labels(entry):
    entry.get_dist2canonical(pi_df=pi_df, alpha_df=alpha_df, three_df=three_df, pp2_df=pp2_df)
    entry.get_labels(dist_min=0.2, pp2_max=0.07)
    entry.detect_hlx(show_plot=False, selected_cols=['curv_norm_smooth', 'tor_norm_smooth', 'wri_smooth'])
    return entry

# ~~! PANDAS APPLY !~~
def load_entries(row):
    entry = PDBx.entry(coord_flpath=row['chain_flspth'],
            xgeo_flpath=row['xgeo_flspth'],
            pdbid=row['pdbid'], chain=row['chain'])
    entry.load_dssp_data(row['dssp_flspth'])
    return entry

# ~~! SS label algorithm !~~~
alpha_df = pck.load(open(CANONICAL_DIR+'alpha_can.p', 'rb'))
pi_df = pck.load(open(CANONICAL_DIR+'pi_can.p', 'rb'))
three_df = pck.load(open(CANONICAL_DIR+'three_can.p', 'rb'))
pp2_df = pck.load(open(CANONICAL_DIR+'pp2_can.p', 'rb'))

def get_labels(entry):
    entry.get_dist2canonical(pi_df=pi_df, alpha_df=alpha_df, three_df=three_df, pp2_df=pp2_df)
    entry.get_labels(dist_min=0.2, pp2_max=0.07)
    entry.detect_hlx(show_plot=False, selected_cols=['curv_norm_smooth', 'tor_norm_smooth', 'wri_smooth'])
    return entry

#-------------------------------------------------------------------------------

### INPUTS #####################################################################
print('----- INPUTS ----------------------------------------------------------')
print(': NCPUS         = ', NCPUS)
print(': lPDB_CSV_PTH  = ', lpdb_csv_pth)
print(': SAVE_DIR      = ', save_dir)
print(': BC_GROUP      = ', bc_group)
print(': CANONICAL_DIR = ', CANONICAL_DIR)
print('-----------------------------------------------------------------------')

# load lpdb data
print('@ loading entries...')
lPDB_df = pd.read_csv(lpdb_csv_pth, index_col=0)
print('   > ', len(lPDB_df), 'total entries on lPDB metadata csv')
# get BC sample
df_filter = get_bc_sample(lPDB_df, bc_group+'_grp')
# generate group
grp = PDBx.group('test')
# load data
grp.entries = df_filter.apply(load_entries, axis=1)
print('>> ', len(grp.entries), 'total loaded entries')

# 1 - Normalization
print('@ normalization...')
workers = mp.Pool(processes=NCPUS)
output = workers.map(norm_entries,grp.entries.values)
workers.close()

# 2 - smoothing
print('@ smoothing...')
workers = mp.Pool(processes=NCPUS)
output = workers.map(smooth_entries,output)
workers.close()

# 3 - standardize C alpha coordinates
print('@ standardize coordinates...')
workers = mp.Pool(processes=NCPUS)
output = workers.map(standardize_coords,output)
workers.close()

# 4 - get labels
print('@ labeling residues...')
workers = mp.Pool(processes=NCPUS)
output = workers.map(get_labels,output)
workers.close()

# store at group object
grp.entries = output

# 5 - loading group dataframes
print('@ loading group dataframes...')
grp.load_grp_dfs()

# 6 - save outputs
print('@ saving group object...')
pck.dump(grp, open(save_dir+'/grp_state_'+bc_group+'.p', 'wb'))
print(' --> at ', save_dir+'/grp_state_'+bc_group+'.p')
# save raw data
print('@ saving residues dataframe...')
pck.dump(grp.grp_df, open(save_dir+'/grp_res_df_'+bc_group+'.p', 'wb'))
print(' --> at ', save_dir+'/grp_res_df_'+bc_group+'.p')
pck.dump(grp.grp_frag_df, open(save_dir+'/grp_frag_df_'+bc_group+'.p', 'wb'))
print(' --> at', save_dir+'/grp_frag_df_'+bc_group+'.p')
print(':: DONE ::')
