#!/usr/bin/env python3
__author__ = "Antonio Marinho da Silva Neto"
__license__ = "MIT License"
__version__ = "1.0_dev"
__maintainer__ = "Antonio Marinho da Silva Neto"
__email__ = "amarinhosn@pm.me"
__status__ = "Alpha"

import PDBx
import pickle as pck
import argparse
import os
import pandas as pd

desc = '''
This script run the Secondary Structure Assignment based on differential
geometry (xgeo) descriptors on a specified entry and save results as csv files.
'''
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('pdb_flpath', type=str,
                    help='protein chain coordinate file path')
parser.add_argument('xgeo_flpath', type=str,
                    help='protein chain xgeo file path')
script_dir = os.path.dirname(os.path.realpath(__file__))
parser.add_argument('-can_dir', type=str, default=script_dir+'/canonical/',
                help='SSE canonical dataframes directory (default=script dir)')
parser.add_argument('-out_dir', type=str, default=os.getcwd(),
                help='Output directory (default=working dir)')
parser.add_argument('-min_dist', type=float, default=0.2,
  help='''
  Minimum distance from canonical Pi/Alpha/3(10) to consider as a label
(default=0.2)
  ''')
parser.add_argument('-pp2_max', type=float, default=0.07,
  help='Maximum distance from canonical PP2 to consider as PP2 (default=0.07)')

parser.add_argument('-prefix', type=str, default='',
  help='Prefix for ouput files (default=None)')

args = parser.parse_args()
parser.print_help()

print('|---------------------------------------------------------------------|')
print('|---- Input Files ----|')
pdb_flpath = args.pdb_flpath
print('| pdb_flpath  = ', pdb_flpath)
xgeo_flpath = args.xgeo_flpath
print('| xgeo_flpath = ', xgeo_flpath)
print('|---- Directories ----|')
can_dir = args.can_dir
print('| -can_dir    = ', can_dir)
out_dir = args.out_dir
print('| -out_dir    = ', out_dir )
prefix = args.prefix
print('| -prefix =', prefix)
print('|--- SSAx parameters ----|')
min_dist = args.min_dist
print('| -min_dist   =', min_dist)
pp2_max = args.pp2_max
print('| -pp2_max    =', pp2_max)
print('|---------------------------------------------------------------------|')

# Load entry data
print('@ loading entry...')
entry = PDBx.entry(pdb_flpath,xgeo_flpath)
print('  > nres          = ', entry.nres)
print('  > is continuous = ', entry.cont)

# plot arrows representation
print('@ ploting arrows...')
entry.plot_arrows(save_fig=True, myDIR=out_dir, show_plot=False)
print('@ running SSAx...')
# Load canonical dataframes
print('  >> loading canonical SSE dataframes...')
alpha_can = pck.load(open(can_dir+'/alpha_can.p','rb'))
pi_can = pck.load(open(can_dir+'/pi_can.p','rb'))
three_can = pck.load(open(can_dir+'/three_can.p','rb'))
pp2_can = pck.load(open(can_dir+'/pp2_can.p','rb'))

print('  >> computing residues to canonical distances...')
# get distances to canonical
entry.get_dist2canonical(alpha_df=alpha_can, pi_df=pi_can, three_df=three_can,
                         pp2_df=pp2_can)
print('  >> seting residues labels...')
# get labels
entry.get_labels(dist_min=0.2, pp2_max=0.07)

# write output files
print('@ writing output files...')
# sort file naming
xdata_flnm = 'xdata_df.csv'
ssax_flnm = 'ssax.csv'

if prefix != '':
    xdata_flnm = prefix+'_'+xdata_flnm
    ssax_flnm = prefix+'_'+ssax_flnm

# write xdata dataframe csv
entry.xdata_df.to_csv(out_dir+xdata_flnm)
print('  >> ', out_dir+xdata_flnm)

# write ssax only data dataframe
cols_to_ssax = ['res', 'D(Alfa)', 'D(Pi)', 'D(3(10))', 'D(PP2)','label']
entry.xdata_df[cols_to_ssax].to_csv(out_dir+ssax_flnm)
print('  >> ', out_dir+ssax_flnm)
print(':: DONE ::')
