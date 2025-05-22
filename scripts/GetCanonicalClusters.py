#import PDBx
#import matplotlib.pyplot as plt
#import pandas as pd
#import random

import pickle as pck
import pandas as pd
import ssaxgeo.MathToolBox as MathToolBox
import seaborn as sns
import matplotlib.pyplot as plt
import hdbscan
import numpy as np
import argparse
import os

### FUNCTIONS ##################################################################

def remove_outliers(grp):
    '''Remove outliers from dataframe'''
    # PS: This was designed for plotting porpose.
    # set tresholds
    curv_max = grp.grp_xdata_df['curv']< 3.0
    tor_min = grp.grp_xdata_df['tor'] > -0.4
    tor_max = grp.grp_xdata_df['tor'] < 0.25
    wri_min = grp.grp_xdata_df['wri'] > -0.3
    wri_max = grp.grp_xdata_df['wri'] < 0.35
    is_cont = grp.grp_xdata_df['is_cont'] == True

    # update attribute
    grp.grp_xdata_df = grp.grp_xdata_df.loc[curv_max & tor_min & tor_max
                                            & wri_min & wri_max & is_cont]

def cluster_grp_frags(frag_df, max_d=1.0, show_plot=False,
    cols = ['c_mean', 'w_mean', 't_mean'], method='kmeans', **cls_kwargs):
    '''Recluster group fragment dataframe.'''
    #TODO add sanity check
    VALID_METHODS = ['kmeans', 'hdbscan', 'hierarc']
    assert(method in VALID_METHODS), '{} not supported'.format(method)

    # get data to cluster
    #vld_frag_df = self.grp_frag_df[cols]
    vld_frag_df = frag_df[cols]

    if method == 'hierarc':
            # TODO improve kwargs handling, instead of using a default value
            #      maybe will be better to let the user have more control
            max_d = cls_kwargs.get('max_d', 1.0)
            # Calculate linkage matrix
            Z = MathToolBox.get_linkage_mtx_of(vld_frag_df)

            # get clusters assignment array
            clusters = fcluster(Z, max_d, criterion='distance')

            # generate plots
            if show_plot is True:
                MathToolBox.fancy_dendrogram(Z,
                    truncate_mode='lastp', p=12,# show only the last p merged clusters
                    leaf_rotation=90.,
                    leaf_font_size=12.,
                    show_contracted=True,
                    annotate_above=0.1,  # useful in small plots so annotations don't overlap
                    max_d=max_d)

    if method == 'kmeans':
            # clustering specific arguments
            err_msg = 'n_clusters = <int> must be provided for kmeans'
            assert('n_clusters' in cls_kwargs), err_msg

            n_clusters = cls_kwargs['n_clusters']
            random_state = cls_kwargs.get('random_state', 0)
            n_jobs= cls_kwargs.get('n_jobs', 1)

            # do the clustering
            clusters = MathToolBox.run_kmeans_at(vld_frag_df,
                        n_clusters, random_state=random_state, n_jobs=n_jobs)

    if method == 'hdbscan':
            # got clustering arguments
            min_cluster_size = cls_kwargs['min_cluster_size']
            do_soft = cls_kwargs.get('do_soft', True)
            use_distMTX = cls_kwargs.get('use_distMTX', False)
            allow_single_cluster=cls_kwargs.get('allow_single_cluster', True)

            # do the clustering
            clusterer = MathToolBox.do_hdbscan_at(vld_frag_df,
                            min_cluster_size,do_soft=do_soft,
                            use_distMTX=use_distMTX,
                            allow_single_cluster=allow_single_cluster)
            clusters = clusterer.labels_

    ### NOTE ###########################################################
    ## for the future, add clusters probabilities for soft clustering
    ####################################################################

    # add clustering results column to group fragments dataframe
    cluster_col = pd.DataFrame(clusters, columns=['grp_frag_clusters'],
                               index=vld_frag_df.index)
    final_df = pd.concat([frag_df, cluster_col], axis=1)
    if method == 'hdbscan' and do_soft == True:
            return final_df, clusterer
    else:
            return final_df

def get_core_labels(df, Pc_LIM=0.60):
    '''
    Check if a residue belongs to the 'core' of the cluster.
    The core points in a cluster are consider each an every point which the
    membership probability is greater than a specified limit (default=0.60).
    '''
    def _check_prob(row):
        Pc_max = np.array(row['membership_vec']).max()
        if Pc_max >= Pc_LIM:
            return True
        else:
            return False

    return df.apply(_check_prob, axis=1)


# ~~~ pandas handy functions ~~~
def load_entries(row):
    entry = PDBx.entry(coord_flpath=row['chain_flspth'],
            xgeo_flpath=row['xgeo_flspth'],
            pdbid=row['pdbid'], chain=row['chain'])
    entry.load_dssp_data(row['dssp_flspth'])
    return entry

# ~~~ functions for ploting~~
def __plot_grp_frag_clusters(df, x_col, y_col):
    for c in df['grp_frag_clusters'].unique():
        clstr_df = df.loc[df['grp_frag_clusters']==c]
        plt.scatter(clstr_df[x_col],clstr_df[y_col],s=1, label=c)
    plt.legend()

#-------------------------------------------------------------------------------

#### INPUT #####################################################################
desc = '''
This script run the bootstrap solution to identify to the canonical regions
based on a representative set of protein structures.
'''
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('grp_flpath', type=str,
                help='pickle file path for a PDBx.group object')
parser.add_argument('-save_dir', type=str, default=os.getcwd()+"/",
                help='directory to save output files (default= working dir)')

script_dir = os.path.dirname(os.path.realpath(__file__))

parser.add_argument('-kmeans_nclstrs', type=int, default=6,
    help='''The number of the cluster to be used at k-means (default=6).''')

parser.add_argument('-pp2_cluster_n', type=str, default=2,
    help='''The number of the cluster which corresponds to PP2(default=2).
    [WARNING: check the plots!]''')

parser.add_argument('-min_cluster_size', type=int, default=15,
  help='''
  Minimum number of points that HDBSCAN should consider for a cluster
  (default=15)
  ''')
parser.add_argument('-pp2_max', type=float, default=0.07,
  help='Maximum distance from canonical PP2 to consider as PP2 (default=0.07)')

parser.add_argument('-Pc_MIN', type=float, default=0.4,
  help='''
  Minimum membership probability for a point to be considered as part of the
'core' (default=0.4)''')
parser.add_argument('-alpha_cluster_n', type=int, default=2,
  help='''cluster number assigned as alpha canonical region''')
parser.add_argument('-pi_cluster_n', type=int, default=0,
  help='''cluster number assigned as pi canonical region''')
parser.add_argument('-three_cluster_n', type=int, default=1,
  help='''cluster number assigned as 3(10) canonical region''')
parser.add_argument('-show_plot', action='store_true', default=False,
  help='''show generated plots at runing (default=False)''')
parser.add_argument('-plot_debug_data', action='store_true', default=False,
  help='''generate plots for intermediate steps for debuging (default=False)''')

args = parser.parse_args()
parser.print_help()

################### INPUT ######################################################
# general
grp_flpath = args.grp_flpath#"/home/users/amarinho/PDBx_data/grp_state_bc-90.p"
save_dir = args.save_dir#"/home/users/amarinho/PDBx_data/canonical/"
show_plot = args.show_plot
plot_debug_data = args.plot_debug_data
#pp2 canonical settings
pp2_kmeans_nclusters = args.kmeans_nclstrs
pp2_cluster_n = args.pp2_cluster_n

# pi, alpha and 3(10) canonical settings
min_cluster_size = args.min_cluster_size # 200 (defautl)
Pc_LIM = args.Pc_MIN #0.50
alpha_cluster_n = args.alpha_cluster_n#2
pi_cluster_n = args.pi_cluster_n#0
three_cluster_n = args.three_cluster_n#1
################################################################################

# 1 - load dataframe (BCXX)
print('-----------------------------------------------------------------------')
print('@ loading group object (', grp_flpath,')')
grp = pck.load(open(grp_flpath, 'rb'))
print('>> ', len(grp.entries), 'total loaded entries')
print(':: DONE ::')

### Bootstrap strategy #########################################################
# 2 - cluster via fragments
#   Here we cluster the fragments by considering the mean values of each
#   fragment. The clusters obtained will be used to define canonical regions.

print('@ Filtering fragments...')
# 2.1 - Filter fragments
# 2.1.1- Remove extreme values
print('  > removing extreme values...')
grp_frag_df = grp.grp_frag_df
t_max_cond = grp_frag_df['t_mean']<0.2
t_min_cond = grp_frag_df['t_mean']>-0.3
w_max_cond = grp_frag_df['w_mean']<0.3
w_min_cond = grp_frag_df['w_mean']>-0.3

grp_fdf_1 = grp_frag_df.loc[t_max_cond & t_min_cond & w_max_cond & w_min_cond]
print('  :: ', len(grp_fdf_1), ' fragments remaining')


# 2.1.2 - Remove fragments with higher standard deviation
#   The idea is to keep only the most 'helical', we can filter by c_std and
# t_std
print('  > removing fragments with high standard deviations...')
c_std_cond = grp_fdf_1['c_std'] <= 0.07
t_std_cond = grp_fdf_1['t_std'] <= 0.01
grp_fdf_2 = grp_fdf_1.loc[c_std_cond & t_std_cond]
print('  :: ', len(grp_fdf_2), ' fragments remaining')

# Generate canonical regions
# 3 - get PP2 canonical
#   PP2 is on the beta region region and currently a different clustering
#   strategy is used. In summary:
print('|---------| PP2 canonical region |---------|')
print('@ clustering |w|<0.1 region (k-means[n_clusters=',pp2_kmeans_nclusters
     ,'])...')
#       1) hierachical clustering of all fragments in which |w|<0.1
wless_cond = ((grp_fdf_2['w_mean'] < 0.1) & (grp_fdf_2['w_mean'] > -0.1))
grp_df_wlt_1 = cluster_grp_frags(grp_fdf_2.loc[wless_cond],
                                 n_clusters=pp2_kmeans_nclusters)

print('@ selecting cluster as pp2 canonical [c =', pp2_cluster_n,']...')
#       2) select the cluster which correspond to PP2 region
pp2_can = grp_df_wlt_1.loc[(grp_df_wlt_1['grp_frag_clusters'] == pp2_cluster_n)]
print('  :: ', len(pp2_can), ' fragments at PP2 canonical region')
print('@ saving pp2 canonical dataframe...')
#       3) save this cluster as canonical dataframe

pck.dump(pp2_can, open(save_dir+'pp2_can.p', 'wb'))
print('  > saved at:',save_dir+'pp2_can.p')

#### WARNING ###################################################################
# CLUSTER NUMBERING MAY CHANGE FROM DEFAULT VALUES, CHECK DATA PLOTS AND SELECT
# THE RIGHT NUMBER ACCORDINGLY.
################################################################################

print('|---------| Alpha/Pi/3(10) canonical region |---------|')
# 4 - generate clusters for Pi, Alpha and 3(10)
#   Pi, Alpha and 3(10) helices residues belongs to the region of |w| > 0.1 and
# those  usually forms globular clusters. The strategy for those SSE is:
print('@ clustering |w| > 0.1...')
print('   > filtering ')
wgrtr_1 = grp_fdf_2['w_mean'] >  0.1
wgrtr_2 = grp_fdf_2['w_mean'] < -0.1
n_data = len(grp_fdf_2.loc[wgrtr_1])
print('  :: ', n_data, ' residues after filtering')

# 4.1) consider only frags with wri mean higher than 0.1 and run hdbscan
print('  > running hdbscan [min_cluster_size=',min_cluster_size,']...')
grp_df_wgrtr, clusterer = cluster_grp_frags(grp_fdf_2.loc[wgrtr_1], #| wgrtr_2],
                cols = ['c_mean', 'w_mean', 't_mean'], method='hdbscan',
                min_cluster_size=min_cluster_size)

# 4.2) Compute the membership vectors for the clustering assignment
print('  > computing membership vectors...')
soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
grp_df_wgrtr['membership_vec'] = tuple(soft_clusters)

# 4.3) Determinate which points belongs to the core of the clusters
print('  > selecting data points belonging to clusters core...')
grp_df_wgrtr['is_at_core'] = get_core_labels(grp_df_wgrtr, Pc_LIM=Pc_LIM)
#print(grp_df_wgrtr[['membership_vec','is_at_core']])
core_df = grp_df_wgrtr.loc[grp_df_wgrtr['is_at_core']==True]
print('  :: ', len(core_df), 'total fragments selected' )
#print(np.unique(core_df['grp_frag_clusters'], return_counts=True))

print('   > saving core dataframes at ', save_dir,'...')
alpha_can = core_df.loc[core_df['grp_frag_clusters']==alpha_cluster_n]
pi_can = core_df.loc[core_df['grp_frag_clusters']==pi_cluster_n]
three_can = core_df.loc[core_df['grp_frag_clusters']==three_cluster_n]

# 5 - save canonical dataframes
pck.dump(alpha_can, open(save_dir+'alpha_can.p', 'wb'))
pck.dump(pi_can, open(save_dir+'pi_can.p', 'wb'))
pck.dump(three_can, open(save_dir+'three_can.p', 'wb'))
pck.dump(pp2_can, open(save_dir+'pp2_can.p', 'wb'))

# 6 - plot final assignment
print('@ generating canonical regions final plot')
plt.scatter(pi_can['c_mean'], pi_can['t_mean'], s=1, label=r'$\pi$')
plt.scatter(alpha_can['c_mean'], alpha_can['t_mean'], s=1, label=r'$\alpha$')
plt.scatter(three_can['c_mean'], three_can['t_mean'], s=1, label=r'$3_{10}$')
plt.scatter(pp2_can['c_mean'], pp2_can['t_mean'], s=.1, label=r'PP2')
plt.title('FINAL CANONICAL ASSIGNMENT')
plt.xlabel(r'$\kappa [\AA^{-1}]$')
plt.ylabel(r'$\tau [\AA^{-1}]$')
plt.legend()
plt.savefig(save_dir+'FINAL_CANONICAL_KT.png', dpi=300)
print('   > ', save_dir+'FINAL_CANONICAL_KT.png')
if show_plot == True:
    plt.show()
plt.close()

plt.scatter(pi_can['w_mean'], pi_can['t_mean'], s=1, label=r'$\pi$')
plt.scatter(alpha_can['w_mean'], alpha_can['t_mean'], s=1, label=r'$\alpha$')
plt.scatter(three_can['w_mean'], three_can['t_mean'], s=1, label=r'$3_{10}$')
plt.scatter(pp2_can['w_mean'], pp2_can['t_mean'], s=.1, label=r'PP2')
plt.title('FINAL CANONICAL ASSIGNMENT')
plt.xlabel(r'$w$')
plt.ylabel(r'$\tau [\AA^{-1}]$')
plt.legend()
plt.savefig(save_dir+'FINAL_CANONICAL_WT.png', dpi=300)
print('   > ', save_dir+'FINAL_CANONICAL_WT.png')
if show_plot == True:
    plt.show()
plt.close()

### PLOT INTERMEDIATES STEPS FOR DEBUGING
#ORGANIZE INTERMEDIATE PLOTS (ONE PLOT FOR EACH STEP =D)
if plot_debug_data == True:
    print('@ generating debug plots...')
    print(' >>> raw data:')
    # plot raw data first
    #print('@ ploting raw data...')
    plt.scatter(grp_frag_df['c_mean'], grp_frag_df['t_mean'],s=1)
    plt.title(r'raw data')
    plt.xlabel(r'$\kappa [\AA^{-1}]$')
    plt.ylabel(r'$\tau [\AA^{-1}]$')
    plt.savefig(save_dir+'DEBUG_non_round_raw_KT.png', dpi=300)
    print('   > ',save_dir+'DEBUG_non_round_raw_KT.png')
    if show_plot==True:
        plt.show()
    plt.close()

    plt.scatter(grp_frag_df['w_mean'], grp_frag_df['t_mean'],s=1)
    plt.title('Non round raw data')
    plt.xlabel(r'$w$')
    plt.ylabel(r'$\tau [\AA^{-1}]$')
    plt.savefig(save_dir+'DEBUG_non_round_fdf_WT.png', dpi=300)
    print('   > ',save_dir+'DEBUG_HDBScan.png')
    plt.close()

    # HDBSCAN
    print(' >>> hdbscan data:')
    # for alpha/pi/3(10) assignment
    n_clusters = len(np.unique(grp_df_wgrtr['grp_frag_clusters']))
    color_palette = sns.color_palette('deep', n_clusters)
    cluster_colors = [sns.desaturate(color_palette[np.argmax(x)], np.max(x))
                    for x in grp_df_wgrtr['membership_vec'].values]
    #cluster_colors = [color_palette[np.argmax(x)] for x in soft_clusters]
    plt.scatter(grp_df_wgrtr['c_mean'], grp_df_wgrtr['t_mean'], s=.1,
                color=cluster_colors)
    #plt.scatter(grp_df_wgrtr['c_mean'], grp_df_wgrtr['t_mean'],
    #            s=.1, linewidth=0, color='grey', alpha=0.55)
    plt.title(r'hdbscan [$|w|>0.1$; min_cluster_size='+str(min_cluster_size)+']')
    plt.savefig(save_dir+'DEBUG_HDBScan.png', dpi=300)
    print('   > ', save_dir+'DEBUG_HDBScan.png')
    if show_plot == True:
        plt.show()
    plt.close()

    print(' >>> core selection:')
    __plot_grp_frag_clusters(core_df, 'c_mean','t_mean')
    plt.savefig(save_dir+'DEBUG_hdbscan_core_selected_KT.png', dpi=300)
    print('   > ', save_dir+'DEBUG_hdbscan_core_selected.png')
    if show_plot == True:
        plt.plot()
    plt.close()

    __plot_grp_frag_clusters(core_df, 'w_mean','t_mean')
    plt.savefig(save_dir+'DEBUG_hdbscan_core_selected_WT.png', dpi=300)
    print('   > ', save_dir+'DEBUG_hdbscan_core_selected_WT.png')
    if show_plot == True:
        plt.plot()
    plt.close()

    # plot clustered data
    print(' >>> KMEANS:')
    __plot_grp_frag_clusters(grp_df_wlt_1,'w_mean','t_mean')
    plt.title('K-MEANS [$\| w \| < 0.1$; n_clusters=6]')
    plt.xlabel(r'$\kappa [\AA^{-1}]$')
    plt.ylabel(r'$\tau [\AA^{-1}]$')
    plt.savefig(save_dir+'DEBUG_KMEANS_KT.png', dpi=300)
    print('   > ', save_dir+'DEBUG_KMEANS_KT.png')
    plt.close()

    __plot_grp_frag_clusters(grp_df_wlt_1,'w_mean','t_mean')
    plt.title('K-MEANS [$\| w \| < 0.1$; n_clusters=6]')
    plt.xlabel(r'$w$')
    plt.ylabel(r'$\tau [\AA^{-1}]$')
    plt.savefig(save_dir+'DEBUG_KMEANS_WT.png', dpi=300)
    print('   > ', save_dir+'DEBUG_KMEANS_WT.png')
    plt.close()

print(':: DONE ::')
