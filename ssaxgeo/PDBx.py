import pandas as pd
import numpy as np
import os
import pickle as pck
# plot stuff
import matplotlib.pyplot as plt
import seaborn as sns
#
import ssaxgeo.MathToolBox as MathToolBox
import ssaxgeo.PDBfileToolBox as PDBfileToolBox
from scipy.cluster.hierarchy import fcluster
# third party data
import ssaxgeo.dssp as dssp

# diff geo
import melodia as mel

'''
        >>>> PDBx <<<<

Here we define classes to deal with PDB entries chains xgeo data, as a single
entry and as a group of entries.

'''
### GLOBAL VAR #################################################################
XGEO_ENGINES = ["melodia", "diffgeo"]

### FUNCTIONS ##################################################################

# ~~~ pandas handy functions ~~~
def intfy_res(row):
    return int(row)

def load_xgeo_df(xgeo_flpath):
    '''Load xgeo dataframe from flexgeo output file'''
    dtypes = {
        "label": "category",
        "ID": "category",
        "res": np.int32,
        "atom": "category",
        "res_name": "category",
        "curv": np.float32,
        "tor": np.float32,
        "wri": np.float32,
        "arc": np.float32,
        "D(Alpha)": np.float32,
        "D(Pi)": np.float32,
        "D(3(10))": np.float32,
        "D(PP2)": np.float32,
        "conf": np.int16,
        "aa_idx": np.int32
    }

    xgeo_df = pd.read_csv(xgeo_flpath, index_col=False, dtype=dtypes)#, header=None)
    xgeo_df.drop(["Unnamed: 0"], inplace=True, axis=1)
    xgeo_df.drop(["res_name"], inplace=True, axis=1)

    if "phi" in xgeo_df.columns:
        xgeo_df.drop(["phi", "psi"], axis=1, inplace=True)
    #xgeo_df.columns=['conf_n', 'res', 'curv', 'tor', 'arc','wri',
    #                     'ca_x', 'ca_y', 'ca_z', 'seq']

    return xgeo_df

def custom_norm(X, min, max, a=0, b=1):
    '''Normalize array using custom min and max values'''
    assert(max>min)
    X_clip = np.clip(X, min, max)
    return pd.DataFrame(a+(X_clip - min)*(b-a) / (max - min))

def load_group_df(entries_lst, frag_df = False):
        '''Get a dataframe with all pdbx xdata.'''
        if frag_df == False:
            print('@ mounting group dataframe...')
        if frag_df == True:
            print('@ mounting group frag dataframe...')

        # get xdata and pdbid chains to mount grp_xdata
        df_list = []
        new_cols_list = []
        total = len(entries_lst)
        for i, pdbx_i in enumerate(entries_lst):
            print("  --> {} out of {}".format(i+1,total), end='\r')
            if pdbx_i == 0:
                continue
            if frag_df == False:
                df = pdbx_i.xdata_df
            if frag_df == True:
                df = pdbx_i.frag_df
            try:
                df_list.append(df)
                new_cols_list.append([[pdbx_i.pdbid, pdbx_i.chain]]*len(df))
            except(TypeError):
                print(df)
                print(pdbx_i.xdata_df)
                print(pdbx_i.frag_df)
                print(f"WARNING: skiping {[pdbx_i.pdbid, pdbx_i.chain]} due to no dataframe available")
                continue
        
        # stack individual dfs
        all_xdata = pd.concat(df_list, axis=0, sort=False)
        all_xdata = all_xdata.reset_index(drop=True)

        # prepare pdbid and chain columns to be added
        flat_list = MathToolBox.flatten(new_cols_list)
        new_cols_df = pd.DataFrame(flat_list)

        return pd.concat([all_xdata, new_cols_df], axis=1)

def dssp_pp2_assgn(row):
    '''
    Get pp2 assignment using rules for DSSP-PPII.

    REFERENCE:
    Chebrek, R., Leonard, S., de Brevern, A. G., & Gelly, J. C. (2014).
    PolyprOnline: polyproline helix II and secondary structure assignment
    database. Database : The Journal of Biological Databases and Curation,
    2014(9), 1–8. https://doi.org/10.1093/database/bau102
    '''
    # set conditions
    phi_pp2 = -75
    psi_pp2 = +145
    e = 29

    # get data
    phi_row = float(row['phi'])
    psi_row = float(row['psi'])
    dssp_row = row['ss']

    # check if is PP2
    cond_Phi = (phi_row >= phi_pp2 - e) and (phi_row <= phi_pp2 + e)
    cond_Psi = (psi_row >= psi_pp2 - e) and (psi_row <= psi_pp2 + e)

    if cond_Phi == True and cond_Psi == True and dssp_row=='':
        return 'P'
    else:
        return dssp_row

# ~~~ single entry handy functions ~~~

def cluster_res_by_wfilter(xdata_df, wri_lim=0.1, r_dist=0.3, nr_dist=0.85,
    min_rsize=5, min_nrsize=4, show_plot=False,
    selected_cols=['curv', 'tor', 'wri']):
    '''
    Cluster residues by 1) apply the wri filter and 2) hierarchical
    clustering. This function return two dataframes with indexes compatible
    with input data.
    '''
    # get filtered dataframe
    wri_col = xdata_df['wri'].astype(float)
    round_df=xdata_df.loc[(abs(wri_col) >= wri_lim)]
    not_round_df = xdata_df.loc[(abs(wri_col) < wri_lim)]

    # separate try helix also assigned correctly by DSSP
    # 1 - get data to cluster
    r_data = round_df[selected_cols].dropna()
    nr_data = not_round_df[selected_cols].dropna()

    # 2 - do the hierarchical clustering
    if len(r_data) >= min_rsize:
        Z_round = MathToolBox.get_linkage_mtx_of(r_data)
        clusters_round = fcluster(Z_round, r_dist, criterion='distance')

    if len(r_data) < min_rsize:
        Z_round = 0
        clusters_round =[]

    if len(nr_data) >= min_nrsize:
        Z_not_round = MathToolBox.get_linkage_mtx_of(nr_data)
        clusters_not_round = fcluster(Z_not_round, nr_dist, criterion='distance')

    if len(nr_data) < min_nrsize:
        Z_not_round = 0
        clusters_not_round =[]

    if show_plot is True:
        if type(Z_round) != int:
            MathToolBox.fancy_dendrogram(Z_round,
            truncate_mode='lastp', p=12,# show only the last p merged clusters
            leaf_rotation=90.,
            leaf_font_size=12.,
            show_contracted=True,
            annotate_above=0.1,  # useful in small plots so annotations don't overlap
            max_d=r_dist)
            plt.title('Hierarchical Clustering "Round"')
            plt.show()
            plt.close()

        if type(Z_not_round) != int:
            MathToolBox.fancy_dendrogram(Z_not_round,
            truncate_mode='lastp', p=12,# show only the last p merged clusters
            leaf_rotation=90.,
            leaf_font_size=12.,
            show_contracted=True,
            annotate_above=0.1,#useful in small plots(annotations don't overlap)
            max_d=nr_dist)
            plt.title('Hierarchical Clustering "Not Round"')
            plt.show()
            plt.close()

    # new columns name
    cls_r_colnm = ['clusters_round']
    cls_nr_colnm = ['clusters_not_round']

    # source indexes
    #src_idxs_r = round_df[selected_cols[0]].dropna().index
    src_idxs_r = r_data.index

    #src_idxs_nr = not_round_df[selected_cols[0]].dropna().index
    src_idxs_nr = nr_data.index

    # generate columns
    cluster_round_col = pd.DataFrame(clusters_round,
                                    columns=cls_r_colnm,
                                    index=src_idxs_r)

    cluster_not_round_col = pd.DataFrame(clusters_not_round,
                                    columns=cls_nr_colnm,
                                    index=src_idxs_nr)

    return cluster_round_col, cluster_not_round_col

def mount_frags_df(xdata_df, pdbid, chain,len_lim=3, label=True):
    '''
    Detect fragments which belongs to the same cluster and mount fragments
    dataframe.
    '''
    # --- LOCAL FUNCTION -------------------------------------------------------
    def __get_pdbx_frag_dict(pdbx_df,frag_len_arr, frag_idx_arr, frag_clstr_arr,
                    list_to_fill, len_lim, kind, pdbid, chain, dssp_pp2=True):
        '''
        Get fragments data, mount a dictionary and append it to a list
        '''
        def get_mean_std_median_of(df_col):
            return df_col.mean(), df_col.median(), df_col.std()
        #----------------------------------------------------------------------#
        # just to assure compatibitility with old code, but should be removed
        # on final version
        if label == True:
            class_col = 'label'
        if label == False:
            class_col = 'cluster'
        #----------------------------------------------------------------------#
        
        for n, i in enumerate(frag_idx_arr):
            # skip fragments bellow lengh limit
            size = frag_len_arr[n]
            if size < len_lim:
                continue
            # get indexes
            start_i = frag_idx_arr[n]
            final_i = start_i + size

            # get hlx dataframe
            hlx_df = pdbx_df.iloc[start_i:start_i+size]
            #get data to store
            res_s = hlx_df['res'].iloc[0]
            res_f = hlx_df['res'].iloc[-1]
            clstr = frag_clstr_arr[n]

            empty_dssp = False
            # handle empty dssp data
            if "ss" not in hlx_df.columns:
                empty_dssp = True
                print(f"WARNING : NO DSSP DATA FOR {pdbid}")

            if empty_dssp == False:
                dssp = hlx_df['ss'].values            
                seq = hlx_df['aa'].values
                if dssp_pp2 == True:
                    dssp_pp2_values = hlx_df['ss_pp2'].values
                # phi psi raw
                phi_raw = hlx_df['phi'].values
                psi_raw = hlx_df['psi'].values

                # h_nho1_en
                ho1en_col = hlx_df['h_nho1_en']
                ho1en_mean, ho1en_median, ho1en_std = get_mean_std_median_of(ho1en_col)
                ho1ens = ho1en_col.values

                # h_nho1_A
                ho1A_col = hlx_df['h_nho1_A']
                ho1A_mean, ho1A_median, ho1A_std = get_mean_std_median_of(ho1A_col)
                ho1As = ho1A_col.values

                # h_nho2_en
                ho2en_col = hlx_df['h_nho2_en']
                ho2en_mean, ho2en_median, ho2en_std = get_mean_std_median_of(ho2en_col)
                ho2ens = ho2en_col.values

                # h_nho2_A
                ho2A_col = hlx_df['h_nho2_A']
                ho2A_mean, ho2A_median, ho2A_std = get_mean_std_median_of(ho2A_col)
                ho2As = ho2A_col.values

                # h_ohn1_en
                hn1en_col = hlx_df['h_ohn1_en']
                hn1en_mean, hn1en_median, hn1en_std = get_mean_std_median_of(hn1en_col)
                hn1ens = hn1en_col.values

                # h_ohn1_A
                hn1A_col = hlx_df['h_ohn1_A']
                hn1A_mean, hn1A_median, hn1A_std = get_mean_std_median_of(hn1A_col)
                hn1As = hn1A_col.values

                # h_ohn2_en
                hn2en_col = hlx_df['h_ohn2_en']
                hn2en_mean, hn2en_median, hn2en_std = get_mean_std_median_of(hn2en_col)
                hn2ens = hn2en_col.values

                # h_ohn2_A
                hn2A_col = hlx_df['h_ohn2_A']
                hn2A_mean, hn2A_median, hn2A_std = get_mean_std_median_of(hn2A_col)
                hn2As = hn2A_col.values
            
            if empty_dssp == True:
                dssp = None
                seq = None
                if dssp_pp2 == True:
                    dssp_pp2_values = None
                # phi psi raw
                phi_raw = None
                psi_raw = None


                # h_nho1_en
                ho1en_col = None
                ho1en_mean, ho1en_median, ho1en_std = None, None, None
                ho1ens = None

                # h_nho1_A
                ho1A_col = None
                ho1A_mean, ho1A_median, ho1A_std = None, None, None
                ho1As = None

                # h_nho2_en
                ho2en_col = None
                ho2en_mean, ho2en_median, ho2en_std = None, None, None
                ho2ens = None

                # h_nho2_A
                ho2A_col = None
                ho2A_mean, ho2A_median, ho2A_std = None, None, None
                ho2As = None

                # h_ohn1_en
                hn1en_col = None
                hn1en_mean, hn1en_median, hn1en_std = None, None, None
                hn1ens = None

                # h_ohn1_A
                hn1A_col = None
                hn1A_mean, hn1A_median, hn1A_std = None, None, None
                hn1As = None

                # h_ohn2_en
                hn2en_col = None
                hn2en_mean, hn2en_median, hn2en_std = None, None, None
                hn2ens = None

                # h_ohn2_A
                hn2A_col = None
                hn2A_mean, hn2A_median, hn2A_std = None, None, None
                hn2As = None
            
            # Xgeo data
            # wri
            wri_col = hlx_df['wri']
            #try:
            w_mean, w_median, w_std = get_mean_std_median_of(wri_col)
            #except(TypeError):
            #    print(f"EITAOH: {wri_col.mean()}")
            #    print(stop)
            w_raw = wri_col.values

            # wri smoothed
            ws_col = hlx_df['wri_smooth']
            ws_mean, ws_median, ws_std = get_mean_std_median_of(ws_col)
            ws_raw = ws_col.values
            # tor
            tor_col = hlx_df['tor']
            t_mean, t_median, t_std = get_mean_std_median_of(tor_col)
            tors = hlx_df['tor'].values
            t_raw = tor_col.values

            # tor smoothed
            tor_col_smth = hlx_df['tor_norm_smooth']
            ts_mean, ts_median, ts_std = get_mean_std_median_of(tor_col_smth)
            ts_raw = tor_col_smth.values

            # curv
            cur_col = hlx_df['curv']
            c_mean, c_median, c_std = get_mean_std_median_of(cur_col)
            c_raw = cur_col.values

            # curv
            cur_col_s = hlx_df['curv_norm_smooth']
            cs_mean, cs_median, cs_std = get_mean_std_median_of(cur_col_s)
            cs_raw = cur_col_s.values



            hlx = {
                # metadata
                'res_start': res_s, 'res_final':res_f, 'kind': kind,
                'size':size, 'label':clstr, 'pdbid':pdbid, 'chain':chain,
                'seq':seq,
                # means and stds
                'c_mean':c_mean,'t_mean': t_mean, 'w_mean':w_mean,
                'c_std':c_std,'t_std':t_std, 'w_std':w_std,
                # raw data
                'c_raw':c_raw, 't_raw':t_raw, 'w_raw':w_raw,
                'cs_raw':cs_raw, 'ts_raw':ts_raw, 'ws_raw':ws_raw,
                # smooth means
                'cs_mean':cs_mean, 'ts_mean': ts_mean,'ws_mean':ws_mean,
                'c_std_s':cs_std,'t_std_s':ts_std, 'ws_std':ws_std,
                # dssp
                'phi':phi_raw, 'psi':psi_raw,
                'ss':dssp, 'ss_pp2':dssp_pp2_values,
                'ho1_en_mean': ho1en_mean, 'ho1_A_mean':ho1A_mean,
                'ho1_en_std': ho1en_std, 'ho1_A_std':ho1A_std,
                'ho1_en_median': ho1en_median, 'ho1_A_median':ho1A_median,
                'ho1_en_vals': ho1ens, 'ho1_A_vals':ho1As,

                'ho2_en_mean': ho2en_mean, 'ho2_A_mean':ho2A_mean,
                'ho2_en_std': ho2en_std, 'ho2_A_std':ho2A_std,
                'ho2_en_median': ho2en_median, 'ho1_A_median':ho2A_median,
                'ho2_en_vals': ho2ens, 'ho2_A_vals':ho2As,

                'hn1_en_mean': hn1en_mean, 'hn1_A_mean':hn1A_mean,
                'hn1_en_std': hn1en_std, 'hn1_A_std':hn1A_std,
                'hn1_en_median': hn1en_median, 'hn1_A_median':hn1A_median,
                'hn1_en_vals': hn1ens, 'hn1_A_vals':hn1As,

                'hn2_en_mean': hn2en_mean, 'hn2_A_mean':hn2A_mean,
                'hn2_en_std': hn2en_std, 'hn2_A_std':hn2A_std,
                'hn2_en_median': hn2en_median, 'hn2_A_median':hn2A_median,
                'hn2_en_vals': hn2ens, 'hn2_A_vals':hn2As,
                }

            list_to_fill.append(hlx)

    # --------------------------------------------------------------------------
    hlxs = []
    pdbx_df = xdata_df

    # sanity check
    if label == False:
        nr_colnm = 'clusters_not_round'
        r_colnm = 'clusters_round'
        assert(nr_colnm in pdbx_df.columns), nr_colnm+" column not found."
        assert(r_colnm in pdbx_df.columns), r_colnm+" column not found."

        nr_arr = pdbx_df['clusters_not_round'].values
        r_arr = pdbx_df['clusters_round'].values

        # detect fragments
        nr_frags_data = MathToolBox.get_repeats_at(nr_arr)
        r_frags_data = MathToolBox.get_repeats_at(r_arr)

        # get data only for fragments higher than the lenght limit
        __get_pdbx_frag_dict(pdbx_df=pdbx_df,
                            frag_len_arr= nr_frags_data[0],
                            frag_idx_arr= nr_frags_data[1],
                            frag_clstr_arr = nr_frags_data[2],
                            len_lim = len_lim, kind='not_round',
                            list_to_fill = hlxs, pdbid=pdbid, chain=chain)

        __get_pdbx_frag_dict(pdbx_df=pdbx_df,
                            frag_len_arr= r_frags_data[0],
                            frag_idx_arr= r_frags_data[1],
                            frag_clstr_arr = r_frags_data[2],
                            len_lim = len_lim, kind='round',
                            list_to_fill = hlxs, pdbid=pdbid, chain=chain)
    if label == True:
        lbl_arr = pdbx_df.label.values
        repeats_data = MathToolBox.get_repeats_at(lbl_arr)
        __get_pdbx_frag_dict(pdbx_df=pdbx_df,#.loc[res_df['label'] != 'p-1'],
                            frag_len_arr= repeats_data[0],
                            frag_idx_arr= repeats_data[1],
                            frag_clstr_arr = repeats_data[2],
                            len_lim = 3, kind='labels',
                            list_to_fill = hlxs, pdbid=pdbid, chain=chain)

    # mount the fragments dataframe
    frag_df = pd.DataFrame(hlxs)
    return frag_df

# ~~~~ Function for paralelization ~~~~ #
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

def detect_hlx(entry):
    entry.detect_hlx(show_plot=False, w_lim=0.08, r_dist=0.323, nr_dist=0.15,
        selected_cols=['curv_norm_smooth', 'tor_norm_smooth', 'wri_smooth'])
    return entry

# ~~~~ Functions for Calfa coordinate operations ~~~~ #
def trans2origin(coord):
    '''
    Compute coordinates translated so the first point is at the origin (0,0)
    '''
    dest=np.array([[0,0,0]])
    p = coord[0,:]
    Ts = dest - p
    return Ts + coord

def get_rotAngl_plane(coord):
    '''get rotation angles'''

    # get second point coordinates
    x = coord[1,0]
    y = coord[1,1]
    z = coord[1,2]

    # compute hypot
    h_xy = np.sqrt((x**2) + (y**2))
    h_zy = np.sqrt((z**2) + (y**2))

    # get rotation theta on XY plane
    cos_theta_xy = x / h_xy
    theta_xy = np.arccos(cos_theta_xy)

    # get rotation theta on ZY plane
    cos_theta_zy = z / h_zy
    theta_zy = np.arccos(cos_theta_zy)

    # phi = 90 - x
    phi_xy = (np.pi/2) - theta_xy
    phi_zy = (np.pi/2) - theta_zy

    return theta_xy,phi_xy, phi_zy, theta_zy

def rotate_on_xy(coord, theta_xy):
    '''rotate coordinates on the xy plane (i. e., around Z axis)'''

    # trigonometric functions
    cos_xy = np.cos(theta_xy)
    sin_xy = np.sin(theta_xy)

    # rotate coordinates
    x = coord[:,0]
    y = coord[:,1]
    z = coord[:,2]

    x_ = (x * cos_xy) - (y * sin_xy)
    y_ = (x * sin_xy) + (y * cos_xy)

    # return new coordinates
    coord_R_xy = np.zeros((len(coord),3))
    coord_R_xy[:,0] = x_
    coord_R_xy[:,1] = y_
    coord_R_xy[:,2] = z

    return coord_R_xy

def rotate_on_zy(coord, theta_zy):
    '''rotate coordinates on the zy plane (i. e. around X axis)'''
    # trigonometric functions
    cos_zy = np.cos(theta_zy)
    sin_zy = np.sin(theta_zy)

    # rotate coordinates
    x = coord[:,0]
    y = coord[:,1]
    z = coord[:,2]

    y_ = (y * cos_zy) - (z * sin_zy)
    z_ = (y * sin_zy) + (z * cos_zy)

    # return new coordinates
    coord_R_zy = np.zeros((len(coord),3))
    coord_R_zy[:,0] = x
    coord_R_zy[:,1] = y_
    coord_R_zy[:,2] = z_

    return coord_R_zy

def __get_single_entry_mtdt(lPDB_df, pdbid, chain):
    '''Get metadata to load a single entry from lPDB dataframe'''
    assert(type(pdbid) == str), 'pdbid must be str, not {}'.format(type(pdbid))
    assert(type(chain) == str), 'chain must be str, not {}'.format(type(chain))
    row = lPDB_df.loc[(lPDB_df['pdbid'] == pdbid) & (lPDB_df['chain']==chain)]
    try:
        assert(len(row)==1)
        return row.index[0]
    except(AssertionError):
        if len(row) > 1:
            print('More than one entry found on the metadata dataframe')
        if len(row) == 0:
            print('None entry found for {} chain {}'.format(pdbid, chain))

#-------------------------------------------------------------------------------
class entry:
    '''
    This class is designed to represent data for a unique chain PDB entry

    Inputs
    ---
    coord_flpath: <str>, filepath of the .pdb file
    xgeo_flpath: <str>, filepath of the .xgeo file
    pdbid: <str> or None, pdb identification [optional]
    chain: <str> or None, chain identification [optional]
    moltype: <str> or None, molecular type [optional]
    ----------------------------------------------------------------------------
                                Attributes
    ----------------------------------------------------------------------------
    coord_flpath = <str>,
        entry coordinate file path

    xgeo_flpath = <str>,
        xgeo file path

    pdbid = <str> or None, (default = None)
        entry identification label
    chain = <str> or None, (default = None)
        chain identification label
    moltype = <str>, (default = None)
        molecular type label
    cont = <bool>,
        True if backbone is continuous, False otherwise
    nres = <int>
        number of residues

    dssp_flpath = <str>
        DSSP output file path (default=None)

    frag_df = <pandas.DataFrame> or None,
        fragments dataframe

    xdata_df = <pandas.DataFrame>,
        dataframe for residues xgeo data

    ----------------------------------------------------------------------------
                                Methods
    ----------------------------------------------------------------------------

    ---| Load data methods |----------------------------------------------------
    self.load_dssp_data,
        load dssp data from dssp output file

    ---| Plot data methods |----------------------------------------------------
    self.plot_arrows,
        generate arrow plots of xgeo data

    self.plot_frag_df,
        generate plot for fragments (mean and std rep)

    ---| Data preparation methods |---------------------------------------------
    self.add_smoothed,
        generate a new column of smoothed values from another one

    self.get_normalize_xdf,
        generate normalized values for curvature and torsion values

    self.standardize_coords,
        translate coordinates so the first Calpha is at the origin (0,0,0) and
        second Calpha at ~(0,1,0)

    ---| Analyses methods |-----------------------------------------------------
    self.detected_hlx,
        run helix detection via wri filter and hard tresholds for hierarchical
        clustering

    self.get_dist2canonical,
        compute distance between entry residues and canonical clusters

    self.get_labels,
        assign labels according to residue distance to canonical groups
    ----------------------------------------------------------------------------
    '''

    def __init__(self, coord_flpath, xgeo_flpath=None, pdbid=None, chain=None,
        moltype=None, name=None, model_i=None):

        # 1 - Set attributes
        self.name = name
        self.coord_flpath = coord_flpath
        self.coord_flname = coord_flpath.split('/')[-1]
        self.pdbid = pdbid
        self.moltype = moltype
        self.chain = chain
        self.xgeo_flpath = xgeo_flpath
        # 2 - load xgeo dataframe
        if self.xgeo_flpath != None:
            self.xdata_df = load_xgeo_df(xgeo_flpath)
        if self.xgeo_flpath == None:
            self.computeDiffGeo(model_i=model_i)
        # 3 - get is count
        self.cont = self.__is_bkb_cont()

        # 4 - get number of residues
        self.nres = len(self.xdata_df['res'])

        # << set placeholders >>
        self.dssp_flpath = None
        self.frag_df = None

    def __is_bkb_cont(self):
        '''
        Check if backbone is continuous. This function just check if res column
        on xgeo data is continuous. If not, than backbone is not continuous.
        '''
        # get consecutive distance and check if there is any which is not one.
        d = np.diff(self.xdata_df['res'])
        if (len(np.unique(d)) == 1) and (np.unique(d)[0] == 1):
            return True
        else:
            return False

    # ~~~ COMPUTE DG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def computeDiffGeo(self, model_i=None, engine="diffgeo"):
        '''
        compute curvature, torsion and writhing number and return a pandas dataframe
        '''
        try:
            assert(engine in XGEO_ENGINES)
        except(AssertionError):
            print(f"ERROR: {engine} is not valid ({XGEO_ENGINES})")
            exit(1)

        if engine == "melodia":
            # compute geometry
            dfi = mel.geometry_from_structure_file(self.coord_flpath)

            # if model is provided, check if it is available
            if model_i != None:
                assert(model_i in dfi.model.unique())

            # if model note provided, choose first model
            if model_i == None:
                model_i = dfi.model.unique()[0]
        
            model = dfi['model'] == model_i
            dfo = dfi[model].copy()
        
            # store xdata as attribute
            # NEED TO BS SURE IT IS CONSISTENT WITH PREVIOUS DF FORMAT!! (check load_xgeo_df)
            dfo.rename(columns={'model':'conf_n', 'order':'res', 'curvature':'curv', 
                                'torsion':'tor', 'arc_length':'arc','writhing':'wri'}, inplace=True)
            self.xdata_df = dfo
        
        if engine == "diffgeo":

            in_path = self.coord_flpath.split('.')[0]
            PDBfileToolBox.run_diffgeo(in_path)
            xgeo_flpath = f"{in_path}.csv"
            assert(os.path.exists(xgeo_flpath))
            # fix columns
            dgo_cols = ["conf","res_name","atom", "res", "curv", "tor", "wri","arc"]
            self.xdata_df = pd.read_csv(xgeo_flpath, names=dgo_cols)
            self.xdata_df["wri"] = self.xdata_df["wri"].astype(float)
            self.xgeo_flpath = xgeo_flpath
            
    # ~~~ LOAD DATA METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def load_dssp_data(self, dssp_flpath, dssp_pp2 = True):
        '''
        Add dssp data to xdata dataframe

        Input:
        ::dssp_flpath: <str>, dssp file path
        ::dssp_pp2: <bool>, add pp2 assignment (default=True)
        '''
        # load dssp dataframe
        dssp_obj = dssp.DSSPData(dssp_flpath)
        dssp_data = dssp_obj.getDSSP_data()
        if len(dssp_data) > 0:
            # rename column for merge and set type as int
            dssp_data.rename(columns={'resnum':'res'}, inplace=True)

            def intfy(row):
                try:
                    return int(row)
                except(ValueError):
                    return None

            dssp_data['res'] = dssp_data['res'].apply(intfy)
            # merge new data =)
            self.xdata_df = pd.merge(self.xdata_df, dssp_data, on='res')
            if "phi_x" in self.xdata_df.columns:
                self.xdata_df.drop(
                    ["phi_x", "psi_x"], axis=True, inplace=True
                    )
                
                self.xdata_df.rename(
                        columns={"phi_y":"phi", "psi_y":"psi"},
                        inplace=True
                    )
            # get pp2 assignment
            if dssp_pp2==True:
                self.xdata_df['ss_pp2'] = self.xdata_df.apply(dssp_pp2_assgn,axis=1)

    # ~~~ PLOT METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def plot_arrows(self, myDIR=None, save_fig=False, show_plot=True,
                     plot_suffix='Arrows',plot_prefix='', highlight_res=None):
        ''' Generate arrow plot for entry xgeo data'''
        def __compute_dvecs(xdata_df):
            '''
            Compute displacement vectors of between pairs of coordinates n and n-1
            '''
            # get displacement vectors
            P_n_m1_arr = xdata_df[['curv', 'tor', 'wri']].iloc[:-1]#0:len(xdata_df)] #-1:2]
            P_n_arr = xdata_df[['curv', 'tor', 'wri']].iloc[1:]#:len(xdata_df)] #:2]
            dP_arr = P_n_arr.values - P_n_m1_arr.values

            #dP_arr = displ_df[['Px_n', 'Py_n', 'Pz_n']].values - P_arr

            # get arrows tail (aka, P(n-1) position)
            X = P_n_m1_arr.values[:,0] #curv
            Y = P_n_m1_arr.values[:,1] #tor
            Z = P_n_m1_arr.values[:,2] #curv

            # get arrows head directions (aka, towards P(n) positions)
            Ux = dP_arr[:,0]
            Uy = dP_arr[:,1]
            Uz = dP_arr[:,2]

            #store results
            dvec_dct = {'X':X,'Y':Y,'Z':Z,'Ux':Ux,'Uy':Uy,'Uz':Uz}
            dvecs_df = pd.DataFrame(dvec_dct)
            return dvecs_df

        def __gen_arrow_plot(X, Y, Ux, Uy, M, plot_flnm, xlabel, ylabel,
                            x_min=None, x_max=None, y_min=None, y_max=None):
            '''generate arrow plot'''
            fig, ax = plt.subplots()
            ax.quiver(X, Y, Ux, Uy, M, units='xy', scale=1, pivot='tail',
                     alpha=0.5, lw=0.1, cmap='coolwarm')
            plt.scatter(X,Y, s=1)
            plt.axis('equal')
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)

            if y_min != None and y_max != None:
                plt.ylim(y_min, y_max)
            if x_min != None and x_max != None:
                plt.xlim(x_min, x_max)
            if save_fig == True:
                plt.savefig(plot_flnm, dpi=300)
            if show_plot == True:
                plt.show()

            plt.close()

        # compute displacement vectors
        if highlight_res == None:
            dvecs_df = __compute_dvecs(self.xdata_df)
        if highlight_res != None:
            err_m = 'highlight_res must be a list of two elements'
            assert(len(highlight_res) == 2), err_m
            res_i = highlight_res[0]
            res_f = highlight_res[1]
            assert(type(res_i) == int), 'res_i must be int'
            assert(type(res_f) == int), 'res_i must be int'
            select = (self.xdata_df['res']>=res_i) & (self.xdata_df['res']<= res_f)
            dvecs_df = __compute_dvecs(self.xdata_df.loc[~select])
            dvecs_df_res = __compute_dvecs(self.xdata_df.loc[select])

        if save_fig is True:
            assert(myDIR != None)
        # get axis to plot
        X = dvecs_df.X.values #curv
        Y = dvecs_df.Y.values # tors
        Z = dvecs_df.Z.values # wri

        # get arros directions
        Ux = dvecs_df.Ux.values
        Uy = dvecs_df.Uy.values
        Uz = dvecs_df.Uz.values

        # get arrows magnitude
        M = np.sqrt(Ux**2 + Uy**2 + Uz**2)#np.hypot(Ux, Uy,Uz)

        # plot stuff
        # get output file paths and names
        plt_flnm1 = None
        plt_flnm2 = None
        plt_flnm3 = None

        if myDIR != None:
            plt_flnm1 = myDIR+'/'+plot_prefix+'wri_tor'+plot_suffix+'.png'
            plt_flnm2 = myDIR+'/'+plot_prefix+'curv_tor'+plot_suffix+'.png'
            plt_flnm3 = myDIR+'/'+plot_prefix+'wri_curv'+plot_suffix+'.png'

        __gen_arrow_plot(X=Z, Y=Y, Ux=Uz, Uy=Uy, M=M, xlabel='Writhing',
                      ylabel=r'$\tau$', plot_flnm=plt_flnm1)
        __gen_arrow_plot(X=X, Y=Y, Ux=Ux, Uy=Uy, M=M, xlabel=r'$\kappa$',
                      ylabel=r'$\tau$', plot_flnm=plt_flnm2)
        __gen_arrow_plot(X=Z, Y=X, Ux=Uz, Uy=Ux, M=M, xlabel='Writhing',
                       ylabel=r'$\kappa$', plot_flnm=plt_flnm3, x_min=-0.3, x_max=0.3 )

    def plot_frag_df(self, show_plot=True, suffix ='', scols=None, label_col='label'):
        '''Plot mean and standard deviation of detected fragments'''
        ## LOCAL FUNCTIONS #####################################################
        def __plot_err_bar(ax, src_df, x_col, y_col, xerr, yerr, lbl):
            ''' '''
            for clstr in src_df[label_col].unique():
                data = src_df.loc[(src_df[label_col] == clstr)]
                if len(data) == 0:
                    continue
                ax.errorbar(data[x_col], data[y_col], xerr=data[xerr],
                            yerr=data[yerr], fmt='.', lw=0.5,
                            label=str(clstr)+lbl)

        def __gen_plot(x_col, y_col, xerr, yerr, xlbl, ylbl):
            # plot cluster curv vs tor
            fig, ax = plt.subplots()
            frag_data = self.frag_df

            __plot_err_bar(ax, frag_data, x_col, y_col, xerr, yerr,suffix)

            plt.xlabel(xlbl)
            plt.ylabel(ylbl)
            ax.legend()
            # TODO add savefig using proper dir
            if show_plot == True:
                plt.show()
            plt.close()
        ########################################################################
        __gen_plot('c_mean', 't_mean', 'c_std', 't_std', r"$\kappa$", r"$\tau$")
        #__gen_plot('cs_mean', 'ts_mean', 'cs_std', 'ts_std', r"$\kappa_{sm}$", r"$\tau_{sm}$")

        __gen_plot('w_mean', 't_mean', 'w_std', 't_std',r"$W$", r"$\tau$")
        #__gen_plot('ws_mean', 'ts_mean', 'ws_std', 'ts_std',r"$W_{sm}$", r"$\tau_{sm}$")

        __gen_plot('w_mean', 'c_mean', 'w_std', 'c_std',r"$W$", r"$\kappa$")
        #__gen_plot('ws_mean', 'cs_mean', 'ws_std', 'cs_std',r"$W_{sm}$", r"$\kappa_{sm}$")

        __gen_plot('ho1_en_mean','w_mean','ho1_en_std','w_std',"ho1_en",
                    "Writhing")
        __gen_plot('hn1_en_mean','w_mean','hn1_en_std','w_std',"hn1_en",
                    "Writhing")
        __gen_plot('hn1_en_median','w_mean','hn1_en_std','w_std',
                    "hn1_en_median","Writhing")

        if 'res_pturn_mean' in self.frag_df.columns:
            __gen_plot('res_pturn_mean', 'c_mean', 'res_pturn_std', 'c_std',
                r'residues per turn', r"$\kappa$" )
            __gen_plot('res_pturn_mean', 't_mean', 'res_pturn_std', 't_std',
                r'residues per turn', r"$\tau$" )

    # ~~~ PRE PROCESSING METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def add_smoothed(self, column_name, method='SGfilter', window_size=3):
        '''
        Smoothed a given column of xdata and add as a new columns to it
        Parameters
        ----------
        column_name : string
            the name of the column of xdata dataframe to be smoothed
        method : string
            the keyword for the smoothing function. Current version supports:
                * Moving Average ('MAvg')
                * Savitzky-Folay filter ('SGfilter')
        window_size : int, non negative
        '''
        # 1 - Input Processing
        try:
            assert(column_name in self.xdata_df.columns)
            data = self.xdata_df[column_name].values.astype(float)

        except(AssertionError):
            raise ValueError("{} is not a columns of xdata".format(column_name))

        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")

        VALID_METHODS = ['MAvg', 'SGfilter']
        if method not in VALID_METHODS:
            raise ValueError("{} is not a valid method keyword.".format(method))

        # 2 - Smoothing
        if method == 'MAvg':
            #smoothed = Analyses.MathToolBox.running_mean(data,window_size)
            smoothed = MathToolBox.running_avg(data, window_size)
        if method == 'SGfilter':
            order = 1
            smoothed = MathToolBox.savitzky_golay(data, window_size,order)

        new_column = pd.DataFrame(smoothed, columns=[column_name+"_smooth"])
        self.xdata_df = pd.concat([self.xdata_df, new_column], axis=1)

    def get_normalize_xdf(self, c_norm=False, **kwargs):
        '''
        Normalize curvature and torsion columns of xgeo data frame.

        The default behaviour is use min and max values of the xdata columns for
        for the normalization procedure.
        Custom values of minimum and maximum values can be set by 'c_norm=True'.
        User must provide 'k_min', 'k_max', 't_min' and 't_max'. This is usefull
        when working with a group of pdbs.

        PARAMETERS
        ==========
        c_norm : bool, DEFAULT=False
                defines if custom values of min and max will be used on the
                normalization.

        **kwargs
        Only used for custom normalization
        k_min : float
                defines curvature minimum value
        k_max : float
                defines curvature maximum value
        t_min : float
                defines torsion minimum values
        t_max : float
                defines curvature minimum values
        '''
        # 0 - Check if there is no normalized data
        #     *avoids the inclusion of a new column on xdata
        assert('curv_norm' not in self.xdata_df.columns)
        assert('tor_norm' not in self.xdata_df.columns)

        # 1 - get columns data
        X_curv = self.xdata_df['curv'].loc[self.xdata_df['curv'] != 0.0]
        X_tor = self.xdata_df['tor'].loc[self.xdata_df['tor'] != 0.0]

        # 2 - get min and max values for normalization
        if c_norm is False:
            k_min = X_curv.min()
            k_max = X_curv.max()

            t_min = X_tor.min()
            t_max = X_tor.max()

        if c_norm is True:
            # process kwargs
            try:
             k_min = kwargs['k_min']
             k_max = kwargs['k_max']

             t_min = kwargs['t_min']
             t_max = kwargs['t_max']

            except(KeyError):
             print("!!ERROR: You must provide the minimum and maximum values")
             print("         to use the custom normalization. The keys to")
             print("         use are 't_min', 't_max', 'k_min', 'k_max'.")
             print("  kwargs provided: {}".format(kwargs))
             exit(1)


        # 3 - Do the normalization
        X_curv_norm = custom_norm(X_curv, min=k_min, max=k_max)
        X_tor_norm = custom_norm(X_tor, min=t_min, max=t_max, a=0, b=1.0)

        # 4 - stick the normalized columns to xdata
        df = pd.concat([X_curv_norm, X_tor_norm], axis=1)
        df.columns = ['curv_norm', 'tor_norm']
        self.xdata_df = pd.concat([self.xdata_df, df], axis=1)

    def standardize_coords(self):
        '''
        Reorient coordinates of Calfa in order that:
         1) the first C alpha sits at the origin
         2) the second C alpha is aligned on the y axis
        '''
        # load coordinates
        coord = self.xdata_df[['ca_x', 'ca_y', 'ca_z']].values

        # translate 1st C alpha to origin
        coord_T = trans2origin(coord)

        # get possible rotational angles
        theta_xy,phi_xy, phi_zy, theta_zy = get_rotAngl_plane(coord_T)

        # rotate on xy
        if coord_T[1,1] >= 0:
            coord_Rxy = rotate_on_xy(coord_T, phi_xy)
        if coord_T[1,1] < 0:
            coord_Rxy = rotate_on_xy(coord_T, theta_xy+(np.pi/2))

        #rotate on zy
        if coord_Rxy[1,2] < 0:
            coord_Rzy = rotate_on_zy(coord_Rxy, phi_zy)
        if coord_Rxy[1,2] >= 0:
            coord_Rzy = rotate_on_zy(coord_Rxy, phi_zy*-1)

        try:
            assert(coord_Rzy[1,0] < 0.1), 'xy rotation error'
            assert(coord_Rzy[1,2] < 0.1), 'zy rotation error'
        except(AssertionError):
            print(self.pdbid, self.chain, ' rotation error!')

        self.xdata_df.loc[:,['ca_x', 'ca_y', 'ca_z']] = coord_Rzy

    # ~~~ ANALYSES METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def detect_hlx(self, show_plot=False, w_lim=0.1, r_dist=0.3, nr_dist=0.85,
                    selected_cols=['curv', 'tor', 'wri'], label=True):
        '''
        Detect helix using writhin number filter, tresholds for 'round' and
        'non-round' and Hierarchical clustering.

        Parameters
        ---
        w_lim = <float>, writhing number trshold module to define the 'round' and
                'non-round'. (default=0.1)
        r_dist = <float>, threshold to cut the 'round' dendogram, which defines
                final cluster solution . (default=0.3)
        nr_dist = <float>, threshold to cut the 'non-round' dendogram, which
                defines cluster solution
        selected_cols = <list>, list of column names to run the algorithm.
                (default=['curv','tor','wri'])
        label = <bool>, Assert if there is already labels assigned, if TRUE
                the fragments dataframe will be generated based on such labels.
                If False, the labels will be assigned.
                (default=True)
        '''

        xdata_df = self.xdata_df
        if label == False:

            c_round_col, c_nround_col = cluster_res_by_wfilter(xdata_df,
                wri_lim=w_lim, min_rsize=3, min_nrsize=3,r_dist=r_dist,
                nr_dist=nr_dist, show_plot=show_plot,
                selected_cols=selected_cols)
            # stick new columns
            xdata_df = pd.concat([xdata_df, c_round_col], axis=1)
            xdata_df = pd.concat([xdata_df, c_nround_col],axis=1)
            xdata_df["wri"] = xdata_df["wri"].astype(float)
            if show_plot is True:
                # plot cluster curv vs tor
                def __plot_XY_clstrs(colx, coly, labelx, labely, r_min=5,
                                    nr_min=5):
                    '''generate scatter plot of clustered residues'''
                    fig, ax = plt.subplots()
                    def __add_scatter(clst_colnm, min_datapts, sfx):
                        err_msg = 'No '+clst_colnm+' column found'
                        assert(clst_colnm in xdata_df.columns), err_msg
                        for clstr in xdata_df[clst_colnm].dropna().unique():
                            data = xdata_df.loc[(xdata_df[clst_colnm] == clstr)]
                            if len(data) <= r_min:
                                continue
                            ax.scatter(data[colx],data[coly],
                                        label=str(clstr)+sfx,
                                        s=.9, alpha=0.8)

                    __add_scatter('clusters_round', r_min, '_r')
                    __add_scatter('clusters_not_round', nr_min, '_nr')

                    plt.xlabel(labelx)
                    plt.ylabel(labely)
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles, labels, loc='upper right',
                              bbox_to_anchor=(1.3,1.0))
                    plt.show()
                    plt.close()

                __plot_XY_clstrs('curv','tor',r'$\kappa$',r'$\tau$',r_min=5,nr_min=5)
                __plot_XY_clstrs('wri','tor',r'$w$',r'$\tau$',r_min=5,nr_min=5)
                __plot_XY_clstrs('wri','curv',r'$w$',r'$\kappa$',r_min=5,nr_min=5)

        if label == True:
            err= "No 'label' column on dataframe."
            assert('label' in self.xdata_df.columns),err
        #try:
        self.frag_df = mount_frags_df(xdata_df, len_lim=3, pdbid=self.pdbid,
                                          chain=self.chain, label=label)
        #except(TypeError):
        #    print(xdata_df["wri"].values)
        #    print(self.xgeo_flpath)
        self.xdata_df = xdata_df

    def get_dist2canonical(self, pi_df, alpha_df, three_df, pp2_df, lalpha_df=None, lthree_df=None):
        ''' '''
        # for pandas apply
        def get_alphad(row):
            '''Compute distance for alpha set'''
            p = row[['curv', 'tor', 'wri']].values
            C = alpha_df[['c_mean', 't_mean', 'w_mean']].values
            d_pA_min = MathToolBox.get_d2_pC_min(p, C)
            return d_pA_min

        def get_pid(row):
            '''Compute distance for alpha set'''
            p = row[['curv', 'tor', 'wri']].values
            C = pi_df[['c_mean', 't_mean', 'w_mean']].values
            d_pPi_min = MathToolBox.get_d2_pC_min(p, C)
            return d_pPi_min

        def get_310d(row):
            '''Compute distance for 3(10) set'''
            p = row[['curv', 'tor', 'wri']].values
            C = three_df[['c_mean', 't_mean', 'w_mean']].values
            d_p310_min = MathToolBox.get_d2_pC_min(p, C)
            return d_p310_min

        def get_pp2d(row):
            '''Compute distance for pp2 set'''
            p = row[['curv', 'tor', 'wri']].values
            C = pp2_df[['c_mean', 't_mean', 'w_mean']].values
            d_pPP2_min = MathToolBox.get_d2_pC_min(p, C)
            return d_pPP2_min

        def get_lalphad(row):
            '''Compute distance for lalpha set'''
            p = row[['curv', 'tor', 'wri']].values
            C = lalpha_df[['c_mean', 't_mean', 'w_mean']].values
            d_plA_min = MathToolBox.get_d2_pC_min(p, C)
            return d_plA_min
        
        def get_lthree(row):
            '''Compute distance for lthree set'''
            p = row[['curv', 'tor', 'wri']].values
            C = lthree_df[['c_mean', 't_mean', 'w_mean']].values
            d_pl310_min = MathToolBox.get_d2_pC_min(p, C)
            return d_pl310_min

        # ----------------------- #
        self.xdata_df['D(Alpha)'] = self.xdata_df.apply(get_alphad, axis=1)
        self.xdata_df['D(Pi)'] = self.xdata_df.apply(get_pid, axis=1)
        self.xdata_df['D(3(10))'] = self.xdata_df.apply(get_310d, axis=1)
        self.xdata_df['D(PP2)'] = self.xdata_df.apply(get_pp2d, axis=1)

        # add left alpha and 3(10) distances if provided
        if lalpha_df is not None:
            self.xdata_df['D(lalpha)'] = self.xdata_df.apply(get_lalphad, axis=1)
        if lthree_df is not None:
            self.xdata_df['D(l3(10))'] = self.xdata_df.apply(get_lthree, axis=1)


    def get_labels(self, dist_min=0.2, pp2_max = 0.07, left_hand=True):
        '''
        Get labels for residues based on distances to canonical groups.
        For alfa, pi and and 3(10), consider only residues with d<0.2, and
        choose the label based on the smallest distance among those three groups
        The PP2 labels need a less permisive max distance in order to avoid
        assign beta-strands as PP2 (because, helix continuum).
        This function adds a label column to entry dataframe.

        PARAMETERS
        ---
        dist_min = <float>, minimum distance for Alfa, pi and 3(10)
        pp2_max = <float>, max distance away from PP2 to still be assigned as such
        '''
        # get distances
        dist_arr = self.xdata_df[['D(Alpha)', 'D(Pi)', 'D(3(10))', 'D(PP2)']].values
        labels_arr = ['Alpha', 'Pi', '3(10)', 'PP2']
        
        if left_hand == True:
            # get left hand distances
            dist_arr = self.xdata_df[['D(Alpha)', 'D(Pi)', 'D(3(10))', 
                                      'D(PP2)','D(lalpha)', 'D(l3(10))']].values
            labels_arr = ['Alpha', 'Pi', '3(10)', 'PP2', 'lAlpha', 'l3(10)']

        # check distances bellow threshold
        row_is, col_is = np.where(dist_arr < dist_min)
        valid_is = np.unique(row_is)
        label_lst = []

        # check if there are any valid distances
        for i in range(0, len(dist_arr)):
            if i not in valid_is:
                label_lst.append('1-p')
                continue
            # if valid, check which one is the closest
            if i in valid_is:
                dist_min = dist_arr[i].min()
                l_min_i = np.where(dist_arr[i] == dist_min)[0][0]
                if labels_arr[l_min_i] == 'PP2' and dist_min >= pp2_max:
                    label_lst.append('1-p')
                    continue
                else:
                    label_lst.append(labels_arr[l_min_i])
                    continue

        # add new column
        label_col = pd.DataFrame(label_lst, columns=['label'],
                                 index=self.xdata_df.index)
        self.xdata_df = pd.concat([self.xdata_df, label_col], axis=1)

def load_single_entry(pdbid,chain,lPDB_df):
    '''load an entry object from lPDB dataframe'''
    q_i = __get_single_entry_mtdt(lPDB_df, pdbid, chain)
    entry_ = entry(coord_flpath=lPDB_df['chain_flspth'].iloc[q_i],
               xgeo_flpath=lPDB_df['xgeo_flspth'].iloc[q_i],
               pdbid=lPDB_df['pdbid'].iloc[q_i],
               chain=lPDB_df['chain'].iloc[q_i])
    return entry_

class group:
    '''
    This class is designed to represent arbitrary groups of entries and
    provide handy methods to create and analyse them.

    INPUT
    ---
    name=<str>,
        group name
    dir_path=<str>,
        path to store output files (default= working dir)[OPTIONAL]

    ----------------------------------------------------------------------------
                                Attributes
    ----------------------------------------------------------------------------

    entries = <list>,
        list to store entries which belong to this group. (default = [])

    grp_df = <pandas.DataFrame> or <None>,
        A single dataframe containing all xdata_df of the entries belonging to
        the group

    grp_frag_df = <pandas.DataFrame> or <None>,
        A single dataframe with all the frag_df of the entries

    ----------------------------------------------------------------------------
                                Methods
    ----------------------------------------------------------------------------

    ---| Load data methods |----------------------------------------------------
    load_from_lpdbmtdt_df,
        load group of entries based on local PDB dataframe obtained by lPDBmngr,
        (check https://github.com/AMarinhoSN/lPDBmngr)

    load_grp_dfs,
        generate a grp_df based on all entries stored at self.entries.

    ---| Plot data methods |----------------------------------------------------
    plot_dssp_data,
        generate KTW plot according to dssp assigned
    '''

    def __init__(self, name, dir_path=None):
        self.name = name
        self.entries = []
        # get directory to store ouput files, if not set as the working dir
        if dir_path == None:
            self.dir_path = os.getcwd()+'/'
        if dir_path != None:
            assert(os.path.isabs(dir_path)), 'dir_path needs to be a vaid path'
            self.dir_path = dir_path
        # << PLACEHOLDERS >>
        self.grp_df = None
        self.grp_frag_df = None

    # ~~~ LOAD DATA METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def load_from_lpdbmtdt_df(self,lPDB_df):
        '''Load group metadata from lPDB dataframe generated by LPDBmgr '''
        # ~~~ define local functions for apply ~~~
        def load_pdbx(row):
            pdbx = entry(xgeo_flpath=row['xgeo_flspth'],
                            coord_flpath=row['chain_flspth'],
                            pdbid=row['pdbid'], chain=row['chain'],
                            moltype=row['mol_type'])
            pdbx.load_dssp_data(row['dssp_flspth'])
            return pdbx

        def get_unit_gdfs(row):
            '''mount unit df for the group dataframe.'''
            ug_df = row.xdata_df[['res', 'curv', 'tor', 'wri', 'ss']]
            nrows = len(ug_df)
            ug_df.loc[:,'pdbid'] = [row.pdbid for i in range(nrows)]
            ug_df.loc[:,'is_cont'] = [row.cont for i in range(nrows)]
            ug_df.loc[:,'chain'] = [row.chain for i in range(nrows)]
            return ug_df
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # load pdbx objects
        pdbxs_lst = lPDB_df.apply(load_pdbx, axis=1)
        # mount unit df for group data
        ugdfs = pdbxs_lst.apply(get_unit_gdfs)
        # concat unit dataframes
        self.grp_xdata_df = pd.concat(ugdfs.values, ignore_index=True)

    def load_grp_dfs(self, frag_df=True, reload=False):
        '''
        load a single dataframe with all entries raw data.
        Parameter
        ---
        frag_df: <bool>, (default = True), mount group fragments dataframe
        reload: <bool>, (default = False), reload group regardless of current
            state
        '''
        # sanity check
        assert(len(self.entries) != 0), 'No entries available'
        if reload == False:
            assert(type(self.grp_df) == type(None)), 'group dataframe already loaded'
        if frag_df == True:
            f_msg = 'fragments group dataframe already loaded'
            assert(self.grp_frag_df == None), f_msg

        # load dataframe
        self.grp_df = load_group_df(self.entries, frag_df=False)
        if frag_df == True:
            self.grp_frag_df = load_group_df(self.entries, frag_df=True)

    # ~~~ PLOT METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def plot_dssp_data(self,only_cont=True,ss_to_plot='All',s_grey=.01,s=.01,
                       mrkr='.', alpha=0.5):
        '''Plot dssp data from group data'''
        # ~~~ Define local functions ~~~
        def __add_plot_colxy(colx,coly):
            # add grey dots
            #plt.scatter(grp_df[colx],grp_df[coly],s=.01, c='grey',alpha=0.8,
            #            marker='o',label='')
            # add dssp states
            for ss in ss_labels:
                df_slice = grp_df.loc[grp_df['ss']==ss]
                if ss == '':
                    plt.scatter(df_slice[colx], df_slice[coly], s=s_grey, label=ss,
                                c='grey', marker=mrkr, alpha=alpha)
                    continue

                if ss not in only_ss:
                    continue

                plt.scatter(df_slice[colx],df_slice[coly],s=s,label=ss,
                            marker=mrkr, alpha=alpha)
            plt.legend()
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Uncontinuous chain creates some unrealistic description on the residues
        # flanking the gap. Usually you want to avoid this or at list remove the
        # residues affected by it.
        print(' -> filtering for continuous fragments')
        grp_df = self.grp_xdata_df.loc[self.grp_xdata_df['is_cont']==only_cont]
        print('    :: ', len(grp_df),' found.')

        dir_path = self.dir_path+'plot/'
        # get labels to be ploted
        ss_labels = list(np.unique(grp_df['ss']))

        if ss_to_plot == 'All':
            only_ss = ss_labels
        else:
            assert(type(ss_to_plot)==list), 'ss_to_plot is not a list.'
            only_ss = ss_to_plot

        # plot wri and tor
        __add_plot_colxy('wri', 'tor')
        plt.xlabel(r'$w$')
        plt.ylabel(r'$\tau [\AA^{-1}]$')
        plt.savefig(dir_path+'WT_dssp.png', dpi=300)
        plt.close()
        # plot wri and curv
        __add_plot_colxy('wri', 'curv')
        plt.xlabel(r'$w$')
        plt.ylabel(r'$\kappa [\AA^{-1}]$')
        plt.savefig(dir_path+'WK_dssp.png', dpi=300)
        plt.close()

        # plot curv and tor
        __add_plot_colxy('curv','tor')
        plt.xlabel(r'$\kappa [\AA^{-1}]$')
        plt.ylabel(r'$\tau [\AA^{-1}]$')
        plt.savefig(dir_path+'KT_dssp.png', dpi=300)
        #plt.show()
        plt.close()
