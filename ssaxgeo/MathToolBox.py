import numpy as np
from math import factorial
from scipy.stats import zscore
import hdbscan
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance
from sklearn.cluster import KMeans


'''
Here handy mathematical functions are defined here
'''
#-----| Smoothing functions |---------------------------------------------------
# Basically, functions to remove spikes. Can be handy noise data analyses.
def running_mean(x, win_size):
    r'''
    Calculate the average mean for a window N
    this code was obtained from:
    https://stackoverflow.com/questions/13728392/moving-average-or-running-mean

    Parameters
    ----------
     x : array_like, shape(N,)
        the 1D function to be smoothed
    win_size : int
        the sliding window size for the averaging procedure

    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed data
    '''
    x = np.nan_to_num(x,copy=False)
    if win_size % 2 != 1 or win_size < 1:
        raise TypeError("win_size size must be a positive odd number")
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[win_size:] - cumsum[:-win_size]) / float(win_size)

def running_avg(x, N):
    '''running average with edge treatment, but less efficient'''
    return np.convolve(x, np.ones((N,))/N, mode='valid')

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
	r"""
    Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
	The Savitzky-Golay filter removes high frequency noise from data.
	It has the advantage of preserving the original shape and
	features of the signal better than other types of filtering
	approaches, such as moving averages techniques.
	Parameters
	----------
	y : array_like, shape (N,)
		the values of the time history of the signal.
	window_size : int
		the length of the window. Must be an odd integer number.
	order : int
		the order of the polynomial used in the filtering.
		Must be less then `window_size` - 1.
	deriv: int
		the order of the derivative to compute (default = 0 means only smoothing)
	Returns
	-------
	ys : ndarray, shape (N)
		the smoothed signal (or it's n-th derivative).
	"""
    # 1 - Input processing
	try:
		window_size = np.abs(int(window_size))
		order = np.abs(int(order))
	except(ValueError):#, msg):
		raise ValueError("window_size and order have to be of type int")
	if window_size % 2 != 1 or window_size < 1:
		raise TypeError("window_size size must be a positive odd number")
	if window_size < order + 2:
		raise TypeError("window_size is too small for the polynomials order")
	order_range = range(order+1)
	half_window = (window_size -1) // 2

    # precompute coefficients
	b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
	m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)

    # pad the signal at the extremes with
	# values taken from the signal itself
	firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
	lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
	y = np.concatenate((firstvals, y, lastvals))
	return np.convolve( m[::-1], y, mode='valid')

#-----| Stats functions |-------------------------------------------------------
def is_normal_dist(x, sig_max=3):
    r'''
    Test if there are values with unlikely Z-score for a normal dist. This test
    is based on the 68-95-99.7 rules.
    Check https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule for more
    information.
    ----------
    Parameters
    ----------
    x : numpy array,
    '''
    # --| RATIONALE |--
    # Values in the range of 6 standard deviations units account for 99.73% of
    # the data IF is a normal distribution. If one or more values with unlikely
    # Z-scores, it is evidence that you are not dealing with normal
    # distribution.
    # -----------------

    # 1 - Calculate z-score for the input distribution
    Zs = zscore(x)

    # 2 - Check for unlikely Z-scores
    if Zs.min() < -1*sig_max or Zs.max() > sig_max:
        return False
    else:
        return True

# ---| Clustering |-------------------------------------------------------------
# -- Hierarchical clustering --
def get_linkage_mtx_of(df, method='ward'):
    '''Get linkage of a pandas dataframe'''
    # calculate distance matrix
    dmtx = scipy.spatial.distance.pdist(df)
    # generate the linkage matrix
    return linkage(dmtx, method)

def fancy_dendrogram(*args, **kwargs):
    '''
    plot dendogram with some nice features
    check:
    https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
    '''
    # set max distance
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

# ---| HANDY NUMPY OPERATIONS |-------------------------------------------------

def get_repeats_at(inarray):
    '''
    Get and count lenght of repetitions on an array
    '''
    ia = np.asarray(inarray)                  # force numpy
    n = len(ia)
    if n == 0:
        return (None, None, None)
    else:
        y = np.array(ia[1:] != ia[:-1])     # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)   # must include last element posi
        z = np.diff(np.append(-1, i))       # run lengths
        p = np.cumsum(np.append(0, z))[:-1] # positions
        return(z, p, ia[i])

def generate_grid_cols(nx, ny, start_x, stop_x, start_y, stop_y,
                       show_plot=False):
    '''
    Generate a 2d grid with nx points for axis x, starting from start_x to
    stop_x, and ny points for axis y starting from start_y to stop_y.
    PARAMETERS
    ==========
    nx: int
        number of points on axis x column of the grid
    ny: int
        number of points on axis x column of the grid
    start_x: float
        number to start grid column x (where the x axis will start counting)
    stop_x: float
        number to stop grid column x (where the x axis will stop counting)
    start_y: float
        number to start grid column y (where the y axis will start counting)
    stop_y: float
        number to stop grid column y (where the y axis will stop counting)

    RETURN
    ======
    Gcol_x = np.array shape (1,m)
        column x of the grid array
    Gcol_y = np.array shape (n,1)
        column y of the grid array
    '''
    # make grid
    x = np.linspace(start=start_x, stop=stop_x, num=nx)
    y = np.linspace(start=start_y, stop=stop_y, num=ny)
    # get x and y axis of grid
    # NOTE: using sparse as TRUE, the meshgrid return the columns withou
    #       generating the matrix array. This can save some time and memory
    #       if computations can be done indexing wise. If not necessary
    #       turn sparse to False and you will get the squared matrix.
    Gcol_x,Gcol_y = np.meshgrid(x,y, sparse=True)

    # generate plot
    if show_plot == True:
        xv_plot, yv_plot = np.meshgrid(x,y)
        plt.scatter(xv_plot, yv_plot, s=0.1)
        plt.show()

    return Gcol_x, Gcol_y

def get_vecs_subgrids_dict(P, dP,Gcol_x,Gcol_y ,show_plot=False):
    '''
    Determinate the subgrid a given set of vector (P, P+dP) is affecting.
    For each vector, the interval on both dimensions considered determinate
    the subgrid where the vactor is acting on.
    PARAMETERS
    ==========
    P: np.array shape(N,2)
        initial position of 2D vectors
    dP: np.array shape(N,2)
        displacement vector
    Gcol_x = np.array shape (1,m)
        column x of the grid array
    Gcol_y = np.array shape (n,1)
        column y of the grid array
    show_plot: boolean, DEFAULT=False
        generate plot with vectors of defined by P+dP
    '''
    # get final position after vectors displacement
    P2 = P+dP

    # get vector subgrid indexes
    vecs_grids_dcts = []
    for n, vec_n in enumerate(P):
        # define first and last point to consider on the line
        # this step avoids problem with negative values on P2
        if dP[n,0] >= 0:
            Pf_x = P2[n,0]
            Pi_x = P[n,0]

        if dP[n,0] < 0:
            Pf_x = P[n,0]
            Pi_x = P2[n,0]

        if dP[n,1] >= 0:
            Pf_y = P2[n,1]
            Pi_y = P[n,1]

        if dP[n,1] <= 0:
            Pf_y = P[n,1]
            Pi_y = P2[n,1]

        # get vald x_i indexes
        i = np.nonzero((Gcol_x >= Pi_x) & (Gcol_x <= Pf_x))[1]
        j = np.nonzero((Gcol_y >= Pi_y) & (Gcol_y <= Pf_y))[0]

        dict = {'is':i, 'js':j, 'vec_n':n}
        vecs_grids_dcts.append(dict)

        if show_plot is True:
            # get vector subgrid
            subgrid_xs = Gcol_x[(Gcol_x >= Pi_x) & (Gcol_x <= Pf_x)]
            subgrid_ys = Gcol_y[(Gcol_y >= Pi_y) & (Gcol_y <= Pf_y)]
            # get mesh to plot
            sub_xv_plot, sub_yv_plot = np.meshgrid(subgrid_xs, subgrid_ys)
            # plot
            plt.scatter(sub_xv_plot, sub_yv_plot, s=2.0)
            plt.quiver(P[:,0],P[:,1], dP[:,0], dP[:,1], units='xy', scale=1)# head_width=0.01, head_length=0.01, lw=0.1)
            plt.axis('equal')
    return vecs_grids_dcts

def get_unq_idxs_arr(arr_):
    '''
    Get unique items indexes interval.

    WARNING: This function assumes unique itens appears in sequential order. So
    the returned indexes can be used to access entries directly and avoid
    DataFrame search on a for loop.

    PARAMETERS
    ==========
    arr_: 1D np.array
          Input array.

    RETURN
    ======
    unq_entries: 1D np.array
            Array of unique entries found on input numpy array
    unq_idxs: 1D np.array
            Index interval os unique entries of input numpy array
    '''
    unq_entries, unq_idxs = np.unique(np.vstack(arr_), axis=0, return_index=True)
    unq_idxs.sort()
    return unq_entries, unq_idxs

def search_sequence_numpy(arr,seq):
    """ Find sequence in an array using NumPy only.

    Parameters
    ----------
    arr    : input 1D array
    seq    : input 1D array

    Output
    ------
    Output : 1D Array of indices in the input array that satisfy the
    matching of input sequence in the input array.
    In case of no match, empty list is returned.

    Check:
    https://stackoverflow.com/questions/36522220/searching-a-sequence-in-a-numpy-array
    """

    # Store sizes of input array and sequence
    Na, Nseq = arr.size, seq.size

    # Range of sequence
    r_seq = np.arange(Nseq)

    # Create 2D array of sliding indices across entire length of input array.
    # Match up with the input sequence & get the matching starting indices.
    M = (arr[np.arange(Na-Nseq+1)[:,None] + r_seq] == seq).all(1)

    return M
    ### BACKUP CODE
    # Get the range of those indices as final output
    #if M.any>0:
    #    return np.where(np.convolve(M,np.ones((Nseq),dtype=int))>0)[0]
    #else:
    #    return []         # No match found

def d2_p_C(p, C):
    '''compute distance between a point, p, and a set, C'''
    d_pC_lst = []
    for i in range(0, len(C)):
        d_pC = np.linalg.norm(p-C[i])
        d_pC_lst.append(d_pC)

    return np.array(d_pC_lst)

def get_d2_pC_min(p, C):
    '''get dmin(p,C)'''
    d_pC_arr = d2_p_C(p, C)
    return np.array(d_pC_arr).min()


def run_kmeans_at(arr_, n_clusters, random_state=0, n_jobs=1, **kwargs):
    '''run kmeans clustering at ND numpy array.'''
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state).fit(arr_)

    return kmeans.labels_

#---| PANDAS HANDY FUNCTIONS |--------------------------------------------------

def get_dist_mtx_of(A_arr, B_arr):
    '''Calculate all pairs euclidean distance of two arrays'''
    m, n = np.meshgrid(A_arr, B_arr)
    distMTX = abs(m-n)
    return distMTX

def do_hdbscan_at(df, min_cluster_size, do_soft=True, use_distMTX=False,
                    allow_single_cluster=True,**kwargs):
    '''Do HDBScan using pandas dataframe'''
    # TODO add sanity check
    if use_distMTX is False:
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size,
                                    prediction_data=do_soft,
                                    allow_single_cluster= allow_single_cluster,
                                    **kwargs)
    if use_distMTX is True:
        clusterer = hdbscan.HDBSCAN(metric='precomputed',
                                    min_cluster_size=min_cluster_size,
                                    allow_single_cluster=allow_single_cluster,
                                    **kwargs)
    return clusterer.fit(df)

def plot_hdbscan_clusters(clusterer, df_x, df_y, s=50, alpha=0.25):
    ''' Plot clustering results from HDBSCAN'''
    # get colors to cycle
    n_clusters = len(np.unique(clusterer.labels_))
    color_palette = sns.color_palette('deep', n_clusters)
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in clusterer.labels_]
    cluster_member_colors = [sns.desaturate(x, p) for x, p in
                             zip(cluster_colors, clusterer.probabilities_)]

    plt.scatter(df_x.values, df_y.values, s=s, linewidth=0,
                c=cluster_member_colors, alpha=alpha)

def plot_2d_df_cols(df, col_x, col_y, **kwargs):
    '''generate a 2D scatter plot using two columns of a pandas dataframe'''
    # plot stuff
    s = kwargs.get('s', .01)
    alpha = kwargs.get('alpha', 0.8)
    save_fig = kwargs.get('save_fig', False)
    fl_name = kwargs.get('fl_name', "Plot_"+col_x+'_'+col_y+'.png')
    dpi = kwargs.get('dpi', 300)
    x_label = kwargs.get('x_label', col_x)
    y_label = kwargs.get('y_label', col_y)
    fig, ax = plt.subplots()
    ax.scatter(df[col_x],df[col_y], s = s,
               alpha = alpha, **kwargs)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if save_fig == True:
        plt.savefig(fl_name, dpi=dpi)
    plt.show()
    plt.close()

def get_stat_for_col(col_df):
    return col_df.mean(), col_df.std(), col_df.median()

# ---| HANDY PYTHON FUNCTIONS

# flat a python list
flatten = lambda l: [item for sublist in l for item in sublist]

# rotation of 3d arrays
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
