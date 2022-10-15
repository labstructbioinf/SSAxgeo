# HelixContinuum

Here all the methods/scripts/notebooks related to Helix Continuum project are stores.

---

## How to? (quick and dirty)

### Install

```{bash}

git clone https://github.com/labstructbioinf/SSAxgeo.git
cd SSAxgeo/
pip install -e ./
```


### 1 - Get dataframe for entries in BCX
To obtain a dataframe with all the descriptors and run the SS_assignment just run the _prepare_BCX_df.py_
```{bash}
>$ python prepare_BCX_df.py /path/to/lPDB_metadata.csv bc-50 -ncpus 10
```

    usage: prepare_BCX_df.py [-h] [-ncpus NCPUS] [-save_dir SAVE_DIR]
                              [-canonical_dir CANONICAL_DIR]
                              lpdb_csv bc_group_col

    This script load precompute data from a lPDB copy, get a sample from a given
    clustered PDB, compute normalized and smoothed values for each entry, standardize
    coordinates and run the SSEx assignment on all residues.

    positional arguments:
      lpdb_csv              path to a localPDB metadata csv
      bc_group_col          column name with cluster number for each
                            entry (ex. 'bc-90' for BC90)

    optional arguments:
      -h, --help            show this help message and exit
      -ncpus NCPUS          number of cpus to use (default=1)
      -save_dir SAVE_DIR    set directory for output files (default=Working dir)
      -canonical_dir CANONICAL_DIR
                            set directory of canonical regions to be loaded
                            (default=canonical folder at this script directory)


----
### 2 - Run SSAx assignment on a single entry
To run the Secondary Structure Assignment based on Differential Geometry (xgeo) on a given pdb file, you can use the _'RunSSax.py'_ script.


```{bash}
>$ ssaxgeo pdbfile.pdb xgeofile.csv
```

    usage: ssaxgeo [-h] [-can_dir CAN_DIR]
                      [-out_dir OUT_DIR] [-min_dist MIN_DIST]
                      [-pp2_max PP2_MAX] [-prefix PREFIX]
                      pdb_flpath xgeo_flpath

    This script run the Secondary Structure Assignment based on differential
    geometry (xgeo) descriptors on a specified entry and save results as csv
    files.

    positional arguments:
      pdb_flpath          protein chain coordinate file path
      xgeo_flpath         protein chain xgeo file path

    optional arguments:
      -h, --help          show this help message and exit
      -can_dir CAN_DIR    SSE canonical dataframes directory (default=script dir)
      -out_dir OUT_DIR    Output directory (default=working dir)
      -min_dist MIN_DIST  Minimum distance from canonical Pi/Alpha/3(10) to consider
                          as a label (default=0.2)
      -pp2_max PP2_MAX    Maximum distance from canonical PP2 to consider as PP2
                          (default=0.07)
      -prefix PREFIX      Prefix for ouput files (default=None)
---
### 3 - Generate canonical groups for SSE
To compute the canonical regions used for SSAx, one can use the _GetCanonicalClusters.py_.
```{bash}
>$ python GetCanonicalClusters.py /path/to/grp_state_BC_file.p
```

    usage: GetCanonicalClusters.py [-h] [-save_dir SAVE_DIR]
                           [-kmeans_nclstrs KMEANS_NCLSTRS]
                           [-pp2_cluster_n PP2_CLUSTER_N]
                           [-min_cluster_size MIN_CLUSTER_SIZE]
                           [-pp2_max PP2_MAX] [-Pc_MIN PC_MIN]
                           [-alpha_cluster_n ALPHA_CLUSTER_N]
                           [-pi_cluster_n PI_CLUSTER_N]
                           [-three_cluster_n THREE_CLUSTER_N] [-show_plot]
                           [-plot_debug_data]
                           grp_flpath

    This script run the bootstrap solution to identify to the canonical regions
    based on a representative set of protein structures.

    positional arguments:
      grp_flpath            pickle file path for a PDBx.group object

    optional arguments:
      -h, --help            show this help message and exit
      -save_dir SAVE_DIR    directory to save output files (default= working dir)
      -kmeans_nclstrs KMEANS_NCLSTRS
                        The number of the cluster to be used at k-means
                        (default=6).
      -pp2_cluster_n PP2_CLUSTER_N
                        The number of the cluster which corresponds to
                        PP2(default=2). [WARNING: check the plots!]
      -min_cluster_size MIN_CLUSTER_SIZE
                        Minimum number of points that HDBSCAN can should
                        consider for a cluster (default=15)
      -pp2_max PP2_MAX      Maximum distance from canonical PP2 to consider as
                            PP2 (default=0.07)
      -Pc_MIN PC_MIN        Minimum membership probability for a point to be considered
                            as part of the 'core' (default=0.4)

      -alpha_cluster_n ALPHA_CLUSTER_N
                        cluster number assigned as alpha canonical region
      -pi_cluster_n PI_CLUSTER_N
                        cluster number assigned as pi canonical region
      -three_cluster_n THREE_CLUSTER_N
                        cluster number assigned as 3(10) canonical region
      -show_plot            show generated plots at runing (default=False)
      -plot_debug_data      generate plots for intermediate steps for debugging
                            (default=False)

### Algorithm description
