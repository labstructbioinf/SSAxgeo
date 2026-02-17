python ssaxgeo_getSampleOfClstrPDB /home/db/localpdb/ -out_dir ../../data_set/ -redundancy 90 -res_lim 2.0 -ncpus 10 -seed 0

# full dataset
python ssaxgeo_computePDBxgeo --mylocalpdb_path ~/tmp/lpdb/ --sampled_clstrd_path ../../data_set/sampled_clust-90.csv --xgeo_output_dir /home/users/sdunin/tmp/lpdb/xgeo_chains --ncpus 10 --out_csv ../../data_set/sampled_clust-90_updated.csv
python ssaxgeo_clusterResidues ../../data_set/sampled_clust-90_updated.csv clust-90 -ncpus 8 --do_res_labeling -canonical_dir ../canonical/

# test dataset
python ssaxgeo_computePDBxgeo --mylocalpdb_path ~/tmp/lpdb/ --sampled_clstrd_path ../../data_set/mini.csv --xgeo_output_dir /home/users/sdunin/tmp/lpdb/xgeo_chains --ncpus 10 --out_csv ../../data_set/mini_upd.csv 
python ssaxgeo_clusterResidues ../../data_set/mini_upd.csv clust-90 -ncpus 8 --do_res_labeling -canonical_dir ../canonical/


