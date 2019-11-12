mkdir data
cd data

# PBMC 10K Healthy V3
#mkdir pbmc_10k_v3
#pushd pbmc_10k_v3
#wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_raw_feature_bc_matrix.h5 -O raw_feature_bc_matrix.h5
#wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5 -O filtered_feature_bc_matrix.h5
#popd

# PBMC 1K Healthy V3
#mkdir pbmc_1k_v3
#pushd pbmc_1k_v3
#wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_raw_feature_bc_matrix.h5 -O raw_feature_bc_matrix.h5
#wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5 -O filtered_feature_bc_matrix.h5
#popd

# PBMC 1K Healthy V2
#mkdir pbmc_1k_v2
#pushd pbmc_1k_v2
#wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5 -O raw_feature_bc_matrix.h5
#wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5 -O filtered_feature_bc_matrix.h5
#popd

# mouse 10k heart cells
mkdir mouse_heart_10k_v3
pushd mouse_heart_10k_v3
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_10k_v3/heart_10k_v3_raw_feature_bc_matrix.h5 -O raw_feature_bc_matrix.h5
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_10k_v3/heart_10k_v3_filtered_feature_bc_matrix.h5 -O filtered_feature_bc_matrix.h5
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_10k_v3/heart_10k_v3_filtered_feature_bc_matrix.tar.gz -O filtered_feature_bc_matrix.tar.gz
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_10k_v3/heart_10k_v3_raw_feature_bc_matrix.tar.gz -O raw_feature_bc_matrix.tar.gz
tar xzvf filtered_feature_bc_matrix.tar.gz
tar xzvf raw_feature_bc_matrix.tar.gz
popd

# mouse 10k heart cells
mkdir mouse_brain_10k_v3
pushd mouse_brain_10k_v3
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_raw_feature_bc_matrix.h5 -O raw_feature_bc_matrix.h5
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_filtered_feature_bc_matrix.h5 -O filtered_feature_bc_matrix.h5
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_filtered_feature_bc_matrix.tar.gz -O filtered_feature_bc_matrix.tar.gz
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_raw_feature_bc_matrix.tar.gz -O raw_feature_bc_matrix.tar.gz
tar xzvf filtered_feature_bc_matrix.tar.gz
tar xzvf raw_feature_bc_matrix.tar.gz
popd

# mouse 9k brain cells v2
mkdir mouse_brain_9k_v2
pushd mouse_brain_9k_v2
wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neuron_9k/neuron_9k_raw_gene_bc_matrices_h5.h5 -O raw_gene_bc_matrices.h5
#wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neuron_9k/neuron_9k_filtered_gene_bc_matrices_h5.h5 -O filtered_gene_bc_matrices.h5
wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neuron_9k/neuron_9k_filtered_gene_bc_matrices.tar.gz -O filtered_gene_bc_matrices.tar.gz
wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neuron_9k/neuron_9k_raw_gene_bc_matrices.tar.gz -O raw_gene_bc_matrices.tar.gz
tar xzvf filtered_gene_bc_matrices.tar.gz
tar xzvf raw_gene_bc_matrices.tar.gz
popd