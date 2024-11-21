# BACT_code_data

Code and Data for reproducing Figures and Table in the manuscript "BACT: nonparametric Bayesian cell typing for single-cell spatial transcriptomics data"



# Remark

Since GitHub's file size limit is 100.00 MB, the raw data files---`BeadLocationsForR.csv` and `MappedDGEForR.csv` for Slide-seq dataset, and `151507.h5ad` for DLPFC section 151507---need to be downloaded by the users themselves, and stored in the "input_data" folder. The two raw data files for Slide-seq are in the barcode file "Puck_180430_1.tar.gz" which can be downloaded from the Broad institute's single-cell repository: 
https://singlecell.broadinstitute.org/single_cell/study/SCP354/slide-seq-study#study-download. 
However, this situation does not hinder users from implementing BACT, because the input data of BACT has already been stored in "coord_and_pc.RData."

Similarly, the processed data `Real_application/Slide_seq/input_data/gene_raw_counts.csv` and `Real_application/DLPFC/input_data/gene_matrix_raw_151507.csv` also exceed GitHub's file size limit of 100.00 MB, so the users need to generate these files via the codes in "code" folder by themselves.
