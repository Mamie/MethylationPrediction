# File name: setup.sh
# Author: Mamie Wang
# Date created: 02/05/2018
# Date last modified: 02/05/2018
# Set up the data folder and download the data

mkdir data
cd data
wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LAML/20160128/gdac.broadinstitute.org_LAML.mRNAseq_Preprocess.Level_3.2016012800.0.0.tar.gz -O geneExpression.tar.gz
mkdir geneExpression
tar xzf geneExpression.tar.gz -C geneExpression --strip-components 1
rm geneExpression.tar.gz

wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LAML/20160128/gdac.broadinstitute.org_LAML.Merge_Clinical.Level_1.2016012800.0.0.tar.gz -O clinical.tar.gz
mkdir clinical
tar xzf clinical.tar.gz -C clinical --strip-components 1
rm clinical.tar.gz

wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LAML/20160128/gdac.broadinstitute.org_LAML.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz -O methylation_hm27.tar.gz
mkdir methylationHM27
tar xzf methylation_hm27.tar.gz -C methylationHM27 --strip-components 1
rm methylation_hm27.tar.gz

wget https://cancergenome.nih.gov/abouttcga/aboutdata/platformdesign/illumina-infinium-methylation-27-adf -O HM27adf.zip
unzip HM27adf.zip
rm HM27adf.zip
rm jhu-usc.edu_TCGA_HumanMethylation27.v2.adf.txt
rm DESCRIPTION.txt

mkdir processed
mkdir processed/Mcluster
mkdir processed/coefs
mkdir processed/GEP
mkdir figures
mkdir figures/moduleVisualization

cd ../R
echo "Preprocessing methylation and mRNA-seq data"
Rscript preprocessing.R
Rscript model.R
Rscript visualization.R

