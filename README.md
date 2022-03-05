# TCGA-ICGCImmuneFractions

The files contained within this repository are the output of ImmuneCellFractions_TCGA_ICGC.R. They include TCGA and ICGC immune fractions. For information on how to interpret these scores, please read "Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology" by Sturm et al.

R files to generate immune cell fractions from TCGA and ICGC cancer samples using previously published bulk RNA-seq deconvolution methods in the "immunedeconv" package (see https://icbi-lab.github.io/immunedeconv/). Can be used to compare immune cell type proportions within and between samples (depending on specific method used). Functions contained within ImmuneCellFractions_TCGA_ICGC.R can be used to transform TCGA and ICGC data into the appropriate format for use by immunedeconv functions as well as to perform deconvolutions and map data to a common cell vocabulary for easy comparison. 

Please feel free to reach out with any questions.
