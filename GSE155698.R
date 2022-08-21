library(patchwork)
library(dplyr)
library(Seurat)
library(hdf5r)

PDAC_T1.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710689_PDAC_TISSUE_1/PDAC_TISSUE_1/filtered_feature_bc_matrix/")
PDAC_T2.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710690_PDAC_TISSUE_2/PDAC_TISSUE_2/filtered_feature_bc_matrix/")
PDAC_T3.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710691_PDAC_TISSUE_3/PDAC_TISSUE_3/filtered_feature_bc_matrix/")
PDAC_T4.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710692_PDAC_TISSUE_4/PDAC_TISSUE_4/filtered_feature_bc_matrix/")
PDAC_T5.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710693_PDAC_TISSUE_5/PDAC_TISSUE_5/filtered_feature_bc_matrix/")
PDAC_T6.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710694_PDAC_TISSUE_6/PDAC_TISSUE_6/filtered_feature_bc_matrix/")
PDAC_T7.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710695_PDAC_TISSUE_7/PDAC_TISSUE_7/filtered_feature_bc_matrix/")
PDAC_T8.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710696_PDAC_TISSUE_8/PDAC_TISSUE_8/filtered_feature_bc_matrix/")
PDAC_T9.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710697_PDAC_TISSUE_9/PDAC_TISSUE_9/filtered_feature_bc_matrix/")
PDAC_T10.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710698_PDAC_TISSUE_10/PDAC_TISSUE_10/filtered_feature_bc_matrix/")
PDAC_T11a.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710699_PDAC_TISSUE_11A/PDAC_TISSUE_11A/filtered_feature_bc_matrix/")
PDAC_T11b.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710700_PDAC_TISSUE_11B/PDAC_TISSUE_11B/filtered_feature_bc_matrix/")
PDAC_T12.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710701_PDAC_TISSUE_12/PDAC_TISSUE_12/filtered_feature_bc_matrix/")
PDAC_T13.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710702_PDAC_TISSUE_13/PDAC_TISSUE_13/filtered_feature_bc_matrix/")
PDAC_T14.data<-Read10X_h5("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710703_PDAC_TISSUE_14/PDAC_TISSUE_14/filtered_feature_bc_matrix.h5")
PDAC_T15.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710704_PDAC_TISSUE_15/PDAC_TISSUE_15/filtered_feature_bc_matrix/")
PDAC_T16.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710705_PDAC_TISSUE_16/PDAC_TISSUE_16/filtered_gene_bc_matrices/")
PDAC_adjn1.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710706_AdjNorm_TISSUE_1/AdjNorm_TISSUE_1/filtered_feature_bc_matrix/")
PDAC_adjn2.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710707_AdjNorm_TISSUE_2/AdjNorm_TISSUE_2/filtered_feature_bc_matrix/")
PDAC_adjn3.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710708_AdjNorm_TISSUE_3/AdjNorm_TISSUE_3/filtered_feature_bc_matrix/")
PDAC_pbmc1.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710709_PDAC_PBMC_1/PDAC_PBMC_1/filtered_feature_bc_matrix/")
PDAC_pbmc2.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710710_PDAC_PBMC_2/PDAC_PBMC_2/filtered_feature_bc_matrix/")
PDAC_pbmc3.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710711_PDAC_PBMC_3/PDAC_PBMC_3/filtered_feature_bc_matrix/")
PDAC_pbmc4.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710712_PDAC_PBMC_4/PDAC_PBMC_4/filtered_feature_bc_matrix/")
PDAC_pbmc5.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710713_PDAC_PBMC_5/PDAC_PBMC_5/filtered_feature_bc_matrix/")
PDAC_pbmc6.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710714_PDAC_PBMC_6/PDAC_PBMC_6/filtered_feature_bc_matrix/")
PDAC_pbmc7.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710715_PDAC_PBMC_7/PDAC_PBMC_7/filtered_feature_bc_matrix/")
PDAC_pbmc8.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710716_PDAC_PBMC_8/PDAC_PBMC_8/filtered_feature_bc_matrix/")
PDAC_pbmc9.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710717_PDAC_PBMC_9/PDAC_PBMC_9/filtered_feature_bc_matrix/")
PDAC_pbmc10A.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710718_PDAC_PBMC_10A/PDAC_PBMC_10A/filtered_feature_bc_matrix/")
PDAC_pbmc10B.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710719_PDAC_PBMC_10B/PDAC_PBMC_10B/filtered_feature_bc_matrix/")
PDAC_pbmc11.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710720_PDAC_PBMC_11/PDAC_PBMC_11/filtered_feature_bc_matrix/")
PDAC_pbmc12.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710721_PDAC_PBMC_12/PDAC_PBMC_12/filtered_feature_bc_matrix/")
PDAC_pbmc13.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710722_PDAC_PBMC_13/PDAC_PBMC_13/filtered_feature_bc_matrix/")
PDAC_pbmc14.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710723_PDAC_PBMC_14/PDAC_PBMC_14/filtered_feature_bc_matrix/")
PDAC_pbmc15.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710724_PDAC_PBMC_15/PDAC_PBMC_15/filtered_feature_bc_matrix/")
PDAC_pbmc16.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710725_PDAC_PBMC_16/PDAC_PBMC_16/filtered_gene_bc_matrix/")
Healthy_pbmc1.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710726_Healthy_PBMC_1/Healthy_PBMC_1/filtered_feature_bc_matrix/")
Healthy_pbmc2.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710727_Healthy_PBMC_2/Healthy_PBMC_2/filtered_feature_bc_matrix/")
Healthy_pbmc3.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710728_Healthy_PBMC_3/Healthy_PBMC_3/filtered_feature_bc_matrix/")
Healthy_pbmc4.data<-Read10X("c:/Users/xjmik/Downloads/GSE155698_RAW/GSM4710729_Healthy_PBMC_4/Healthy_PBMC_4/filtered_feature_bc_matrix/")

PDAC_T1.data<-as.data.frame(PDAC_T1.data)
PDAC_T2.data<-as.data.frame(PDAC_T2.data)
PDAC_T3.data<-as.data.frame(PDAC_T3.data)
PDAC_T4.data<-as.data.frame(PDAC_T4.data)
PDAC_T5.data<-as.data.frame(PDAC_T5.data)
PDAC_T6.data<-as.data.frame(PDAC_T6.data)
PDAC_T7.data<-as.data.frame(PDAC_T7.data)
PDAC_T8.data<-as.data.frame(PDAC_T8.data)
PDAC_T9.data<-as.data.frame(PDAC_T9.data)
PDAC_T10.data<-as.data.frame(PDAC_T10.data)
PDAC_T11a.data<-as.data.frame(PDAC_T11a.data)
PDAC_T11b.data<-as.data.frame(PDAC_T11b.data)
PDAC_T12.data<-as.data.frame(PDAC_T12.data)
PDAC_T13.data<-as.data.frame(PDAC_T13.data)
PDAC_T14.data<-as.data.frame(PDAC_T14.data)
PDAC_T15.data<-as.data.frame(PDAC_T15.data)
PDAC_T16.data<-as.data.frame(PDAC_T16.data)
PDAC_adjn1.data<-as.data.frame(PDAC_adjn1.data)
PDAC_adjn2.data<-as.data.frame(PDAC_adjn2.data)
PDAC_adjn3.data<-as.data.frame(PDAC_adjn3.data)
PDAC_pbmc1.data<-as.data.frame(PDAC_pbmc1.data)
PDAC_pbmc2.data<-as.data.frame(PDAC_pbmc2.data)
PDAC_pbmc3.data<-as.data.frame(PDAC_pbmc3.data)
PDAC_pbmc4.data<-as.data.frame(PDAC_pbmc4.data)
PDAC_pbmc5.data<-as.data.frame(PDAC_pbmc5.data)
PDAC_pbmc6.data<-as.data.frame(PDAC_pbmc6.data)
PDAC_pbmc7.data<-as.data.frame(PDAC_pbmc7.data)
PDAC_pbmc8.data<-as.data.frame(PDAC_pbmc8.data)
PDAC_pbmc9.data<-as.data.frame(PDAC_pbmc9.data)
PDAC_pbmc10A.data<-as.data.frame(PDAC_pbmc10A.data)
PDAC_pbmc10B.data<-as.data.frame(PDAC_pbmc10B.data)
PDAC_pbmc11.data<-as.data.frame(PDAC_pbmc11.data)
PDAC_pbmc12.data<-as.data.frame(PDAC_pbmc12.data)
PDAC_pbmc13.data<-as.data.frame(PDAC_pbmc13.data)
PDAC_pbmc14.data<-as.data.frame(PDAC_pbmc14.data)
PDAC_pbmc15.data<-as.data.frame(PDAC_pbmc15.data)
PDAC_pbmc16.data<-as.data.frame(PDAC_pbmc16.data)
Healthy_pbmc1.data<-as.data.frame(Healthy_pbmc1.data)
Healthy_pbmc2.data<-as.data.frame(Healthy_pbmc2.data)
Healthy_pbmc3.data<-as.data.frame(Healthy_pbmc3.data)
Healthy_pbmc4.data<-as.data.frame(Healthy_pbmc4.data)

for (i in 1:length(colnames(PDAC_T1.data))) {
  colnames(PDAC_T1.data)[i] <- paste(colnames(PDAC_T1.data)[i],"PDAC_T1",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T2.data))) {
  colnames(PDAC_T2.data)[i] <- paste(colnames(PDAC_T2.data)[i],"PDAC_T2",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T3.data))) {
  colnames(PDAC_T3.data)[i] <- paste(colnames(PDAC_T3.data)[i],"PDAC_T3",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T4.data))) {
  colnames(PDAC_T4.data)[i] <- paste(colnames(PDAC_T4.data)[i],"PDAC_T4",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T5.data))) {
  colnames(PDAC_T5.data)[i] <- paste(colnames(PDAC_T5.data)[i],"PDAC_T5",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T6.data))) {
  colnames(PDAC_T6.data)[i] <- paste(colnames(PDAC_T6.data)[i],"PDAC_T6",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T7.data))) {
  colnames(PDAC_T7.data)[i] <- paste(colnames(PDAC_T7.data)[i],"PDAC_T7",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T8.data))) {
  colnames(PDAC_T8.data)[i] <- paste(colnames(PDAC_T8.data)[i],"PDAC_T8",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T9.data))) {
  colnames(PDAC_T9.data)[i] <- paste(colnames(PDAC_T9.data)[i],"PDAC_T9",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T10.data))) {
  colnames(PDAC_T10.data)[i] <- paste(colnames(PDAC_T10.data)[i],"PDAC_T10",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T11a.data))) {
  colnames(PDAC_T11a.data)[i] <- paste(colnames(PDAC_T11a.data)[i],"PDAC_T11a",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T11b.data))) {
  colnames(PDAC_T11b.data)[i] <- paste(colnames(PDAC_T11b.data)[i],"PDAC_T11b",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T12.data))) {
  colnames(PDAC_T12.data)[i] <- paste(colnames(PDAC_T12.data)[i],"PDAC_T12",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T13.data))) {
  colnames(PDAC_T13.data)[i] <- paste(colnames(PDAC_T3.data)[i],"PDAC_T13",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T14.data))) {
  colnames(PDAC_T14.data)[i] <- paste(colnames(PDAC_T14.data)[i],"PDAC_T14",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T15.data))) {
  colnames(PDAC_T15.data)[i] <- paste(colnames(PDAC_T15.data)[i],"PDAC_T15",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_T16.data))) {
  colnames(PDAC_T16.data)[i] <- paste(colnames(PDAC_T16.data)[i],"PDAC_T16",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_adjn1.data))) {
  colnames(PDAC_adjn1.data)[i] <- paste(colnames(PDAC_adjn1.data)[i],"PDAC_adjn1",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_adjn2.data))) {
  colnames(PDAC_adjn2.data)[i] <- paste(colnames(PDAC_adjn2.data)[i],"PDAC_adjn2",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_adjn3.data))) {
  colnames(PDAC_adjn3.data)[i] <- paste(colnames(PDAC_adjn3.data)[i],"PDAC_adjn3",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc1.data))) {
  colnames(PDAC_pbmc1.data)[i] <- paste(colnames(PDAC_pbmc1.data)[i],"PDAC_pbmc1",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc2.data))) {
  colnames(PDAC_pbmc2.data)[i] <- paste(colnames(PDAC_pbmc2.data)[i],"PDAC_pbmc2",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc3.data))) {
  colnames(PDAC_pbmc3.data)[i] <- paste(colnames(PDAC_pbmc3.data)[i],"PDAC_pbmc3",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc4.data))) {
  colnames(PDAC_pbmc4.data)[i] <- paste(colnames(PDAC_pbmc4.data)[i],"PDAC_pbmc4",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc5.data))) {
  colnames(PDAC_pbmc5.data)[i] <- paste(colnames(PDAC_pbmc5.data)[i],"PDAC_pbmc5",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc6.data))) {
  colnames(PDAC_pbmc6.data)[i] <- paste(colnames(PDAC_pbmc6.data)[i],"PDAC_pbmc6",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc7.data))) {
  colnames(PDAC_pbmc7.data)[i] <- paste(colnames(PDAC_pbmc7.data)[i],"PDAC_pbmc7",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc8.data))) {
  colnames(PDAC_pbmc8.data)[i] <- paste(colnames(PDAC_pbmc8.data)[i],"PDAC_pbmc8",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc9.data))) {
  colnames(PDAC_pbmc9.data)[i] <- paste(colnames(PDAC_pbmc9.data)[i],"PDAC_pbmc9",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc10A.data))) {
  colnames(PDAC_pbmc10A.data)[i] <- paste(colnames(PDAC_pbmc10A.data)[i],"PDAC_pbmc10A",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc10B.data))) {
  colnames(PDAC_pbmc10B.data)[i] <- paste(colnames(PDAC_pbmc10B.data)[i],"PDAC_pbmc10B",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc11.data))) {
  colnames(PDAC_pbmc11.data)[i] <- paste(colnames(PDAC_pbmc11.data)[i],"PDAC_pbmc11",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc12.data))) {
  colnames(PDAC_pbmc12.data)[i] <- paste(colnames(PDAC_pbmc12.data)[i],"PDAC_pbmc12",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc13.data))) {
  colnames(PDAC_pbmc13.data)[i] <- paste(colnames(PDAC_pbmc13.data)[i],"PDAC_pbmc13",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc14.data))) {
  colnames(PDAC_pbmc14.data)[i] <- paste(colnames(PDAC_pbmc14.data)[i],"PDAC_pbmc14",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc15.data))) {
  colnames(PDAC_pbmc15.data)[i] <- paste(colnames(PDAC_pbmc15.data)[i],"PDAC_pbmc15",i,sep = "-")  
}

for (i in 1:length(colnames(PDAC_pbmc16.data))) {
  colnames(PDAC_pbmc16.data)[i] <- paste(colnames(PDAC_pbmc16.data)[i],"PDAC_pbmc16",i,sep = "-")  
}

for (i in 1:length(colnames(Healthy_pbmc1.data))) {
  colnames(Healthy_pbmc1.data)[i] <- paste(colnames(Healthy_pbmc1.data)[i],"Healthy_pbmc1",i,sep = "-")  
}

for (i in 1:length(colnames(Healthy_pbmc2.data))) {
  colnames(Healthy_pbmc2.data)[i] <- paste(colnames(Healthy_pbmc2.data)[i],"Healthy_pbmc2",i,sep = "-")  
}

for (i in 1:length(colnames(Healthy_pbmc3.data))) {
  colnames(Healthy_pbmc3.data)[i] <- paste(colnames(Healthy_pbmc3.data)[i],"Healthy_pbmc3",i,sep = "-")  
}

for (i in 1:length(colnames(Healthy_pbmc4.data))) {
  colnames(Healthy_pbmc4.data)[i] <- paste(colnames(Healthy_pbmc4.data)[i],"Healthy_pbmc4",i,sep = "-")  
}


PDAC_T1.metadata<-data.frame(colnames(PDAC_T1.data),"T1","T")
PDAC_T2.metadata<-data.frame(colnames(PDAC_T2.data),"T2","T")
PDAC_T3.metadata<-data.frame(colnames(PDAC_T3.data),"T3","T")
PDAC_T4.metadata<-data.frame(colnames(PDAC_T4.data),"T4","T")
PDAC_T5.metadata<-data.frame(colnames(PDAC_T5.data),"T5","T")
PDAC_T6.metadata<-data.frame(colnames(PDAC_T6.data),"T6","T")
PDAC_T7.metadata<-data.frame(colnames(PDAC_T7.data),"T7","T")
PDAC_T8.metadata<-data.frame(colnames(PDAC_T8.data),"T8","T")
PDAC_T9.metadata<-data.frame(colnames(PDAC_T9.data),"T9","T")
PDAC_T10.metadata<-data.frame(colnames(PDAC_T10.data),"T10","T")
PDAC_T11a.metadata<-data.frame(colnames(PDAC_T11a.data),"T11a","T")
PDAC_T11b.metadata<-data.frame(colnames(PDAC_T11b.data),"T11b","T")
PDAC_T12.metadata<-data.frame(colnames(PDAC_T12.data),"T12","T")
PDAC_T13.metadata<-data.frame(colnames(PDAC_T13.data),"T13","T")
PDAC_T14.metadata<-data.frame(colnames(PDAC_T14.data),"T14","T")
PDAC_T15.metadata<-data.frame(colnames(PDAC_T15.data),"T15","T")
PDAC_T16.metadata<-data.frame(colnames(PDAC_T16.data),"T16","T")
PDAC_adjn1.metadata<-data.frame(colnames(PDAC_adjn1.data),"adjn1","adj")
PDAC_adjn2.metadata<-data.frame(colnames(PDAC_adjn2.data),"adjn2","adj")
PDAC_adjn3.metadata<-data.frame(colnames(PDAC_adjn3.data),"adjn3","adj")
PDAC_pbmc1.metadata<-data.frame(colnames(PDAC_pbmc1.data),"P1","P")
PDAC_pbmc2.metadata<-data.frame(colnames(PDAC_pbmc2.data),"P2","P")
PDAC_pbmc3.metadata<-data.frame(colnames(PDAC_pbmc3.data),"P3","P")
PDAC_pbmc4.metadata<-data.frame(colnames(PDAC_pbmc4.data),"P4","P")
PDAC_pbmc5.metadata<-data.frame(colnames(PDAC_pbmc5.data),"P5","P")
PDAC_pbmc6.metadata<-data.frame(colnames(PDAC_pbmc6.data),"P6","P")
PDAC_pbmc7.metadata<-data.frame(colnames(PDAC_pbmc7.data),"P7","P")
PDAC_pbmc8.metadata<-data.frame(colnames(PDAC_pbmc8.data),"P8","P")
PDAC_pbmc9.metadata<-data.frame(colnames(PDAC_pbmc9.data),"P9","P")
PDAC_pbmc10A.metadata<-data.frame(colnames(PDAC_pbmc10A.data),"P10A","P")
PDAC_pbmc10B.metadata<-data.frame(colnames(PDAC_pbmc10B.data),"P10B","P")
PDAC_pbmc11.metadata<-data.frame(colnames(PDAC_pbmc11.data),"P11","P")
PDAC_pbmc12.metadata<-data.frame(colnames(PDAC_pbmc12.data),"P12","P")
PDAC_pbmc13.metadata<-data.frame(colnames(PDAC_pbmc13.data),"P13","P")
PDAC_pbmc14.metadata<-data.frame(colnames(PDAC_pbmc14.data),"P14","P")
PDAC_pbmc15.metadata<-data.frame(colnames(PDAC_pbmc15.data),"P15","P")
PDAC_pbmc16.metadata<-data.frame(colnames(PDAC_pbmc16.data),"P16","P")
Healthy_pbmc1.metadata<-data.frame(colnames(Healthy_pbmc1.data),"H1","H")
Healthy_pbmc2.metadata<-data.frame(colnames(Healthy_pbmc2.data),"H2","H")
Healthy_pbmc3.metadata<-data.frame(colnames(Healthy_pbmc3.data),"H3","H")
Healthy_pbmc4.metadata<-data.frame(colnames(Healthy_pbmc4.data),"H4","H")

rownames(PDAC_T1.metadata)<-PDAC_T1.metadata[,1]
rownames(PDAC_T2.metadata)<-PDAC_T2.metadata[,1]
rownames(PDAC_T3.metadata)<-PDAC_T3.metadata[,1]
rownames(PDAC_T4.metadata)<-PDAC_T4.metadata[,1]
rownames(PDAC_T5.metadata)<-PDAC_T5.metadata[,1]
rownames(PDAC_T6.metadata)<-PDAC_T6.metadata[,1]
rownames(PDAC_T7.metadata)<-PDAC_T7.metadata[,1]
rownames(PDAC_T8.metadata)<-PDAC_T8.metadata[,1]
rownames(PDAC_T9.metadata)<-PDAC_T9.metadata[,1]
rownames(PDAC_T10.metadata)<-PDAC_T10.metadata[,1]
rownames(PDAC_T11a.metadata)<-PDAC_T11a.metadata[,1]
rownames(PDAC_T11b.metadata)<-PDAC_T11b.metadata[,1]
rownames(PDAC_T12.metadata)<-PDAC_T12.metadata[,1]
rownames(PDAC_T13.metadata)<-PDAC_T13.metadata[,1]
rownames(PDAC_T14.metadata)<-PDAC_T14.metadata[,1]
rownames(PDAC_T15.metadata)<-PDAC_T15.metadata[,1]
rownames(PDAC_T16.metadata)<-PDAC_T16.metadata[,1]
rownames(PDAC_adjn1.metadata)<-PDAC_adjn1.metadata[,1]
rownames(PDAC_adjn2.metadata)<-PDAC_adjn2.metadata[,1]
rownames(PDAC_adjn3.metadata)<-PDAC_adjn3.metadata[,1]
rownames(PDAC_pbmc1.metadata)<-PDAC_pbmc1.metadata[,1]
rownames(PDAC_pbmc2.metadata)<-PDAC_pbmc2.metadata[,1]
rownames(PDAC_pbmc3.metadata)<-PDAC_pbmc3.metadata[,1]
rownames(PDAC_pbmc4.metadata)<-PDAC_pbmc4.metadata[,1]
rownames(PDAC_pbmc5.metadata)<-PDAC_pbmc5.metadata[,1]
rownames(PDAC_pbmc6.metadata)<-PDAC_pbmc6.metadata[,1]
rownames(PDAC_pbmc7.metadata)<-PDAC_pbmc7.metadata[,1]
rownames(PDAC_pbmc8.metadata)<-PDAC_pbmc8.metadata[,1]
rownames(PDAC_pbmc9.metadata)<-PDAC_pbmc9.metadata[,1]
rownames(PDAC_pbmc10A.metadata)<-PDAC_pbmc10A.metadata[,1]
rownames(PDAC_pbmc10B.metadata)<-PDAC_pbmc10B.metadata[,1]
rownames(PDAC_pbmc11.metadata)<-PDAC_pbmc11.metadata[,1]
rownames(PDAC_pbmc12.metadata)<-PDAC_pbmc12.metadata[,1]
rownames(PDAC_pbmc13.metadata)<-PDAC_pbmc13.metadata[,1]
rownames(PDAC_pbmc14.metadata)<-PDAC_pbmc14.metadata[,1]
rownames(PDAC_pbmc15.metadata)<-PDAC_pbmc15.metadata[,1]
rownames(PDAC_pbmc16.metadata)<-PDAC_pbmc16.metadata[,1]
rownames(Healthy_pbmc1.metadata)<-Healthy_pbmc1.metadata[,1]
rownames(Healthy_pbmc2.metadata)<-Healthy_pbmc2.metadata[,1]
rownames(Healthy_pbmc3.metadata)<-Healthy_pbmc3.metadata[,1]
rownames(Healthy_pbmc4.metadata)<-Healthy_pbmc4.metadata[,1]

colnames(PDAC_T1.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T2.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T3.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T4.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T5.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T6.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T7.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T8.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T9.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T10.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T11a.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T11b.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T12.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T13.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T14.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T15.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_T16.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_adjn1.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_adjn2.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_adjn3.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc1.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc2.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc3.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc4.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc5.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc6.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc7.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc8.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc9.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc10A.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc10B.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc11.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc12.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc13.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc14.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc15.metadata)<-c("Barcode","group1","group2")
colnames(PDAC_pbmc16.metadata)<-c("Barcode","group1","group2")
colnames(Healthy_pbmc1.metadata)<-c("Barcode","group1","group2")
colnames(Healthy_pbmc2.metadata)<-c("Barcode","group1","group2")
colnames(Healthy_pbmc3.metadata)<-c("Barcode","group1","group2")
colnames(Healthy_pbmc4.metadata)<-c("Barcode","group1","group2")

pbmc.data<-cbind(PDAC_adjn1.data,PDAC_adjn2.data,PDAC_adjn3.data,PDAC_pbmc1.data,PDAC_pbmc10A.data,PDAC_pbmc10B.data,PDAC_pbmc11.data,PDAC_pbmc12.data,PDAC_pbmc13.data,PDAC_pbmc14.data,PDAC_pbmc15.data,PDAC_pbmc16.data,PDAC_pbmc2.data,PDAC_pbmc3.data,PDAC_pbmc4.data,PDAC_pbmc5.data,PDAC_pbmc6.data,PDAC_pbmc7.data,PDAC_pbmc8.data,PDAC_pbmc9.data,
                 PDAC_T1.data,PDAC_T10.data,PDAC_T11a.data,PDAC_T11b.data,PDAC_T12.data,PDAC_T13.data,PDAC_T14.data,PDAC_T15.data,PDAC_T16.data,PDAC_T2.data,PDAC_T3.data,PDAC_T4.data,PDAC_T5.data,PDAC_T6.data,PDAC_T7.data,PDAC_T8.data,PDAC_T9.data,
                 Healthy_pbmc1.data,Healthy_pbmc2.data,Healthy_pbmc3.data,Healthy_pbmc4.data)
pbmc.metadata<-rbind(PDAC_adjn1.metadata,PDAC_adjn2.metadata,PDAC_adjn3.metadata,PDAC_pbmc1.metadata,PDAC_pbmc10A.metadata,PDAC_pbmc10B.metadata,PDAC_pbmc11.metadata,PDAC_pbmc12.metadata,PDAC_pbmc13.metadata,PDAC_pbmc14.metadata,PDAC_pbmc15.metadata,PDAC_pbmc16.metadata,PDAC_pbmc2.metadata,PDAC_pbmc3.metadata,PDAC_pbmc4.metadata,PDAC_pbmc5.metadata,PDAC_pbmc6.metadata,PDAC_pbmc7.metadata,PDAC_pbmc8.metadata,PDAC_pbmc9.metadata,
                     PDAC_T1.metadata,PDAC_T10.metadata,PDAC_T11a.metadata,PDAC_T11b.metadata,PDAC_T12.metadata,PDAC_T13.metadata,PDAC_T14.metadata,PDAC_T15.metadata,PDAC_T16.metadata,PDAC_T2.metadata,PDAC_T3.metadata,PDAC_T4.metadata,PDAC_T5.metadata,PDAC_T6.metadata,PDAC_T7.metadata,PDAC_T8.metadata,PDAC_T9.metadata,
                     Healthy_pbmc1.metadata,Healthy_pbmc2.metadata,Healthy_pbmc3.metadata,Healthy_pbmc4.metadata)
remove(PDAC_adjn1.data,PDAC_adjn2.data,PDAC_adjn3.data,PDAC_pbmc1.data,PDAC_pbmc10A.data,PDAC_pbmc10B.data,PDAC_pbmc11.data,PDAC_pbmc12.data,PDAC_pbmc13.data,PDAC_pbmc14.data,PDAC_pbmc15.data,PDAC_pbmc16.data,PDAC_pbmc2.data,PDAC_pbmc3.data,PDAC_pbmc4.data,PDAC_pbmc5.data,PDAC_pbmc6.data,PDAC_pbmc7.data,PDAC_pbmc8.data,PDAC_pbmc9.data,
       PDAC_T1.data,PDAC_T10.data,PDAC_T11a.data,PDAC_T11b.data,PDAC_T12.data,PDAC_T13.data,PDAC_T14.data,PDAC_T15.data,PDAC_T16.data,PDAC_T2.data,PDAC_T3.data,PDAC_T4.data,PDAC_T5.data,PDAC_T6.data,PDAC_T7.data,PDAC_T8.data,PDAC_T9.data,
       Healthy_pbmc1.data,Healthy_pbmc2.data,Healthy_pbmc3.data,Healthy_pbmc4.data,PDAC_adjn1.metadata,PDAC_adjn2.metadata,PDAC_adjn3.metadata,PDAC_pbmc1.metadata,PDAC_pbmc10A.metadata,PDAC_pbmc10B.metadata,PDAC_pbmc11.metadata,PDAC_pbmc12.metadata,PDAC_pbmc13.metadata,PDAC_pbmc14.metadata,PDAC_pbmc15.metadata,PDAC_pbmc16.metadata,PDAC_pbmc2.metadata,PDAC_pbmc3.metadata,PDAC_pbmc4.metadata,PDAC_pbmc5.metadata,PDAC_pbmc6.metadata,PDAC_pbmc7.metadata,PDAC_pbmc8.metadata,PDAC_pbmc9.metadata,
       PDAC_T1.metadata,PDAC_T10.metadata,PDAC_T11a.metadata,PDAC_T11b.metadata,PDAC_T12.metadata,PDAC_T13.metadata,PDAC_T14.metadata,PDAC_T15.metadata,PDAC_T16.metadata,PDAC_T2.metadata,PDAC_T3.metadata,PDAC_T4.metadata,PDAC_T5.metadata,PDAC_T6.metadata,PDAC_T7.metadata,PDAC_T8.metadata,PDAC_T9.metadata,
       Healthy_pbmc1.metadata,Healthy_pbmc2.metadata,Healthy_pbmc3.metadata,Healthy_pbmc4.metadata)
remove(i)
pbmc<-CreateSeuratObject(counts = pbmc.data,meta.data = pbmc.metadata,min.cells = 3,min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 2)
pbmc <- RunUMAP(pbmc, dims = 1:30)
DimPlot(pbmc)
FeaturePlot(pbmc,features = "FOXP3")
VlnPlot(pbmc,features = "FOXP3")
Idents(pbmc)<-pbmc@meta.data$group2
DimPlot(pbmc)
pbmc_T_P<-subset(pbmc,idents = c("T","P"))
pbmc_T_P <- NormalizeData(pbmc_T_P)
pbmc_T_P <- FindVariableFeatures(pbmc_T_P, selection.method = "vst", nfeatures = 2000)
pbmc_T_P <- ScaleData(pbmc_T_P, features = VariableFeatures(object = pbmc_T_P))
pbmc_T_P <- RunPCA(pbmc_T_P, features = VariableFeatures(object = pbmc_T_P))
ElbowPlot(pbmc_T_P,ndims = 50)
pbmc_T_P <- FindNeighbors(pbmc_T_P, dims = 1:22)
pbmc_T_P <- FindClusters(pbmc_T_P, resolution = 2)
pbmc_T_P <- RunUMAP(pbmc_T_P, dims = 1:22)
DimPlot(pbmc_T_P)
FeaturePlot(pbmc_T_P,features = "FOXP3")
VlnPlot(pbmc_T_P,features = "FOXP3")
Treg<-subset(pbmc_T_P,idents = "18")
Treg <- NormalizeData(Treg)
Treg <- FindVariableFeatures(Treg, selection.method = "vst", nfeatures = 2000)
Treg <- ScaleData(Treg, features = VariableFeatures(object = Treg))
Treg <- RunPCA(Treg, features = VariableFeatures(object = Treg))
ElbowPlot(Treg,ndims = 50)
Treg <- FindNeighbors(Treg, dims = 1:14)
Treg <- FindClusters(Treg, resolution = 1)
Treg <- RunUMAP(Treg, dims = 1:14)
DimPlot(Treg)
FeaturePlot(Treg,features = "FOXP3")
VlnPlot(Treg,features = "FOXP3",sort = TRUE)
Treg_new<-subset(Treg,idents = c("0","3","4","8","9"))
Idents(Treg_new)<-Treg_new@meta.data$group2
VlnPlot(Treg_new,features = "JMJD1C",sort = TRUE,pt.size = 0)
Idents(Treg)<-Treg@meta.data$group2
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)
FeaturePlot(Treg,features = "JMJD1C",split.by = "group2")
FeaturePlot(Treg_new,features = "JMJD1C",split.by = "group2",order = TRUE)
