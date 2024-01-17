## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

####################################################################################################################
DNARepairGMTFile <- "/Volumes/Expansion/CCBR_XWANGLAB10_original/T88_LeukemiaTREM1_Wei_2023/06b_CorrelationPlot_DNARepair/HALLMARK_DNA_REPAIR.v2023.2.Hs.gmt"

CSCCExpFile <- "/Volumes/Expansion/CCBR_XWANGLAB10_original/T88_LeukemiaTREM1_Wei_2023/06a_cBioportal_MulticancerDataset/cesc_tcga_PancanAtlas_2018_CervicalSquamousCellCarcinoma/data_mrna_seq_v2_rsem_zscores_ref_all_samples_CSCC.txt"
AMLExpFile <- "/Volumes/Expansion/CCBR_XWANGLAB10_original/T88_LeukemiaTREM1_Wei_2023/06a_cBioportal_MulticancerDataset/laml_tcga_PancanAtlas_2018_AcuteMyeloidLeukemia/data_mrna_seq_v2_rsem_zscores_ref_all_samples_AML.txt"
LHCExpFile <- "/Volumes/Expansion/CCBR_XWANGLAB10_original/T88_LeukemiaTREM1_Wei_2023/06a_cBioportal_MulticancerDataset/lihc_tcga_PancanAtlas_2018_LiverHepatocellularCarcinoma/data_mrna_seq_v2_rsem_zscores_ref_all_samples_LHC.txt"
LSCCExpFile <- "/Volumes/Expansion/CCBR_XWANGLAB10_original/T88_LeukemiaTREM1_Wei_2023/06a_cBioportal_MulticancerDataset/lusc_tcga_Nature2012_LungSquamouseCellCarcinoma/data_mrna_seq_rpkm_zscores_ref_all_samples_LSCC.txt"

####################################################################################################################

##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### =====
## Step1. Read DNA repair GMT file.
##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### =====
## Read GMT file - GOBP, GOCC
DNARepairGMTDataList <- read_concepts(DNARepairGMTFile, min=2); class(DNARepairGMTDataList); length(DNARepairGMTDataList); print(DNARepairGMTDataList[[1]]); length(DNARepairGMTDataList[[1]]) # 150

##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### =====
## Step2. Read TCGA Gene Exp data
##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### ===== ##### =====
## ======== Cervical Squamous Cell Carcinoma
CSCEExpData <- fread(CSCCExpFile, header=TRUE, stringsAsFactors=FALSE); dim(CSCEExpData); CSCEExpData[1:3,1:3] # 20531   296
## ======== AcuteMyeloidLeukemia
AMLExpData <- fread(AMLExpFile, header=TRUE, stringsAsFactors=FALSE); dim(AMLExpData); AMLExpData[1:3,1:3] # 20531   175
## ======== LiverHepatocellularCarcinoma
LHCExpData <- fread(LHCExpFile, header=TRUE, stringsAsFactors=FALSE); dim(LHCExpData); LHCExpData[1:3,1:3] # 20531   368
## ======== LungSquamousCellCarcinoma
LSCCExpData <- fread(LSCCExpFile, header=TRUE, stringsAsFactors=FALSE); dim(LSCCExpData); LSCCExpData[1:3,1:3] # 19451   180


usethis::use_data(DNARepairGMTDataList, internal=TRUE)
