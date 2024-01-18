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


usethis::use_data(DNARepairGMTDataList, compress="xz", overwrite=TRUE)
usethis::use_data(CSCEExpData, compress="xz", overwrite=TRUE)
usethis::use_data(AMLExpData, compress="xz", overwrite=TRUE)
usethis::use_data(LHCExpData, compress="xz", overwrite=TRUE)
usethis::use_data(LSCCExpData, compress="xz", overwrite=TRUE)

DNARepairGene="AAAS";XAxisGene="TREM1"
usethis::use_data(DNARepairGene, compress="xz", overwrite=TRUE)
usethis::use_data(XAxisGene, compress="xz", overwrite=TRUE)

# CSCE_ProcExp[[DNARepairGene]]
# usethis::use_data(XAxisGene, compress="xz", overwrite=TRUE)
#

CSCE_ProcExp <- purrr::map(DNARepairGene[1], ~Func_CorrPlotInputDF(DNARepairGene[1], XAxisGene, CSCEExpData)); names(CSCE_ProcExp) <- DNARepairGene[1]
AML_ProcExp <- purrr::map(DNARepairGene[1], ~Func_CorrPlotInputDF(DNARepairGene[1], XAxisGene, AMLExpData)); names(AML_ProcExp) <- DNARepairGene[1]
LHC_ProcExp <- purrr::map(DNARepairGene[1], ~Func_CorrPlotInputDF(DNARepairGene[1], XAxisGene, LHCExpData)); names(LHC_ProcExp) <- DNARepairGene[1]
LSCC_ProcExp <- purrr::map(DNARepairGene[1], ~Func_CorrPlotInputDF(DNARepairGene[1], XAxisGene, LSCCExpData)); names(LSCC_ProcExp) <- DNARepairGene[1]
usethis::use_data(CSCE_ProcExp, compress="xz", overwrite=TRUE)
usethis::use_data(AML_ProcExp, compress="xz", overwrite=TRUE)
usethis::use_data(LHC_ProcExp, compress="xz", overwrite=TRUE)
usethis::use_data(LSCC_ProcExp, compress="xz", overwrite=TRUE)



## inside "Func_CorrPlotInputDF" function,
MyExpData=CSCEExpData

colnames(MyExpData) <- gsub("-","_",colnames(MyExpData)); MyExpData[1:3,1:5]
table(grepl("\\S",MyExpData$Hugo_Symbol));  MyExpData$Hugo_Symbol[!grepl("\\S",MyExpData$Hugo_Symbol)] # 13 20518  [1] ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""

## Remove 'Entrez_Gene_Id' column and  rows of empty gene
MyExpData_RmvEmptyRow <- MyExpData %>% dplyr::select(-Entrez_Gene_Id) %>% dplyr::filter(grepl("\\S" , MyExpData$Hugo_Symbol)); dim(MyExpData_RmvEmptyRow); MyExpData_RmvEmptyRow[1:2,1:3] #  20518   295
## Remove rows of NA and duplicated gene symbol
MyExpData_RmvNARow <- MyExpData_RmvEmptyRow %>% dplyr::filter(!is.na(rowSums(MyExpData_RmvEmptyRow[, 2:ncol(MyExpData_RmvEmptyRow)])),  !duplicated(MyExpData_RmvEmptyRow$Hugo_Symbol)); dim(MyExpData_RmvNARow) #  19979   295

## zscore
MyExpData_RowHugoSymb <- MyExpData_RmvNARow %>% tibble::column_to_rownames("Hugo_Symbol")
# MyExpData_RowHugoSymb_Zscore <- t(apply(MyExpData_RowHugoSymb, 1, cal_z_score)); dim(MyExpData_RowHugoSymb_Zscore) # 5000 33

usethis::use_data(MyExpData_RowHugoSymb, compress="xz", overwrite=TRUE)

NumberVector<- c(1:10)
usethis::use_data(NumberVector, compress="xz", overwrite=TRUE)




## EmptyDF
Xaxis<-c(3:22); Yaxis<-rnorm(20,1,0.5);
EmptyDF<-data.frame(Xaxis, Yaxis)
usethis::use_data(EmptyDF, compress="xz", overwrite=TRUE)

MyColor <- rep("white", 8)
usethis::use_data(MyColor, compress="xz", overwrite=TRUE)



UserInterface <- fluidPage(
  titlePanel("DR.GEXCOR - DNA Repair Gene EXpression COR in CSCE, AML, LHC, and LSCC - "),
  sidebarLayout(
    sidebarPanel(
      helpText("DR.GEXCOR is to see gene expression correlation between DNA repair genes and your interest gene in different cancer types."),
      helpText(h3("Expand the ShinyApp windows, and scroll right to see more plots")),
      selectizeInput("DNARepairGene", label="Gene Query:", choices=DNARepairGMTDataList[[1]],
                     options=list(placeholder="You can choose/type just one gene of DNA repair", maxOptions=30000)), # multiple=TRUE means that a user can input multiple choices
    ),
    mainPanel(
      plotOutput("CorrelationDotplot"),   plotOutput("NoDNARepair")
      # plotOutput("PAM50_Boxplot"),   plotOutput("Histology_Boxplot"), ## br() between plotOutput is to give space between two plots
      # plotOutput("HistologyPAM50_Boxplot"), br(),br(),br(),br(),br(),br(),br(),br(),
      # plotOutput("HistologyER_Boxplot"),
    )
  )
)

usethis::use_data(UserInterface, compress="xz", overwrite=TRUE)
