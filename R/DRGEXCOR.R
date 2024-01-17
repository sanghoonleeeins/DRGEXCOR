#' z-score calibratin
#'
#' @param x data.frame input
#'
#' @return z-score data.frame
#' @export
#'
#' @examples cal_z_score(x)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}


#' Title To make data.frame to be input for plot
#'
#' @param DNARepairGene My interest of DNA repair gene
#' @param XAxisGene TREM1
#' @param MyExpData Expression data.frame
#'
#' @return Processed data.frame to be input for plot
#' @export
#'
#' @examples Func_CorrPlotInputDF(DNARepairGene="AAAS",XAxisGene="TREM1", MyExpData=CSCEExpData)
Func_CorrPlotInputDF <- function(DNARepairGene="AAAS",XAxisGene="TREM1", MyExpData=CSCEExpData) {
  colnames(MyExpData) <- gsub("-","_",colnames(MyExpData)); MyExpData[1:3,1:5]
  table(grepl("\\S",MyExpData$Hugo_Symbol));  MyExpData$Hugo_Symbol[!grepl("\\S",MyExpData$Hugo_Symbol)] # 13 20518  [1] ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""

  ## Remove 'Entrez_Gene_Id' column and  rows of empty gene
  MyExpData_RmvEmptyRow <- MyExpData %>% dplyr::select(-Entrez_Gene_Id) %>% dplyr::filter(grepl("\\S" , MyExpData$Hugo_Symbol)); dim(MyExpData_RmvEmptyRow); MyExpData_RmvEmptyRow[1:2,1:3] #  20518   295
  ## Remove rows of NA and duplicated gene symbol
  MyExpData_RmvNARow <- MyExpData_RmvEmptyRow %>% dplyr::filter(!is.na(rowSums(MyExpData_RmvEmptyRow[, 2:ncol(MyExpData_RmvEmptyRow)])),  !duplicated(MyExpData_RmvEmptyRow$Hugo_Symbol)); dim(MyExpData_RmvNARow) #  19979   295

  ## zscore
  MyExpData_RowHugoSymb <- MyExpData_RmvNARow %>% tibble::column_to_rownames("Hugo_Symbol")
  MyExpData_RowHugoSymb_Zscore <- t(apply(MyExpData_RowHugoSymb, 1, cal_z_score)); dim(MyExpData_RowHugoSymb_Zscore) # 5000 33

  ## Transpose
  MyExpData_RmvNARow_Tp <- MyExpData_RmvNARow %>% tibble::column_to_rownames("Hugo_Symbol") %>% t %>% data.frame; dim(MyExpData_RmvNARow_Tp); MyExpData_RmvNARow_Tp[1:2,1:3] # 294 19979
  MyExpData_RmvNARow_Tp_DNARepairGene <- MyExpData_RmvNARow_Tp %>% dplyr::select_if(colnames(MyExpData_RmvNARow_Tp) %in% c(XAxisGene, DNARepairGene))

  return(MyExpData_RmvNARow_Tp_DNARepairGene)
}









