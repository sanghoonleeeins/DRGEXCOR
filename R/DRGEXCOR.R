#' z-score calibratin
#'
#' @param x data.frame input
#'
#' @return z-score data.frame
#' @export
#'
#' @examples cal_z_score(MyExpData_RowHugoSymb)
cal_z_score <- function(MyExpData_RowHugoSymb){
  (data.frame(MyExpData_RowHugoSymb) - mean(data.frame(MyExpData_RowHugoSymb))) / sd(data.frame(MyExpData_RowHugoSymb))
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
  require(dplyr)
  colnames(MyExpData) <- gsub("-","_",colnames(MyExpData)); MyExpData[1:3,1:5]
  table(grepl("\\S",MyExpData$Hugo_Symbol));  MyExpData$Hugo_Symbol[!grepl("\\S",MyExpData$Hugo_Symbol)] # 13 20518  [1] ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""  ""

  ## Remove 'Entrez_Gene_Id' column and  rows of empty gene
  MyExpData_RmvEmptyRow <- MyExpData %>% dplyr::select(-Entrez_Gene_Id) %>% dplyr::filter(grepl("\\S" , MyExpData$Hugo_Symbol)); dim(MyExpData_RmvEmptyRow); MyExpData_RmvEmptyRow[1:2,1:3] #  20518   295
  ## Remove rows of NA and duplicated gene symbol
  MyExpData_RmvNARow <- MyExpData_RmvEmptyRow %>% dplyr::filter(!is.na(rowSums(MyExpData_RmvEmptyRow[, 2:ncol(MyExpData_RmvEmptyRow)])),  !duplicated(MyExpData_RmvEmptyRow$Hugo_Symbol)); dim(MyExpData_RmvNARow) #  19979   295

  ## zscore
  MyExpData_RowHugoSymb <- MyExpData_RmvNARow %>% tibble::column_to_rownames("Hugo_Symbol")
  MyExpData_RowHugoSymb_Zscore <- t(apply(MyExpData_RowHugoSymb, 1, DRGEXCOR::cal_z_score)); dim(MyExpData_RowHugoSymb_Zscore) # 5000 33

  ## Transpose
  MyExpData_RmvNARow_Tp <- MyExpData_RmvNARow %>% tibble::column_to_rownames("Hugo_Symbol") %>% t %>% data.frame; dim(MyExpData_RmvNARow_Tp); MyExpData_RmvNARow_Tp[1:2,1:3] # 294 19979
  MyExpData_RmvNARow_Tp_DNARepairGene <- MyExpData_RmvNARow_Tp %>% data.frame %>% dplyr::select_if(colnames(MyExpData_RmvNARow_Tp) %in% c(XAxisGene, DNARepairGene))

  return(MyExpData_RmvNARow_Tp_DNARepairGene)
}


#' Title Func_CorPlot_v2
#'
#' @param ExpData Gene exp data frame.
#' @param YAxisGene DNA repair
#' @param XAxisGene TREM1
#' @param CancerType
#'
#' @return a plot
#' @export
#'
#' @examples Func_CorPlot_v2(ExpData=CSCE_ProcExp[[DNARepairGene]],  YAxisGene=DNARepairGene, XAxisGene="TREM1", CancerType="CSCE" )
Func_CorPlot_v2 <- function(ExpData=CSCE_ProcExp[[DNARepairGene]],  YAxisGene=DNARepairGene, XAxisGene="TREM1", CancerType="CSCE" ) {
  ## remove rows that have NA
  ExpData <- ExpData[!(is.na(ExpData[, colnames(ExpData)==XAxisGene]) | is.na(ExpData[, colnames(ExpData)==YAxisGene])), ]; dim(ExpData)
  CorrTest <- cor.test(as.numeric(ExpData[, colnames(ExpData)==XAxisGene]), as.numeric(ExpData[, colnames(ExpData)==YAxisGene]), methods=spearman)
  round(CorrTest$estimate, 5); round(CorrTest$p.value, 5) #

  XMin <- min(ExpData[, colnames(ExpData)==XAxisGene])-0.2; XMax <- max(ExpData[, colnames(ExpData)==XAxisGene])+0.2; print(XMin); print(XMax) # -2.85  3.53
  YMin <- min(ExpData[, colnames(ExpData)==YAxisGene])-0.2; YMax <- max(ExpData[, colnames(ExpData)==YAxisGene])+0.2; print(YMin); print(YMax) # -2.85  3.53

  ## Change TREM1 in column name to "Gene_TREM1"
  colnames(ExpData)[colnames(ExpData) == XAxisGene] <- "Gene_TREM1" #
  ## Change AAAS in column name to "DNARepairGene"
  colnames(ExpData)[colnames(ExpData) == YAxisGene] <- "DNARepairGene" #  "AAAS"

  MyDotPlot <- ggplot2::ggplot(ExpData, ggplot2::aes(Gene_TREM1, DNARepairGene))+ ggplot2::geom_point( alpha=1, size=1, show.legend=T) +  # increase 'alpha' to make the color stronger. 1 is maximum.
    ggplot2::xlim(c(XMin, XMax)) + ggplot2::ylim(c(YMin, YMax)) +
    ggplot2::labs(x=paste0("TREM1 expression"), y=paste0(paste0(YAxisGene, " expression")),
         title=paste0("Spearman Rho=", print(round(CorrTest$estimate,5)), "\npval=", print(round(CorrTest$p.value,5)) ) ) + ggplot2::theme_bw() +  # title="ILC ER+/HER2- (n=13) vs IDC ER+/HER2+ (n=27)"
    ggplot2::theme(#legend.title="MacroFraction", #element_blank(),  #legend.position="none",
      axis.text.x=ggplot2::element_text(colour="Black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  #face="plain" or "bold"
      axis.text.y=ggplot2::element_text(colour="Black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
      plot.title = ggplot2::element_text(size=14),
      axis.title.x=ggplot2::element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
      axis.title.y=ggplot2::element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
      axis.line.y=ggplot2::element_line(color="black", linewidth=0.7),
      axis.line.x=ggplot2::element_line(color="black", linewidth =0.7),
      axis.ticks=ggplot2::element_line(colour="black",linewidth=1),
      axis.ticks.length=ggplot2::unit(.22, "cm"), text=ggplot2::element_text(size=22),
      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    # geom_abline(intercept=intercept, slope=slope, color="black",linetype="dashed", size=1.5)
    ggplot2::stat_smooth(method = "glm", formula = y ~ x, geom = "smooth", color="grey", se=FALSE)

  ## Return the column name to TREM1
  colnames(ExpData)[colnames(ExpData) == "Gene_TREM1"] <- XAxisGene#
  ## Return the column name to the original DNA reapir gene
  colnames(ExpData)[colnames(ExpData) == "DNARepairGene"] <- YAxisGene #  "AAAS"

  return(MyDotPlot)
}


#' Title Func_EmptyPlot
#'
#' @param EmptyDF When no DNA repaire gene
#' @param Xaxis TREM1
#' @param Yaxis DNA repair gene
#' @param YAxisGene
#'
#' @return an empty plot
#' @export
#'
#' @examples Func_EmptyPlot(EmptyDF, Xaxis,Yaxis, YAxisGene)
Func_EmptyPlot <- function(EmptyDF, Xaxis,Yaxis, YAxisGene) {
  MyDotPlot <- ggplot2::ggplot(EmptyDF, ggplot2::aes(Xaxis, Yaxis)) + ggplot2::geom_point( alpha=1, size=1, show.legend=T) +  # increase 'alpha' to make the color stronger. 1 is maximum.
    ggplot2::xlim(c(-1, 1)) + ggplot2::ylim(c(-1, 1 )) +
    ggplot2::ggtitle("<span style='font-size: 30pt;'> This DNA reapir gene doesn't exist in the expression data </font>") + ggplot2::theme(plot.title = ggtext::element_markdown())  + ggplot2::geom_blank()  # + theme_bw()
}



#' Title ShinySever
#'
#' @param input DNA repaire gene
#' @param output a plot
#'
#' @return shiny app
#' @export
#'
#' @examples ShinySever(input, output)
ShinySever <- function(input, output) {
  shiny::observe({
    # MyDNARepairGene <- DNARepairGMTDataList[[1]][1]; print(MyDNARepairGene) ## AAAS
    XAxisGene="TREM1";

    Xaxis<-c(3:22); Yaxis<-rnorm(20,1,0.5); EmptyDF<-data.frame(Xaxis, Yaxis)
    EmptyPlot <- DRGEXCOR::Func_EmptyPlot(EmptyDF, Xaxis,Yaxis, input$DNARepairGene)

    if ( req(input$DNARepairGene) %in% CSCEExpData$Hugo_Symbol ) {
      ## 1) Cervical Squamous cell carcinoma
      ## Process CSCE expression data
      CSCE_ProcExp <- purrr::map(input$DNARepairGene, ~DRGEXCOR::Func_CorrPlotInputDF(input$DNARepairGene, XAxisGene, CSCEExpData)); names(CSCE_ProcExp) <- input$DNARepairGene
      ## Make corr dot plot.
      CSCE_CorrDotPlot <- purrr::map(input$DNARepairGene, ~DRGEXCOR::Func_CorPlot_v2(CSCE_ProcExp[[input$DNARepairGene]], input$DNARepairGene, XAxisGene, "CSCE"))
      CSCE_CorrDotPlot_Arrange <- ggpubr::ggarrange(plotlist = CSCE_CorrDotPlot, ncol = 1)

      ## 2) Acute Myeloid Leukemia
      ## Process AML expression data
      AML_ProcExp <- purrr::map(input$DNARepairGene, ~DRGEXCOR::Func_CorrPlotInputDF(input$DNARepairGene, XAxisGene, AMLExpData)); names(AML_ProcExp) <- input$DNARepairGene
      ## Make corr dot plot.
      AML_CorrDotPlot <- purrr::map(input$DNARepairGene, ~DRGEXCOR::Func_CorPlot_v2(AML_ProcExp[[input$DNARepairGene]], input$DNARepairGene, XAxisGene, "AML"))
      AML_CorrDotPlot_Arrange <- ggpubr::ggarrange(plotlist = AML_CorrDotPlot, ncol = 1)

      ## 3) Liver Hepatocellular Carcinoma
      ## Process LHC expression data
      LHC_ProcExp <- purrr::map(input$DNARepairGene, ~DRGEXCOR::Func_CorrPlotInputDF(input$DNARepairGene, XAxisGene, LHCExpData)); names(LHC_ProcExp) <- input$DNARepairGene
      ## Make corr dot plot.
      LHC_CorrDotPlot <- purrr::map(input$DNARepairGene, ~DRGEXCOR::Func_CorPlot_v2(LHC_ProcExp[[input$DNARepairGene]], input$DNARepairGene, XAxisGene, "LHC"))
      LHC_CorrDotPlot_Arrange <- ggpubr::ggarrange(plotlist = LHC_CorrDotPlot, ncol = 1)

      ## 4) Lung Squamous Cell carcinoma
      ## Process LSCC expression data
      LSCC_ProcExp <- purrr::map(input$DNARepairGene, ~DRGEXCOR::Func_CorrPlotInputDF(input$DNARepairGene, XAxisGene, LSCCExpData)); names(LSCC_ProcExp) <- input$DNARepairGene
      ## Make corr dot plot.
      LSCC_CorrDotPlot <- purrr::map(input$DNARepairGene, ~DRGEXCOR::Func_CorPlot_v2(LSCC_ProcExp[[input$DNARepairGene]], input$DNARepairGene, XAxisGene, "LSCC"))
      LSCC_CorrDotPlot_Arrange <- ggpubr::ggarrange(plotlist = LSCC_CorrDotPlot, ncol = 1)

      PlotArrange_FourCancer <- ggpubr::ggarrange(plotlist = list(CSCE_CorrDotPlot_Arrange, AML_CorrDotPlot_Arrange, LHC_CorrDotPlot_Arrange, LSCC_CorrDotPlot_Arrange),
                                                  ncol=4, nrow=1, labels = c("CSCE","AML","LHC","LSCC"),
                                                  font.label = list(size =16, color="black", face="bold", family = NULL), align = c("h"),
                                                  label.x = 0, label.y = 1, hjust = 0.05, vjust = 0.15) + ggplot2::theme(plot.margin = ggplot2::margin(0.5,0.1,0.1,0.2, "cm"))
      output$CorrelationDotplot <- renderPlot({ PlotArrange_FourCancer }, width=1000, height=260)
      output$NoDNARepair <- NULL

    } else { ## AGO4
      output$CorrelationDotplot <- NULL
      ErrorMessage <- print("This DNA repair gene doesn't exist in the expression data")
      output$NoDNARepair <- shiny::renderPlot({ EmptyPlot })
      output$CorrelationDotplot <- NULL

    } # end of if-else phrase
  }) # close of observe
}  # close of Server function.



# Run the app
# shinyApp(ui=UserInterface, server=ShinySever)





