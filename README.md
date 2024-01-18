---
title: "DrGEXCOR"
author: "SanghoonLee"
date: "2024-01-17"
output: 
  html_document: default
  pdf_document: default
---
<span style="color:red">**D**</span>NA <span style="color:red">**r**</span>epair <span style="color:red">**G**</span>ene <span style="color:red">**EX**</span>pression <span style="color:red">**COR**</span>relation (DrGEXCOR) was built on R version 4.3.1 (2023-06-16) -- "Beagle Scouts"  This means you need to install R version >= 4.3.1
DrGEXCOR was tested on both Windows and Mac environment.


## A. Instruction 
DrGEXCOR is a tool to visualize the correlation of DNA repair gene and your interst gene expression in multiple cancer types 

This is an instruction to explain how to use DrGEXCOR. Considering your computer programming skills, I made a R package to minimize your computer work and coding. Here are brief overview to use DrGEXCOR.

## B. Basic requirements

### Step1. Download R and Rstudio, and install them. 
You can just google for **"R download and install"** and **"Rstudio download"** to complete this step. If you have some experience already,

You can install R from CRAN: https://cran.r-project.org/

You can install a user-friendly interface, R-Studio, from here: https://www.rstudio.com/products/rstudio/download/

If you want to learn more in structure or if you got stuck, before asking me, visit https://hsls.libguides.com/video_intro2R There are a video lecture, ppt slides, and whole basic explanation.

### Step2. Start Rstudio
What does 'start Rstudio' mean? Visit the lecture I introduced in Step1.

### Step3. In the R console, install necessary R packages and install DrGEXCOR. 
What are 'R console' and 'R package'? The lecture video I introduced in Step1 will explain everything.

Simply, run these code lines in your R console. Just copy the lines, paste them to console, and hit enter key.

```{r}
# Install and load multiple packages at a time.  This will take 10~15 min if you are installing any packages for the first time.
nstall.packages(c("devtools","usethis","pacman","data.table","dplyr","tidyr","ggplot2","ggpubr","purrr","shiny"), repos="http://cran.us.r-project.org")
pacman::p_load(devtools,usethis,pacman,data.table,dplyr,tidyr,ggplot2,ggpubr,purrr,shiny)
```

![Installing required packages](/Users/sanghoonlee/Library/CloudStorage/OneDrive-UniversityofPittsburgh/H45_ShinyApp_METABRICSCANB_TROP2/06b_Rpacakge_GEXPLOER/InstalldevtoolsPackage.png)

### Step4. Install DrGEXCOR R package

```{r}
# Install DRGEXCOR R package. This takes 1~2 minutes. If you see a message like a screenshot blow, type "1" and hit Enter. 
devtools::install_github("sanghoonleeeins/DRGEXCOR")   

# Load your DrGEXCOR package. Note you DON'T need quotation. 
library(DRGEXCOR)
```

![Install Dr.GEXCOR package](/Volumes/Expansion/CCBR_XWANGLAB10_original/T88_LeukemiaTREM1_Wei_2023/06d_ShinyApp_Rpackage_Github/InstallDRGEXCOR.png)

</br>

Once you install DrGEXCOR package one time, you don't need to install it again when you restarted your R session. Just you need to load DrGEXCOR 

> library(DRGEXCOR)


## C. Play with DrGEXCOR

Copy the code line below and run it in your R console. 

```{r}
## DrGEXCOR. It takes about 10 seconds to start a ShinyApp depending on your computer spec.
shinyApp(ui=UserInterface, server=ShinySever)
```

For now, you should choose a gene of DNA repair in the Gene Query box. You can't type in your gene of interest. 

![Caption for the picture.](/Volumes/Expansion/CCBR_XWANGLAB10_original/T88_LeukemiaTREM1_Wei_2023/06d_ShinyApp_Rpackage_Github/RunDRGEXCORShinyApp.png)


