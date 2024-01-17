# DR_GEXCOR
Correlation of Gene Expression between TREM1 and DNA repair genes 
---
title: "DR_GEXCOR"
author: "SanghoonLee"
date: "2024-01-17"
output: 
  html_document: default
  pdf_document: default
---
<span style="color:red">**D**</span>NA <span style="color:red">**R**</span>pair <span style="color:red">**G**</span>ene <span style="color:red">**Ex**</span>pression <span style="color:red">**COR**</span>relation (DR_GEXCOR) was built on R version 4.3.1 (2023-06-16) -- "Beagle Scouts"  This means you need to install R version >= 4.3.1
DR_GEXCOR was tested on both Windows and Mac environment.


## A. Instruction 
DR_GEXCOR is a tool to visualize the correlation of DNA repair gene and your interst gene expression in multiple cancer types 

This is an instruction to explain how to use DR_GEXCOR. Considering your computer programming skills, I made a R package to minimize your computer work and coding. Here are brief overview to use DR_GEXCOR.

## B. Basic requirements

### Step1. Download R and Rstudio, and install them. 
You can just google for **"R download and install"** and **"Rstudio download"** to complete this step. If you have some experience already,

You can install R from CRAN: https://cran.r-project.org/

You can install a user-friendly interface, R-Studio, from here: https://www.rstudio.com/products/rstudio/download/

If you want to learn more in structure or if you got stuck, before asking me, visit https://hsls.libguides.com/video_intro2R There are a video lecture, ppt slides, and whole basic explanation.

### Step2. Start Rstudio
What does 'start Rstudio' mean? Visit the lecture I introduced in Step1.

### Step3. In the R console, install necessary R packages and install DR_GEXCOR. 
What are 'R console' and 'R package'? The lecture I introduced in Step1 will explain everything.

Simply, run these code lines in your R console. Just copy the lines, paste them to console, and hit enter key.

```{r}
# Install and load multiple packages at a time.  This will take 10~15 min if you are installing any packages for the first time.
# install.packages(c("devtools","usethis","pacman","data.table","dplyr","tidyr","ggplot2","ggpubr","purrr","shiny"), repos="http://cran.us.r-project.org")
# pacman::p_load(devtools,usethis,pacman,data.table,dplyr,tidyr,ggplot2,ggpubr,purrr,shiny)
```

![Installing required packages](/Users/sanghoonlee/Library/CloudStorage/OneDrive-UniversityofPittsburgh/H45_ShinyApp_METABRICSCANB_TROP2/06b_Rpacakge_GEXPLOER/InstalldevtoolsPackage.png)

### Step4. Download DR_GEXCOR R package file from UPMC OneDrive
I couldn't store DR_GEXCOR R package to the lab Github because METABRIC and SCAN-B data file sizes are too big, 218 Mb and 375 Mb. Github doesn't allow to upload a big size file > 25 Mb.

Therefore, I made R package .zip file and stored it at the lab Pitt OneDrive. Click [here](https://pitt-my.sharepoint.com/:f:/r/personal/avleelab_pitt_edu/Documents/Lee-Oesterreich%20Lab/Lee-Oesterreich%20General%20Lab%20Items/Bioinformatics_Basics_And_Scripts/07_DR_GEXCOR?csf=1&web=1&e=gMN2lr) or go to 

**Lee-Oesterreich lab > Lee-Oesterreich General Lab Items > Bioinformatics_Basics_And_Scripts > 07_DR_GEXCOR**

Download "DR_GEXCOR_0.0.1.1.tar.gz" (650 Mb). It takes long due to the file size. Store the .tar.gz file to your local computer, but you don't need to untar the file. For example, I stored it at "..../H45_ShinyApp_METABRICSCANB_TROP2" (screenshot below)

**Mac users:** Mouse right click on "DR_GEXCOR_0.0.1.1.tar.gz" file and click "Get Info." Drag the folder address and click Command+C (screenshot below)
Window users: You can copy the address in Windows Explorer.

![Download DR_GEXCOR_0.0.1.1.tar.gz file](/Users/sanghoonlee/Library/CloudStorage/OneDrive-UniversityofPittsburgh/H45_ShinyApp_METABRICSCANB_TROP2/06b_Rpacakge_GEXPLOER/GEXPLORERZipDownload.png)
### Step5. Install DR_GEXCOR R package. 
In the R console, let's change your working directory to the folder you stored "DR_GEXCOR_0.0.1.1.tar.gz." The code line should be like this. Look at the screenshot below. 

setwd("/DirectoryAddress/YouCopied/FromGetInfo)    # setwd means set working directory. 

**Windows users:** When you copy and paste your address in your R code, the folder split is backslash like this, "\user\MyFolder\" You should revise it to slash. R doesn't understand backslash. e.g. setwd(/user/MyFolder/")

```{r}
## Let's change your working directory to the folder you stored "DR_GEXCOR_0.0.11.tar.gz" file. Note you need quotation inside parenthesis, like setwd("YourFolderAddress")
# setwd("/Volumes/Expansion/CCBR_XWANGLAB10_original/T88_LeukemiaTREM1_Wei_2023/06c_ShinyAppRpackage_CorrPlot_DNARepair/DR_GEXCOR")

## Install DR_GEXCOR. This takes 2~5 min depending on your computer spec.  Note you need quotation outside the file name. 
# install.packages("DR_GEXCOR_0.0.1.1.tar.gz")

## Load your DR_GEXCOR package. Note you DON'T need quotation. The package name is DR_GEXCOR, not DR_GEXCOR_0.0.1.1.tar.gz
# library(DR_GEXCOR) 
```

Once you install DR_GEXCOR package one time, you don't need to install it again when you restarted your R session. Just you need to load DR_GEXCOR 
> library(DR_GEXCOR)

![Install DR_GEXCOR package](/Users/sanghoonlee/Library/CloudStorage/OneDrive-UniversityofPittsburgh/H45_ShinyApp_METABRICSCANB_TROP2/06b_Rpacakge_GEXPLOER/InstallGEXPLORERZip.png)


## C. Play with DR_GEXCOR

Copy the code line below and run it in your R console. 

```{r}
## DR_GEXCOR_CellLine. It takes about 10 seconds to start a ShinyApp, depending on your computer spec.
# shiny::shinyApp(ui=ui_CellLine, server=server_CellLine)
```

```{r}
## DR_GEXCOR-TMS. It takes about 10 seconds to start a ShinyApp, depending on your computer spec. 
# shiny::shinyApp(ui=ui_TMS, server=server_TMS)
```

In the Gene Query box, type in your gene of interest. You can choose the dataset and clinical subtype of your interest below. 


![Caption for the picture.](/Users/sanghoonlee/Library/CloudStorage/OneDrive-UniversityofPittsburgh/H45_ShinyApp_METABRICSCANB_TROP2/06b_Rpacakge_GEXPLOER/GEXPLORER_TMS.png)

### FAQ
1) While playing with Shiny App, the ShinyApp windows is suddenly closed and I can't see a cursor in R condole like a screenshot below. What should I do?
Answer: Click "Stop" button and restart the ShinyApp (screenshot below)

# ![Caption for the picture.](/Users/sanghoonlee/Library/CloudStorage/OneDrive-UniversityofPittsburgh/H45_ShinyApp_METABRICSCANB_TROP2/06b_Rpacakge_GEXPLOER/ShinyAppShutdown.png)

