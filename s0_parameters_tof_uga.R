####################################################################################################
####################################################################################################
## Set environment variables
## Contact remi.dannunzio@fao.org adapted to Ugandan context
## 2018/05/04 AND 2018/12/26
####################################################################################################
####################################################################################################

####################################################################################################
options(stringsAsFactors = FALSE)

packages <- function(x){
  x <- as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}
packages(gfcanalysis)
packages(Hmisc)

### Load necessary packages
library(raster)
library(rgeos)
library(ggplot2)
library(rgdal)
library(stringr)

## Set the working directory
rootdir       <- "~/tof/uganda_draft_201812/"
gfcstore_dir  <- "~/downloads/gfc_2016/tof/"
gadm_dir <- "~/downloads/gfc_2016/"
the_country   <- "UGA"

setwd(rootdir)
rootdir <- paste0(getwd(),"/")

scriptdir<- paste0(rootdir,"scripts/")
data_dir <- paste0(rootdir,"data/")
gfc_dir  <- paste0(rootdir,"data/gfc/")
dd_dir   <- paste0(rootdir,"data/dd_map/")
wmz_dir  <- paste0(rootdir,"data/aoi/")

dir.create(data_dir,showWarnings = F)
dir.create(gadm_dir,showWarnings = F)
dir.create(gfc_dir,showWarnings = F)
dir.create(dd_dir,showWarnings = F)
dir.create(gfcstore_dir,showWarnings = F)

#################### GFC PRODUCTS ADAPTED TO FOREST DEFINITION IN UGANDA (1 HA AND 30%) AND FREL 00_15
gfc_threshold <- 30
beg_year <- 2000
end_year <- 2015
mmu <- 11

#################### PRODUCTS AT THE THRESHOLD
gfc_tc       <- paste0(gfc_dir,"gfc_th",gfc_threshold,"_tc.tif")
gfc_ly       <- paste0(gfc_dir,"gfc_th",gfc_threshold,"_ly.tif")
gfc_gn       <- paste0(gfc_dir,"gfc_gain.tif")
gfc_15       <- paste0(gfc_dir,"gfc_th",gfc_threshold,"_F_",end_year,".tif")
gfc_00       <- paste0(gfc_dir,"gfc_th",gfc_threshold,"_F_",beg_year,".tif")

