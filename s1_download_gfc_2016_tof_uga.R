####################################################################################################
####################################################################################################
## DOWNLOAD GFC DATA UGANDA
## Contact remi.dannunzio@fao.org 
## 2017/11/02 AND 2018/12/26
####################################################################################################
####################################################################################################

options(stringsAsFactors = FALSE)
#install.packages("gfcanalysis")
### Load necessary packages
library(gfcanalysis)
library(rgeos)
library(ggplot2)
library(rgdal)
library(raster)


## Country code
countrycode <- c("UGA")

#######################################################################
############  PART I: Check GFC data availability - download if needed
#######################################################################
### Make vector layer of tiles that cover the country
aoi             <- getData('GADM',
                           path=gfcstore_dir, 
                           country= countrycode, 
                           level=0)
tiles           <- calc_gfc_tiles(aoi)
bb <- extent(aoi)

proj4string(tiles) <- proj4string(aoi)
tiles <- tiles[aoi,]

### Find the suffix of the associated GFC data for each tile
tmp         <- data.frame(1:length(tiles),rep("nd",length(tiles)))
names(tmp)  <- c("tile_id","gfc_suffix")

for (n in 1:length(tiles)) {
  gfc_tile <- tiles[n, ]
  min_x <- bbox(gfc_tile)[1, 1]
  max_y <- bbox(gfc_tile)[2, 2]
  if (min_x < 0) {min_x <- paste0(sprintf("%03i", abs(min_x)), "W")}
  else {min_x <- paste0(sprintf("%03i", min_x), "E")}
  if (max_y < 0) {max_y <- paste0(sprintf("%02i", abs(max_y)), "S")}
  else {max_y <- paste0(sprintf("%02i", max_y), "N")}
  tmp[n,2] <- paste0("_", max_y, "_", min_x, ".tif")
}

### Store the information into a SpatialPolygonDF
df_tiles <- SpatialPolygonsDataFrame(tiles,tmp,match.ID = F)
rm(tmp)

### Display the tiles and area of interest (country) to check in R
plot(df_tiles)
plot(aoi,add=T)


download_gfc_2016 <- function (tiles, output_folder, images = c("treecover2000", "loss", 
                                                                "gain", "lossyear", "datamask")) 
{ 
  base <- "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2016-v1.4/"
  
  stopifnot(all(images %in% c("treecover2000", "loss", "gain", 
                              "lossyear", "datamask", "first", "last")))
  if (!file_test("-d", output_folder)) {
    stop("output_folder does not exist")
  }
  message(paste(length(tiles), "tiles to download/check."))
  successes <- 0
  failures <- 0
  skips <- 0
  for (n in 1:length(tiles)) {
    gfc_tile <- tiles[n, ]
    min_x <- bbox(gfc_tile)[1, 1]
    max_y <- bbox(gfc_tile)[2, 2]
    if (min_x < 0) {
      min_x <- paste0(sprintf("%03i", abs(min_x)), "W")
    }
    else {
      min_x <- paste0(sprintf("%03i", min_x), "E")
    }
    if (max_y < 0) {
      max_y <- paste0(sprintf("%02i", abs(max_y)), "S")
    }
    else {
      max_y <- paste0(sprintf("%02i", max_y), "N")
    }
    file_root   <- "Hansen_GFC-2016-v1.4_"
    file_suffix <- paste0("_", max_y, "_", min_x, ".tif")
    filenames   <- paste0(file_root, images, file_suffix)
    tile_urls   <- paste0(paste0(base, filenames))
    local_paths <- file.path(output_folder, filenames)
    for (i in 1:length(filenames)) {
      tile_url <- tile_urls[i]
      local_path <- local_paths[i]
      if (file.exists(local_path)) {
        print(paste0("skipping ",local_path))
        skips <- skips + 1
        next
      }
      # system(sprintf("wget %s -O %s",
      #                tile_url,
      #                local_path))
      download.file(tile_url,local_path,method="auto")
      if (file.exists(local_path)) {
        successes <- successes + 1
      }
      else {
        failures <- failures + 1
      }
    }
  }
  message(paste(successes, "file(s) succeeded,", skips, "file(s) skipped,", 
                failures, "file(s) failed."))
}

### Check if tiles are available and download otherwise : download can take some time
beginCluster()
download_gfc_2016(tiles,
                  gfcstore_dir,
                  images = c("treecover2000","lossyear","gain","datamask"))
endCluster()
df_tiles@data

### MERGE THE TILES TOGETHER, FOR EACH LAYER SEPARATELY
prefix <- "Hansen_GFC-2016-v1.4_"
tiles  <- unlist(paste0(prefix,"datamask",df_tiles@data$gfc_suffix))
tilesx <- substr(tiles,31,38)

types <- c("treecover2000","lossyear","gain","datamask")

for(type in types){
  print(type)
  to_merge <- paste(prefix,type,"_",tilesx,".tif",sep = "")
  system(sprintf("gdal_merge.py -o %s -v -co COMPRESS=LZW %s",
                 paste0(gfc_dir,"tmp_merge_",type,".tif"),
                 paste0(gfcstore_dir,to_merge,collapse = " ")
  ))
  
  system(sprintf("gdal_translate -ot Byte -projwin %s %s %s %s -co COMPRESS=LZW %s %s",
                 floor(bb@xmin),
                 ceiling(bb@ymax),
                 ceiling(bb@xmax),
                 floor(bb@ymin),
                 paste0(gfc_dir,"tmp_merge_",type,".tif"),
                 paste0(gfc_dir,"gfc_",type,".tif")
  ))
  
  system(sprintf("rm %s",
                 paste0(gfc_dir,"tmp_merge_",type,".tif")
  ))
  print(to_merge)
}

