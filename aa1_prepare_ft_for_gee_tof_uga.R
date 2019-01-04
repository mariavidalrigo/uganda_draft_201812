####################################################################################################
####################################################################################################
## Prepare data for time series clipping for TOF in Uganda: create the KMl of the CE points and the grids 
## in Landsat and Sentinel to include all of them
## Contact remi.dannunzio@fao.org
## 2018/08/24 and 2019/01/03
####################################################################################################
####################################################################################################
options(stringsAsFactors=FALSE)

library(Hmisc)
library(sp)
library(rgdal)
library(raster)
library(plyr)
library(foreign)


#######################################################################
##############################     SETUP YOUR DATA 
#######################################################################
## Set the working directory
rootdir       <- "~/tof/uganda_draft_201812/"

## Go to the root directory
setwd(rootdir)

rootdir <- paste0(getwd(),"/")

dd_dir <- paste0(rootdir,"data/dd_map/")
the_map <- paste0(dd_dir,"dd_map_0015_gt30_20181226.tif")

sae_dir  <- paste0(dirname(the_map),"/",
                   "sae_design_",
                   substr(basename(the_map),
                          1,
                          nchar(basename(the_map))-4
                          ),
                   "/"
                   )

## Read the datafile and setup the correct names for the variables
point_file <- list.files(sae_dir,glob2rx("pts_CE*.csv"))
pts <- read.csv(paste0(sae_dir,point_file))

## Check that names match
names(pts)
map_code <- "map_class"
point_id <- "id"
xcoord   <- "XCoordinate"
ycoord   <- "YCoordinate"


##############################     Perform some checks
names(pts)
head(pts)
summary(pts)

###    Are the point ID unique indeed ?
nrow(pts)==length(unique(pts[,point_id]))


#######################################################################
### PART I: Create square boxes around the sampling points (2km) 
### each of the points is identified as square boxes that will be used as KML in GEE to identify the landsat and sentinel tiles to be downloaded
#######################################################################

############### Set the size of the boxes for the points of sampling (2km = 2*1000m)
ysize <- 1000/111321


#################### Loop through points to create a box around each point
lp<-list()

for(i in 1:nrow(pts)){
  ymin <- pts[i,ycoord]-ysize
  ymax <- pts[i,ycoord]+ysize
  xmin <- pts[i,xcoord]-ysize*cos(pts[1,ycoord]*pi/180)
  xmax <- pts[i,xcoord]+ysize*cos(pts[1,ycoord]*pi/180)
  
  p  <- Polygon(cbind(c(xmin,xmin,xmax,xmax,xmin),c(ymin,ymax,ymax,ymin,ymin)))
  ps <- Polygons(list(p), pts[i,1])
  lp <- append(lp,list(ps))
}

#################### Transform the list of polygons into a SPDF
spdf<-SpatialPolygonsDataFrame(
  SpatialPolygons(lp,1:nrow(pts)), 
  pts[,c(map_code,point_id)],#,xcoord,ycoord)], 
  match.ID = F
)

proj4string(spdf)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

head(spdf)

#################### Export as shapefile or KML
writeOGR(obj=spdf,dsn=paste0(sae_dir,"pts_2km_boxes.kml"),layer="pts_2km_boxes",driver="KML",overwrite_layer = T)
#writeOGR(obj=spdf,dsn="pts_2km_boxes.shp",layer="pts_2km_boxes",driver="ESRI Shapefile",overwrite_layer = T)


#######################################################################
### PART II: Create rectangular grid for export in GEE-API LANDSAT grid of tiles to be created and downloaded
#######################################################################

### What grid size do we need ? (in degrees)
grid_size = 2           ## for Landsat  @ 30m spatial resolution

### Create a set of regular SpatialPoints on the extent of the created polygons  
sqr <- SpatialPoints(makegrid(spdf,offset=c(0,0),cellsize = grid_size))

### Convert points to a square grid
grid <- points2grid(sqr)

### Convert the grid to SpatialPolygonDataFrame
SpP_grd <- as.SpatialPolygons.GridTopology(grid)

sqr_df <- SpatialPolygonsDataFrame(
  Sr=SpP_grd,
  data=data.frame(rep(1,length(SpP_grd))),
  match.ID=F)

### Assign the right projection
proj4string(sqr_df) <- proj4string(spdf)
sqr_df_selected <- sqr_df[spdf,]

### Plot the results
plot(sqr_df_selected)
plot(spdf,add=T)

### Give the output a decent name, with unique ID
names(sqr_df_selected) <- "tileID"
sqr_df_selected@data$tileID <-row(sqr_df_selected@data)[,1]

### Check how many tiles will be created
nrow(sqr_df_selected@data)


#######################################################################
### PART III: Export as KML
#######################################################################
base_sqr <- "download_area_grid_lsat"
writeOGR(obj=sqr_df_selected,dsn=paste0(sae_dir,base_sqr,".kml"),layer=base_sqr,driver = "KML",overwrite_layer = T)


#######################################################################
### PART IV: Create rectangular grid for export in GEE-API SENTINEL idem sentinel
#######################################################################

### What grid size do we need ? (in degrees)
grid_size = 10000/11200 ## for Sentinel @ 10m spatial resolution

### Create a set of regular SpatialPoints on the extent of the created polygons  
sqr <- SpatialPoints(makegrid(spdf,offset=c(0,0),cellsize = grid_size))

### Convert points to a square grid
grid <- points2grid(sqr)

### Convert the grid to SpatialPolygonDataFrame
SpP_grd <- as.SpatialPolygons.GridTopology(grid)

sqr_df <- SpatialPolygonsDataFrame(
  Sr=SpP_grd,
  data=data.frame(rep(1,length(SpP_grd))),
  match.ID=F)

### Assign the right projection
proj4string(sqr_df) <- proj4string(spdf)
sqr_df_selected <- sqr_df[spdf,]

### Plot the results
plot(sqr_df_selected)
plot(spdf,add=T)

### Give the output a decent name, with unique ID
names(sqr_df_selected) <- "tileID"
sqr_df_selected@data$tileID <-row(sqr_df_selected@data)[,1]

### Check how many tiles will be created
nrow(sqr_df_selected@data)


#######################################################################
### PART V: Export as KML
#######################################################################
base_sqr <- "download_area_grid_stnl"
writeOGR(obj=sqr_df_selected,dsn=paste0(sae_dir,base_sqr,".kml"),layer=base_sqr,driver = "KML",overwrite_layer = T)

