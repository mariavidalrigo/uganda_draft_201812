########################################################################################## 
# Contact: remi.dannunzio@fao.org
# Last update: 2018-11-28 ADAPTED TO UGANDAN CONTEXT 2018/12/26
##########################################################################################

time_start  <- Sys.time()

aoi <- paste0(wmz_dir,"WMZ_latlong.shp")
aoi_field <- "WMZ_code"

####################################################################################
####### COMBINE GFC LAYERS
####################################################################################

#################### CREATE GFC TREE COVER MAP IN 2000 AT THRESHOLD 
### diferenciate pixels in TC2000 higher than 30% as forests in 2000 (A>threshold) in 2000 and losses are ALL included as losses
### as we consider 2000 to 2015 (no losses from gfc_lossyear are incorporated as non forest as in NIgeria)

system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(gfc_dir,"gfc_treecover2000.tif"),
               paste0(gfc_dir,"gfc_lossyear.tif"),
               paste0(dd_dir,"tmp_gfc_tc_start.tif"),
               paste0("(A>",gfc_threshold,")*((B==0)+(B>1))")
))

#################### SIEVE TO THE MMU
### probably this is to separate the polygons of minimum 1 ha (11 pixels) and the ones smaller. The smaller ones, if they have minimum 30% can be considered TOF
### ask remi how the aggreggation has been done
system(sprintf("gdal_sieve.py -st %s %s %s ",
               mmu,
               paste0(dd_dir,"tmp_gfc_tc_start.tif"),
               paste0(dd_dir,"tmp_gfc_tc_start_fsieve.tif")
))

#################### FIX THE HOLES
### ??
system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(dd_dir,"tmp_gfc_tc_start.tif"),
               paste0(dd_dir,"tmp_gfc_tc_start_fsieve.tif"),
               paste0(dd_dir,"tmp_gfc_tc_start_sieve.tif"),
               paste0("(A>0)*(B>0)*B")
))

#################### DIFFERENCE BETWEEN SIEVED AND ORIGINAL
### to separate the smaller polygons (than 1 ha) that contain at least 30%??
system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(dd_dir,"tmp_gfc_tc_start.tif"),
               paste0(dd_dir,"tmp_gfc_tc_start_sieve.tif"),
               paste0(dd_dir,"tmp_gfc_tc_start_inf.tif"),
               paste0("(A>0)*(A-B)+(A==0)*(B==1)*0")
))


#################### CREATE GFC LOSS MAP AT THRESHOLD between 2000 and 2015
### all the losses from 2000 to 2015 (Uganda FREL)
system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(gfc_dir,"gfc_treecover2000.tif"),
               paste0(gfc_dir,"gfc_lossyear.tif"),
               paste0(dd_dir,"tmp_gfc_loss.tif"),
               paste0("(A>",gfc_threshold,")*(B>0)*(B<15)")
))

#################### SIEVE TO THE MMU
### to separate the polygons of losses minimum or bigger than 1 ha and smaller than 1 ha (the potential losses to consider and combine with the TOF identified in the previous step for identify degradation in TOF)
system(sprintf("gdal_sieve.py -st %s %s %s ",
               mmu,
               paste0(dd_dir,"tmp_gfc_loss.tif"),
               paste0(dd_dir,"tmp_gfc_loss_fsieve.tif")
))

#################### FIX THE HOLES
### ??
system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(dd_dir,"tmp_gfc_loss.tif"),
               paste0(dd_dir,"tmp_gfc_loss_fsieve.tif"),
               paste0(dd_dir,"tmp_gfc_loss_sieve.tif"),
               paste0("(A>0)*(B>0)*B")
))

#################### DIFFERENCE BETWEEN SIEVED AND ORIGINAL
### to separate the losses on polygons bigger than 1 ha from polygons smaller than 1 ha??
system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(dd_dir,"tmp_gfc_loss.tif"),
               paste0(dd_dir,"tmp_gfc_loss_sieve.tif"),
               paste0(dd_dir,"tmp_gfc_loss_inf.tif"),
               paste0("(A>0)*(A-B)+(A==0)*(B==1)*0")
))


#################### CREATE GFC TREE COVER MASK IN 2015 AT THRESHOLD
### to obtain forests in 2000 (from tmp_gfc_tc_start or A>0) and no losses (B==0) and losses after 2015 (B>15)
system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(dd_dir,"tmp_gfc_tc_start.tif"),
               paste0(gfc_dir,"gfc_lossyear.tif"),
               paste0(dd_dir,"tmp_gfc_tc_end.tif"),
               paste0("(A>0)*((B>=15)+(B==0))")
))


#################### SIEVE TO THE MMU
### minimum area
system(sprintf("gdal_sieve.py -st %s %s %s ",
               mmu,
               paste0(dd_dir,"tmp_gfc_tc_end.tif"),
               paste0(dd_dir,"tmp_gfc_tc_end_fsieve.tif")
))

#################### FIX THE HOLES
###??
system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(dd_dir,"tmp_gfc_tc_end.tif"),
               paste0(dd_dir,"tmp_gfc_tc_end_fsieve.tif"),
               paste0(dd_dir,"tmp_gfc_tc_end_sieve.tif"),
               paste0("(A>0)*(B>0)*B")
))

#################### DIFFERENCE BETWEEN SIEVED AND ORIGINAL
### to identify forest in 2015 bigger or smaller than 1 ha?? is it not the same than the first step?
system(sprintf("gdal_calc.py -A %s -B %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(dd_dir,"tmp_gfc_tc_end.tif"),
               paste0(dd_dir,"tmp_gfc_tc_end_sieve.tif"),
               paste0(dd_dir,"tmp_gfc_tc_end_inf.tif"),
               paste0("(A>0)*(A-B)+(A==0)*(B==1)*0")
))

#################### COMBINATION INTO DD MAP (1==Forest, 2==NonForest, 3==gain IT DOES NOT EXIST, 4==deforestation, 5==Degradation 6==ToF, 7==Dg_TOF)
system(sprintf("gdal_calc.py -A %s -B %s -C %s -D %s -E %s -F %s -G %s --co COMPRESS=LZW --outfile=%s --calc=\"%s\"",
               paste0(dd_dir,"tmp_gfc_tc_start.tif"),
               paste0(dd_dir,"tmp_gfc_loss_sieve.tif"),
               paste0(dd_dir,"tmp_gfc_loss_inf.tif"),
               paste0(dd_dir,"tmp_gfc_tc_end_sieve.tif"),
               paste0(dd_dir,"tmp_gfc_tc_end_inf.tif"),
               paste0(dd_dir,"tmp_gfc_tc_start_fsieve.tif"),
               paste0(dd_dir,"tmp_gfc_tc_start_inf.tif"),
               paste0(dd_dir,"tmp_dd_map.tif"),
               paste0("(A==0)*2+",
                      "(A>0)*((B==0)*(C==0)*((D>0)*1+(E>0)*6)+",
                             "(B>0)*4+",
                             "(C>0)*((F>0)*5+(G>0)*7))")
))

#############################################################
### CROP TO COUNTRY BOUNDARIES
system(sprintf("python %s/oft-cutline_crop.py -v %s -i %s -o %s -a %s",
               scriptdir,
               aoi,
               paste0(dd_dir,"tmp_dd_map.tif"),
               paste0(dd_dir,"tmp_dd_map_aoi.tif"),
               aoi_field
))

#################### CREATE A COLOR TABLE FOR THE OUTPUT MAP
my_classes <- c(0,2,1,4,5,6,7)
my_colors  <- col2rgb(c("black","grey","darkgreen","red","orange","lightgreen","yellow"))

pct <- data.frame(cbind(my_classes,
                        my_colors[1,],
                        my_colors[2,],
                        my_colors[3,]))

write.table(pct,paste0(dd_dir,"color_table.txt"),row.names = F,col.names = F,quote = F)




################################################################################
#################### Add pseudo color table to result
################################################################################
system(sprintf("(echo %s) | oft-addpct.py %s %s",
               paste0(dd_dir,"color_table.txt"),
               paste0(dd_dir,"tmp_dd_map_aoi.tif"),
               paste0(dd_dir,"tmp_dd_map_aoi_pct.tif")
))

################################################################################
#################### COMPRESS
################################################################################
system(sprintf("gdal_translate -ot Byte -co COMPRESS=LZW %s %s",
               paste0(dd_dir,"tmp_dd_map_aoi_pct.tif"),
               paste0(dd_dir,"tmp_dd_map_aoi_pct_geo.tif")
))

################################################################################
#################### REPROJECT in UTM 32636 UGANDAN UTM COORDINATE SYSTEM
################################################################################
system(sprintf("gdalwarp -t_srs EPSG:32636 -overwrite -co COMPRESS=LZW %s %s",
               paste0(dd_dir,"tmp_dd_map_aoi_pct_geo.tif"),
               paste0(dd_dir,"dd_map_0015_gt",gfc_threshold,"_20181226.tif")
))

################################################################################
#################### COMPUTE AREAS
################################################################################
system(sprintf("oft-stat %s %s %s",
               paste0(dd_dir,"dd_map_0015_gt",gfc_threshold,"_20181226.tif"),
               paste0(dd_dir,"dd_map_0015_gt",gfc_threshold,"_20181226.tif"),
               paste0(dd_dir,"stats.txt")
))

df <- read.table(paste0(dd_dir,"stats.txt"))[,1:2]
names(df) <- c("class","pixels")
res <- res(raster(paste0(dd_dir,"dd_map_0015_gt",gfc_threshold,"_20181226.tif")))[1]

df$area <- df$pixels * res * res /10000
df
write.csv(df,paste0(dd_dir,"map_areas.csv"),row.names = F)

system(sprintf("rm %s",
               paste0(dd_dir,"tmp*.tif")
))

Sys.time() - time_start

