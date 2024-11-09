library(terra)
library(raster)
library(sf)
library(mapview)
library(dplyr)
library(purrr)
library(landscapemetrics)


# Read data and create buffers around cameras -----------------------------

cams <- st_read("cams2021_Lambert.shp") # cameras
cams <- st_read("cams2017_Lambert.shp")

lc <- rast("OH_nlcd_Lambert.tif") # land cover NLCD data
plot(lc)

# create buffers around cameras
buf_500 <- st_buffer(cams, dist=1640) # 500 m
buf_1000 <- st_buffer(cams, dist=3280) # 1000 m
buf_1500 <- st_buffer(cams, dist=4920) # 1500 m

camID <- buf_500$Camera
camID <- buf_500$SiteCam # Heidi data

mapview::mapview(buf_1500)


# Extract land cover types from NLCD (as number of cells within a  --------

lc_500 <- extract(x=lc, y=buf_500,  fun='table')
str(lc_500)

lc_500 <- lc_500 |>
  mutate(camID = camID) |>
  rename("nodata" = "0",
             "OpenWater500" = "1",
              "DevOpenSpace500" = "2",
              "DevLowInt500" = "3",
              "DevMedInt500" = "4",
              "DevHighInt500" = "5",
              "Barren500" = "6",
              "DecForest500" = "7",
              "ConForest500" = "8",
              "mixedForest500" = "9",
              "Shrub500" = "10",
              "GrassHerb500" = "11",
              "Pasture500" = "12",
              "Crops500" = "13",
              "WoodyWetland500" = "14",
              "EmergWetland500" = "15")



lc_1000 <- extract(x=lc, y=buf_1000,  fun='table')
str(lc_1000)

lc_1000 <- lc_1000 |>
  mutate(camID = camID) |>
  rename("nodata" = "0",
         "OpenWater1000" = "1",
         "DevOpenSpace1000" = "2",
         "DevLowInt1000" = "3",
         "DevMedInt1000" = "4",
         "DevHighInt1000" = "5",
         "Barren1000" = "6",
         "DecForest1000" = "7",
         "ConForest1000" = "8",
         "mixedForest1000" = "9",
         "Shrub1000" = "10",
         "GrassHerb1000" = "11",
         "Pasture1000" = "12",
         "Crops1000" = "13",
         "WoodyWetland1000" = "14",
         "EmergWetland1000" = "15")



lc_1500 <- extract(x=lc, y=buf_1500,  fun='table')
str(lc_1500)

lc_1500 <- lc_1500 |>
  mutate(camID = camID) |>
  rename("nodata" = "0",
         "OpenWater1500" = "1",
         "DevOpenSpace1500" = "2",
         "DevLowInt1500" = "3",
         "DevMedInt1500" = "4",
         "DevHighInt1500" = "5",
         "Barren1500" = "6",
         "DecForest1500" = "7",
         "ConForest1500" = "8",
         "mixedForest1500" = "9",
         "Shrub1500" = "10",
         "GrassHerb1500" = "11",
         "Pasture1500" = "12",
         "Crops1500" = "13",
         "WoodyWetland1500" = "14",
         "EmergWetland1500" = "15")


# combine land cover for all scales into a single object
all_lc <- purrr::reduce(list(lc_1500, lc_1000, lc_500), 
                         dplyr::left_join, by = 'camID')
str(all_lc)



# Calculate and extract Shannon Diversity Index for habitat hetero --------

# THIS TAKES A LONG TIME
# create a larger buffer round cams and use it to calculate Shannon Div Index

buf_10000 <- st_buffer(cams, dist=10000) # 3000 m
buf_10000_1 <- vect(buf_10000)
small_lc <- crop(lc, buf_10000_1, mask=T)
plot(small_lc)

mapview(small_lc) +
  mapview(cams)

# create moving 500 x 500 m moving window and calculate SHDI
window500 <- matrix(1, 17, 17) # cell is 30 x 30 m
shannon500 <- window_lsm(small_lc, window = window500, what = "lsm_l_shdi")
str(shannon500)
plot(shannon500$layer_1$lsm_l_shdi)
writeRaster(shannon500$layer_1$lsm_l_shdi, "ShannonDiv500_heidi.tif")

# extract Shannon Div Index for each camera location
SHDI_500 <- extract(x=shannon500$layer_1$lsm_l_shdi, y=cams)
SHDI_500 <- SHDI_500 |>
  mutate(camID = camID) |>
  rename("SHDI_500" = "OID")
str(SHDI_500)

# create moving 1000 x 1000 m moving window and calculate SHDI
window1000 <- matrix(1, 33, 33) # cell is 30 x 30 m
shannon1000 <- window_lsm(small_lc, window = window1000, what = "lsm_l_shdi")
plot(shannon1000$layer_1$lsm_l_shdi)
writeRaster(shannon1000$layer_1$lsm_l_shdi, "ShannonDiv1000_heidi.tif")

# extract Shannon Div Index for each camera location
SHDI_1000 <- extract(x=shannon1000$layer_1$lsm_l_shdi, y=cams)
SHDI_1000 <- SHDI_1000 |>
  mutate(camID = camID) |>
  rename("SHDI_1000" = "OID")
str(SHDI_1000)


window1500 <- matrix(1, 51, 51) # cell is 30 x 30 m
shannon1500 <- window_lsm(small_lc, window = window1500, what = "lsm_l_shdi")
plot(shannon1500$layer_1$lsm_l_shdi)
writeRaster(shannon1500$layer_1$lsm_l_shdi, "ShannonDiv1500_heidi.tif")

# extract Shannon Div Index for each camera location
SHDI_1500 <- extract(x=shannon1500$layer_1$lsm_l_shdi, y=cams)
SHDI_1500 <- SHDI_1500 |>
  mutate(camID = camID) |>
  rename("SHDI_1500" = "OID")
str(SHDI_1500)

# combine all SHDI in an object
all_SHDI <- purrr::reduce(list(SHDI_500, SHDI_1000, SHDI_1500), 
                        dplyr::left_join, by = 'camID')
str(all_SHDI)


# Extract distance to nearest main road and distance to streams -----------

# we already had these distance rasters calculated from a previous project

# distance to main roads
distroads <- rast("distmainrds.tif")
dist.to.rds <- extract(distroads, cams)
dist.to.rds$camID <- camID
str(dist.to.rds)

# distance to streams
diststreams <- rast("DistanceStreams.tif")
dist.to.streams <- extract(diststreams, cams)
dist.to.streams$camID <- camID
str(dist.to.streams)



# Calculate and extract road density (length within a given buffer) -------

# read in layer of main (high traffic) roads
main_roads <- st_read("roads_IR_US_SR.shp")

# intersect the road layer with the 500 m buffer and calculate the length of road segments
ints500 <- st_intersection(main_roads, buf_500)
str(ints500)
rds500 <- tapply(st_length(ints500), ints500$Camera, sum) |> #Site
  as.data.frame() |>
  rename("MainRds500" = "tapply(st_length(ints500), ints500$Camera, sum)")
str(rds500)
rds500$camID <- rownames(rds500)

# intersect the road layer with the 1000 m buffer and calculate the length of road segments
ints1000 <- st_intersection(main_roads, buf_1000)
rds1000<- as.data.frame(tapply(st_length(ints1000), ints1000$Camera, sum)) |>
  as.data.frame() |>
  rename("MainRds1000" = "tapply(st_length(ints1000), ints1000$Camera, sum)")
str(rds1000)
rds1000$camID <- rownames(rds1000)

# intersect the road layer with the 1500 m buffer and calculate the length of road segments
ints1500 <- st_intersection(main_roads, buf_1500)
rds1500 <- as.data.frame(tapply(st_length(ints1500), ints1500$Camera, sum)) |>
  as.data.frame() |>
  rename("MainRds1500" = "tapply(st_length(ints1500), ints1500$Camera, sum)")
str(rds1500)
rds1500$camID <- rownames(rds1500)

# combine main roads in a single object
main_rds <- purrr::reduce(list(rds1500, rds1000, rds500), 
                          dplyr::left_join, by = 'camID')


# read in layer of ALL roads
all_roads <- st_read("roads_all.shp")

# intersect the road layer with the 500 m buffer and calculate the length of road segments
ints500 <- st_intersection(all_roads, buf_500)
str(ints500)
rds500 <- tapply(st_length(ints500), ints500$Camera, sum) |>
  as.data.frame() |>
  rename("AllRds500" = "tapply(st_length(ints500), ints500$Camera, sum)")
str(rds500)
rds500$camID <- rownames(rds500)

# intersect the road layer with the 1000 m buffer and calculate the length of road segments
ints1000 <- st_intersection(all_roads, buf_1000)
rds1000<- as.data.frame(tapply(st_length(ints1000), ints1000$Camera, sum)) |>
  as.data.frame() |>
  rename("AllRds1000" = "tapply(st_length(ints1000), ints1000$Camera, sum)")
str(rds1000)
rds1000$camID <- rownames(rds1000)

# intersect the road layer with the 1500 m buffer and calculate the length of road segments
ints1500 <- st_intersection(all_roads, buf_1500)
rds1500 <- as.data.frame(tapply(st_length(ints1500), ints1500$Camera, sum)) |>
  as.data.frame() |>
  rename("AllRds1500" = "tapply(st_length(ints1500), ints1500$Camera, sum)")
str(rds1500)
rds1500$camID <- rownames(rds1500)

# combine all the road information into a single file
all_rds <- purrr::reduce(list(rds1500, rds1000, rds500), 
                          dplyr::left_join, by = 'camID')

roads <- left_join(all_rds, main_rds, by = 'camID')
str(roads)



# Combine ALL data (land cover, Shannon, roads, proximity to roads --------

# ADD Shannon Div Index when ready

all_data <- purrr::reduce(list(all_lc, all_SHDI, roads, dist.to.rds, dist.to.streams), 
                         dplyr::left_join, by = 'camID')

write.csv(all_data, "cam2021_envdata.csv")
