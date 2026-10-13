##Code for Workshop##
##Author: Jon Rodemann
##Date started: 10/12/2026

####YAPS####

#packages
install.packages(c('tidyverse', 'devtools', 'TMB'))
library(devtools)
library(tidyverse)
library(TMB)
devtools::install_github("baktoft/yaps")
library(yaps)
options(digits.secs = 3)
#test to make sure YAPS will work
testYaps()

#load data
vue <- data.table::fread('data/detections.csv', fill=TRUE)
hydros <- data.table::fread('data/hydros.csv', fill=TRUE)
vue$`Date and Time (UTC)` <- as.character(vue$`Date and Time (UTC)`)
detections <- prepDetections(raw_dat = vue, type='vemco_vue')
detections <- detections %>% dplyr::filter(!is.na(frac))
#now we have our 2 datasets needed!

#sync hydrophones
#set parameters
max_epo_diff <- 140 #Should be a little less than half of your shortest sync tag delay (ours are 6 - 10 minute)
min_hydros <- 4 #minimum amount of hydrophones needed for each sync tag ping to be heard on. Start high, then go low
time_keeper_idx <- 13 #Receiver with the "best clock"
fixed_hydros_idx <- c(3, 7:8, 11, 13:17, 20:21) #Receivers that have the best known positions. The more, the better
n_offset_day <- 2 #1 or 2
n_ss_day <- 2 #number of times you want the model to calculate speed of sound a day. 2 is good, best to use water temp if you have it
keep_rate <- 75 #number of sync tag detections for each hydro X offset_idx combinations to keep. Keeps the model doable

#compile list of dataframes
RYAPS <- list(hydros, detections)
names(RYAPS) <- c('hydros', 'detections')
RYAPS$hydros$sync_tag <- as.double(RYAPS$hydros$sync_tag)
RYAPS$detections$tag <- as.integer(RYAPS$detections$tag)
RYAPS$detections$epo <- as.integer(RYAPS$detections$epo)
RYAPS$detections$serial <- as.integer(RYAPS$detections$serial)
#compile data for syncing
inp_sync <- getInpSync(sync_dat = RYAPS, max_epo_diff, min_hydros, time_keeper_idx,
                       fixed_hydros_idx, n_offset_day, n_ss_day, keep_rate)

#check if it is good
getSyncCoverage(inp_sync, plot=T)


#we are good. Let's try syncing
sync_model <- getSyncModel(inp_sync, silent=TRUE)

#check quality of synchronization
plotSyncModelResids(sync_model, by='overall')
plotSyncModelResids(sync_model, by='sync_tag')
plotSyncModelResids(sync_model, by='hydro')

plotSyncModelCheck(sync_model, by="sync_bin_sync")
plotSyncModelCheck(sync_model, by="sync_bin_hydro")
plotSyncModelCheck(sync_model, by="sync_tag")
plotSyncModelCheck(sync_model, by="hydro")

#success! there is only a problem of 6 using SHYAPS1. I think it is fine. save it
save(sync_model, file = 'sync_model_4')

#apply sync model to detection data
detections_synced <- applySync(toa=SHYAPS$detections, hydros=SHYAPS$hydros, sync_model)

#Now that you have the sync model, run the tracks
#use the estimated positions of hydros
hydros_yaps <- data.table::data.table(sync_model$pl$TRUE_H)
colnames(hydros_yaps) <- c('hx','hy','hz')

#specify focal tag
focal_tag <- 34507
rbi_min <- 20
rbi_max <- 40

#extract relevant tag data from detections
synced_dat_34507 <- detections_synced[tag == focal_tag]


#compile TOA matrix
toa_YAPS_34507 <- getToaYaps(synced_dat_34507, hydros_yaps, rbi_min, rbi_max, pingType='rbi')

hydros_yaps

#compile all input data needed
inp_SHYAPS <- getInp(hydros_yaps, toa_YAPS_34507, E_dist="Mixture", n_ss=3, pingType="rbi",
                     sdInits=1, rbi_min=rbi_min, rbi_max=rbi_max, ss_data_what="est", ss_data=0)



inp_SHYAPS$datTmb$rbi_min
inp_SHYAPS$params$top
min(diff(inp_SHYAPS$params$top))

checkInp(inp_SHYAPS)


#run YAPS!
yaps_out_34507 <- runYaps(inp_SHYAPS, silent=F)

plotYaps(yaps_out_34507)

###So now we have our positions. Time to filter and set them up for the HMMs