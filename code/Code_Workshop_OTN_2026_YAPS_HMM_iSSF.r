##Code for Workshop##
##Authors: Jon Rodemann, Sophia Costa, Jennifer Rehage
##Source Data - Rodemann et al., 2025 Sci. Rep., Costa et al., 2026 (In Prep)

####YAPS####

##This is the general workflow for syncing receivers and running tracks in YAPS. We won't be running this code in this workshop
##and in fact, this code does not work because the syncing fails. In reality, you need to play around with your sync data 
##(including splitting up the time periods if there are "blackouts", possibly removing receivers, and other adjustments) for the sync model
##to be successful. We will pick up the workshop with the positions that were calculated by the developer of YAPS (WOW!), Henrik Baktoft.

#Install and load required libraries
# install.packages(c('tidyverse', 'devtools', 'TMB'))
# library(devtools)
# library(tidyverse)
# library(TMB)
# devtools::install_github("baktoft/yaps")
# library(yaps)
# options(digits.secs = 3)
# #test to make sure YAPS will work
# testYaps()
# 
# #load data
# vue <- data.table::fread('data/detections.csv', fill=TRUE)
# hydros <- data.table::fread('data/hydros.csv', fill=TRUE)
# vue$`Date and Time (UTC)` <- as.character(vue$`Date and Time (UTC)`)
# detections <- prepDetections(raw_dat = vue, type='vemco_vue')
# detections <- detections %>% dplyr::filter(!is.na(frac))
# #now we have our 2 datasets needed!
# 
# #sync hydrophones
# #set parameters
# max_epo_diff <- 140 #Should be a little less than half of your shortest sync tag delay (ours are 6 - 10 minute)
# min_hydros <- 4 #minimum amount of hydrophones needed for each sync tag ping to be heard on. Start high, then go low
# time_keeper_idx <- 13 #Receiver with the "best clock"
# fixed_hydros_idx <- c(3, 7:8, 11, 13:17, 20:21) #Receivers that have the best known positions. The more, the better
# n_offset_day <- 2 #1 or 2
# n_ss_day <- 2 #number of times you want the model to calculate speed of sound a day. 2 is good, best to use water temp if you have it
# keep_rate <- 75 #number of sync tag detections for each hydro X offset_idx combinations to keep. Keeps the model doable
# 
# #compile list of dataframes
# RYAPS <- list(hydros, detections)
# names(RYAPS) <- c('hydros', 'detections')
# RYAPS$hydros$sync_tag <- as.double(RYAPS$hydros$sync_tag)
# RYAPS$detections$tag <- as.integer(RYAPS$detections$tag)
# RYAPS$detections$epo <- as.integer(RYAPS$detections$epo)
# RYAPS$detections$serial <- as.integer(RYAPS$detections$serial)
# #compile data for syncing
# inp_sync <- getInpSync(sync_dat = RYAPS, max_epo_diff, min_hydros, time_keeper_idx,
#                        fixed_hydros_idx, n_offset_day, n_ss_day, keep_rate)
# 
# #check if it is good
# getSyncCoverage(inp_sync, plot=T)
# 
# 
# #we are good. Let's try syncing
# sync_model <- getSyncModel(inp_sync, silent=TRUE)
# 
# #check quality of synchronization
# plotSyncModelResids(sync_model, by='overall')
# plotSyncModelResids(sync_model, by='sync_tag')
# plotSyncModelResids(sync_model, by='hydro')
# 
# plotSyncModelCheck(sync_model, by="sync_bin_sync")
# plotSyncModelCheck(sync_model, by="sync_bin_hydro")
# plotSyncModelCheck(sync_model, by="sync_tag")
# plotSyncModelCheck(sync_model, by="hydro")
# 
# #success! there is only a problem of 6 using SHYAPS1. I think it is fine. save it
# save(sync_model, file = 'sync_model_4')
# 
# #apply sync model to detection data
# detections_synced <- applySync(toa=SHYAPS$detections, hydros=SHYAPS$hydros, sync_model)
# 
# #Now that you have the sync model, run the tracks
# #use the estimated positions of hydros
# hydros_yaps <- data.table::data.table(sync_model$pl$TRUE_H)
# colnames(hydros_yaps) <- c('hx','hy','hz')
# 
# #specify focal tag
# focal_tag <- 34507
# rbi_min <- 20
# rbi_max <- 40
# 
# #extract relevant tag data from detections
# synced_dat_34507 <- detections_synced[tag == focal_tag]
# 
# 
# #compile TOA matrix
# toa_YAPS_34507 <- getToaYaps(synced_dat_34507, hydros_yaps, rbi_min, rbi_max, pingType='rbi')
# 
# hydros_yaps
# 
# #compile all input data needed
# inp_SHYAPS <- getInp(hydros_yaps, toa_YAPS_34507, E_dist="Mixture", n_ss=3, pingType="rbi",
#                      sdInits=1, rbi_min=rbi_min, rbi_max=rbi_max, ss_data_what="est", ss_data=0)
# 
# 
# 
# inp_SHYAPS$datTmb$rbi_min
# inp_SHYAPS$params$top
# min(diff(inp_SHYAPS$params$top))
# 
# checkInp(inp_SHYAPS)
# 
# 
# #run YAPS!
# yaps_out_34507 <- runYaps(inp_SHYAPS, silent=F)
# 
# plotYaps(yaps_out_34507)

##This is where we will leave it for the YAPS run. We encourage all of you if you are interested to explore the package, as 
##it can even run simulations! Now onto the data exploration and filtering

####Data filtering, exploration and graphing####
##In this section, we will look at the data, filter it based on the positional error from the YAPS algorithm, explore the data, and then
##graph it as an abacus plot (yay!) as well as on a map. Then we will prep it for the HMM analysis.

#Load in libraries
library(splitstackshape)
library(chron)
library(scales)
library(ggmap)
library(sf)

data <- sfheadersdata <- read_csv('data/fish_positions.csv')
head(data)

#The columns of this data include tag (last digits of tag number), x (longitude in UTM), y (latitude in UTM), 
#x_sd (standard deviation (error) of longitude in meters), y_sd (standard deviation (error) of latitude in meters),
#z_sd (standard deviation (error) of depth in meters - not applicable), nobs (number of receivers detected on), and k (positions that represent a track)

#make a column with total error
data <- data %>%
  mutate(E = sqrt(x_sd^2 + y_sd^2))

#filter by error. We are going to filter out anything with an error greater than 10 m. This will have to match your environmental variables later
dat_10 <- data %>% filter(E < 10)

#summary stats
#number of positions per fish
dat_10 %>% group_by(tag) %>% summarize(n = n())

#abacus plot for positions
dat_10$ts <- as.POSIXct(dat_10$ts, format = '%Y-%m-%dT%H:%M:%OS')

# Break up the new datetime column into  a date, and a time
data <- dat_10

data$Date_Time2 <- data$ts

data <- cSplit(data, "Date_Time2", sep = " ", type.convert = F)

head(data)


# Rename columns so they make sense, it'll just be your last 2 column numbers, in this case the 15th and 16th column s
colnames(data)[12:14]<-c('E', "Date", "Time")
head(data)


# Then I repeat this and parse the date into year, month, and day (then hour, minute, second), so I can easily subset by year
# Copy date to parse out, maintain original date
data$Date2 <- data$Date


data<-cSplit(data, "Date2", sep = "-", type.convert = FALSE)
head(data)
colnames(data)[15:17]<-c("Year", "Month", "Day")
head(data)

data$Time2 <- data$Time

data <-cSplit(data, "Time2", sep = ":", type.convert = FALSE)
head(data)
colnames(data)[18:20]<-c("Hour", "Minute", "Second")
head(data)

# Assign classes so data plays nice for plotting and summaries
data$Date <- as.Date(data$Date)
data$Time <- as.times(data$Time)
data$Year <- as.numeric(data$Year)
data$Month <- as.numeric(data$Month)
data$Day <- as.numeric(data$Day)
data$Hour <- as.numeric(data$Hour)
data$Minute <- as.numeric(data$Minute)
data$Second <- as.numeric(data$Second)
data$tag <- as.factor(data$tag)


data$ts <- as.POSIXct(data$ts, format='%Y-%m-%d %H:%M:%S', tz='America/New_York')

head(data)

fish_det <- ggplot(data,aes(ts, tag)) + 
  geom_point() +
  scale_x_datetime(labels = date_format("%Y-%m-%d"),
                   date_breaks = "2 weeks", limits = as.POSIXct(c('2021-12-15 00:00:00', '2022-05-10 00:00:00')))+
  labs(x = "Date",  y = "Transmitter") +
  theme_bw()+
  theme(axis.text.x=element_text(angle= 50, size=10, vjust = 1, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y=element_text(size=8))+
  labs(title = "Trout Detections YAPS")
fish_det

#We need to make the data into a spatial dataframe to graph onto a map
sfdat <- data %>% 
  st_as_sf(coords = c('x', 'y')) %>% 
  st_set_crs(2958) %>% 
  st_transform(4326) %>% 
  mutate(
    long = st_coordinates(.)[, 1],
    lat  = st_coordinates(.)[, 2])

#Fok, let's plot this in ggmap to get our bearings
# Note,a API key from Google Maps is required to load the map and run the below code. Visit https://rpubs.com/ktcross/1154003 for instructions on  how to install and run ggamps in R
# For now, we will use the already saved maps but the code to generate the maps if you have an existing API code is below as well.

# FLmap_zoom_out <- get_googlemap(center = c(lon=mean(sfdat$long), lat=mean(sfdat$lat)),
#                         zoom = 6,
#                         maptype = c("satellite"))

load("data/FLmap_zoom_out.RData")

ggmap(FLmap_zoom_out, extent='normal') + 
  geom_point(data = sfdat,
             aes(x = long, y = lat), col ="yellow") 

# FLmap_zoom_in <- get_googlemap(center = c(lon=mean(sfdat$long), lat=mean(sfdat$lat)),
#                                zoom = 16,
#                                maptype = c("satellite"))
load("data/FLmap_zoom_in.RData")

ggmap(FLmap_zoom_in, extent='normal') + 
  geom_point(data = sfdat,
             aes(x = long, y = lat,), col ="yellow") +
  facet_wrap(~tag)

#So now we know where we are in space. Let's now set up the tracks for the HMM analysis
#We are going to filter it even more so that each position we keep is calculated by detections heard by 2 or more receivers
dat20 <- data %>% dplyr::filter(nobs >= 2)

#Now we are going to create a column that allows us to identify the start of a new track
#for this exercise, we decided that if the fish was not detected for 5 minutes, a new track would be created
dat20t <- dat20 %>% group_by(tag) %>% mutate(trackstart = ifelse(as.numeric(unclass(ts)-unclass(lag(ts))) > 300, 'Y', 'N')) %>% ungroup()
dat20t[is.na(dat20t)] <- 'Y'

str(dat20t)

#create a column that we can put track numbers in
dat20tt <- dat20t %>% group_by(tag) %>% mutate(track = ifelse(row_number() == 1, 1, NA)) %>% ungroup()

#create a function that numbers the tracks and inputs what track each row is in
tracknum1 <- function(x){
  for (i in 2:nrow(x)){
    if (x$trackstart[i] == 'Y'){
      x$track[i] <- x$track[i-1]+1
    } else {
      x$track[i] <- x$track[i-1]
    }
  }
  x
}

#apply across the data
dat_HMM <- dat20tt %>% group_by(tag) %>% group_modify(~tracknum1(.x)) %>% ungroup()

write.csv(dat_HMM, file = 'data/dat_HMM.csv')

####HMMs#####

#Load in libraries
library(tidyverse)
library(lubridate)
library(moveHMM)
library(purrr)
library(readr)
library(parallel)
library(circular)

dat <- read_csv('data/dat_HMM.csv')

# Parse timestamp
dat <- dat %>%
  mutate(ts = as.POSIXct(ts,
                         format = "%Y-%m-%d %H:%M:%OS",
                         tz = "UTC"))

# The model expects at least these columns:
#   tag          : individual identifier
#   track        : track number/ID for that individual's track segment
#   ts           : timestamp string 
#   x, y         : coordinates (assumed UTM/metres). If they're lat/long, convert to projected meters first.

glimpse(dat)

ts_col <- "ts"
id_col <- "tag"
track_col <- "track"
x_col <- "x"
y_col <- "y"

##Hidden Markov Models require the data to have a regular time step. Because these positions are calculated from random burst
##tags, we need to interpolate the tracks to a regular interval. This is the third big decision - how do you parse out the time
##interval to interpolate the tracks? Our suggestion is to start at twice the average burst interval of the tag to account for 
##possible missed detections.

#Create a function that interpolates the tracks from the positions. For this function to work, you need to group the data before running it into unique ID X track combinations

interpolate_track <- function(df_group, dt_seconds = 180) { #The default time of the interpolated tracks is 3 minutes, but you can change that
  # df must contain columns: ts, x, y, plus id/track identifiers
  if (nrow(df_group) < 2) return(NULL)  # can't interpolate single-point tracks
  
  tmin <- min(df_group$ts, na.rm = TRUE)
  tmax <- max(df_group$ts, na.rm = TRUE)
  
  new_times <- seq(from = tmin, to = tmax, by = dt_seconds) #sequence the track out to regular intervals based on your start and end times and the interval you identify in the function
  if (length(new_times) < 2) return(NULL)
  
  # numeric times for approx - makes math easier
  orig_t_num <- as.numeric(df_group$ts)
  new_t_num <- as.numeric(new_times)
  
  # ensure sorted by time - earliest to latest
  ord <- order(orig_t_num)
  orig_t_num <- orig_t_num[ord]
  x_orig <- df_group[[x_col]][ord]
  y_orig <- df_group[[y_col]][ord]
  
  # linear interpolation (NA handling: approx requires finite values) - takes actual points and times, then interpolates the point between the actual points based on actual times
  ok_x <- is.finite(x_orig) & is.finite(orig_t_num)
  ok_y <- is.finite(y_orig) & is.finite(orig_t_num)
  if (!any(ok_x) || !any(ok_y)) return(NULL)
  
  x_interp <- approx(x = orig_t_num[ok_x], y = x_orig[ok_x], xout = new_t_num, rule = 2)$y #The y extracts only the interpolated points
  y_interp <- approx(x = orig_t_num[ok_y], y = y_orig[ok_y], xout = new_t_num, rule = 2)$y

  tibble(
    ts = new_times, #interpolated times
    !!id_col := df_group[[id_col]][1],
    !!track_col := df_group[[track_col]][1],
    x = x_interp,
    y = y_interp
  )
}

#Now let's run it on our data!
interp_list <- dat %>%
  group_by(.data[[id_col]], .data[[track_col]]) %>%
  group_map(~ interpolate_track(.x, dt_seconds = 60), .keep = TRUE)

# Remove NULLs (tracks that couldn't be interpolated)
interp_list <- keep(interp_list, ~ !is.null(.x) && nrow(.x) > 0)

interp_df <- bind_rows(interp_list) #make dataframe
interp_df <- interp_df %>% mutate(ID = paste(tag, track, sep = "_")) #make ID column

##We have our interpolated tracks! As you can see, we didn't actually lose a lot of data (went from 8436 to 8254 observations).
##Just to check that these interpolated tracks are similar to the original, let's do a comparison on a few of them! real will
##be dark green, interpolated purple

#34501_2
int_34501_2 <- interp_df %>% dplyr::filter(ID == '34501_2')
real_34501_2 <- dat %>% dplyr::filter(tag == '34501') %>% dplyr::filter(track == 2)

ggplot()+
  geom_point(data = int_34501_2, aes(x = x, y = y), color = 'purple')+
  geom_point(data = real_34501_2, aes(x = x, y = y), color = 'darkgreen', inherit.aes = F)+
  theme_classic()

#34507_31
int_34507_31 <- interp_df %>% dplyr::filter(ID == '34507_31')
real_34507_31 <- dat %>% dplyr::filter(tag == '34507') %>% dplyr::filter(track == 31)

ggplot()+
  geom_point(data = int_34507_31, aes(x = x, y = y), color = 'purple')+
  geom_point(data = real_34507_31, aes(x = x, y = y), color = 'darkgreen', inherit.aes = F)+
  theme_classic()

#34518_12
int_34518_12 <- interp_df %>% dplyr::filter(ID == '34518_12')
real_34518_12 <- dat %>% dplyr::filter(tag == '34518') %>% dplyr::filter(track == 12)

ggplot()+
  geom_point(data = int_34518_12, aes(x = x, y = y), color = 'purple')+
  geom_point(data = real_34518_12, aes(x = x, y = y), color = 'darkgreen', inherit.aes = F)+
  theme_classic()

##They match pretty well! Now comes the fourth big decision. What tracks are worth keeping in the array? Is a 5 minute
##track worth it? Probably not because it will only include 5 positions. But, you don't want to get rid of too much data.
##For this study, we decided to cut it off at 30 minutes or longer in the array for each track. This represents a good amount
##of time where we might see two types of behaviors. So let's filter some more!

id20 <- interp_df %>% group_by(ID) %>% mutate(tott = max(ts)-min(ts)) %>% ungroup()
id20$tottime <- as.numeric(id20$tott)
id20_30 <- id20 %>% dplyr::filter(tottime >= 1800)
nrow(as.data.frame(unique(id20_30$ID))) #number of unique tracks kept through all of the filtering - 59

##We are now able to fit some HMMs! To do this, we need to first prep the data and then remove missing data

interp_df_30 <- id20_30 %>% arrange(ID, ts)
# prepData: coordinate names are "x" and "y" in our frame
hmm_data_30 <- prepData(interp_df_30, type = "UTM", coordNames = c("x", "y"))

# moveHMM uses a column 'step' and 'angle' (radians) created by prepData
# Remove rows with zero or NA step (stay/missing movement) because gamma requires positive steps
hmm_data_30 <- hmm_data_30 %>% filter(!is.na(step) & step > 0) %>% filter(!is.na(angle))

##We have prepped our data. However, it is still not just plug and play. We need to investigate our data to fit 
##general curves to our data to inform the model. We don't need to do this, but it creates a more informed model based
##on the original data. We are trying to fit these data into two behavioral states, so we will want to inform the model with 2 curves

hist(hmm_data_30$step, breaks = 20, freq = F, main = "Step lengths")
curve(dgamma(x, shape = 1.5, rate=0.9), n=200, add=TRUE, col="royalblue", lwd=4) #shape is how right-skewed the gamma curve is
curve(dgamma(x, shape = 1.9, rate=0.6), n=200, add=TRUE, col="cadetblue", lwd=4) #rate is how tall the curve apex is
curve(dgamma(x, shape = 2.2, rate=0.4), n=200, add=TRUE, col="navyblue", lwd=4)

##this is a little subjective, but there is a good starting out point. The first curve encompasses the lower step lengths well.
##The third curve does not really do a good job because it is spread out too far. So let's go with the first and second curves

#assign initial step length distributions
sl_init_mean <- c(1.9/0.6, 1.5/0.9) #initial means should be the shape/rate for both
sl_init_sd <- c(sqrt(1.9)/0.6, sqrt(1.5)/0.9) #initial standard deviations should be square root of shape/rate

##Now we do the same thing for the turn angles. However, we want the model to try to find states with opposite turn angles (turn 180 degrees vs straight)
##Therefore, we will play with the numbers a little.
hist(hmm_data_30$angle, breaks = 30, freq = F, main = "Turn angles")
curve(circular::dvonmises(x, mu = pi, kappa = 0.01), n=200, add=TRUE, col="royalblue", lwd=4) #mu is the mean direction of the angle
curve(circular::dvonmises(x, mu = -0.1, kappa = 7), n=200, add=TRUE, col="cadetblue", lwd=4) #kappa is the concentration (how tall it is)

##The first curve is almost a straight line. But we want the model to try to fit 180 degree turn and straight. The second curve is 
##straight ahead and the first is the turn around. It needs to be as close to flat as possible to fit the data the best.

#assign initial angle distributions
ta_init_mean <- c(-0.1, pi) #mean is from mu
ta_init_con <- c(7, 0.01) #concentration is kappa

##Once you have decided on the initial conditions, you can run the HMMs!

mod <- fitHMM(data = hmm_data_30,
               nbStates = 2,
               stepPar0 = c(sl_init_mean, sl_init_sd),
               anglePar0 = c(ta_init_mean, ta_init_con),
               formula = ~1)

mod

mod %>% plotPR() #gives the goodness of fit. Looks pretty good!

##Looking at the step length and angle distributions, you can see a few things. State 1 is categorized by longer step lengths
##(average of 4.13 meters) and straight turn angles (0.007 radians, basically zero). State 2, on the other hand, has a much shorter
##step length (average of 0.9 meters) but a similar turn angle (-0.02 radians). However, the concentration is a lot lower, so state 2
##encompasses the larger turn angles. These results also give you a transistion probability to transition between states. Pretty low for both,
##but it does happen. If you run the following code, you can visualize the distributions for each state as well as which state each step on all the
##tracks is! Since that is 59 tracks, we won't run it online but we will go through it in the workshop

statecol <- c('lightblue', 'purple1')
plot(mod, col = statecol)

##Now we can assign states directly to the dataframe

states <- viterbi(mod)

state_30 <- mod$data %>%
  mutate(state = states)

# Convert angle to degrees for interpretability (angles are in radians)
state_30 <- state_30 %>%
  mutate(angle_deg = angle * 180 / pi,
         abs_angle_deg = abs(angle_deg))

output_HMM_df <- state_30 %>%
  dplyr::select(ID, track, tag, ts, x, y, step, angle, angle_deg, state)

##We have our final HMM dataset! Now we can move onto looking at what environmental variables we can use to model against the
##HMMs and for the iSSFs!

####Environmental Variables####

# This section goes over what environmental variables we are using for our iSSFs. In this example, we will be using a two class (dense and sparse/bare) map of 
#submerged aquatic vegetation (SAV) mapped using aerial imagery flown over our study site. The map was made using photointerpreted points to train a Random Forest algorithm on the spectral information within the aerial imagery. 
#We will be using this map to also calculate spatial pattern metrics using the dense class of SAV, which are measures of composition and configuration of the seascape. 
#Included in this section will also be the setup of the data for iSSFs, including calculating the "false" steps that the iSSFs use to determine habitat selection.

# We will end up with 6 environmental variables: class of SAV, PLAND (Percentage of Landscape; the percentage that the dense SAV class makes up), ED (Edge Density; the density of edges of the dense SAV class), 
#PD (Patch Density; the density of patches of the dense SAV class), Division (a measure of likelihood that two randomly selected cells are not the same patch), and Shape_mn (Mean Shape; a measure of how complex the patches of dense SAV are). 
#This is where you can substitute any environmental data you may have based on your questions.

#Load in libraries
library(amt) #Package that iSSFs are in
library(tidyverse)
library(survival)
library(terra)
library(sf)
library(furrr) #for parallel processing for landscape metrics. Makes it much quicker
library(purrr)
library(ggeffects)
library(landscapemetrics) #package for calculating spatial pattern metrics

rdat <- output_HMM_df
rmap <- rast('maps/rmap_dense.tif')

plot(rmap)

# Make fish tracks with the amt functionality
fish_trackr <- rdat %>%
  mutate(ID = as.factor(ID),
         ts = as.POSIXct(ts, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")) %>% 
  #group_by(ID) %>% 
  make_track(
    .x = x, # X coordinate column
    .y = y, # Y coordinate column
    .t = ts, # Timestamp column
    id = ID, # Individual ID column
    crs = st_crs(rmap),
    state = state) 

# Convert the track to observed steps, grouped by track
# This calculates step length (sl_) and turn angle (ta_) for the observed steps. You have to nest by ID or else 
#it thinks it is one long track
steps_obsr <- fish_trackr %>%
  # Isolate each track into a list-column called 'data'
  nest(data = -id) %>%
  # Apply steps() to each nested track independently
  mutate(
    step_data = map(data, ~ steps(.x, keep_cols = "end"))
  ) %>%
  # Keep only the ID and the new step data, then expand back to a normal dataframe
  select(id, step_data) %>%
  unnest(cols = step_data)

##Now that we have our observed data set up for environmental extraction and iSSF calculation, we need to calculate the random steps
##that the iSSFs need to compare against our observed step. In this workshop, we are just going to do a 10 - 1 ratio, but
##in the actual paper we are using a 100 to 1 ratio. This will lead to lower error in the model, but take A LOT longer.

# Generate random steps (10 random steps for every 1 observed step)
set.seed(19)
steps_randr <- steps_obsr %>%
  random_steps(n_control = 10)

#We have our random steps, now we need to extract the habitat variables at the end of each step - The HMM state is preserved across each step
head(steps_randr)

sr <- steps_randr %>% drop_na() %>% dplyr::mutate(ID = 1:n()) #drop NAs and create an ID column to calculate spatial pattern metrics across

##We have our data in the format we need to both extract classes and run spatial pattern metrics on them. However, do you remember when we
##told you that the environmental data needs to match the resolution of your environmental data? Well, here is where it comes into play.
##For the extraction of SAV classes, we are going to use a kernel to calculate the class most present within a 10m x 10m square around each
##cell in our map. This mimics the 10m error we filtered by for our positions. We will then extract the class at each position using this new map.

#create a kernel for the map to calculate the majority class with
kernel <- matrix(1, nrow = 31, ncol = 31)

#Apply the kernel across all of the cells in the original map and plot it
rmap_kernel <- focal(rmap, w = kernel, fun = "modal", na.rm = TRUE)

plot(rmap_kernel) #it reduces the complexity, but is better for fitting the data. We will use the complexity of the original map to calculate spatial pattern metrics

#create a spatialpointsdataframe and extract from map
sr2 <- st_as_sf(sr, coords = c("x2_", "y2_"), crs = st_crs(rmap_kernel)) %>% terra::extract(x = rmap_kernel, y = ., ID = F) %>% bind_cols(sr) %>% 
  mutate(class = if_else(focal_modal == 1,2,1)) %>% drop_na() %>% dplyr::select(-focal_modal)

##If the class is 2, it is dense SAV and if it is 1, it is sparse SAV/bare sediment. Now that we have our classes extracted, we are going to
##calculate some spatial pattern metrics. To do this, we are going to calculate the 5 variables brought up in the intro: PLAND, ED, PD, Divisionm and Shape_mn.

##This next code block runs the spatial pattern metrics. However, it takes about 6.5 minutes to run this on a fast computer. For the sake of time, we are going to
##provide the resulting dataset. Unlike the YAPS position calculations, this code WILL work if you decide to run it. This will be the last step
##before we run the iSSFs!

#Just keep coordinates for the landscapemetrics package to recognize it
rcoord <- as.data.frame(cbind(sr2$x2_, sr2$y2_))
rcoord_sf <- st_as_sf(rcoord, coords = c("V1", "V2"), crs = st_crs(rmap))


plan(multisession, workers = max(1, parallel::detectCores() - 1))  # Use multisession to run on multiple cores


result_metr <- 1:nrow(rcoord_sf) %>% future_map(
  function(i) {
    rmap_r <- rast("maps/rmap_dense.tif")
    sample_lsm(
      rmap_r,
      y = st_geometry(rcoord_sf[i, ]),  # Ensure geometry is passed
      size = 15, #15 is for a 5m perception radius and a 10m positional error
      what = c(
        "lsm_c_pd", "lsm_c_ed",
        "lsm_c_division",
        "lsm_c_pland", "lsm_c_shape_mn"
      ),
      shape = "circle"
    )
  }
)
plan(sequential)

rt <-lapply(result_metr, as.data.frame)

rt<- lapply(rt, function(x) if (is.null(x)) data.frame() else x)

all_columns <- unique(unlist(lapply(rt, colnames)))
rt <- lapply(rt, function(df) {
  df[setdiff(all_columns, colnames(df))] <- NA  # Add missing columns with NA
  df <- df[, all_columns]  # Reorder columns
  df
})
metr <- bind_rows(rt, .id = 'source')

unique_sources <- unique(metr$source)
unique_classes <- unique(metr$class)
unique_metric <- unique(metr$metric)



# Create a complete grid of all combinations of `source` and `class`
full_combinations <- expand.grid(source = unique_sources, class = unique_classes, metric = unique_metric)

# Perform a left join with the original dataset to add NAs for missing combinations
mettry <- full_combinations %>%
  left_join(metr, by = c("source", "class", 'metric')) %>% distinct()

metr1 <- mettry %>% pivot_wider(names_from = metric, values_from = value) %>% dplyr::filter(class == 1)
names(metr1)

strue <- metr1 %>% dplyr::select(c(division, ed, pd, pland, shape_mn)) %>% bind_cols(sr2) %>% dplyr::select(-c(ID)) %>% drop_na()

##Now we are ready to run some iSSFs!

####iSSFs####

#Load in libraries
library(amt) #Package that iSSFs are in
library(tidyverse)
library(survival)
library(terra)
library(sf)
library(ggeffects)

strue <- read_csv('data/data_environmental.csv')

# Standardize continuous habitat variables (z-score transformation)
# This helps with model convergence and interpretation of interaction terms.
steps_finalr <- strue %>%
  mutate(
    # Standardize continuous habitat covariates
    pland_z = scale(pland)[, 1],
    ed_z = scale(ed)[, 1],
    pd_z = scale(pd)[, 1],
    shape_mn_z = scale(shape_mn)[, 1],
    division_z = scale(division)[,1],
    
    # Convert state to a factor for proper modeling interaction
    state_f = factor(state)
  ) %>% drop_na()

head(steps_finalr)
corrplot(cor(steps_finalr[,2:6]))

#PLAND and division are highly correlated. Division is also correlated with pd. Let's remove division from our analysis

# Fit the iSSF using conditional logistic regression (clogit).
issf_modelr <- steps_finalr %>%
  drop_na() %>% 
  fit_clogit(case_ ~
               # Movement terms 
               sl_ + log(sl_) +
               cos(ta_) +
               # Habitat selection with interaction with HMM state
               pland_z * state_f + 
               ed_z * state_f +
               pd_z * state_f +
               shape_mn_z * state_f +
               class * state_f +
               # The required stratification term for iSSF
               strata(step_id_), model = T
  )

summary(issf_modelr)

##There are no significant effects except for step length and turn angle. But because this is a workshop, we are going to
##graph things anyway! The most single significant effect is the complexity of the shape of the dense SAV patches. There are no significant interaction effects, but pd and shape are close 
##So we will go over how to graph these interactions and then we can answer any questions! 

mean(steps_finalr$sl_)
max(steps_finalr$sl_)
min(steps_finalr$sl_)
mean(steps_finalr$ta_)
sl_mean <- 2.222 #mean step length (meters)
ang_mean <- 0.01 #mean turn angle (radians)
z_mean <- 0 #mean environmental variable (z-score)
sl_range <- seq(0.1, 22.56, length.out = 100) #range of step lengths (meters)
ang_range <- seq(-3.14, 3.14, length.out = 100) #range of turn angles (radians)
z_range <- seq(-2.1, 2.1, length.out = 100) #range of environmental variables (z-score)
n_points <- length(sl_range) # 100
hab_vars <- c("pland_z", "ed_z", "pd_z", "shape_mn_z")


#create the dataset to test across for step length, keeping everything else constant
tsl <- data.frame(
  sl_ = sl_range, 
  ta_ = rep(ang_mean, n_points),             
  pland_z = rep(0, n_points), ed_z = rep(0, n_points), pd_z = rep(0, n_points), shape_mn_z = rep(0, n_points), #all of the environmental variables are centered around 0 because they are z-scored
  state_f = factor(rep("1", n_points), levels = c("1", "2")), 
  class = factor(rep("1", n_points), levels = c("1", "2"))  
)

#Then we need a 1-row dataset representing a single location. This will be the same as the previous one with the sl set to the mean

msl <- data.frame(
  sl_ = sl_mean, 
  ta_ = ang_mean, 
  pland_z = 0, ed_z = 0, pd_z = 0, shape_mn_z = 0,
  state_f = factor("1", levels = c("1", "2")),
  class = factor("1", levels = c("1", "2"))
)

#Calculate the log relative selection strength (log-RSS) for step length
rss_sl <- log_rss(issf_modelr, x1 = tsl, x2 = msl, ci = 'se')

#We can plot the step length log-RSS curve
plot_sl <- ggplot(rss_sl$df, aes(x = sl_x1, y = log_rss)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkblue") +
  theme_bw() +
  labs(
    title = "Log-RSS for Step Length (zero turn angle)",
    x = "Step Length (meters)",
    y = "Log Relative Selection Strength"
  )
plot_sl


##Looks like our fish prefer really short steps  compared to the random steps created. The random steps were
##created using distributions from the actual data. 

#Now we can plot turn angle! 

#create the dataset to test across for turn angle, keeping everything else constant
tang <- data.frame(
  sl_ = rep(sl_mean, n_points), 
  ta_ = ang_range,             
  pland_z = rep(0, n_points), ed_z = rep(0, n_points), pd_z = rep(0, n_points), shape_mn_z = rep(0, n_points), #all of the environmental variables are centered around 0 because they are z-scored
  state_f = factor(rep("1", n_points), levels = c("1", "2")), 
  class = factor(rep("1", n_points), levels = c("1", "2"))  
)

#1-row at mean turn angle
mang <- data.frame(
  sl_ = sl_mean, 
  ta_ = ang_mean, 
  pland_z = 0, ed_z = 0, pd_z = 0, shape_mn_z = 0,
  state_f = factor("1", levels = c("1", "2")),
  class = factor("1", levels = c("1", "2"))
)

#Calculate the log relative selection strength (log-RSS) for turn angle
rss_ang <- log_rss(issf_modelr, x1 = tang, x2 = mang, ci = 'se')



plot_ang <- ggplot(rss_ang$df, aes(x = ta_x1, y = log_rss)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkblue") +
  theme_bw() +
  labs(
    title = "Log-RSS for Turn Angle (mean step length)",
    x = "Turn angle (radians)",
    y = "Log Relative Selection Strength"
  )
plot_ang

##Looks like our fish liked to go less straight than the random steps. We can also do the same thing
##for the environmental variables. Since shape is the most significant one, we will only show that one.

#create the dataset to test across for shape, keeping everything else constant
tshape <- data.frame(
  sl_ = rep(sl_mean, n_points), 
  ta_ = rep(ang_mean, n_points),             
  pland_z = rep(0, n_points), ed_z = rep(0, n_points), pd_z = rep(0, n_points), shape_mn_z = z_range,  #all of the environmental variables are centered around 0 because they are z-scored
  state_f = factor(rep("1", n_points), levels = c("1", "2")), 
  class = factor(rep("1", n_points), levels = c("1", "2"))  
)

#1-row
mshape <- data.frame(
  sl_ = sl_mean, 
  ta_ = ang_mean, 
  pland_z = 0, ed_z = 0, pd_z = 0, shape_mn_z = 0,
  state_f = factor("1", levels = c("1", "2")),
  class = factor("1", levels = c("1", "2"))
)

#Calculate the log relative selection strength (log-RSS) for shape
rss_shape <- log_rss(issf_modelr, x1 = tshape, x2 = mshape, ci = 'se')


plot_shape <- ggplot(rss_shape$df, aes(x = shape_mn_z_x1, y = log_rss)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkblue") +
  theme_bw() +
  labs(
    title = "Log-RSS for Mean Shape",
    x = "Shape_mn (z-score)",
    y = "Log Relative Selection Strength"
  )
plot_shape

##Looks like our fish prefer seascapes with less complex patches! Now we can look at the interactions between state and environmental variables.
##If there was a significant result, that would mean that the fish preferred certain environmental variables when in a certain state.
##The most significant interaction was state and shape, but pd was also close. We will graph both.

#shape first

#create the dataset to test across for both state and shape, keeping everything else constant
tshapestate <- data.frame(
  sl_ = rep(sl_mean, n_points*2), 
  ta_ = rep(ang_mean, n_points*2),             
  pland_z = rep(0, n_points*2), ed_z = rep(0, n_points*2), pd_z = rep(0, n_points*2), shape_mn_z = rep(z_range, 2),#all of the environmental variables are centered around 0 because they are z-scored
  state_f = factor(rep(c("1", "2"), each = n_points)), 
  class = factor(rep("1", n_points), levels = c("1", "2"))  
)

#1-row
mshapestate <- data.frame(
  sl_ = sl_mean, 
  ta_ = ang_mean, 
  pland_z = 0, ed_z = 0, pd_z = 0, shape_mn_z = 0, 
  state_f = factor("1", levels = c("1", "2")),
  class = factor("1", levels = c("1", "2"))
)

#Calculate the log relative selection strength (log-RSS) for both states relative to baseline
rss_shapestate <- log_rss(issf_modelr, x1 = tshapestate, x2 = mshapestate, ci = 'se')

#separate out based on state
df_state1 <- rss_shapestate$df %>% filter(state_f_x1 == "1") %>% arrange(get("shape_mn_z_x1"))
df_state2 <- rss_shapestate$df %>% filter(state_f_x1 == "2") %>% arrange(get("shape_mn_z_x1"))

#calculate SE of the difference (using a conservative approximation, ignoring covariance)
#SE_diff = sqrt(SE_1^2 + SE_2^2)
se_diff <- sqrt(((df_state1$upr-df_state1$lwr)/3.92)^2 + ((df_state2$upr-df_state2$lwr)/3.92)^2)

#create a new data frame for the difference plot
df_diff <- data.frame(
  z_score = df_state1[["shape_mn_z_x1"]],
  log_rss_diff = df_state2$log_rss - df_state1$log_rss,
  # Calculate 95% CI bounds (Z-score 1.96 for 95%)
  lwr_diff = (df_state2$log_rss - df_state1$log_rss) - 1.96 * se_diff,
  upr_diff = (df_state2$log_rss - df_state1$log_rss) + 1.96 * se_diff
)

#plot the difference with confidence interval
plot_shapestate <- ggplot(df_diff, aes(x = z_score, y = log_rss_diff)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_line(linewidth = 1.2, color = "darkblue") +
  # Add the confidence interval ribbon
  geom_ribbon(aes(ymin = lwr_diff, ymax = upr_diff), alpha = 0.2, fill = "darkblue") +
  theme_bw() +
  labs(
    title = "State-Dependent Effect: Difference in Mean Shape Selection",
    x = "Shape_mn (z-score)",
    y = "log-RSS Difference (State 2 - State 1)"
  )

plot_shapestate

##Not significant, but shows that fish are more likely to be in state 2 ("foraging state") over areas with more complex shapes
##So how about PD?

#create the dataset to test across for both state and pd, keeping everything else constant
tpdstate <- data.frame(
  sl_ = rep(sl_mean, n_points*2), 
  ta_ = rep(ang_mean, n_points*2),             
  pland_z = rep(0, n_points*2), ed_z = rep(0, n_points*2), pd_z = rep(z_range, 2), shape_mn_z = rep(0, n_points*2), division_z = rep(0, n_points*2), #all of the environmental variables are centered around 0 because they are z-scored
  state_f = factor(rep(c("1", "2"), each = n_points)), 
  class = factor(rep("1", n_points), levels = c("1", "2"))  
)

#1-row
mpdstate <- data.frame(
  sl_ = sl_mean, 
  ta_ = ang_mean, 
  pland_z = 0, ed_z = 0, pd_z = 0, shape_mn_z = 0, division_z = 0,
  state_f = factor("1", levels = c("1", "2")),
  class = factor("1", levels = c("1", "2"))
)

#Calculate the log relative selection strength (log-RSS) for both states relative to baseline
rss_pdstate <- log_rss(issf_modelr, x1 = tpdstate, x2 = mpdstate, ci = 'se')

#separate out based on state
df_state1_ps <- rss_pdstate$df %>% filter(state_f_x1 == "1") %>% arrange(get("pd_z_x1"))
df_state2_ps <- rss_pdstate$df %>% filter(state_f_x1 == "2") %>% arrange(get("pd_z_x1"))

#calculate SE of the difference (using a conservative approximation, ignoring covariance)
#SE_diff = sqrt(SE_1^2 + SE_2^2)
se_diff_ps <- sqrt(((df_state1_ps$upr-df_state1_ps$lwr)/3.92)^2 + ((df_state2_ps$upr-df_state2_ps$lwr)/3.92)^2)

#create a new data frame for the difference plot
df_diff_ps <- data.frame(
  z_score = df_state1_ps[["pd_z_x1"]],
  log_rss_diff = df_state2_ps$log_rss - df_state1_ps$log_rss,
  # Calculate 95% CI bounds (Z-score 1.96 for 95%)
  lwr_diff = (df_state2_ps$log_rss - df_state1_ps$log_rss) - 1.96 * se_diff_ps,
  upr_diff = (df_state2_ps$log_rss - df_state1_ps$log_rss) + 1.96 * se_diff_ps
)

#plot the difference with confidence interval
plot_pdstate <- ggplot(df_diff_ps, aes(x = z_score, y = log_rss_diff)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_line(linewidth = 1.2, color = "darkblue") +
  # Add the confidence interval ribbon
  geom_ribbon(aes(ymin = lwr_diff, ymax = upr_diff), alpha = 0.2, fill = "darkblue") +
  theme_bw() +
  labs(
    title = "State-Dependent Effect: Difference in PD Selection",
    x = "PD (z-score)",
    y = "log-RSS Difference (State 2 - State 1)"
  )

plot_pdstate

##Looks like our fish like areas with more patches while "foraging"!

##And thus concludes our workshop on going from raw acoustic telemetry data through analyses to results. Yes, yes, yes. 
##In the title, we said "Downloads to Discussions." So where is the discussion section? Well friends, that is up to you to interpret your data.
##For this example, you will just have to read the paper when it is out! Thanks for joining us!

