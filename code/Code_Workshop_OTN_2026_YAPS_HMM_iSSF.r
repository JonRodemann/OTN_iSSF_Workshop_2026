##Code for Workshop##
##Authors: Jon Rodemann, Sophia Costa, Jennifer Rehage
##Source Data - Rodemann et al., 2025 Sci. Rep., Costa et al., 2026 (In Prep)

####YAPS####

##This is the general workflow for syncing receivers and running tracks in YAPS. We won't be running this code in this workshop
##and in fact, this code does not work because the syncing fails. In reality, you need to play around with your sync data 
##(including splitting up the time periods if there are "blackouts", possibly removing receivers, and other adjustments) for the sync model
##to be successful. We will pick up the workshop with the positions that were calculated by the developer of YAPS (WOW!), Henrik Baktoft.

#Install and load required libraries
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

##Now we have the data we need for the HMMs!

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

##We are not able to fit some HMMs! To do this, we need to first prep the data and then remove missing data

interp_df_30 <- id20_30 %>% arrange(ID, ts)
# prepData: coordinate names are "x" and "y" in our frame
hmm_data_30 <- prepData(interp_df_30, type = "UTM", coordNames = c("x", "y"))

# moveHMM uses a column 'step' and 'angle' (radians) created by prepData
# Remove rows with zero or NA step (stay/missing movement) because gamma requires positive steps
hmm_data_30 <- hmm_data_30 %>% filter(!is.na(step) & step > 0) %>% filter(!is.na(angle))

hmm_data_30 <- hmm_data_30 %>%
  filter(is.finite(step), is.finite(angle), step > 0) %>%
  mutate(
    step = pmax(step, quantile(step, 0.01, na.rm = TRUE)),
    step = pmin(step, quantile(step, 0.99, na.rm = TRUE))
  )

##We have prepped our data. However, it is still not just plug and play. We need to investigate our data to fit 
##general curves to our data to inform the model. We don't need to do this, but it creates a more informed model based
##on the original data. We are trying to fit these data into two behavioral states, so we will want to inform the model with 2 curves

hist(hmm_data_30$step, breaks = 20, freq = F, main = "Step lengths")
curve(dgamma(x, shape = 1.5, rate=0.9), n=200, add=TRUE, col="royalblue", lwd=4) #shape is how right-skewed the gamma curve is
curve(dgamma(x, shape = 1.9, rate=0.7), n=200, add=TRUE, col="cadetblue", lwd=4) #rate is how tall the curve apex is
curve(dgamma(x, shape = 2.2, rate=0.4), n=200, add=TRUE, col="navyblue", lwd=4)

##this is a little subjective, but there is a good starting out point. The first curve encompasses the lower step lengths well.
##The second curve does not really add any more information either than including the step lengths that are longer. However,
##the third curve does this but also includes the higher amount of longer steps. So let's go with the first and third curves

#assign initial step length distributions
sl_init_mean <- c(2.2/0.4, 1.5/0.9) #initial means should be the shape/rate for both
sl_init_sd <- c(sqrt(2.2)/0.7, sqrt(1.5)/0.9) #initial standard deviations should be square root of shape/rate

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
##(average of 4 meters) and straight turn angles (0.007 radians, basically zero). State 2, on the other hand, has a much shorter
##step length (average of 0.9 meters) but a similar turn angle (-0.02 radians). However, the concentration is a lot lower, so state 2
##encompasses the larger turn angles. These results also give you a transistion probability to transition between states. Pretty low for both,
##but it does happen. If you run the following code, you can visualize the distributions for each state as well as which state each step on all the
##tracks is! Since that is 59 tracks, we won't run it online but we will go through it in the workshop

# statecol <- c('lightblue', 'purple1')
# plot(mod, col = statecol)

##Now we can assign states directly to the dataframe

states <- viterbi(mod)

mod$data
state_30 <- mod$data %>%
  mutate(state = states)

# Convert angle to degrees for interpretability (angles are in radians)
state_30 <- state_30 %>%
  mutate(angle_deg = angle * 180 / pi,
         abs_angle_deg = abs(angle_deg))

output_HMM_df <- state_30 %>%
  dplyr::select(ID, track, tag, ts, x, y, step, angle, angle_deg, state)

write_csv(output_HMM_df, "data/interpolated_hmm_states_30min.csv")

##We have our final HMM dataset! Now we can move onto looking at what environmental variables we can use to model against the
##HMMs and for the iSSFs!

