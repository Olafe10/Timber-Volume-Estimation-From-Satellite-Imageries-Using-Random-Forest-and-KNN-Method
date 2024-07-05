
#========================================================================
# TASK 2. ESTIMATION OF TIMBER VOLUME WITH THE RANDOM FOREST (RF) METHOD
#========================================================================

rm(list = ls()) #Clean R's memory

# When running the code first time, these additional libraries need to be installed.
# If you didn't install them yet, decomment the lines below, run to install, and 
# comment again. 

# Note that this lab uses an older geospatial library (sp) than the lidR lab (sf). 
# Some of the geoprocessing operation are therefore different.

# install.packages("Rcpp") 
# install.packages("raster")
# install.packages("sp")
# install.packages("rgeos")
# install.packages("rgdal")
# install.packages("randomForest")
# install.packages("yaImpute")
# install.packages("gower")

# Load the installed libraries
require(Rcpp) 
require(raster)
require(sp)
require(rgeos)
require(rgdal)
require(randomForest)
require(yaImpute)
require(gower)


# Set your own working directory
setwd("C:/Users/Olamidayo Fasalejo/OneDrive - University of Eastern Finland/Documents/Remote sensing/Lab2")

# Define function for computing RMSE and RMSE-%
rmse <- function(x1,x2) { sqrt(sum((x1-x2)^2) / length(x1)) }
relrmse <- function(x1,x2) { sqrt(sum((x1-x2)^2) / length(x1)) / mean(x1) }

#=============================================================
# Extract Sentinel-2 for Juupajoki plots 
#=============================================================

# Read in field plots
juu_plots <- read.csv("Juupajoki_field_data.csv")
tail(juu_plots)

# Set coordinate system
coordinates(juu_plots) <- c("xcoordinate","ycoordinate")
proj4string(juu_plots) <- CRS("+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs") # TM35FIN parameters

# Read Sentinel-2
sentinel2 <- brick("S2_Juupajoki.tif")

# Plot satellite image and add the field plot locations
plotRGB(sentinel2, r=7, g=3, b=2, stretch="lin") 
plot(juu_plots, col="yellow", pch=16, cex=.75, add=T)

# Set plot radius based on sample plot type
juu_plots$rad <- 9
juu_plots$rad[juu_plots$sampleplottype==2] <- 5.64
juu_plots$rad[juu_plots$sampleplottype==4] <- 12.62

# Convert plot centers to polygons based on the set radius
buffered <- gBuffer(juu_plots, width = juu_plots$rad, byid=T)

# Extract as data frame, get a normalized weight for each pixel within a polygon
image_df <- extract(sentinel2, buffered, df=T, weights=T, normalizeWeights=T) 

# Rename columns
names(image_df) <- c("ID", "blue", "green", "red", "re1", "re2", "re3", "nir", "nirnar", "swir1", "swir2", "weight")
tail(image_df)

# Function for computing the weighted mean by band
wmean <- function(band){ 
  wsum <- image_df[,band] * image_df$weight
  tapply(wsum, image_df$ID, sum)
}

# Apply function with the named columns
juu_s2 <- as.data.frame(cbind(juu_plots$sampleplotid, wmean("blue"), wmean("green"), wmean("red"), wmean("re1"), wmean("re2"),
                              wmean("re3"), wmean("nir"), wmean("nirnar"), wmean("swir1"), wmean("swir2")))

# Rename columns
names(juu_s2) <- c("sampleplotid", "blue", "green", "red", "re1", "re2", "re3", "nir", "nirnar", "swir1", "swir2")

# View last rows
tail(juu_s2)


#=============================================================
# Extract PALSAR-2 for Juupajoki plots 
#=============================================================

# Read PALSAR-2 bands
hh <- raster("PALSAR-2_HH_Juupajoki.tif")
hv <- raster("PALSAR-2_HV_Juupajoki.tif")

# Merge bands into a single raster brick
palsar2 <- stack(hh,hv)

# Plot HH band
plotRGB(palsar2, r=1, g=1, b=1, stretch="lin") 
plot(juu_plots, col="yellow", pch=16, cex=.75, add=T)
?plotRGB
# Extract radar data from a fixed 20 m radius around each plot to decrease noise
juu_plots$rad <- 20
buffered <- gBuffer(juu_plots, width = juu_plots$rad, byid=T)
image_df <- extract(palsar2, buffered, df=T, weights=T, normalizeWeights=T) 

# Rename columns
names(image_df) <- c("ID", "hh", "hv", "weight")
tail(image_df)

# Function for computing the weighted mean by band
wmean <- function(band){ 
  wsum <- image_df[ , band] * image_df$weight
  tapply(wsum, image_df$ID, sum)
}

#Apply function with the named columns
juu_p2 <- as.data.frame(cbind(juu_plots$sampleplotid, wmean("hh"), wmean("hv")))

names(juu_p2) <- c("sampleplotid", "hh", "hv")
tail(juu_p2)


#=============================================================
# Merge data sets
#=============================================================

# Merge Juupajoki data frames 
juu_merged <- cbind(as.data.frame(juu_plots), juu_s2, juu_p2)

# Remove repeated "sampleplotid" columns
juu_merged <- juu_merged[ , -c(34,45)]
dim(juu_merged)

# Read plots from the remaining Orivesi area
ori_merged <- read.csv("Orivesi_field+RS_data.csv")
dim(ori_merged)

# Merge Juupajoki and Orivesi plots into new data frame d
d <- rbind(juu_merged, ori_merged)


#=============================================================
# Construct random forest models
#=============================================================

dev.off() # Close previous plots

d# Sentinel-2 + PALSAR-2 model
s2p2 <- randomForest(volume~ blue+green+red+re1+re2+re3+nir+nirnar+swir1+swir2+hh+hv, data=d)    # Write model here
predict(s2p2)   #report the best
summary(s2p2)

varImpPlot(s2p2, type=2)    #report the best
rmse(d$volume, predict(s2p2))  #85.45 #report
relrmse(d$volume, predict(s2p2)) #0.46 #report

# Scatter plot base code
plot(d$volume, predict(s2p2), xlim=c(0,800), ylim=c(0,800)); abline(0,1)  #report


# Sentinel-2 model

s2 <- randomForest(volume~ blue+green+red+re1+re2+re3+nir+nirnar+swir1+swir2, data = d)
plot(d$volume, predict(s2), xlim=c(0,800), ylim=c(0,800)); abline(0,1)    #report
varImpPlot(s2, type=2)
rmse(d$volume, predict(s2))  #86.17 #report
relrmse(d$volume, predict(s2)) #0.46 #report
#=============================================================
# Processing the map data
#=============================================================

# Resample PALSAR-2 to same resolution as Sentinel-2
palsar2 <- resample(palsar2, sentinel2, method="bilinear") 

# Stack S2 and  layers into a single raster
sentinel2 <- stack(sentinel2, palsar2) 

plotRGB(sentinel2, r=11, g=4, b=3, stretch="lin")     #report

# Convert raster stack into data frame
imagedataframe <- as.data.frame(sentinel2) 

# Add column names for bands
names(imagedataframe) <- c("blue", "green", "red", "re1", "re2", "re3", "nir", "nirnar", "swir1", "swir2", "hh", "hv")

# Read NFI volume raster
nfivol <- raster("Juupajoki_NFI_Volume.tif")

# Set up a forest mask: areas that are not nodata are forest
formask <- !is.na(nfivol)

# Plot the mask
plot(formask)

# Convert also forest mask into data frame
maskdataframe <- as.data.frame(formask)


#=============================================================
# Predict with the earlier model object s2p2
#=============================================================

# Apply the previous rf model with also radar bands
pred <- predict(s2p2, imagedataframe)

# Replace the predicted values for non-forest areas with nodata
pred[! maskdataframe$layer] <- NA

# Construct a new map
# Make a copy of the first image band
volmap <- sentinel2[[1]] 

# Replace original values by model predictions
values(volmap) <- pred 

# Show the  map
plot(volmap) 

# Write the map as a tif file
writeRaster(volmap,"RF_Volume.tif",overwrite=T)


#=============================================================
# Comparison
#=============================================================

# Initialize plotting
dev.off()
par(mfrow=c(1,2))

v <- pred[!is.na(pred)] # Volume vector without NA's

# Compute mean volume and make a histogram from your map
mean(v)  #report mean value 152.64  
hist(v)  #report this

# Compute mean volume and make a histogram from the NFI map
nfivol_df<- as.data.frame(nfivol)
mean(nfivol@data$Juupajoki_NFI_Volume)

#mean 145.6553   #report this
mean(nfivol_df$Juupajoki_NFI_Volume, na.rm=T)


hist(nfivol_df$Juupajoki_NFI_Volume)

#===============================================================================
# TASK 3: SPECIES-SPECIFIC VOLUME ESTIMATION WITH THE K NEAREST NEIGHBOR METHOD 
#===============================================================================

# The code continues from where the previous task finished

# Define a vector with response variable names
yvar <- c("pinevol", "sprucevol", "decidvol") 

# Define a weight for each response variable
weights <- c(1,1,1) 

# Define the column range that contains the predictor variables 
colmin <- 34 # 34 = predictors start from 34th column (blue)
colmax <- 43 # 43 = predictors end at 43th column (swir2)

# KNN parameters

prednum <- 5        # How many predictors to search for
KNNK <- 5           # How many nearest neighbors for KNN
met <- "msn"        # Distance metric for knn. Type ?yai for more info.
wm <- "dstWeighted" # Weighting of neighbors = inverse distance .

# Simulated annealing parameters

t_ini <- 0.2        # Initial temperature at simulated annealing
n_iter <- 10000     # How many iterations in optimization; larger is better but takes longer

# Extract response variables into their own data frame "ytrain" 
ytrain<-as.data.frame(d[,yvar]) 

# Rename columns
names(ytrain)<-yvar 

#=============================================================
# Variable selection functions
#=============================================================

# A function that generates alternative predictor variable combinations

PickNewSolution <- function(s,kk,n_iterr){
  
  # s = current solution, k=number of iterations
  # Initially change 2/3 of the variable
  # When k/n >=.8, change only one variable
  s <- sort(s)
  keepmax <- prednum-1          # Maximum number of predictors to keep = prednum-1
  keepmin <- round(1/3*prednum) # Minimum number of predictors to keep = 1/3*length
  
  kperc <- kk/n_iterr           # Iteration progress
  
  # Current number of predictors to keep
  keepcount <- min(length(s)-1, round((keepmax-keepmin) / 0.8*kperc + keepmin)); keepcount
  
  idxx <- colmin:colmax                                 # List of variables to select from
  keep <- sort(sample(s,keepcount,replace=F))           # List of untouchable variables
  ssel <- idxx[! idxx %in% s]                           # List of selectable variables - current sample not allowed
  newvars <- sample(ssel, prednum-keepcount, replace=F) # Sampling from selectables
  sol<-sort(c(newvars,keep))                            # Merge
  sol                                                  # Return new variable list
}

# yaisel-function fits a KNN model and returns the weighted mean RMSE

yaisel<-function(x){ 
  
  xtrain <- as.data.frame(d[,x])                         # List X-variables
  knn <- yai(y=ytrain, x=xtrain, method = met, k=KNNK);  # Fit KNN model
  pred <- impute(knn,k=KNNK,method=wm)                   # Get fitted values
  
  # Initialize rmse vector
  rvec <- -999 
  
  # Calculate mean relative rmse for each response variable
  for(i in 1:length(yvar)) rvec[i] <- relrmse(pred[,i], pred[ , i+length(yvar)])
  
  # Returned value = mean of weighted RMSEs
  mean(rvec*weights) 
}


#=============================================================
# Variable selection by simulated annealing
#=============================================================

t <- t_ini                                     # Initial temperature
s <- sample(colmin:colmax, prednum, replace=F) # Initial variables

e <- yaisel(s) # Run KNN  for initial variables

ebest <- e     # Save initial mean rmse
sbest <- s     # Save initial variable combination

k <- 0         # Initialize iteration counter

while(k < n_iter){ 
  
  # sdot  = new experimental solution
  # s     = current solution to be improved,
  # sbest = best solution ever found
  
  sdot <- PickNewSolution(s, k, n_iter) # New candidate variables
  edot <- yaisel(sdot)                  # KNN result for the new candidate variables
  
  # Implement the simulated annealing algorithm
  if(exp((-(edot-e))/t) > runif(1)){
    e <- edot
    s <- sdot
  }
  if(edot < ebest){
    ebest <- edot
    sbest <- sdot
  }
  t <- max(0, -.2/.8*k/n_iter+.2)      # Cool temperature
  k <- k+1
  
}


names(d)[sbest] # Print selected variable names   # selected variable are green, nir, nirnar, swir1, swir2
ebest           # Print the mean rmse with these variables     #the mean rmse is 1.07

#=============================================================
# Fitting a KNN model with the selected variables
#=============================================================

# Type a vector with the selected variable combination here
xvar <- c("green", "nir", "nirnar", "swir1", "swir2")

xtrain <- d[, xvar] # Extract x-variables into their own raster

tail(xtrain) # View last rows

# Rename columns to avoid problems later on
row.names(ytrain) <- 1:nrow(ytrain)
row.names(xtrain) <- 1:nrow(xtrain)

# Train the model
knn <- yai(y=ytrain, x=xtrain, method = met, k=KNNK); 

# Get and view predicted values
pred <- impute(knn, k=KNNK, method=wm) 

tail(pred)  # Variables with .o = observed, without .o = predicted

# Calculate relative RMSEs for all species!
relrmse(pred$pinevol.o, pred$pinevol)     #report rmse 1.14
relrmse(pred$sprucevol.o, pred$sprucevol)  #report rmse 0.88
relrmse(pred$decidvol.o, pred$decidvol)    #report rmse 1.25

# Draw scatterplots with estimated volume on x and observed volume on y axis
# for all species! Remember to add an 1:1 line.
dev.off()
plot(pred$pinevol.o, pred$pinevol); abline(1,1)    #report
plot(pred$sprucevol.o, pred$sprucevol); abline(1,1)  #report
plot(pred$decidvol.o, pred$decidvol); abline(1,1)  #report

#=============================================================
# Species mapping
#=============================================================

# Use data extraced previously from the images
df <- imagedataframe[, xvar]

# Change the row names to start where the training data row names ended
row.names(df) <- as.numeric(row.names(df)) + nrow(d)

# Finding the nn references for the test data based on the previous model
knn2 <- newtargets(knn, newdata=df)

# Prediction
pred  <- impute(knn2, vars=yvars(knn2), method=wm, k=KNNK)

tail(pred)

# Replace the predicted values for non-forest areas with nodata
pred$pinevol[! maskdataframe$layer] <- NA
pred$sprucevol[! maskdataframe$layer] <- NA
pred$decidvol[! maskdataframe$layer] <- NA

# Construct a new map

# Make copies of an image band
pinevmap <- sentinel2[[1]] 
sprucevmap <- sentinel2[[1]]
decidvmap <- sentinel2[[1]] 

# Replace values
values(pinevmap) <- pred$pinevol 
values(sprucevmap) <- pred$sprucevol 
values(decidvmap) <- pred$decidvol 

# Stack species layers
s <- stack(pinevmap, sprucevmap, decidvmap)

# Write stacked raster as .tif
writeRaster(s, "Species_stack.tif", overwrite=T) 
a