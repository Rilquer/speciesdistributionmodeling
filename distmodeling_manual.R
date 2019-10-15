##Rilquer Mascarenhas
##Last updates: Aug 27 2019

#This script is designed to model the distribution of species using the Biomod package, which can implement up to
#ten algorithms and create the ensembled model from them.
#This script requires a list of species to be modeled. Then, it searches gbif database for species coordinates (using
#package rgbif), filters gbif erros (using package CoordinateCleaner), filters localities geographically too close, by
#a minimum distance set by the user (using custom-made functions) and then proceeds to modeling with package biomod2.
#The script also requires athe path to a folder containing climatic variables for the current period (for calibration),
#as well as the number of pseudo-absence points to be used and the strategy to select those pseudo-absences. Finally,
#it requires files maxent.jar and maxent.bat in the same directory the script is run (these files are provided along
#with the script in the github page. More specific info on biomod modeling can be seen in the Biomod help files (specif-
#ically acessed by using ??BIOMOD_FormatingData in R).
#More information is described in the commented sections below."

#####################################################
######### Preparing packages and functions ##########
#####################################################

##Checking if packages are present; if not, installing them
if (!(require(raster))) {
  install.packages('raster')
}
if (!(require(rgdal))) {
  install.packages('rgdal')
}
if (!(require(rgeos))) {
  install.packages('rgeos')
}
if (!(require(biomod2))) {
  install.packages('biomod2')
}
if (!(require('rgbif'))) {
  devtools::install_github("ropensci/rgbif")
}
if (!(require('CoordinateCleaner'))) {
  devtools::install_github("ropensci/CoordinateCleaner")
}
if (!(require('countrycode'))) {
  devtools::install_github("vincentarelbundock/countrycode")
}

##Loading packages
library(raster)
library(rgdal)
library(rgeos)
library(biomod2)
library(rgbif)
library(CoordinateCleaner)
library(countrycode)

##Setting distcalc and geodel functions to geographic filtering
#Function to caculate distance matrix
distcalc <- function(x,m) { 
  coordinates(x) <- c('long', 'lat')
  proj4string(x) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  locs <- spTransform(x, CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  dist<-gWithinDistance(locs, dist = m, byid = TRUE)
  return(dist)
}

#Function that deal with pairs from previous distance matrix
geodel <- function(x,m) {
  points <- data.frame(x)
  nloci <- nrow(points)
  message("Calculating distance matrix among ",nloci," points and removing points...")
  
  #1. Calculating distance matrix. This section calculates a distance matrix and
  #then sorts out which points are mostly closer to others than the minimum distance
  #set. The point that is close to the highest number of neighbouring points is removed
  #and the matrix is calculated again. This process is repeated to remove those points
  #are close to many others, until each point is close only to one other point.
  nremoved=0
  repeat {
    distmat <- distcalc(points,m)
    diag(distmat) <- NA
    nloc <- ncol(distmat)
    quant<-c()
    for (i in 1:nloc) {
      q <- c(i,(nrow(subset(distmat,distmat[,i]==TRUE))))
      quant <- data.frame(rbind(quant,q))
    }
    colnames(quant) <- c("loc","quant")
    quant<-quant[order(-quant$quant),] ##Sorting points from highest to lowest
    ##This condition tests if all values in column quant is equal or less than 1.
    #If so, the loop breaks. If not, the locality with the highest number of
    #neighbouring points is removed from the original dataset and the distance
    #matrix is calculated again.
    if (all(quant$quant<=1)) {
      break  
    } else { 
      remove<-quant[1,1]
      points<-points[-remove,]
      nremoved<-nremoved+1
      cat(".")
    }
  }
  
  message("\n")
  message("Removing final points...")
  
  #2. Look for close pairs and removes the first close pair found (in the for loop). Then stops the loop
  #and calculates the distance matrix again. The repeat loop will only break if the pair object is created
  #but no value is stored in it during the for loop (meaning no close localities were found in the dataset).
  repeat {
    pair <- c()
    distmat[lower.tri(distmat, diag=TRUE)]<-NA
    nloc<-nrow(points)
    b <- FALSE
    for (i in 2:nloc) {
      for (j in 1:(i-1)) {
        if (distmat[j,i]==TRUE) {
          pair <- c(i,j)
          locremove<-pair[sample(c(1:2),1)]
          points<-points[-locremove,]
          nremoved=nremoved+1
          cat(".")
          b <- TRUE
        }
      }
      if (b==TRUE) {
        break
      }
    }
    if (b==TRUE) {
      distmat <- distcalc(points,m)
      diag(distmat) <- NA
    } else {
      break
    }
  }
  return(points)
  message("\n")
  message("Dataset successfully filtered. ",nremoved," points were removed.")
}

#####################################################
#########       Setting user options       ##########
#####################################################

##This is where the user will set parameters for the analysis, into the
##list object "opt".
opt <- list(species=NULL,file=NULL,biomes=NULL,fdist=50000,minpoints=15,dir=NULL,format='tif',projection=NULL,panumber=1000,strategy='random',min=NULL,max=NULL,sre.proportion=NULL)

#On the line above, all settings were created, and now they need to be set. For each setting below, the previous
#line is commented and explains what the setting is about. Pay attention on mandatory settings, that need to be
#set for the script to run properly.

#Inform path to textfile with list of species below. File must have one species name per line.
opt$species <- ''

#Alternatively, inform path to file containing the coordinates in decimal format, in Long-lat order, no header
##(i.e., first line is first coordinate) and tab separated.
opt$file <- ''

##Getting biomes shapefile
##"Complete path to directory containing shapefile to be used when plotting species coordinates. This plot is used to show the final set of coordinates used for modeling each species. MUST contain the name of the layer (name before .shp extension). Mandatory."

opt$biomes <- ''
  
#Minimum distance for geographical filtering. Only localities away from each other more than the distance set here
#will be kept for modeling. Distance must be in meter. Type 0 to deactivate it. Type: numeric. Default: 50000 (50 km).
opt$fdist

#Set here the minimum number of occurrence points to be used for modeling (only species with a number of points equal to
#or higher than the number set here will be modeled). Type: numeric. Default: 15
opt$minpoints

#Path to the directory with the climatic variables from the current period (i.e., present climatic conditions)
#to be used in model calibration.
#Must be typed between quotes.
opt$dir

#Format of variables files (either asc, grd or tif). Default is tif.
#Must be typed between quotes.
opt$format

#Path to the directory with climatic variables for future or past projections (optional, leave blank if
#you do not wish to project the model into other periods.). MUST NOT have the final slash '/' in the end of path.
#Climatic variables for different periods must be separated in different subfolders inside this main folder.
#The name of the subfolders will be the name used for the projections outputed by Biomod for the respective
#period. For example, to model mid-holocene, LGM and LIG, create a folder named, e.g., "Past" and indicate
#the path to it in this setting. Within the "Past" folder, create three sub-folders named mid-holocene, LGM and LIG,
#and place the climatic layers for each period within the respective folder.
opt$projection

#Number of pseudo-asbences to be selected. Default = 1000. This value depends on the number of presence points
#you have. Note that higher values (proportionally to the number of presence points) will tend to create overfitted
#model. This value can be varied in different runs to test the effect of this setting on the final result.
opt$panumber

#Strategy for selection of pseudo-absence points. Either 'random', 'disk' or 'sre'. Random strategy samples pseudo-
#absence from all the background environment in the climatic variables provided to the script. Disk strategy requires
#you to determine a minimum and maximum distance around occurence points within which pseudo-absence points will
#be sampled (this is useful to avoid sampling in cells to close to occurence points, which can result in overfit, and
#to secure inclusion of similar abiotic conditions out of the range of the focal species). Sre strategy first creates
#a bioclimatic envelope of the occurence points and samples pseudo-absence points from areas which differs from this
#envelope at least by a proportion set by the user. As the biomod help informs, this can lead to overoptimistic models.
#Default='random'
opt$strategy

#Minimum distance for strategy 'disk', in meters.
opt$min

#Maximum distance for strategy 'disk', in meters.
opt$max

#Similarity proportion for strategy 'sre' (between 0.25 and 0.5).
opt$sre.proportion

#####################################################
#########       Running the models         ##########
#####################################################

##Setting important variables
coords <- list()
notmodeledsp <- c()
if (!is.null(opt$species)) {
  species <- scan(opt$species,what='character',sep = '\n')  
} else {
  species <- strsplit(opt$file,'/')[[1]][[length(strsplit(opt$file,'/')[[1]])]]
  coords[[1]] <- unique(data.frame(read.table(opt$file,sep='\t')))
}
##ALTERNATIVELY, you can inform the list of species as a character vector in the line below (uncomment line before running)
#species <- c('')

message('\n')
message('List of species that will be modeled: ')
species
message('\n')

##Creating directory for species plots
dir.create('coordinates_plots')
bsplit <- strsplit(opt$biomes,'/')[[1]]
if (bsplit[1]=='') {
  biomes <- readOGR(dsn = paste0('/',paste0(bsplit[2:length(bsplit)-1],'/',collapse = '')),layer = bsplit[length(bsplit)])
} else {
  biomes <- readOGR(dsn = paste0(bsplit[2:length(bsplit)-1],'/',collapse = ''),layer = bsplit[length(bsplit)])
}

##Getting all climatic data
##Current for calibration and projection
message('Reading climatic variables for current period...')
calibperiod <- list.files(paste0(opt$dir,'/'),pattern=paste0('.',opt$format,'$'))
current <- stack(paste0(opt$dir,'/',calibperiod))
layers <- paste0('var',1:NROW(calibperiod))
names(current)<-layers
message('Done!')
message('')

##Climatic data for projection
if (!is.null(opt$projection)) {
  projperiods <- list.dirs(path=opt$projection)
  projperiods <- projperiods[2:NROW(projperiods)]
  for (i in 1:NROW(projperiods)) {
    projvec <- strsplit(projperiods[i],"/")[[1]]
    projperiods[i] <- projvec[length(projvec)]
  }
  projvar <- list()
  for (i in 1:length(projperiods)) {
    message('Reading climatic variables for period ',projperiods[i],'...')
    files <- list.files(paste0(opt$projection,'/',projperiods[i],'/'),pattern=paste0('.',opt$format,'$'))
    projvar[[i]] <- stack(paste0(opt$projection,'/',projperiods[i],'/',files))
    names(projvar[[i]]) <- layers
    message('Done!')
    message('')
  }
}

for (n in 1:length(species)) {

  if (!is.null(opt$species)) {
    ##Getting coordinates and filtering data
    message('Downloading coordinates for ',species[n])
    coords[[n]] <- data.frame(na.omit(occ_search(scientificName = species[n], fields = c('species','decimalLongitude','decimalLatitude','countryCode','institutionCode'))$data))
    message('\n')
    message('Looking for errors in coordinates dataset ',species[n])
    coords[[n]]$countryCode <- countrycode(coords[[n]]$countryCode, origin =  'iso2c', destination = 'iso3c')
    flags <- clean_coordinates(x = coords[[n]], lon = "decimalLongitude", lat = "decimalLatitude",
                               countries = "countryCode", 
                               species = "species",
                               tests = c("capitals","centroids", "equal","gbif", "institutions",
                                         "zeros", "countries"))
  
    coords[[n]] <- coords[[n]][flags$.summary,which(colnames(coords[[n]]) %in% c('decimalLongitude',
                                                                                                         'decimalLatitude'))]
    colnames(coords[[n]]) <- c('long','lat')
  }
  
  if (opt$fdist!=0) {  
    message('\n')
    message('Removing spatial autocorrelation for ',species[n],' (filtering localities less than ',opt$fdist/1000,' km apart).')
    coords[[n]] <- geodel(coords[[n]],opt$fdist)
  }
  
  if (nrow(coords[[n]])<opt$minpoints) {
    message('\n')
    message('Not enough points left to model species ',species[n],' (only ', nrow(coords[[n]]),' points available. Skipping!')
    notmodeledsp<-c(notmodeledsp,species[n])
  } else {
    message('\n')
    message('Plotting coordinates for ',species[n])
    pdf(file = paste0('coordinates_plots/',species[n],'.pdf'))
    plot(biomes,col='seashell')
    plot(SpatialPoints(coords[[n]]), cex = 0.8, pch = 20, col = 'darkgreen', add = T)
    mtext(species[n],side=1,at = c(-40))
    dev.off()
    tiff(file = paste0('coordinates_plots/',species[n],'.tif'))
    plot(biomes,col='seashell')
    plot(SpatialPoints(coords[[n]]), cex = 0.8, pch = 20, col = 'darkgreen', add = T)
    mtext(species[n],side=1,at = c(-40))
    dev.off()
    
    ##BIOMOD_FormatingData
    message('Reading coordinates...')
    coordinates <- SpatialPoints(coords[[n]])
    message('Done!')
    message('')
    
    if (opt$strategy=='random') {
      dados <- BIOMOD_FormatingData(resp.var = coordinates, expl.var = current, resp.name = species[n], PA.nb.rep = 10,
                                    PA.nb.absences = as.numeric(opt$panumber), PA.strategy = opt$strategy)
    } else if (opt$strategy=='disk') {
      dados <- BIOMOD_FormatingData(resp.var = coordinates, expl.var = current, resp.name = species[n], PA.nb.rep = 10,
                                    PA.nb.absences = as.numeric(opt$panumber), PA.strategy = opt$strategy,
                                    PA.dist.min=as.numeric(opt$min),PA.dist.max = as.numeric(opt$max))
    } else {
      dados <- BIOMOD_FormatingData(resp.var = coordinates, expl.var = current, resp.name = species[n], PA.nb.rep = 10,
                                    PA.nb.absences = as.numeric(opt$panumber), PA.strategy = opt$strategy,
                                    PA.sre.quant=as.numeric(opt$sre-proportion))
    }
    
    options <- BIOMOD_ModelingOptions()
    
    ##Calibrando o modelo
    message('Running model calibration for ',species[n])
    message('')
    
    model <- BIOMOD_Modeling(data = dados, models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips'),
                             models.options = options, DataSplit = 75, models.eval.meth = 'TSS', SaveObj = TRUE,
                             do.full.models= TRUE, modeling.id= species[n])
    dir.create(paste0('rastertemp_',species[n]))
    rasterOptions(tmpdir = paste0('rastertemp_',species[n]))
    
    current.projection <- BIOMOD_Projection(modeling.output=model, new.env = current,
                                            proj.name='current.projection', selected.models="all",
                                            binary.meth= 'TSS', output.format=".grd", do.stack=FALSE, build.clamping.mask=TRUE)
    
    if (!is.null(opt$projection)) {
      projmodels <- list()
      for (i in 1:length(projvar)) {
        message('Running projection of model for directory ',projperiods[i],' (species ',species[n],')...')
        message('')
        projmodels[[i]] <- BIOMOD_Projection(modeling.output=model, new.env = projvar[[i]],
                                             proj.name=paste0('period_',projperiods[i],'.projection'), selected.models="all",
                                             binary.meth= 'TSS', output.format=".grd", do.stack=FALSE, build.clamping.mask=TRUE)
      }
    }
    
    ##BIOMOD_EnsembleModeling (biomod2)
    message('Running calibration of ensembled model (species ',species[n],')...')
    message('')
    ensemble.model <- BIOMOD_EnsembleModeling(model, chosen.models = 'all', em.by = 'all',
                                              eval.metric = 'all', eval.metric.quality.threshold = 0.9,
                                              models.eval.meth = 'TSS', prob.mean = TRUE,prob.cv=TRUE,prob.mean.weight=TRUE)
    
    ##BIOMOD_EnsembleForecasting (biomod2)
    current.ensembled <- BIOMOD_EnsembleForecasting(ensemble.model, projection.output = current.projection,
                                                    selected.models = 'all', proj.name='current.ensembled', binary.meth= 'TSS')
    
    if (!is.null(opt$projection)) {
      projensemble=list()
      for (i in 1:length(projmodels)) {
        message('Running projection of ensembled model for directory ',projperiods[i],' (species ',species[n],')...')
        message('')
        projensemble[[i]] <- BIOMOD_EnsembleForecasting(ensemble.model, projection.output = projmodels[[i]],
                                                        selected.models = 'all', proj.name=paste('period_',projperiods[i],'.ensembled',sep=''),
                                                        binary.meth= 'TSS')
      }
    }
    
    unlink(paste0('rastertemp_',species[n]), recursive = TRUE)
  }
}
