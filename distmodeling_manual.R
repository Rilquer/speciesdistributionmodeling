##Rilquer Mascarenhas
##Last updates: Jun 05 2019

#This script is designed to model the distribution of species using the Biomod package,
#which can implement up to ten algorithms and create the ensembled model from them.
#This script requires certain R packages (which will be installed when the script is run,
#if not already present). This script requires you provide a file with coordinates, folder
#paths containing climatic variables and a name for the analysis. It also requires you to
#inform the number of pseudo-absence points to be used and the strategy to select those
#pseudo-absences. Finally, it requires files maxent.jar and maxent.bat in the same directory
#the script is run. More info can be seen in the Biomod help files (specifically acessed by
#using ??BIOMOD_FormatingData in R). This script can be directly called (i.e., Rscript ./script_modeling)
#when in the folder it is located, or copied to the bin folder to be called as an executable from any location.

#This line creates a list of all setting necessary for the script. Don't modify this!!
opt <- list(title=NULL,file=NULL,dir=NULL,format='tif',projection=NULL,panumber=1000,strategy='random',min=NULL,max=NULL,
         sre.proportion=NULL,threads=1)

#On the line above, all settings were created, and now they need to be set. For each setting below, the previous
#line is commented and explains what the setting is about. Pay attention on mandatory settings, that need to be
#set for the script to run properly.

#Name for this modeling analysis. Must be typed between quotes.
opt$title

#Path to file with coordinates. The file MUST be in long-lat format, tab separated and without headers.
#Must be typed between quotes.
opt$file

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

#Number of maximum threads to be used by biomod2. Default=1
opt$threads

#Checking arguments
if (is.null(opt$title)) {
  stop("You must provide the name for the analysis.", call.=FALSE)
}
if (is.null(opt$file)) {
  stop("You must provide a coordinates file (long-lat format, tab delimited and no header).", call.=FALSE)
}
if (is.null(opt$dir)) {
  stop("You must provide a path for the folder containing current climatic variables for model calibration.", call.=FALSE)
}

##Checking necessary packages and installing those not installed
if (!require(sp)) {
  install.packages('sp',dependencies=TRUE,repos='http://cran.us.r-project.org')
  library(sp)
}
if (!require(raster)) {
  install.packages('raster',dependencies=TRUE,repos='http://cran.us.r-project.org')
  library(raster)
}
if (!require(biomod2)) {
  install.packages('biomod2',dependencies=TRUE,repos='http://cran.us.r-project.org')
  library(biomod2)
}
if (!require(doParallel)) {
  install.packages('doParallel',dependencies=TRUE,repos='http://cran.us.r-project.org')
  library(doParallel)       
}

##BIOMOD_FormatingData
message('Reading coordinates...')
coordinates <- SpatialPoints(read.table(opt$file, sep='\t'))
message('Done!')
message('')

message('Reading climatic variables for current period...')
calibperiod <- list.files(paste0(opt$dir,'/'),pattern=c('.asc','.grd','.tif'))
current <- stack(paste0(opt$dir,'/',calibperiod))
layers <- paste0('var',1:NROW(calibperiod))
names(current)<-layers
message('Done!')
message('')

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
    files <- list.files(paste0(opt$projection,'/',projperiods[i],'/'),pattern=c('.asc','.grd','.tif'))
    projvar[[i]] <- stack(paste0(opt$projection,'/',projperiods[i],'/',files))
    names(projvar[[i]]) <- layers
    message('Done!')
    message('')
  }
}

if (opt$strategy=='random') {
  dados <- BIOMOD_FormatingData(resp.var = coordinates, expl.var = current, resp.name = opt$title, PA.nb.rep = 10,
                             PA.nb.absences = as.numeric(opt$panumber), PA.strategy = opt$strategy)
} else if (opt$strategy=='disk') {
  dados <- BIOMOD_FormatingData(resp.var = coordinates, expl.var = current, resp.name = opt$title, PA.nb.rep = 10,
                            PA.nb.absences = as.numeric(opt$panumber), PA.strategy = opt$strategy,
                            PA.dist.min=as.numeric(opt$min),PA.dist.max = as.numeric(opt$max))
} else {
  dados <- BIOMOD_FormatingData(resp.var = coordinates, expl.var = current, resp.name = opt$title, PA.nb.rep = 10,
                             PA.nb.absences = as.numeric(opt$panumber), PA.strategy = opt$strategy,
                             PA.sre.quant=as.numeric(opt$sre-proportion))
}

options <- BIOMOD_ModelingOptions()

##Calibrando o modelo
message('Running model calibration')
message('')

model <- BIOMOD_Modeling(data = dados, models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips'),
                      models.options = options, DataSplit = 75, models.eval.meth = 'TSS', SaveObj = TRUE,
                      do.full.models= TRUE, modeling.id= opt$title)
dir.create(paste0('rastertemp_',opt$title))
rasterOptions(tmpdir = paste0('rastertemp_',opt$title))

current.projection <- BIOMOD_Projection(modeling.output=model, new.env = current,
                                     proj.name='current.projection', selected.models="all",
                                     binary.meth= 'TSS', output.format=".grd", do.stack=FALSE, build.clamping.mask=TRUE)

if (!is.null(opt$projection)) {
  projmodels <- list()
  for (i in 1:length(projvar)) {
    message('Running projection of model for directory ',projperiods[i])
    message('')
    projmodels[[i]] <- BIOMOD_Projection(modeling.output=model, new.env = projvar[[i]],
                                         proj.name=paste0('period_',projperiods[i],'.projection'), selected.models="all",
                                         binary.meth= 'TSS', output.format=".grd", do.stack=FALSE, build.clamping.mask=TRUE)
  }
}

##BIOMOD_EnsembleModeling (biomod2)
message('Running calibration of ensembled model')
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
    message('Running projection of ensembled model for directory ',projperiods[i])
    message('')
    projensemble[[i]] <- BIOMOD_EnsembleForecasting(ensemble.model, projection.output = projmodels[[i]],
                                                    selected.models = 'all', proj.name=paste('period_',projperiods[i],'.ensembled',sep=''),
                                                    binary.meth= 'TSS')
  }
}

unlink(paste0('rastertemp_',opt$title), recursive = TRUE)
