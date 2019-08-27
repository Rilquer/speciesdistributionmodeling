# Species Distribution Modeling

distmodeling.R: Species Distribution Modeling using BIOMOD package

This is a script to perform SDM using as input a list of species to be modeled, the folder to climatic variables for calibration and a few other specific parameters. This script uses the packaged BIOMOD2 (Thuiller et al., 2009) and implements all ten algorithms available in the package.

This script requires certain R packages (which will be installed when the script is run, if not already present).

The main input a list of species to be modeled. Using this list, the script searches gbif database for species coordinates (using
R package rgbif), filters gbif errors (using R package CoordinateCleaner), filters localities geographically too close, by
a minimum distance set by the user (using custom-made functions) and then proceeds to modeling with R package biomod2.

The script also requires a path to the folder containing climatic variables for the current period (for calibration), as well as the number of pseudo-absence points to be used and the strategy to select those pseudo-absences. Finally, it requires files maxent.jar and maxent.bat in the same directory the script is run (these files are provided along with the script in the github pag).

More info can be seen in the Biomod help files (specifically acessed by using ??BIOMOD_FormatingData in R).

Bioclim variables for current and past periods (mid-Holocene, Last Glaciam Maximum and Last Interglacial) can be found in the link https://www.dropbox.com/s/jbu19octwma3l27/bioclim.zip?dl=0. The three AOGCMs CCSM4, MIROC-ESM and MPI-ESM-P were combined into a single bioclimatic dataset, calculating the mean across these three models for each variable, for mid-Holocene and Last Glacial Maximum periods, using the R package _raster_.

Another source of bioclimatic variables data for species modeling are http://paleoclim.org/ and http://ecoclimate.org/. Note that variables in PaleoClim are already in tif format, whereas variables in ecoClimate need to first be rasterized and saved in tif format to be used in this script.

Two versions of the script are available. _distmodeling.R_ can be called directly (i.e., Rscript distmodeling.R) in the folder it is located, or copied to the bin folder to be called as an executable from any location. In this case, the options below need to be provided. The script _distmodeling\_manual_ is proved to run in R GUI, allowing manual editing of all the inner options in the script. In this case, the options need to be provided within the script (check comments within the script).

	[OPTIONS]

	-i CHARACTER, --species=CHARACTER
		Path to file with the list of species to be analyzed. File must have one species name per line. Mandatory.

	-b CHARACTER, --biomes=CHARACTER
		Complete path to directory containing shapefile to be used when plotting species coordinates. This plot is used to show the final set of coordinates used for modeling each species. MUST contain the name of the layer (name before .shp extension). Mandatory.

	-f CHARACTER, --fdist=CHARACTER
		Minimum distance for filtering geographically close localities. Only localities away from each other more than the distance set here will be retained for analyses. Distance must be provided in meters. [default= 50000].

	-t CHARACTER, --minpoints=CHARACTER
		Minimum number of geographic points to be used for modeling. Only species with a number of localities equal to or higher than the number set here (after error and geographic distance filtering) will be modeled. [default= 15].

	-d CHARACTER, --dir=CHARACTER
		Path to the directory with climatic variables from current period, to be used in model calibration.  Mandatory.

	-m CHARACTER, --format=CHARACTER
		Format of variables files (either asc, grd or tif). [default= tif]

	-p CHARACTER, --projection=CHARACTER
		Path to the directory with climatic variables for future or past projections (optional). MUST NOT have final slash '/' in the end of path. Climatic variables for different periods must be separated in different subfolders inside the folder. The name of the subfolders will be the name used for the projections outputed by Biomod for the respective period.

	-n CHARACTER, --panumber=CHARACTER
		Number of pseudo-asbences to be selected. [default= 1000]

	-s CHARACTER, --strategy=CHARACTER
		Strategy for selection of pseudo-absence points. Either 'random', 'disk' or 'sre'. Random strategy samples pseudo-absence from all the background environment
              in the climatic variables provided to the script. Disk strategy requires you to determine a minimum and maximum distance around occurence points within which pseudo-absence points will
              be sampled (this is useful to avoid sampling in cells to close to occurence points, which can result in overfit, and to secure inclusion of similar abiotic conditions out of the range of
              the focal species). Sre strategy first creates a bioclimatic envelope of the occurence points and samples pseudo-absence points from areas which differs from this envelope at least by
              a proportion set by the user. As the biomod help informs, this can lead to overoptimistic models. [default= random]

	--min=CHARACTER
		Minimum distance for strategy 'disk', in meters

	--max=CHARACTER
		Maximum distance for strategy 'disk', in meters

	--sre-proportion=CHARACTER
		Similarity proportion for strategy 'sre' (between 0.25 and 0.5). [default= 0.5]

	-h, --help
		Show this help message and exit
		
Usage example:
Modeling two bird species, which names are in the file "species.txt". The shapefile to plot coordinates is in the path '/Users/rilquer/shapefile/'. Variables layers are in asc format. Variables for current period are in the folder "current". Past projections will be made for Last Glacial Maximum and Last Interglacial; layers for each period are separated in the folder "LGM" and "LIG", which are in turn inside the folder "past". A total of 1000 pseudo-absence points will be sampled, using a random strategy.

	Rscript distmodeling.R -i species -b /Users/rilquer/shapefile/ -m asc -d current -p past -n 1000

### References
Thuiller, W., Lafourcade, B., Engler, R., & Araújo, M. B. (2009). BIOMOD–a platform for ensemble forecasting of species distributions. Ecography, 32(3), 369-373.
