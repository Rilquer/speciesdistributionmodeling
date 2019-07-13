# Species Distribution Modeling

distmodeling: Species Distribution Modeling using BIOMOD package

This is a script to perform SDM using as input a coordinates file, the folder to climatic variables and a few other specific settings. This script uses the packaged BIOMOD2 (Thuiller et al., 2009) and implements all ten algorithms available in the package.

This script requires certain R packages (which will be installed when the script is run, if not already present). Additionallty, you must inform the number of pseudo-absence points to be used and the strategy to select those pseudo-absences. Finally, it requires files maxent.jar and maxent.bat in the same directory the script is run. More info can be seen in the Biomod help files (specifically acessed by using ??BIOMOD_FormatingData in R).

Bioclim variables for current and past periods (mid-Holocene, Last Glaciam Maximum and Last Interglacial) can be found in the link https://www.dropbox.com/s/jbu19octwma3l27/bioclim.zip?dl=0. The three AOGCMs CCSM4, MIROC-ESM and MPI-ESM-P were combined into a single bioclimatic dataset, calculating the mean across these three models for each variable, for mid-Holocene and Last Glacial Maximum periods, using the R package _raster_.

Two versions of the scripr are available. _distmodeling_ can be called directly (i.e., Rscript distmodeling) in the folder it is located, or copied to the bin folder to be called as an executable from any location. In this case, the options below need to be provided. The script _distmodeling\_manual_ is proved to allow manual editing of all the inner options of biomod2. In this case, the options need to be provided within the script (check comments in script).

	[OPTIONS]

	-t CHARACTER, --title=CHARACTER
		Name for this modeling analysis.

	-f CHARACTER, --file=CHARACTER
		File with coordinates in long-lat format, tab separated and without header

	-d CHARACTER, --dir=CHARACTER
		Path to the directory with climatic variables from current period, to be used in model calibration. Accepted formats are .asc, .grd or .tif

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
		Similarity proportion for strategy 'sre' (between 0.25 and 0.5)

	-h, --help
		Show this help message and exit
		
Usage example:
Modeling a bird species, which coordinates are in the file "cyanocorax_cyanopogon.txt". Variables layers for current conditions are in the folder "current". Past projections will be made for Last Glacial Maximum and Last Interglacial; layers for each period are separated in the folder "LGM" and "LIG", which are in turn inside the folder "past". 1000 pseudo-absence points will be sampled, using a random strategy.

	Rscript distmodeling -t crow_modeling -f cyanocorax_cyanopogon.txt -d current/ -p past -n 1000

### References
Thuiller, W., Lafourcade, B., Engler, R., & Araújo, M. B. (2009). BIOMOD–a platform for ensemble forecasting of species distributions. Ecography, 32(3), 369-373.
