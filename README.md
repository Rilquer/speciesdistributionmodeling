# Species Distribution Modeling with BIOMOD

This is a script to perform SDM using as input a coordinates file, the folder to climatic variables and a few other specific settings. This script uses the packaged BIOMOD2 (Thuiller et al., 2009) and implements all ten algorithms available in the package.

This script requires certain R packages (which will be installed when the script is run, if not already present). Additionallty, you must inform the number of pseudo-absence points to be used and the strategy to select those pseudo-absences. Finally, it requires files maxent.jar and maxent.bat in the same directory the script is run. More info can be seen in the Biomod help files (specifically acessed by using ??BIOMOD_FormatingData in R).

The script can be directly called (i.e., Rscript ./script_modeling) in the folder it is located, or copied to the bin folder to be called as an executable from any location.

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

### References
Thuiller, W., Lafourcade, B., Engler, R., & Araújo, M. B. (2009). BIOMOD–a platform for ensemble forecasting of species distributions. Ecography, 32(3), 369-373.
