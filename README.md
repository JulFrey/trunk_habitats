# trunk_habitats

by Julian Frey
2025-01-31

R scripts to characterize tree trunks according to their habitat structures for lichens and mosses.
This collection of functions can be used on already segmented SfM, or LiDAR point clouds of tree stems
to calculate the surface area, volume, and geometric features, such as the average curvature.  

If you want to use it, you need a compiler first. 
* On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
* On Mac, install Xcode from the app store.
* On Linux, `sudo apt-get install r-base-dev` or similar.

Second you'll need to install the required packages:
`install.packages('concaveman', 'lidR','conicfit','sf','Rcpp')`

The script ends with an example of running the functions on multiple point clouds.
All functions are documented in roxygen style.
