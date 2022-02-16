##ACS code master file
##Installs necessary packages

##All directories are defined off of the current working directory
##If another directory is desired, change here.
##Assumes under the parent folder that there is a data, results, and packages folder
source(config.R)

##Installing packages into the packages dir
install.packages(c("data.table","MCMCpack","plyr","MASS","ggplot2","synthpop","rlang","tidyr"),lib = packages.dir,repos='http://cran.us.r-project.org')


