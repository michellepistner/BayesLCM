##Sets up the appropriate directories
##Change base.dir to change major directory. Everything else is relative
##All directories are defined off of the current working directory
##If another directory is desired, change here.
##Assumes under the parent folder that there is a data, results, and packages folder

base.dir = getwd()

paste("The base directory is ", base.dir ,sep="")

data.dir = paste(base.dir,"/data",sep="")
paste("The data directory is ", base.dir ,sep="")

results.dir = paste(base.dir,"/results",sep="")
paste("The results directory is ", base.dir ,sep="")

packages.dir = paste(base.dir,"/packages",sep="")
paste("The packages directory is ", base.dir ,sep="")
