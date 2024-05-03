##this script will visualise the movement tracks for you on a 2D seascape
#and run some basic analysis of settlement sites

#required libraries
library(rgl)
library(igraph)
library(raster)

#set your working directory
setwd("C:/Users/Rey/Documents/MerMADEsoftware/Pipeline/") #!!! this is the filepath to the Pipeline folder, so change this to wherever you have that folder saved

#seascape parameters (needed for mapping)
n_layers=9 #how many layer files do you have?
min_depth=0 #surface=0, otherwise which minimum depth to plot from
dint=10 #depth interval in metres
res=1500 #resolution in metres
min_x=428430 #this is the minimum x coordinate IN METRES!
min_y=6123500 #this is the minimum y coordinate IN METRES!
nr=245 #number of rows
nc=198 #number of columns
dir= "Sims/Test" #you don't need to change this since it's relative to being in the Pipeline folder
test_num=1 #!!! this is the simulation that you are wanting to visualise, so this would be concatenated with dir to become Sims/Test1, for example

##################################################
#plot the seascape on a 2D map
flat_hab<- read.table(paste(dir, test_num,"/flat_map.txt", sep="")) #every simulation folder has this file, it's just the habitat map in 2D
max_y=min_y+(res*nr)
max_x=min_x+(res*nc)
hab_ras<- raster(as.matrix(flat_hab))
extent(hab_ras)<- c(min_x, max_x, min_y, max_y)
plot(hab_ras)

#plot the tracks
track_file<- paste(dir, test_num, "/Outputs/Indiv_Mov_rep0.txt", sep="")
no_col<- max(count.fields(track_file, sep="\t"))
indiv_mov<- read.table(track_file, skip=1, sep="\t", fill=T, col.names=1:no_col) #this can take a minute to read in because the filesize gets big

nindivs=(unique(indiv_mov[,4]))#creates a list of the unique IDs
length(nindivs) #look at how many individuals you have, then decide what proportion of them you want to plot
indiv_prop=0.05 #!!! if you plotted all, it takes forever and the many tracks obscure the patterns
#for the simulations where you've released individuals from the entire map, 5% still plots over 1k indivs
#for the simulations where you've released individuals from an enhancement site, you might all because they're not that big!

for(i in 1:length(nindivs)){
  to_plot=runif(1)
  if(to_plot<indiv_prop){
    indiv_data<- indiv_mov[indiv_mov[,4]==nindivs[i],]
    
    x<- indiv_data[1,-(1:7)][!is.na(indiv_data[1,-(1:7)])]
    
    y<- indiv_data[2,-(1:7)][!is.na(indiv_data[2,-(1:7)])]
    
    lines(x, y, col="red") #!!! you can change the colour as you like
    
  }
}

###################################################
#plot sites of settlement of COMPETENT individuals
plot(hab_ras) #plot a clean map
comp_data<- indiv_mov[indiv_mov[,7]==1,] #competent individuals
stat_data<- comp_data[comp_data[,6]==3,] #this is the status of settled individuals
nindivs=(unique(stat_data[,4]))
length(nindivs) #check how many settled and then decide your proportion to plot
indiv_prop=0.5 #!!! this will plot half the tracks
for(i in 1:length(nindivs)){
  to_plot=runif(1)
  if(to_plot<indiv_prop){
  indiv_data<- stat_data[stat_data[,4]==nindivs[i],]
  if(indiv_data[1,6]!=3){next}
  
  x<- indiv_data[1,-(1:7)][!is.na(indiv_data[1,-(1:7)])]
  
  y<- indiv_data[2,-(1:7)][!is.na(indiv_data[2,-(1:7)])]
  
  
  points(x[length(x)-1], y[length(y)-1], col="red", pch=16) 
  
  }
}
