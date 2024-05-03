####
### written by Dr. Rebekka Allgayer, University of Aberdeen, April 2024
#This script is for generating habitat layers for MerMADE inputs

#libraries
library(raster)
library(rgl)
library(igraph)

#You need:
#A .txt file for habitat types in your seascapes(this file is 2D)
#If you are running a patch-based model, you need to have your patch designations
  #in a .txt file as well
#You need the .txt file of seafloor depths
#NB: these files should have the SAME DIMENSIONS! same number of rows and columns

#set your working directory
setwd("C:/Users/Rey/Documents/MerMADE/Pipeline/") #!!! this is the filepath to the Pipeline folder, so change this to wherever you have that folder saved
#to use this on your own files, change it to where you have your Simulation folders saved
#!!! Change this to FALSE if you are NOT running a patch-based model
patch_based=TRUE

#read in the .txt file for habitat types in your seascape
#this is the file with 1,2,3 etc for each cell of your seascape
#remember that habitat type 0 is always water!
habs<- read.table("Data/hab_types.txt") #!!!if your file is called something else, change this

#IF YOU ARE RUNNING A PATCH BASED MODEL
#read in the .txt file for patch designations in your seascape 
#there should be -9999 in the cells that are not part of a patch
if(patch_based==T){
  patches<- read.table("Data/patches.txt") #!!! if your patch designation file is called something else, change this
}

#read in the seafloor depth data
sfd<- read.table("Data/sf_depth.txt") #!!! if your seafloor depth file is called something else, change this

#check that patches and habtypes match
#there should be no cell that has a patch designation but no habitat type
#also, there should be no cell designated as open water that has a patch number

suitable_habs<-c(1) #!!! you need to add all the habitat types that are suitable and could be part of patches
                    #so this might look like c(1,2,5) for example
if(patch_based==T){
 for(r in 1:nrow(habs)){
  for(c in 1:ncol(habs)){
    if(patches[r,c]!=-9999 && is.na(match(habs[r,c], suitable_habs))){ #if there's a patch designation but it's not suitable
      print(paste("there is a mismatch at ", r, ", ", c, sep=""))
    }
   }
  } 
}


#visualise your landscape in 2D
#seascape parameters (needed for mapping)
#!!!
res=1500 #resolution in metres
min_x=428430 #this is the minimum x coordinate IN METRES!
min_y=6123500 #this is the minimum y coordinate IN METRES!
nr=245 #number of rows
nc=198 #number of columns
#plot the seascape on a 2D map
max_y=min_y+(res*nr)
max_x=min_x+(res*nc)
hab_ras<- raster(as.matrix(habs))
extent(hab_ras)<- c(min_x, max_x, min_y, max_y)
nhabs <- c(0,unique(hab_ras)[-1])
plot(hab_ras, breaks=nhabs, col=terrain.colors(length(nhabs)+1),legend=F)

#depending on your depth intervals, the depth measures will not line up properly
#so we want to split the depths into layers

#figure out sfd's layers first
depth_layers=sfd #this will hold info on which layer it sits in
depth_interval=10 #!!! change this to the interval you are using in metres!

#because of the method of interpolation of hydrodynamic layers in MATLAB,
#there might be depth measurements even where there is land
land_type<- 2 #!!! change this value to whichever habitat type you have assigned to land!

for( r in 1:nrow(sfd)){
  for(c in 1:ncol(sfd)){
    if(habs[r,c]==land_type){ #if that cell is actually land
      depth_layers[r,c]=0 #make the depth 0
      next
    }
    
    la=floor(sfd[r,c]/depth_interval)
    
    #i want to assign habitat based on nearest to depth interval, so if seafloor depth is 29m, i dont want to block out the entire 20-30m depth layer
    half_depth=(la*depth_interval)+depth_interval/2; #min depth+half the interval
    if(sfd[r,c]>= half_depth){ #if the seafloor depth is deeper than half-depth
      #for example, if sfd=26, i want the 20-30interval to be water
      depth_layers[r,c]=la+1 #the value will belong to the next layer
    }
    else{
      #ie sfd=22, then i want the 20-30 interval to be solid
      depth_layers[r,c]=la
    }
  }
}


seafloor=2; #!!! designate a hab type for seafloor/land (this can be the same type as unsuitable if you dont want your indivs to settle there)
max_depth=90 #!!! change this to whatever maximum depth you want to model to!
max_lay = (max_depth/depth_interval)-1

for(la in 0:max_lay){
  copy_h=habs
  copy_p=patches
  for(r in 1:nrow(sfd)){
    for(c in 1:ncol(sfd)){
      if(la >= depth_layers[r,c]){ #if it should be solid
        if(habs[r,c]==-9999){#but there is no designation
          copy_h[r,c]=seafloor;
          copy_p[r,c]=-9999; #unsuitable hab, outside of patch
        }
        else{#otherwise make it that habitat/patch
          #since copy_h and copy_p are already copies, they already hold the hab and patch values
          #dont need to do anything here
        }
      }
      
      else{
        copy_h[r,c]=0; #otherwise, give it water because it's above the seafloor
        copy_p[r,c]=-9999; #since it's water, there is no patch designation
      }
    }
  }
  #!!! change the filenames to wherever  you want them to save. Replace "Sims/Test1/" 
  #with whichever folder you want 
  #Keep the "Inputs/hab_" part!
  write.table(copy_h, file=paste("Sims/Test1/Inputs/hab_", la, ".txt", sep=""), row.names=F, col.names=F)
  write.table(copy_p, file=paste("Sims/Test1/Inputs/patch_", la, ".txt", sep=""), row.names=F, col.names=F)
  print(paste("done with layer ", la, sep=""))
}
# 


###see the seascape in 3D
#!!! change this information to match your seascape
n_layers=9 #how many layer files do you have?
min_depth=0 #surface=0, otherwise which minimum depth to plot from
dint=10 #depth interval in metres
res=1500 #resolution in metres
min_x=428430 #this is the minimum x coordinate IN METRES!
min_y=6123500 #this is the minimum y coordinate IN METRES!
nr=245 #number of rows
nc=198 #number of columns
dir= "C:/Users/Rey/Documents/MerMADE/Pipeline/Sims/Test" #this is where you have the data saved, in the Inputs folder of a simulation
test_num=1 #if you have set up your sims as Test1, Test2 etc, this lets you loop over the sim folders


usrm<- par3d()$userMatrix
plot_results<- function(dir, test_num, n_layers, min_depth, dint, res, min_x, min_y, nr, nc, usrmatrix){
  
  max_depth=-1* (min_depth + n_layers*dint)
  
  layer_names<-unlist(read.table(paste(dir, test_num, "/Inputs/hab_layers.txt", sep="")), use.names=F)
  
  
  layers<- read.table(paste(dir, test_num, "/Inputs/", layer_names[1], sep=""))
  layers<- layers[1:nr, 1:nc]
  for(r in 1:nrow(layers)){
    for(c in 1:ncol(layers)){
      if(layers[r,c]>0){layers[r,c]=1}
    }
  }
  for(la in 2:(n_layers)){ 
    layer<- read.table(paste(dir, test_num, "/Inputs/", layer_names[la], sep=""))
    layer<- layer[1:nr, 1:nc]
    for(r in 1:nrow(layer)){
      for(c in 1:ncol(layer)){
        if(layer[r,c]>0){layer[r,c]=1}
      }
    }
    layers<-layers+layer
  }
  surface_test<- matrix(nrow=nrow(layers), ncol=ncol(layers))
  for(r in 1:nrow(layers)){
    for(c in 1: ncol(layers)){
      surface_test[r,c]=max_depth+(layers[r,c]*dint) #it's + because everything is on a negative axis
      
    }
  }
  
  surface_test<- as.matrix(surface_test)
  
  mx_vals<- seq(min_x+res,(min_x+res*nc),res)
  mx<- matrix(nrow=nrow(surface_test), ncol=ncol(surface_test))
  
  my_vals<- seq((min_y+nr*res),min_y+res,(-1*res))
  my<-matrix(nrow=nrow(surface_test), ncol=ncol(surface_test))
  
  for(i in 1:length(my_vals)){
    mx[i,]<- mx_vals
  }
  for(i in 1:length(mx_vals)){
    my[,i]<- my_vals
  }
  open3d()
  plot3d(0, 0, 0, type="l", xlab="x", ylab="y", zlab="z", xlim=c(min_x,min_x+res*nr), ylim=c(min_y,min_y+res*nr), zlim=c(min_depth,max_depth))
  
  z=surface_test
  surface_layer=surface_test
  
  colorlut <- terrain.colors(max(z)-min(z) +1) # height color lookup table
  
  col <- colorlut[ z - min(z) + 1 ] # assign colors to heights for each point
  
  surface3d(x=mx, y=my, z=surface_test, color=col)
  aspect3d(1,1,.08)
  view3d(userMatrix=usrm, zoom=.5) #zoom 0.5 is a good value for seeing it not from above
  
}

#plot your 3D landscape
plot_results(dir, test_num, n_layers, min_depth, dint, res, min_x, min_y, nr, nc, usrmatrix)
#this step might take a while, just be patient

#it will open a new little plotting window with a 3D rendering of your seascape
#you can zoom in, drag it around and change the perspective on it. Have a play!

