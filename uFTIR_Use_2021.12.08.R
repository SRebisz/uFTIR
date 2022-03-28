# https://github.com/fcorra/uFTIR


# remove.packages("uFTIR") #remove the uFTIR package
# remove.packages('Rcpp') #remove the rcpp package
.libPaths()
install.packages('Rcpp', dependencies = TRUE)

# I.  Installation (only for the 1st time) ####

  # * Instal Rtools 
  # https://cran.r-project.org/bin/windows/Rtools/rtools40.html 
    
  # *	Install "devtools" and "remotes"
    if(!require(devtools)){install.packages("devtools")}
    if(!require(remotes)){install.packages("remotes")}
    library(remotes)
    library(devtools)
   
  # * Install "uFTIR"
    remotes::install_github("fcorra/uFTIR") #, ref=dev ref=mod_plot
    remotes::install_github("fcorra/GDALPolygonize")
    library(uFTIR)
    library(GDALPolygonize) #package exporting shapefile
    packageVersion("uFTIR") 
    packageDescription("uFTIR")
    
  # *	Install other packages
    install.packages("tidyverse")
    install.packages( 'rgdal')
    install.packages("rgeos")
    install.packages("prospectr")
    install.packages("parallel")
    library(parallel)
    library(rgdal)
    library(rgeos)
    library(dplyr)
    library(prospectr)
    
# II. uFTIR Mosaic processing ####
  # *	Load packages 
    library(uFTIR)
    library(GDALPolygonize)
    library(rgdal)
    library(rgeos)
    library(dplyr)
    library(prospectr)
    
    rm(list=ls()) # cleaning console
    graphics.off() # cleaning plots
  
  # * Load reference spectra Library
    Library_ref=data.frame(Substance=primpke@substances,
                           Cluster_name=primpke@clusternames[primpke@clusterlist]  )
    
  # * Set your working directory: 
    wd="C:/Users/berio001/Documents/uFTIR/Ines" 
    setwd(wd)
  # * Set a sample name 
    File="35"
    src <- paste(wd,File,sep="/") 
    setwd(src)
    list.files(src)[10] #check the folder
  
  # * Load files .dmd   
    Files.dmd=list.files(src)
    Files.dmd=Files.dmd[grep(".dmd",Files.dmd)]
  
  # * mosaic_info .dmt
  x.mosa <- mosaic_info(dmtfile =paste(File,".dmt", sep = "")) #paste("2_", File,"_pm.dmt", sep = "") basic information about the image, does not load the measured spectra
  
  # *  spectral angle mapper algorithm SAM ####
  #Spectral angle mapper algorithm chunk by chunk for the mosaics files
  x.mosa@path <- system.time(mosaic_sam(x.mosa, primpke, derivative = NULL, 
                                   base_corr = TRUE,
                                   n_cores = NULL, temporal = TRUE,
                                   FUN = function(x){prospectr::savitzkyGolay(x, 1, 3, 7)}) ) #https://rdrr.io/cran/prospectr/man/savitzkyGolay.html
  
  # * Compose ####
  system.time(x.sam <-mosaic_compose(x.mosa@path, clusterlist = primpke@clusterlist, drop_raw=TRUE, nslices =2 )) # Load the SAM results in R
  
  # * Smooth ####
  system.time(x.smooth <- smooth_sam(x.sam ,nclusters =  32L, window = 3, nslices = 1 )) #Smooth SAM results  window = 3px min
  

  # * Clip #### 
    #clip the picture, needed to use th highlight
               # toClip(rad = ,segments =   , centre = c(x,y))  
    par(mfrow=c(1,1))
    plot(x.smooth, legend = FALSE)
  
    clip_mask <- toClip(500, 20, c(450,560) ) # FPE1
    
    polygon(clip_mask @xycoords)
    x.clip <- clipper(x.smooth, centre = clip_mask @centre, rad = clip_mask @rad)
    
    par(mfrow = c(1,3))
    plot(x.sam, legend = FALSE)
    plot(x.smooth, legend = FALSE)
    plot(x.clip, legend = FALSE)
  
  
  # * Export summary ####
    x.sum_clip <- summary_sam(x.smooth, mask = clip_mask , clusternames = primpke@clusternames, slice=1, smooth = FALSE, temporal = TRUE)
    
    x.sum_red=subset(x.sum_clip, x.sum_clip$area>10)
    B=x.sum_red %>% group_by(clname) %>% 
      summarise(c_number = mean(cluster), total_area = sum(area), n_particles = n())
    B
    
    setwd(wd)
    #write.csv(B, paste(File,"_Sumary2.csv",sep=""))
    #write.csv(x.sum_clip, paste(File,"_Results2.csv",sep=""))
  
  # * Check spectra ####
  # load the shape file in R
    
    #we used temporal=TRUE in the examples, so we can load back the files by
    x.shape <- rgdal::readOGR(paste(src, "/shape_out", sep = ""), "clusters" , verbose = F)
    
    #the you can do whatever you want with them. Here we subset and plot.
    x.polyethylene <- x.shape[x.shape@data[,1] == 30, ]
    x.polypropylene <- x.shape[x.shape@data[,1] == 4, ]
    x.plant.fibres  <- x.shape[x.shape@data[,1] == 14, ]
    
    par(mfrow = c(1,3))
    sp::plot(x.shape)
    sp::plot(x.plant.fibres, xlim = x.shape@bbox[1,], ylim = x.shape@bbox[2,])
    
    # you can overlay them too!
    sp::plot(x.shape)
    sp::plot(x.polyethylene, xlim = x.shape@bbox[1,], ylim = x.shape@bbox[2,], col ='red', add = T)
    sp::plot(x.polypropylene, xlim = x.shape@bbox[1,], ylim = x.shape@bbox[2,], col ='blue', add = T)# get the spectra with x.profile
    sp::plot(x.plant.fibres, xlim = x.shape@bbox[1,], ylim = x.shape@bbox[2,], col ='blue', add = T)
  
    # Plot the spectra
    
    x.prof <- get_profile_sinfo(x.smooth, 
                                   where = list("info" = x.mosa, "dmdfile" =  "G1_0003_0002.dmd"),  
                                   dst_cluster =4, plotpol = FALSE, plotpt = TRUE, cluster = TRUE, 
                                   clusternames = primpke@clusternames)
    
    highlight_substance(x.clip, 7) #highlight particules of cluster 
    
    
    Files.dmd
    
    par(mfrow = c(1,2), las = 1)
    cluster=11 #4 polypropylene "G1_0000_0000.dmd" # 11 polyester "G1_0003_0002.dmd"
    x.prof<- get_profile(x.smooth, 
                where = list('info' = x.mosa, 'dmdfile' = "G1_0003_0002.dmd"),
                dst_cluster = cluster, plotpol = TRUE, axes = F, box = F, legend = F) 
    
    plot(rnorm(1), type = "n", xlim = range(x.mosa@wavenumbers), 
         ylim = c(0, max(x.prof)), xlab = 'Wavenumbers', ylab = "Absorbance (%)")
    lines(colMeans(x.prof) ~ x.mosa@wavenumbers, col = "red", lty = 2)
    lines(colMeans(
      primpke@Spectra[primpke@clusterlist == cluster,]) ~ primpke@wavenumbers, 
      col = 'blue' )
    