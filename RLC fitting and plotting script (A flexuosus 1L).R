
# Script to fit rapid light curves using R package "phytotools" and the model by Platt et al. 1980, 
# to extract the coefficients/parameters alpha, beta, ETRmax, Ek and ps from the fitted curve, 
# The script follows the example provided by the authors of the phytotools package and ideas from
# Antti Takolander (2019)

# Written by Suleiman Dauda January, 2020.


rm (list=ls()) # To clear the workspace
# ctrl+l or cmd+l (in Mac) to clear the console


# Step 1 ##########################################################################
## To load the RLC data - having 3 columns (par, etr, id)

rlc.data <- read.csv("/Users/suleiman/My Drive/Phd research/RESULTS/1 Litre experiments/A. flexuosus/RLC A.flexuosus Day 3 R format.csv")

## read.table(file.choose(),header=TRUE,sep=",")


# If loaded data has replicates - continue to step 2, otherwise proceed to step 3


# Step 2 ##########################################################################

library(dplyr)

r <- rlc.data %>%                ## To get the means of each PAR value for every treatment/species
  group_by(par, id) %>% 
  summarise(etr = mean(etr, na.rm = TRUE)) %>% 
  arrange(id)
str(r)
View(r)

r$id <- as.factor(r$id) ### It is important to convert the Id character vector to factor;
                        ## without doing this, an error will be gotten when plotting multiple plots, due to the pch argument
str(r)
## Proceed to step 4!!!


# Step 3 ##########################################################################

rlc.data <- r
## Proceed to step 4!!!


# Step 4 ##########################################################################

library(phytotools)  ## To load the phytotools package

ncurves <- length(unique(r$id)) # number of unique ids in the data 
ids <- unique(r$id) # store the unique ids 


# create a data frame to store the extracted curve parameters after model fitting
rlc.parameters <- data.frame(
  id = ids, 
  alpha = 0, 
  beta = 0, 
  ETRmax = 0, 
  Ek = 0, 
  ps = 0
)



# Step 5 ##########################################################################
# Fitting the RLC parameters
for (i in 1:ncurves){
  
  temp.id = ids[i] # extract the id of the curve to be fitted
  
  temp.rlc.data <- r[r$id==temp.id,] # extract the the data of a single curve into a temporary variable
  PAR = temp.rlc.data$par 
  ETR = temp.rlc.data$etr
  
  fit = fitPGH(PAR, ETR, fitmethod = "Port") # curve fitting using the Platt et al. 1980 model
  
  # store the fitted RLC values into temporary variables
  alpha.rlc = fit$alpha[1]
  beta.rlc = fit$beta[1]
  ps.rlc = fit$ps[1]
  
  # store the parameters
  rlc.parameters$id[i] <- temp.id
  rlc.parameters$alpha[i] <- alpha.rlc
  rlc.parameters$beta[i] <- beta.rlc
  rlc.parameters$ps[i] <- ps.rlc
  
  # calculate ETRmax and Ek for the PGH model (see e.g.Ralph & Gademann 2005 Aquatic Botany 82 (3): 222 - 237). 
  
  
  ETRmax = ps.rlc*(alpha.rlc/(alpha.rlc + beta.rlc))*(beta.rlc/(alpha.rlc+beta.rlc))^(beta.rlc/alpha.rlc)
  Ek = ETRmax/alpha.rlc 
  
  # store the variables
  rlc.parameters$ETRmax[i] <- ETRmax
  rlc.parameters$Ek[i] <- Ek
}


# Step 6 ##########################################################################
# To get the fitted RLC parameters
rlc.parameters

# To save the fitted RLC parameters as a csv file in the working directory
write.csv(rlc.parameters, "RLC parameters (A.flexuosus 1L Exp. (Day 3))(03.06.2020).csv",row.names=FALSE) # Always change the no and date
                                                                            # inorder not to overwrite previous result     

# Step 7  ##########################################################################
range(r$etr) ## Run to know the limit of the y-axis

# Step 8 ###########################################################################
# To plot individual a .tiff graphs for each treatment/species and save in the working directory

for (i in 1:ncurves){
  
  temp.id = ids[i] 
  
  print(paste("Now fitting curve ", as.character(temp.id))) # to keep track what's happening if 
                                                            # the data has many curves
  
  temp.rlc.data <- r[r$id==temp.id,] 
  PAR = temp.rlc.data$par 
  ETR = temp.rlc.data$etr
  
  fit = fitPGH(PAR, ETR, fitmethod = "Port") 
  alpha.rlc = fit$alpha[1]
  beta.rlc = fit$beta[1]
  ps.rlc = fit$ps[1]
  
  tiff(file=paste0(temp.id, ".tiff"), width = 2000, height = 1500, pointsize= 15, # To save plot, edit "file name.tiff" to match
       restoreConsole = TRUE,units = "px", res = 300,  compression = "lzw")       # every new plot!
  
  par(xpd = T, mar = par()$mar + c(0,0.5,-1,1),tck = -0.021)
  
  # plot the data, 
  plot(x=PAR, y=ETR, main=temp.id,
       xlim = c(0,1200), 
       ylim = c(0,150),     ## set axis limit from the range of the etr obtained in step 
       xaxt = "n", yaxt = "n", cex.axis = 0.8, cex.lab = 1.2,
       ylab = (expression(paste("rETR"," ","(",mu,"mol"," ", "electrons"," ",m^-2," ",s^-1,")" ))),
       xlab = (expression(paste("PAR"," ","(",mu,"mol"," ", "photons"," ",m^-2," ",s^-1,")" ))))
  
  axis(side=1, seq(0, 1200, by = 200), las = 1) ## x-axis range
  axis(side=2, seq(0, 150, by =25), las =1) ## y-axis range
  
  # plot the model fit
  with(fit, {
    P <- ps.rlc*(1-exp(-1*alpha.rlc*PAR/ps.rlc))*exp(-1*beta.rlc*PAR/ps.rlc) # the PGH model equation
    lines(PAR,P)
  }
  ) # end of with
  dev.off() #close the plotting devide. if this is not done, the next run of the loop will override the plot. 
  
}



###################################################################################################################
### Multiple plots on one sheet

ncurves <- length(unique(r$id))

ids <- unique(r$id)  

# Step 
range (r$etr) ## Run to know the limit of the y-axis

## TIFF plot, 600 DPI or PPI

tiff(file = "RLC A. flexuosus 1L.tiff", width = 6500, height = 4500, pointsize= 25,       # To save plot, edit "file name.tiff" to match
     restoreConsole = TRUE,units = "px", res = 600,  compression = "lzw")  # every new plot!

par(xpd = T, bty="o", mar = par(mgp=c(3,0.7,0),mar=c(5,4,4,2)+0.1)$mar + c(-1.5,-0.5,-3,2),tck = -0.021)

#Par(mar=c(5,4,4,2))

plot (0,type = "n",xlim = c(0,1200), ## Make an empty plot
      ann = FALSE,
      ylim = c(0,180), 
      xaxt = "n", yaxt = "n", cex.axis = 1, cex.lab = 1,
      # ylab = (expression(paste("rETR"," ","(",mu,"mol"," ", "electrons"," ",m^-2," ",s^-1,")" ))),
      # xlab = (expression(paste("PAR"," ","(",mu,"mol"," ", "photons"," ",m^-2," ",s^-1,")" ))))
)
axis(side=1, cex.axis = 1, seq(0, 1200, by = 200), las = 1) ## x-axis range
axis(side=2, cex.axis = 1, seq(0, 180, by =30), las =1) ## y-axis range

mtext(side = 1, cex =1, (expression(paste("PAR"," ","(",mu,"mol"," ", "photons"," ",m^-2," ",s^-1,")" ))), line = 1.9)
mtext(side = 2, cex =1, text = (expression(paste("rETR"," ","(",mu,"mol"," ", "electrons"," ",m^-2," ",s^-1,")" ))), line = 1.9)
## To plot the data points and curve fits
for (i in 1:ncurves){
  
  emp.id = ids[i] # extract the id of the curve to be fitted
  
  print(paste("Now fitting curve ", as.character(emp.id))) # to keep track what's happening if the data has many curves
  
  emp.rlc.data <- r[r$id==emp.id,] # extract the the data of a single curve into a temporary variable
  PAR = emp.rlc.data$par 
  ETR = emp.rlc.data$etr
  
  bbfit = fitPGH(PAR, ETR, fitmethod = "Port") # curve fitting using the Platt et al. 1980 model
  
  # store the fitted RLC values into temporary variables
  alpha.rlc = bbfit$alpha[1]
  beta.rlc = bbfit$beta[1]
  ps.rlc = bbfit$ps[1]
  points(r$par, r$etr,col =  c("black"), bg = c("black","black","black","black","black","black","black"),
         pch = c(1,0,2,15,8,23,25)[(r$id)]) # edit the numbers (0-25) in the parenthesis to match. Needs the id to be factor or numeric
  # the number of curves (treatments/species) to plot.
  
  with(bbfit, {
    P <- ps.rlc*(1-exp(-1*alpha.rlc*PAR/ps.rlc))*exp(-1*beta.rlc*PAR/ps.rlc) # the PGH model equation
    lines(PAR,P)
  })
}

#legend(1260,100, c("A","B","C","D","E","F","G"),cex=0.7,pch = c(16,6, 8, 15,5,17,4), # Type legend labels corresponding to
#col= "black", lty = c(0,0),text.font = 1,title=expression(bold("Treatment")), # each curve (treatment/species)
#box.lty = 0, bg = "transparent") 

legend(1250,160, c("1.7","2.5","3.4","3.8","7.4","21.4","589.0"),cex=1,
       pt.bg =  c('black', "black", "black", "black", 'black', "black", 'black'),
       pt.cex = 1, # Expansion factor for the point only
       pch = c(1,0,2,15,8,23,25), # Type legend labels corresponding to
       lty = c(0,0),text.font = 1.3,
       title=expression(paste("Cu"^'2+'," ","(","nM",")")), # each curve (treatment/species)
       box.lty = 0, bg = "transparent")   


text(50, 170, expression(bold('a.')))
text(310, 170, expression(italic("A. flexuosus")))


dev.off()  




####################################################################################################
################################# REPLICATES CALCULATION ############################################
####################################################################################################

# Step 2 ##########################################################################

library(dplyr)

r <- rlc.data %>%                ## To get the means of each PAR value for every treatment/species
  group_by(par, id2) %>% 
  summarise(etr = mean(etr, na.rm = TRUE)) %>% 
  arrange(id2)
str(r)
View(r)

r$id2 <- as.factor(r$id2) ### It is important to convert the id2 character vector to factor;
## without doing this, an error will be gotten when plotting multiple plots, due to the pch argument
str(r)
## Proceed to step 4!!!


# Step 3 ##########################################################################

rlc.data <- r
## Proceed to step 4!!!


# Step 4 ##########################################################################

library(phytotools)  ## To load the phytotools package

ncurves <- length(unique(r$id2)) # number of unique ids in the data 
ids <- unique(r$id2) # store the unique ids 


# create a data frame to store the extracted curve parameters after model fitting
rlc.parameters <- data.frame(
  id2 = ids, 
  alpha = 0, 
  beta = 0, 
  ETRmax = 0, 
  Ek = 0, 
  ps = 0
)



# Step 5 ##########################################################################
# Fitting the RLC parameters
for (i in 1:ncurves){
  
  temp.id2 = ids[i] # extract the id2 of the curve to be fitted
  
  temp.rlc.data <- r[r$id2==temp.id2,] # extract the the data of a single curve into a temporary variable
  PAR = temp.rlc.data$par 
  ETR = temp.rlc.data$etr
  
  fit = fitPGH(PAR, ETR, fitmethod = "Port") # curve fitting using the Platt et al. 1980 model
  
  # store the fitted RLC values into temporary variables
  alpha.rlc = fit$alpha[1]
  beta.rlc = fit$beta[1]
  ps.rlc = fit$ps[1]
  
  # store the parameters
  rlc.parameters$id2[i] <- temp.id2
  rlc.parameters$alpha[i] <- alpha.rlc
  rlc.parameters$beta[i] <- beta.rlc
  rlc.parameters$ps[i] <- ps.rlc
  
  # calculate ETRmax and Ek for the PGH model (see e.g.Ralph & Gademann 2005 Aquatic Botany 82 (3): 222 - 237). 
  
  
  ETRmax = ps.rlc*(alpha.rlc/(alpha.rlc + beta.rlc))*(beta.rlc/(alpha.rlc+beta.rlc))^(beta.rlc/alpha.rlc)
  Ek = ETRmax/alpha.rlc 
  
  # store the variables
  rlc.parameters$ETRmax[i] <- ETRmax
  rlc.parameters$Ek[i] <- Ek
}


# Step 6 ##########################################################################
# To get the fitted RLC parameters
rlc.parameters

# To save the fitted RLC parameters as a csv file in the working directory
write.csv(rlc.parameters, "RLC parameters REPLICATE (A.flexuosus 1L Exp. (Day 3))(03.06.2020).csv",row.names=FALSE) # Always change the no and date
# inorder not to overwrite previous result     

