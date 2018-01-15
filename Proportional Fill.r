#### Background geographic sampling occupancy marine invertebrates ####
# Calculate the total number of possible 1 degree by 1 degree paleogeographic grid
# squares for marine invertebrates through time and also largest possible linear
# distance and geographic area
# iteration is on each stage
library(geosphere)
library(functional)

Stages <- read.csv('GTS 2012 ages Grossman ch.csv', header = TRUE)
Stage <- Stages$Stage


NumLocations <- vector()
RoundLocations <- vector()
MaxDist <- vector()
MaxArea <- vector()

allInvertTotal <- read.table('https://paleobiodb.org/data1.2/occs/list.csv?interval=phanerozoic&envtype=marine&show=paleoloc',
                             header=TRUE,
                             sep=',')

for(i in 1:(length(Stage))){
  
  link <- paste('https://paleobiodb.org/data1.2/occs/list.csv?interval=', Stage[i], '&envtype=marine&show=paleoloc', sep = "")
  
  all.invert<-read.table(link,
                           header=TRUE,
                           sep=',')
  
  PLatRound <- round(all.invert$paleolat)
  PLonRound <- round(all.invert$paleolng)
  PGeoComb <- as.factor(paste(PLatRound, PLonRound, sep = "_"))
  
  NumLocations <- append(NumLocations, length(levels(PGeoComb)))
  
  PGeoComb <- as.factor(paste(levels(PGeoComb), rep(Stage[i], times = length(levels(PGeoComb))), sep = "_"))
  RoundLocations <- append(RoundLocations, levels(PGeoComb))
  
  locations <- as.matrix(cbind(all.invert$paleolng,all.invert$paleolat))
  locations <- locations[apply(locations, 1, Compose(is.finite, all)),]
  #MaxDist<-append(MaxDist,(max(distm(locations,locations,distCosine))/1000))
  
  MaxArea<-append(MaxArea,(areaPolygon(locations)/1000000))

}

plot(Stages$Base.Age[1:length(MaxArea)], MaxArea, type = 'o', pch = 15)
plot(Stages$Base.Age[1:length(MaxArea)], NumLocations, type = 'o', pch = 15)
plot(MaxArea, NumLocations)

AllLocations <- read.table(RoundLocations, sep = '_')

ArmAge <- seq(from = 495, to = 10, by = -10)
hist(locations)
plot(MidAge, locations)
plot(duration, locations)
Baseline <- as.data.frame(cbind(MidAge, ArmAge, duration,locations))
#Baseline <- cbind(Baseline, Stages)
Baseline <- Baseline[order(ArmAge),]
plot(Baseline$ArmAge, Baseline$locations, type = 'o', pch = 15, ylab = 'Number of marine 1x1 in PBDB', xlab = 'Age (Ma)')
abline(v = 252, col = 'red')
abline(v = 66, col = 'green')
Baseline$Stages[Baseline$locations==1]

Totalbaseline <- Baseline

plot(Totalbaseline$MidAge, Totalbaseline$locations, type = 'o', pch = 15, ylab = 'Number of marine 1x1 in PBDB', xlab = 'Age (Ma)',
     ylim = c(0, 800))
points(MidAge, locations, pch = 15, col = 'purple', type = 'o')

plot(MidAge, ((locations/Totalbaseline$locations[Totalbaseline$MidAge<395&Totalbaseline$MidAge>60])*100),
     ylab = '% Ammonoid coverage')

#### Ammonoid location calculations ####

AmmoniteMorpTable <- read.csv('Classificationmod.csv', header = TRUE)
AmmNumLocations <- vector()
AmmRoundLocations <- vector()
AmmMaxDist <- vector()
AmmMaxArea <- vector()

AmmStages <- Stages[Stages$Base.Age>64&Stages$Base.Age<412,]
AmmStage <- as.character(AmmStages$Stage)

amm.total<-read.table('https://paleobiodb.org/data1.2/occs/list.csv?base_name=Ammonoidea&envtype=marine&show=paleoloc',
                      header=TRUE,
                      sep=',')

for(i in 1:(length(AmmStage))){
  
  all.amm <- amm.total[amm.total$early_interval==AmmStage[i],]
  
  PLatRound <- round(all.amm$paleolat)
  PLonRound <- round(all.amm$paleolng)
  PGeoComb <- as.factor(paste(PLatRound, PLonRound, sep = "_"))
  
  AmmNumLocations <- append(AmmNumLocations, length(levels(PGeoComb)))
  
  PGeoComb <- as.factor(paste(levels(PGeoComb), rep(Stage[i], times = length(levels(PGeoComb))), sep = "_"))
  AmmRoundLocations <- append(AmmRoundLocations, levels(PGeoComb))
  
  locations <- as.matrix(cbind(all.amm$paleolng,all.amm$paleolat))
  locations <- locations[apply(locations, 1, Compose(is.finite, all)),]
  #MaxDist<-append(MaxDist,(max(distm(locations,locations,distCosine))/1000))
  
  AmmMaxArea<-append(AmmMaxArea,(areaPolygon(locations)/1000000))
  
}

plot(AmmStages$Base.Age, AmmMaxArea)
plot(Stages$Base.Age[1:length(MaxArea)], MaxArea, type = 'o', pch = 15, ylim = c(0,5*10^8))
plot(AmmStages$Base.Age, AmmMaxArea, type = 'o', pch = 17, col = 'red')

plot(Stages$Base.Age[1:length(NumLocations)], 
     NumLocations, 
     type = 'o', 
     pch = 15, 
     ylim = c(0,650), 
     xlim = c(550,0),
     xlab = 'Age (Ma)')

points(AmmStages$Base.Age, AmmNumLocations[1:53], type = 'o', pch = 17, col = 'red')

NormalizedLocations <- AmmNumLocations[1:53]/NumLocations[Stages$Base.Age>64&Stages$Base.Age<412]

plot(AmmStages$Base.Age,NormalizedLocations*100,
     type = 'o',
     pch = 15,
     xlim = c(550,0),
     xlab = 'Age (Ma)')

abline(v=c(252, 200, 66), lwd = 2, col = 'red')
abline(h=mean(NormalizedLocations*100))

plot(AmmStages$Base.Age,AmmMaxArea,
     type = 'o',
     pch = 15,
     xlim = c(550,0),
     xlab = 'Age (Ma)',
     ylab = 'Area Km^2')
