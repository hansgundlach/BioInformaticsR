#Heatmaps for Colon
library("ggplot2")
library(ggplot2)
library(reshape2)
library(plotly)
install.packages("ddplyr")
#x axis is age

AllRace_Male_CRC_location_SEER9_1975.2013 <- read.csv("~/AllRace_Male_CRC_location_SEER9_1975-2013.csv", header=FALSE)
colset <- AllRace_Male_CRC_location_SEER9_1975.2013
male85 <- colset[349:359,2:87]
male95 <- colset[362:367,2:87]
male05 <- colset[375:380,2:87]
male75 <- colset[338:343,2:87]
male <- male05
male <- as.matrix(male)
male <- male[,40:86]
#heta is cases per hundred thousand
r = c("Cecum","Ascending","Transverse","Descending","Sigmoid","Rectum")
plot_ly(z=male,type="heatmap",x = seq(from= 40, to= 86),y = r ,zmin= 0,zmax = 200 ) # %>% layout(title = "Cancer Incidance by Colon Location")



heatmap(as.numeric(male))

bettdf <- grid.to.xyz(male)


















#melted_cormat <- melt(male)
#head(melted_cormat)
#ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
#  geom_tile()




#ggplot(data = male, aes(x=Var1, y=Var2, fill=value)) + 
#geom_tile()



#p <- ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = rescale),
 #               +     colour = "white") + scale_fill_gradient(low = "white",
                 #                                                                                  +     high = "steelblue"))