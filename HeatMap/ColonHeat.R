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
#male75[1,] = male75[1,]+50
male <- male75
male <- as.matrix(male)
male <- male[,50:86]
#heta is cases per hundred thousand
male[1,] = 199


#r = c("Cecum","Ascending","Transverse","Descending","Sigmoid","Rectum")
plot_ly(z=male,type="heatmap",x = seq(from= 50, to= 86) ,zmin= 0,zmax = 200 ) # %>% layout(title = "Cancer Incidance by Colon Location")



plotCont <- function(input){

#more advanced contour plot
best_fit <- Tps(PBCpos,PBCz,cost= 1.4)
grid.list<- list( x= seq(0,4,.1), y=seq(0,12,.1))   #Note - these ranges will need to be changed for your needs
xg<- make.surface.grid(grid.list)
dev.off()
dev.new()
f1<- predict( best_fit, xg)
out.p<- as.surface( xg, f1)
best_plot = plot.surface(out.p,zlim=c(14,40),xlab='Age',ylab='Part Of Colon',main='Colon',cex.main=1.3,cex.lab=1.3,cex.axis = .75, xaxt = 'n',col = myPal)
axis(side=2,at = c(1,2,3,4,5,6),labels = r)
}







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