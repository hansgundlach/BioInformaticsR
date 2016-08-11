mset <- read.csv("~/mset.csv", na.strings= c("",NA))
library("RColorBrewer")
library(plotly)
library(stringr)
library("fields")
library("ggplot2")
library("lattice")



#code only retrieves High Grade rows
#mset <- mset[grep("HGD",mset$indication),]

#only retrieve non-highgrade rows
#mset <- mset[-c(grep("HGD", mset$indication)),]

#only retrieve rows with tissue age
mset <- mset[c(grep("[0-9]",mset$tissue_age)),]

#removes samples that are normal from data frame or any other bad string
#mset <- mset[-c(grep("normal",mset$Sample.Clinical.Pathology)),]
better_dat <- mset$Tissue_location
#better dat contains all the values with length or centimeter data
better_dat <- mset[grep("cm", better_dat),]
clocker <- lapply(better_dat$Tissue_location,FUN = function(x) gsub(".+cm","",x))


#caculate distance fromGEJ junction
GEJ <- better_dat[19]
GEJ <- c(as.numeric(unlist(junction)))

#distone and disttwo are both 
dist <- lapply(better_dat$Tissue_location, FUN = function(x) gsub("cm.+","",x))
dist <- c(as.numeric(lapply(dist, FUN = function(x) gsub("[^0-9]","",x))))
dif <- GEJ - dist


clock <- c(unlist(clocker))
#replace 5:00 with 6:00 and organizes data into data frame
clock <- gsub("5:00", "6:00",clock)
dist <- c(dif)
onset <- mset$tissue_age
age <- mset$Age
dwell <- age -onset


df <- data.frame(clock,dist,dwell)
names(df) <- c("clock","dist","age")



#all the rows with clock info
#get methylation info

#clock_rows picks out all the indices with clock location
#X_coord and y_coord are justin the cases where there is clock data
clock_rows = grep(":", df$clock)

#puts clock levels in right order ie 12:00 -> 4 9:00 -> 3 etc
check <- factor(df[clock_rows, ]$clock)
levels(check) <- c(4,1,2,3)
x_coord = as.numeric(levels(check))[check]

#these are the y coordinates at least from R console
y_coord = as.numeric(df[clock_rows, ]$dist)
mether = as.numeric(df[clock_rows, ]$age)
esoph = matrix(rep(0,len=80), nrow = 20, ncol= 4)
weight = matrix(rep(0,len=80), nrow=20,ncol=4)


#adds position and heat data to esoph diagram 
for(vale in 1:length(x_coord)){
  
  #print(mether[vale])
  esoph[y_coord[vale],x_coord[vale]] = esoph[y_coord[vale],x_coord[vale]]+ mether[vale]
  weight[y_coord[vale],x_coord[vale]] =  weight[y_coord[vale],x_coord[vale]]+1
}
plot_ly(z=esoph, type="heatmap")

heatmap(weight)

# for loops adds data for biopsises without clock data
y_coord_no = as.numeric(df[-c(grep(":", df$clock)), ]$dist  )
meth_value = as.numeric(df[-c(grep(":", df$clock)), ]$age)
for(index in 1:length(y_coord_no)){
  # add .25 times val
  print(index)
  esoph[y_coord_no[index],]=  esoph[y_coord_no[index],] + meth_value[index]*.25
  weight[y_coord_no[index],]= weight[y_coord_no[index],] +.25
}

#computes weighted averages by dividing the sum of ages by the sum of weights in the weight matrix 
best_matrix <- esoph*weight^-1
best_matrix
#makes 2D heatplot 
plot_ly(z=best_matrix,x= c("3:00","6:00","9:00","12:00"), type="heatmap")

heatmap(na.omit(best_matrix))

#Part 2 -- smooth graph of data
#retireves the coordinates and z values of the parts of the matrix without NAs

#making copies of position to make peridoic boundary conditions
pos = which(!is.na(best_matrix),TRUE)
left  <-  pos
left[,2] = left[,2]-4 
right <-  pos
right[,2] = right[,2] +4
PBCpos = rbind(pos,right,left)
PBCpos = PBCpos[,c(2,1)]



pre_val = as.vector(best_matrix)
z_index = which(!is.na(pre_val), TRUE)
z_val <- pre_val[z_index]
rorm_mat = matrix(z_val,ncol = 1)
PBCz = rbind(rorm_mat,rorm_mat,rorm_mat)
best_fit <- Tps(PBCpos,PBCz,cost= 1.4)
grid.list<- list( x= seq(0,4,.1), y=seq(0,12,.1))   #Note - these ranges will need to be changed for your needs
xg<- make.surface.grid(grid.list)

#plotting function - modify as needed, including changing zlim to cover the range of y-values 
dev.off()
dev.new()
f1<- predict( best_fit, xg)
out.p<- as.surface( xg, f1)
best_plot = plot.surface(out.p,zlim=c(14,40),xlab='clock position',ylab='length in (cm) from GEJ',main='Esophagus',cex.main=1.3,cex.lab=1.3,cex.axis = .75, xaxt = 'n',col = myPal)
axis(side=1,at = c(1,2,3,4),labels = c("3:00","6:00","9:00","12:00"))
#the code below deals with plotting points on the graph

#high_clock_index is index of high grade biopsies with clock values
p_clock <-  mset[grep(":",mset$Tissue_location),]
high_clock_index <- grep("HGD",p_clock$indication)
# index of nondisplastic tissue with clock
nondist_index <- setdiff(seq(1:length(high_clock_index)), high_clock_index)      

#plot points that don't have clocks and are nondisplastic
#no_clock <- mset_2[grep("[^:]*", mset_2$Tissue_location),]
clock_index <- grep(":",  mset$Tissue_location)
#no_clock <- mset_2[setdiff(seq(1:length(mset_2$Tissue_location)),clock_index),]
no_clock <- mset[-c(clock_index),]

low_clock_index = grep("[^HGD]", no_clock$indication)
high_no_index <- grep("HGD",no_clock$indication)
high_no_clock <- y_coord_no[high_no_index]
low_no_clock <- y_coord_no[low_clock_index]


#tissue age for displastic tissue with clock
highclock_tissue <- p_clock[high_clock_index,]$Age -  p_clock[high_clock_index,]$tissue_age

#create colors 
#color is based on displastic tissue age
myPal <- colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "darkorange", "red","darkred"))(200) 
color_val <- myPal[findInterval( highclock_tissue, seq( 14, 35, length.out= 200), rightmost.closed= T ) ]

color_val <- myPal(200)[as.numeric(highclock_tissue)]
#points for displastic without clock
points(runif(length(high_no_clock),0,4), high_no_clock, pch = "X",col = "blue")

#points for nondisplastic without clock 
points(runif(length(low_no_clock),0,4),low_no_clock, pch = "X", col = "white")


#plot the nondisplastic points with clock values
points(x_coord[nondist_index],y_coord[nondist_index],pch="X",cex = 1.3)

#plot the displastic points with clock values
points(x_coord[high_clock_index], y_coord[high_clock_index],pch=20,col = color_val ,cex = 3)
points(x_coord[high_clock_index], y_coord[high_clock_index],pch=1,cex = 2,lwd = 2)

#more plots of data for fun 
#wireframe(best_matrix ,colorkey = TRUE)
#plot_ly(z=best_matrix, type="surface",x= c("3:00","6:00","9:00","12:00"))
#plot_ly(z=best_matrix, type = "contour",x= c("3:00","6:00","9:00","12:00"))


#retrieve samples by person
#all the rows with a particular id ie 498
#x<- split(new_pdat_072816, new_pdat_072816$Patient_ID)
#typeof(x["498"])
#x["498"][[1]]$Sample_Name


#try <- lapply(better_dat$Tissue_location, FUN = function(x) gsub("\\d\\d(?=cm)","",x,invert=TRUE))
#m <- regexpr("[0-9][0-9](?=cm)","20cm dog are good",perl=T)
#frog <- regmatches("20cm dog are good",m)
