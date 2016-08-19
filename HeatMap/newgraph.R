#load mset and


mset <- read.csv("~/mset.csv", na.strings= c("",NA))
library("RColorBrewer")
library(plotly)
library(stringr)
library("fields")
library("ggplot2")
library("lattice")

#data cleaning
#get it in the proper form for the rest of the model




#only retrieve rows with tissue age
mset_1 <- mset[c(grep("[0-9]",mset$tissue_age)),]

#removes samples that are normal from data frame
mset_2 <- mset_1[-c(grep("normal",mset$Sample.Clinical.Pathology)),]

clean <- mset_2$Tissue_location

#better dat contains all the values with length or centimeter data
better_dat <- mset_2[grep("cm", clean),]
clocker <- lapply(better_dat$Tissue_location,FUN = function(x) gsub(".+cm","",x))


#age
junction <- better_dat[19]
GEJ <- c(as.numeric(unlist(junction)))

#distone and disttwo are both 
distone <- lapply(better_dat$Tissue_location, FUN = function(x) gsub("cm.+","",x))
disttwo <- c(as.numeric(lapply(distone, FUN = function(x) gsub("[^0-9]","",x))))
#I'm making it absolute value for now the values should not be negative
dif <- GEJ - disttwo


clock <- c(unlist(clocker))
#replace 5:00 with 6:00 and organizes data into data frame
clock[12] = " 6:00"

clock <- gsub("5:00", "6:00",clock)
dist <- c(dif)
#age is age of tissue
#age <- c(floor(rexp(length(dif))*50))
#age <- c(seq(1:length(dif)))
onset <- mset_2$tissue_age
age <- mset_2$Age
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
y_coord_no = as.numeric(df[grep("[^:]", df$clock), ]$dist  )
meth_value = as.numeric(df[grep("[^:]", df$clock), ]$age)
for(index in 1:length(y_coord_no)){
  # add .25 times val
  esoph[y_coord_no[index],]=  esoph[y_coord_no[index],] + meth_value[index]*.25
  weight[y_coord_no[index],]= weight[y_coord_no[index],] +.25
}

#computes weighted averages by dividing the sum of ages by the sum of weights in the weight matrix 
best_matrix <- esoph*weight^-1
best_matrix
#makes 2D heatplot 
plot_ly(z=best_matrix,x= c("3:00","6:00","9:00","12:00"), y= c("1","2","3","4","5","6","7","8","9","10","11","12","13"), type="heatmap")
heatmap(na.omit(best_matrix))

#retrieve samples by person
#all the rows with a particular id ie 498
x<- split(new_pdat_072816, new_pdat_072816$Patient_ID)
typeof(x["498"])
x["498"][[1]]$Sample_Name


#secondary plot to make smooth grap of data
#retireves the coordinates and z values of the parts of the matrix without NAs

#making copies for boundary

pos = which(!is.na(best_matrix),TRUE)

left  = pos
left[,2] = left[,2]-4 

right = pos
right[,2] = right[,2] +4

#remove above
PBCpos = rbind(pos, left,right)
PBCpos_2 = pos[,c(2,1)]


#last ditch
x <- pos[,1]
y <- pos[,2]

pos_new = cbind(y,x)

pre_val = as.vector(best_matrix)
z_index = which(!is.na(pre_val), TRUE)

z_val <- pre_val[z_index]

rorm_mat = matrix(z_val,ncol = 1)

PBCz = rbind(rorm_mat,rorm_mat)

best_fit <- Tps(pos_new,z_val)

grid.list<- list( x= seq(0,4,.1), y=seq(0,15,.1))   #Note - these ranges will need to be changed for your needs
xg<- make.surface.grid(grid.list)

#plotting function - modify as needed, including changing zlim to cover the range of y-values 
dev.off()
dev.new()
f1<- predict( best_fit, xg)
out.p<- as.surface( xg, f1)
best_plot = plot.surface(out.p,zlim=c(15,35),xlab='clock position',ylab='length in (cm) from GEJ',main='Esophagus',cex.main=1.3,cex.lab=1.3,xaxt = 'n')
axis(side=1,at = c(1,2,3,4),labels = c("3:00","6:00","9:00","12:00"))
#the code below deals with plotting points on the graph

#give an x for all the values 

#nondisplastic or not high grade
#need points specifically with clock
p_clock <-  mset_2[grep(":",mset_2$Tissue_location),]
high_clock <- grep("yes",p_clock$high_grade)
# index of nondisplastic tissue with clock
nondist_index <- setdiff(seq(1:length(high_clock)), high_clock)      #which(is.na(p_clock$high_grade))

highclock_tissue <- p_clock[high_clock]$tissue_age



#plot points that don't have clocks and are nondisplastic
#no_clock <- mset_2[grep("[^:]*", mset_2$Tissue_location),]
clock_index <- grep(":",  mset_2$Tissue_location)
no_clock <- mset_2[setdiff(seq(1:length(mset_2$Tissue_location)),clock_index),] 
low_clock_index = which(is.na(no_clock$high_grade))
#lowclock_index <- grep("[^y]",no_clock$high_grade) 
high_no_clock <- grep("y",no_clock$high_grade)
high_coord <- y_coord_no[high_no_clock]
low_coord <- y_coord_no[nonclock_index]


#tissue age for displastic tissue with clock
highclock_tissue <- p_clock[high_clock,]$tissue_age

#create colors 
#color is based on displastic tissue age
myPal <- colorRampPalette(c('red','blue'))
color_val <- myPal(5)[as.numeric(cut(highclock_tissue,breaks = 5))]

#points for displastic without clock
#points(high_coord, runif(length(high_coord),0,4), pch = "X")
#points for nondisplastic without clock 
points(high_clock,runif(length(high_clock),0,4), pch = "X")


#plot the nondisplastic points with clock values
points(x_coord[nondist_index],y_coord[nondist_index],pch="X")

#plot the displastic points with clock values
points(x_coord[high_clock], y_coord[high_clock],pch=20,col = color_val,cex.main = 4)




#more plots of data for fun 
wireframe(best_matrix ,colorkey = TRUE)
plot_ly(z=best_matrix, type="surface",x= c("3:00","6:00","9:00","12:00"))
plot_ly(z=best_matrix, type = "contour",x= c("3:00","6:00","9:00","12:00"))



# heatmap of data without high grade







