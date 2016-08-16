mset <- read.csv("~/msetHans16.csv", na.strings= c("",NA))
library("RColorBrewer")
library(plotly)
library(stringr)
library("fields")
library("ggplot2")
library("lattice")
mset <- mset[grep("yes", mset$BE),]
#only retrieve rows with tissue age
mset <- mset[c(grep("[0-9]",mset$tissue_age)),]

mset<- mset[grep("cm", mset$Tissue_location),]
clock <- lapply(mset$Tissue_location,FUN = function(x) gsub(".+cm","",x))


#caculate distance from GEJ junction
GEJ <- mset$gejunction_cm 
GEJ <- c(as.numeric(unlist(GEJ)))

dist <- lapply(mset$Tissue_location, FUN = function(x) gsub("cm.+","",x))
dist <- c(as.numeric(lapply(dist, FUN = function(x) gsub("[^0-9]","",x))))
dif <- GEJ - dist


clock <- c(unlist(clock))
#replace 5:00 with 6:00 and organizes data into data frame
clock <- gsub("5:00", "6:00",clock)
dif <- c(dif)
onset <- mset$tissue_age
age <- mset$Age
dwell <- age -onset

df <- data.frame(clock,dif,dwell)
names(df) <- c("clock","dif","age")

#bifuricate the data into two groups displastic and non displastic 
circles <- df[grep("yes",mset$HGD.LGD),]

df <- df[-c(grep("yes", mset$HGD.LGD)),]
clock_rows = grep(":", df$clock)
#clock_rows picks out all the indices with clock location
#X_coord and y_coord are justt in the cases where there is clock data

#puts clock levels in right order ie 12:00 -> 4 9:00 -> 3 etc
check <- factor(df[clock_rows, ]$clock)
levels(check) <- c(4,1,2,3)
x_coord = as.numeric(levels(check))[check]

#these are the y coordinates at least from R console
y_coord = as.numeric(df[clock_rows, ]$dif)
mether = as.numeric(df[clock_rows, ]$age)
esoph = matrix(rep(0,len=80), nrow = 20, ncol= 4)
weight = matrix(rep(0,len=80), nrow=20,ncol=4)



#I don't think this is correctly filling 
#this shifts up the y values
y_coord  = y_coord + 1

#adds position and heat data to esoph diagram 
for(vale in 1:length(x_coord)){
  
  #print(mether[vale])
  esoph[y_coord[vale],x_coord[vale]] = esoph[y_coord[vale],x_coord[vale]]+ mether[vale]
  weight[y_coord[vale],x_coord[vale]] =  weight[y_coord[vale],x_coord[vale]]+1
}
plot_ly(z=esoph, type="heatmap")

heatmap(weight)

# for loops adds data for biopsises without clock data

#do these map to the same numeric value 
y_coord_no = as.numeric(df[-c(grep(":", df$clock)), ]$dif )
y_coord_no = y_coord_no +1
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
#makes 2D discrete heatmap
plot_ly(z=best_matrix,x= c("3:00","6:00","9:00","12:00"), type="heatmap")


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
PBCpos[,2] = PBCpos[,2] - 1



z_index = as.vector(best_matrix)
z_index = which(!is.na(z_index), TRUE)
z_index <- pre_val[z_index]
rorm_mat = matrix(z_index,ncol = 1)
PBCz = rbind(rorm_mat,rorm_mat,rorm_mat)


best_fit <- Tps(PBCpos,PBCz,cost= 1.2)

grid.list<- list( x= seq(0,4,.1), y=seq(0,12,.1))   #Note - these ranges will need to be changed for your needs
xg<- make.surface.grid(grid.list)

#plotting function - modify as needed, including changing zlim to cover the range of y-values 
dev.off()
dev.new()
f1<- predict( best_fit, xg)
out.p<- as.surface( xg, f1)
myPal <- colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "orange", "darkorange", "red","darkred"))
best_plot = plot.surface(out.p,zlim=c(15,35),xlab='clock position',ylab='length in (cm) from GEJ',main='Esophagus',cex.main=1.3,cex.lab=1.3,cex.axis = 1.5, xaxt = 'n',col = myPal(200))
axis(side=1,at = c(1,2,3,4),labels = c("3:00","6:00","9:00","12:00"),cex.lab = 2)

#the code below deals with plotting points on the graph


#plot the nondisplastic without clock
points(runif(length(y_coord_no),0,4),y_coord_no-1, pch = "X", col = "white")
#plot the nondisplastic points with clock values
points((x_coord+runif(length(x_coord),-.5,.5))%%4,y_coord-1, pch = "X")

#p_lock contains  all the displastic tissue with clock values
p_clock <-  circles[grep(":",circles$clock),]
#color is based on displastic tissue age
color_val <- myPal(200)[findInterval( p_clock$age, seq( 15,35, length.out= 200), rightmost.closed= T ) ]


#convert <- function(set){
  #check_1 <- factor(set)
  #print check_1
  #levels(set) <- c(4,1,2,3)
  #fixedcoord = as.numeric(levels(set))[set] 
  #fixedcoord
 # }

convert(p_clock$clock)
#convert the pclock clock values to numerals
check_1 <- factor(p_clock$clock)
levels(check_1) <- c(4,1,2,3)
fixedcoord = as.numeric(levels(check_1))[check_1]
#plot the displastic points with clock values
jitter = runif(length(p_clock$dif),-.5,.5)
points((fixedcoord+jitter)%%4, p_clock$dif,pch=20,col = color_val ,cex = 3)
points((fixedcoord+jitter)%%4, p_clock$dif,pch=1,cex= 2,lwd = 2)

#plot the displastic points without clock values

no_clock <-  circles[-c(grep(":",circles$clock)),]
no_color_val <- myPal(200)[findInterval( no_clock$age, seq( 15,35, length.out= 200), rightmost.closed= T ) ]
clock_jit = runif(length(no_clock),0,4)
points(clock_jit, no_clock$dif,pch=20,col = no_color_val ,cex = 3)
points(clock_jit, no_clock$dif,pch=1,cex= 2,lwd = 2)


















