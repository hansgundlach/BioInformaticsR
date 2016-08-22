library(ggplot2)
#modifed set
mset$length_dif <- dif
#split by length after junction
splitted <- split(mset,mset$length_dif)
#graph percent displastic 
changed <- lapply(splitted, FUN = function(x) length(na.omit(x$HGD.LGD))/ length(x$HGD.LGD))
changed <- as.numeric(changed)
xdirect = qplot(y = changed,x = seq(from = 0,to = length(changed)-1)) + labs(title="Plot of Displastic Percentage") +
labs(x="level", y="% of BE Displastic") + stat_smooth(se = FALSE) + geom_point() + geom_bar(stat = "identity") + ylim(0,1)
xdirect + coord_flip()


