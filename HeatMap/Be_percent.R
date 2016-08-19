library(ggplot2)

#split by user id
splitted <- split(mset,mset$length_dif)
#graph percent displastic 
changed <- lapply(splitted, FUN = function(x) length(na.omit(x$HGD.LGD))/ length(x$HGD.LGD))
changed <- as.numeric(changed)
xdirect = qplot(y = changed,x = seq(length(changed))) + labs(title="Plot of Displastic Percentage") +
labs(x="level", y="% yaaaaaaof BE Displastic") + stat_smooth() + geom_point() + geom_bar(stat = "identity") + ylim(0,1)
#xdirect + coord_flip()


