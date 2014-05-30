qplot(data[["distance"]],data[["Ka_pos"]],
      geom=("boxplot"),
      outlier.color=NA, outlier.shape = NA,
      fill=data[["distance"]]
) + 
  xlab("Distance") + 
  ylab("Kz+") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  ) + 
  scale_x_discrete(labels=c("Low","Medium","High","Super")) +   
 # scale_x_discrete(labels=c("Active","Both","Inactive")) + 
   #scale_x_discrete(labels=c("Low","High")) +
  coord_cartesian(ylim=c(-0.03, 0.07)
  )

ggsave(file="comparisons/Kz+distanceWObreadth.svg")
mysvg <- grid.export("comparisons/Kz+distanceWOsize.svg")

qplot(data[["breadth"]],data[["Ka_pos"]],
      geom=("boxplot"),
      outlier.color=NA, outlier.shape = NA,
      fill=data[["breadth"]]
) + 
  xlab("Expression bias") + 
  ylab("Kz+") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  ) + 
  scale_x_discrete(labels=c("Low","Medium","High", "Super")) +   
  #scale_x_discrete(labels=c("Active","Both","Inactive")) + 
  #scale_x_discrete(labels=c("Low","High")) +
  coord_cartesian(ylim=c(-0.035, 0.075)
  )

ggsave(file="comparisons/Kz+breadthWOdistance.svg")
mysvg <- grid.export("comparisons/Kz+sizeWOdistance.svg")

qplot(data[["distance"]],data[["Ka_neg"]],
      geom=("boxplot"),
      outlier.color=NA, outlier.shape = NA,
      fill=data[["distance"]]
) + 
  xlab("Distance") + 
  ylab("Kz-") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  ) + 
  scale_x_discrete(labels=c("Low","Medium","High")) +   
  #scale_x_discrete(labels=c("Active","Both","Inactive")) + 
  #scale_x_discrete(labels=c("Low","High")) +
  coord_cartesian(ylim=c(0, 0.09)
  )

ggsave(file="comparisons/Kz-distanceWObreadth.svg")
mysvg <- grid.export("comparisons/Kz-breadthWOsize.svg")

qplot(data[["breadth"]],data[["Ka_neg"]],
      geom=("boxplot"),
      outlier.color=NA, outlier.shape = NA,
      fill=data[["breadth"]]
) + 
  xlab("Expression bias") + 
  ylab("Kz-") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  ) + 
  scale_x_discrete(labels=c("Low","Medium","High", "Super")) +   
  #scale_x_discrete(labels=c("Active","Both","Inactive")) +   
  #scale_x_discrete(labels=c("Low","High")) +
  coord_cartesian(ylim=c(0, 0.07)
  )

ggsave(file="comparisons/Kz-breadthWOdistance.svg")
mysvg <- grid.export("comparisons/Kz-sizeWObreadth.svg")