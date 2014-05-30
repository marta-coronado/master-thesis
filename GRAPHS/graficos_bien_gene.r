qplot(data[["mcomplexity"]],data[["Ka_pos"]],
      geom=("boxplot"),
      outlier.color=NA, outlier.shape = NA,
      fill=data[["mcomplexity"]]
) + 
  xlab("Messenger complexity") + 
  ylab("Kz+") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(plot.title=element_text(face="bold",size="20",color="black"),
        legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  ) + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(-0.04, 0.07)
  )

ggsave(file="comparisons/Kz+mcomplexityWObreadth.svg")

qplot(data[["transcripts"]],data[["Ka_neg"]],
      geom=("boxplot"),
      outlier.color=NA, outlier.shape = NA,
      fill=data[["transcripts"]]
) + 
  xlab("Number of transcripts") + 
  ylab("Kz-") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(plot.title=element_text(face="bold",size="20",color="black"),
        legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  ) + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  #scale_x_discrete(labels=c("Active","Both","Inactive")) + 
  coord_cartesian(ylim=c(0, 0.06)
  )

ggsave(file="comparisons/Kz-tranWObreadth.svg")
