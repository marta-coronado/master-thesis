mod2=lm(kins ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

mod2=lm(kins ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

mod2=lm(kins ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")



## CHR_DOMAIN

#Expression

qplot(Chr_domain, Breadth, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Expression bias") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 1)
  ) 

ggsave(file="factor/bias_dom.svg")


#Expression

qplot(Chr_domain, Expression, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Expression") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 3)
  ) 

ggsave(file="factor/exp_dom.svg")

#Kz

qplot(Chr_domain, AAsubstitutions, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 0.14)
) 

ggsave(file="factor/aas_dom.svg")

# Recombination

qplot(Chr_domain, Recombination, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Recombination") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 10)
  ) 

ggsave(file="factor/rec_dom.svg")

# Order

qplot(Chr_domain, Order, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Order") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 18)
  ) 

ggsave(file="factor/order_dom.svg")

# Distance

qplot(Chr_domain, Distance, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Distance") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 14000)
  ) 

ggsave(file="factor/dist_dom.svg")

# Num_Exons

qplot(Chr_domain, Num_Exons, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Number of exons per exon") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 29)
  ) 

ggsave(file="factor/exons_dom.svg")

# Erep

qplot(Chr_domain, Erep, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Inclusion levels") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 1)
  ) 

ggsave(file="factor/erep_dom.svg")

# Mutation 4

qplot(Chr_domain, Mutation, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("K4") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 0.5)
  ) 

ggsave(file="factor/mutk4_dom.svg")

# Transcripts

qplot(Chr_domain, Transcripts, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Number of transcripts per exon") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 8.5)
  ) 

ggsave(file="factor/tran_dom.svg")

# Size

qplot(Chr_domain, m, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Size") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 1200)
  ) 

ggsave(file="factor/m_dom.svg")

# Mutation ins

qplot(Chr_domain, kins, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Ki") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 0.5)
  )

ggsave(file="factor/kins_dom.svg")

## CHR_STATE

#Kz
qplot(Chr_state, Breadth, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Expression bias") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 1)
  ) 

ggsave(file="factor/bias_state.svg")



#Kz
qplot(Chr_state, Expression, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Expression") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 3)
  ) 

ggsave(file="factor/exp_state.svg")


#Kz
qplot(Chr_state, AAsubstitutions, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.17)
  ) 

ggsave(file="factor/aas_state.svg")

# Recombination

qplot(Chr_state, Recombination, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Recombination") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 7.8)
  ) 

ggsave(file="factor/rec_state.svg")

# Order

qplot(Chr_state, Order, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Order") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 25)
  ) 

ggsave(file="factor/order_state.svg")

# Distance

qplot(Chr_state, Distance, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Distance") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 15000)
  ) 

ggsave(file="factor/dist_state.svg")

# Num_Exons

qplot(Chr_state, Num_Exons, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Number of exons per exon") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 41)
  ) 

ggsave(file="factor/exons_state.svg")

# Erep

qplot(Chr_state, Erep, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Erep") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 1)
  ) 

ggsave(file="factor/erep_state.svg")

# Mutation 4

qplot(Chr_state, Mutation, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("K4") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.47)
  ) 

ggsave(file="factor/mutk4_state.svg")

### Recombination ###

qplot(Chr_state, Recombination, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Recombination rate") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 8)
  ) 

ggsave(file="factor/rec_state.svg")

# Transcripts

qplot(Chr_state, Transcripts, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Number of transcripts per exon") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 9)
  ) 

ggsave(file="factor/tran_state.svg")


# Size

qplot(Chr_state, m, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Size") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 1100)
  ) 

ggsave(file="factor/m_state.svg")


# Mutation ins

qplot(Chr_state, kins, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Ki") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.48)
  )

ggsave(file="factor/kins_state.svg")


## CHROMOSOME

#Kz
qplot(Chromosome, Breadth, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Expression bias") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 1)
  ) 

ggsave(file="factor/bias_chr.svg")
#Kz
qplot(Chromosome, Expression, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Expression") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 3.5)
  ) 

ggsave(file="factor/exp_chr.svg")


#Kz
qplot(Chromosome, AAsubstitutions, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.12)
  ) 

ggsave(file="factor/aas_chr.svg")

# Recombination

qplot(Chromosome, Recombination, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Recombination") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 7.8)
  ) 

ggsave(file="factor/rec_chr.svg")

# Order

qplot(Chromosome, Order, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Order") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 15)
  ) 

ggsave(file="factor/order_chr.svg")

# Distance

qplot(Chromosome, Distance, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Distance") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 4300)
  ) 

ggsave(file="factor/dist_chr.svg")

# Num_Exons

qplot(Chromosome, Num_Exons, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Number of exons per exon") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 25)
  ) 

ggsave(file="factor/exons_chr.svg")

# Erep

qplot(Chromosome, Erep, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Erep") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 1)
  ) 

ggsave(file="factor/erep_chr.svg")

# Mutation 4

qplot(Chromosome, Mutation, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("K4") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.44)
  ) 

ggsave(file="factor/mutk4_chr.svg")

# Transcripts

qplot(Chromosome, Transcripts, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Number of transcripts per exon") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 6.8)
  ) 

ggsave(file="factor/tran_chr.svg")


# Size

qplot(Chromosome, m, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Size") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 1100)
  ) 

ggsave(file="factor/m_chr.svg")

# Mutation ins

qplot(Chromosome, kins, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Ki") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.5)
  )

ggsave(file="factor/kins_chr.svg")
