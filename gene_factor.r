library(car)
library(lattice)
library(ppcor)

qplot(Chr_domain,Recombination, geom=("boxplot"),fill=Chr_domain, binwidth=3) + xlab("Chromatin domain") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right") 
qplot(Chr_domain,Recombination, geom=("boxplot"),fill=Chromosome, binwidth=3) + xlab("Chromatin domain") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right") 
qplot(Chromosome,Mutation, geom=("boxplot"), binwidth=3,fill=Chromosome) + xlab("Chromosome") + ylab("Mutation") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE)

# para saber si son significativos: 	kruskal.test(recombination~chromosome)
# para saber si son significativos mejor usar LM con anova:
mod2=lm(Recombination ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

# cambiar fuentes.... mirar graficos_bien.r

######### GENES #############

## CHR_DOMAIN
#Kz
qplot(Chr_domain,AAsubstitutions, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("AAsubstitutions") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right") + coord_cartesian(ylim = quantile(AAsubstitutions, c(0, 0.99)))
mod2=lm(AAsubstitutions ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

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
  coord_cartesian(ylim=c(0.0, 0.23)
) 

ggsave(file="factor/aas_dom.svg")

# Recombination

qplot(Chr_domain,Recombination, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Recombination, c(0, 0.995)))
mod2=lm(Recombination ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

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
  coord_cartesian(ylim=c(0.0, 9.5)
  ) 

ggsave(file="factor/rec_dom.svg")

# Breadth

qplot(Chr_domain,Breadth, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Breadth") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Breadth, c(0, 1)))
mod2=lm(Breadth ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

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

ggsave(file="factor/br_dom.svg")

# Distance

qplot(Chr_domain,Distance, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Distance") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Distance, c(0, 0.9925)))
mod2=lm(Distance ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

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
  coord_cartesian(ylim=c(0.0, 33000)
  ) 

ggsave(file="factor/dist_dom.svg")

# Exons

qplot(Chr_domain,Exons, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Exons") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Exons, c(0, 0.9885)))
mod2=lm(Exons ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain, Exons, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Number of exons per gene") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 20)
  ) 

ggsave(file="factor/exons_dom.svg")

# Expression

qplot(Chr_domain,Expression, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Expression") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Expression, c(0, 0.99)))
mod2=lm(Expression ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

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
  coord_cartesian(ylim=c(0.0, 3.8)
  ) 

ggsave(file="factor/exp_dom.svg")

# Mutation 4

qplot(Chr_domain,Mutation, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Mutation") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Mutation, c(0, 0.995)))
mod2=lm(Mutation ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

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
  coord_cartesian(ylim=c(0.0, 0.4)
  ) 

ggsave(file="factor/mutk4_dom.svg")

# Transcripts

qplot(Chr_domain,Transcripts, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Transcripts") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Transcripts, c(0.1, 0.96)))
mod2=lm(Transcripts ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain, Transcripts, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Number of transcripts per gene") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 7)
  ) 

ggsave(file="factor/tran_dom.svg")

# Size

qplot(Chr_domain,m, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Size") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(m, c(0.1, 0.96)))
mod2=lm(m ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

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
  coord_cartesian(ylim=c(0.0, 6000)
  ) 

ggsave(file="factor/m_dom.svg")

# Watterson

qplot(Chr_domain, watterson, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromatin domain") + 
  ylab("Watterson estimator") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Active","Both","Inactive")) +  
  coord_cartesian(ylim=c(0.0, 0.055)
  ) 

ggsave(file="factor/theta_dom.svg")

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
  coord_cartesian(ylim=c(0.0, 0.52)
  )

ggsave(file="factor/kins_dom.svg")

## CHR_STATE

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
  coord_cartesian(ylim=c(0.0, 0.23)
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

# Breadth

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

ggsave(file="factor/br_state.svg")

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
  coord_cartesian(ylim=c(0.0, 27000)
  ) 

ggsave(file="factor/dist_state.svg")

# Exons

qplot(Chr_state, Exons, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Number of exons per gene") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 16)
  ) 

ggsave(file="factor/exons_state.svg")

# Expression

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
  coord_cartesian(ylim=c(0.0, 3.8)
  ) 

ggsave(file="factor/exp_state.svg")

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
  coord_cartesian(ylim=c(0.0, 0.4)
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
  ylab("Number of transcripts per gene") + 
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
  coord_cartesian(ylim=c(0.0, 6200)
  ) 

ggsave(file="factor/m_state.svg")


# Watterson

qplot(Chr_state, watterson, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chr_state
) +
  xlab("Chromatin state") + 
  ylab("Watterson estimator") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.045)
  ) 

ggsave(file="factor/theta_state.svg")

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
  coord_cartesian(ylim=c(0.0, 0.52)
  )

ggsave(file="factor/kins_state.svg")

qplot(Chr_state,AAsubstitutions, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("AAsubstitutions") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(AAsubstitutions, c(0, 0.99)))
mod2=lm(AAsubstitutions ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Recombination, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Recombination, c(0, 0.995)))
mod2=lm(Recombination ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Breadth, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Breadth") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Breadth, c(0, 1)))
mod2=lm(Breadth ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Distance, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Distance") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Distance, c(0, 0.9925)))
mod2=lm(Distance ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Exons, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Exons") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Exons, c(0, 0.9885)))
mod2=lm(Exons ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Expression, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Expression") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Expression, c(0, 0.99)))
mod2=lm(Expression ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Mutation, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Mutation") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Mutation, c(0, 0.995)))
mod2=lm(Mutation ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Transcripts, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Transcripts") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Transcripts, c(0.1, 0.96)))
mod2=lm(Transcripts ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,m, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Size") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(m, c(0, 0.97)))
mod2=lm(m ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

## CHROMOSOME

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
  coord_cartesian(ylim=c(0.0, 0.13)
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

# Breadth

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

ggsave(file="factor/br_chr.svg")

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
  coord_cartesian(ylim=c(0.0, 10000)
  ) 

ggsave(file="factor/dist_chr.svg")

# Exons

qplot(Chromosome, Exons, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Number of exons per gene") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 13)
  ) 

ggsave(file="factor/exons_chr.svg")

# Expression

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
  coord_cartesian(ylim=c(0.0, 0.4)
  ) 

ggsave(file="factor/mutk4_chr.svg")

# Transcripts

qplot(Chromosome, Transcripts, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Number of transcripts per gene") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 6.2)
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
  coord_cartesian(ylim=c(0.0, 4600)
  ) 

ggsave(file="factor/m_chr.svg")


# Watterson

qplot(Chromosome, watterson, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Chromosome
) +
  xlab("Chromosome") + 
  ylab("Watterson estimator") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.045)
  ) 

ggsave(file="factor/theta_chr.svg")

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
  coord_cartesian(ylim=c(0.0, 0.52)
  )

ggsave(file="factor/kins_chr.svg")



qplot(Chromosome,AAsubstitutions, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("AAsubstitutions") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(AAsubstitutions, c(0, 0.95)))
mod2=lm(AAsubstitutions ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Mutation, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Mutation") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Mutation, c(0, 0.99)))
mod2=lm(Mutation ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Breadth, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Breadth") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Breadth, c(0, 1)))
mod2=lm(Breadth ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Distance, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Distance") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Distance, c(0, 0.93)))
mod2=lm(Distance ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Exons, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Exons") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Exons, c(0, 0.95)))
mod2=lm(Exons ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Expression, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Expression") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Expression, c(0, 0.99)))
mod2=lm(Expression ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Recombination, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Recombination, c(0, 0.99)))
mod2=lm(Recombination ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Transcripts, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Transcripts") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Transcripts, c(0, 0.96)))
mod2=lm(Transcripts ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")


qplot(Chromosome,m, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Size") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(m, c(0, 0.95)))
mod2=lm(m ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

######### EXONS #############
## CHR_DOMAIN
qplot(Chr_domain,AAsubstitutions, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("AAsubstitutions") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(AAsubstitutions, c(0, 0.96)))
AAsubstitutions[which(is.nan(AAsubstitutions))] = NA
AAsubstitutions[which(AAsubstitutions==Inf)] = NA
mod2=lm(AAsubstitutions ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain,Recombination, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Recombination, c(0, 0.998)))
mod2=lm(Recombination ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain,Erep, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Inclusion levels") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Erep, c(0, 1)))
mod2=lm(Erep ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain,Distance, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Distance") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Distance, c(0, 0.975)))
mod2=lm(Distance ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain,Num_Exons, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Number of exons") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Num_Exons, c(0, 0.965)))
mod2=lm(Num_Exons ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain,Order, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Order") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Order, c(0, 0.95)))
mod2=lm(Order ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

Mutation2=na.omit(Mutation)
qplot(Chr_domain,Mutation, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Mutation") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Mutation2, c(0, 0.99)))
mod2=lm(Mutation ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain,Transcripts, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Transcripts") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Transcripts, c(0, 0.965)))
mod2=lm(Transcripts ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_domain,m, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin domain") + ylab("Size") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(m, c(0.1, 0.95)))
mod2=lm(m ~ Chr_domain)
summary(mod2)
Anova(mod2,type="2")


## CHR_STATE
qplot(Chr_state,AAsubstitutions, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("AAsubstitutions") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(AAsubstitutions, c(0, 0.96)))
AAsubstitutions[which(is.nan(AAsubstitutions))] = NA
AAsubstitutions[which(AAsubstitutions==Inf)] = NA
mod2=lm(AAsubstitutions ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Recombination, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Recombination, c(0, 0.996)))
mod2=lm(Recombination ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Erep, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Inclusion levels") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Erep, c(0, 1)))
mod2=lm(Erep ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Distance, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Distance") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Distance, c(0, 0.99)))
mod2=lm(Distance ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Num_Exons, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Number of exons") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Num_Exons, c(0, 0.995)))
mod2=lm(Num_Exons ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Order, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Order") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Order, c(0, 0.99)))
mod2=lm(Order ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

Mutation2=na.omit(Mutation)
qplot(Chr_state,Mutation, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Mutation") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Mutation2, c(0, 0.995)))
mod2=lm(Mutation ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,Transcripts, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Transcripts") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(Transcripts, c(0, 0.9987)))
mod2=lm(Transcripts ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")

qplot(Chr_state,m, geom=("boxplot"),fill=Chromosome, binwidth=3,outlier.color=NA,outlier.shape = NA) + xlab("Chromatin state") + ylab("Size") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right")  + coord_cartesian(ylim = quantile(m, c(0.1, 0.95)))
mod2=lm(m ~ Chr_state)
summary(mod2)
Anova(mod2,type="2")


## CHROMOSOME
qplot(Chromosome,AAsubstitutions, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("AAsubstitutions") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(AAsubstitutions, c(0, 0.92)))
AAsubstitutions[which(is.nan(AAsubstitutions))] = NA
AAsubstitutions[which(AAsubstitutions==Inf)] = NA
mod2=lm(AAsubstitutions ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

Mutation2=na.omit(Mutation)
qplot(Chromosome,Mutation, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Mutation") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Mutation2, c(0, 0.99)))
mod2=lm(Mutation ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Erep, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Inclusion levels") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Erep, c(0, 1)))
mod2=lm(Erep ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Distance, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Distance") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Distance, c(0, 0.9)))
mod2=lm(Distance ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Num_Exons, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Exons") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Num_Exons, c(0, 0.97)))
mod2=lm(Num_Exons ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Order, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Order") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Order, c(0, 0.95)))
mod2=lm(Order ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Recombination, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Recombination, c(0, 0.99)))
mod2=lm(Recombination ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,Transcripts, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Transcripts") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(Transcripts, c(0, 0.94)))
mod2=lm(Transcripts ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

qplot(Chromosome,m, geom=("boxplot"), binwidth=3,fill=Chromosome,outlier.color=NA,outlier.shape = NA) + xlab("Chromosome") + ylab("Size") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"))  +  guides(fill=FALSE) + coord_cartesian(ylim = quantile(m, c(0, 0.94)))
mod2=lm(m ~ Chromosome)
summary(mod2)
Anova(mod2,type="2")

ggsave(file="graph.svg")







