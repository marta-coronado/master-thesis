rm(list = ls())

data<-read.table(file="EXON",header=TRUE,sep="\t")
data<-read.table(file="EXON.SHORTINTRON.65_8-30",header=TRUE,sep="\t")

data<-na.omit(data) 
data <- subset(data,comeron_100kb > 0)
nrow(data) 

m <- data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]]
k4f <- (data[["div_4f"]])/(data[["mdyak_4f"]])
k0f <- (data[["div_0f"]]+1)/(data[["mdyak_0f"]])
kins <- (data[["div_ins"]])/(data[["mdyak_ins"]])

### Size ###

m_factor = cut(m,c(0,149,328,14000))

qplot(m_factor, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=m_factor
) +
  xlab("Size") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(0.0, 0.13)
) 

ggsave(file="graphs/size.svg")

### Transcripts ###

num_transcripts_exon_factor = cut(data[["num_transcripts_exon"]],c(0,1,2,75))

qplot(num_transcripts_exon_factor, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=num_transcripts_exon_factor
) +
  xlab("Number of transcripts per exon") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(0.0, 0.13)
) 

ggsave(file="graphs/transcripts.svg")

### Inclusion Levels ### 

Erep <-(data[["num_transcripts_exon"]]/data[["num_transcripts_gene"]])
alternative = cut(Erep,c(0,0.99,1))

qplot(alternative, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=alternative
) +
  xlab("Inclusion levels") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(0.0, 0.11)
) 

ggsave(file="graphs/erep.svg")

### Order ####

order = cut(data[["order"]],c(0,1,2,4,112))

qplot(order, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=order
) +
  xlab("Order") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(0.0, 0.16)
) 

ggsave(file="graphs/order.svg")

### Chromosome ###

qplot(data[["chr"]], k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=data[["chr"]]
) +
  xlab("Chromosome") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.12)
  ) 

ggsave(file="graphs/chr.svg")

### Chromatin State ###

qplot(data[["state"]], k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=data[["state"]]
) +
  xlab("Chromatin state") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.175)
  ) 

ggsave(file="graphs/state.svg")

### Chromatin Domain ###

qplot(data[["domain"]], k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=data[["domain"]]
) +
  xlab("Chromatin domain") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  coord_cartesian(ylim=c(0.0, 0.12)
  ) 

ggsave(file="graphs/domain.svg")

### Recombination ###

recombination_factor = cut(data[["comeron_100kb"]],c(0,1.42,2.9,15))

qplot(recombination_factor, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=recombination_factor
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
  scale_x_discrete(labels=c("Low","Medium","High")) +  
  coord_cartesian(ylim=c(0.0, 0.11)
  ) 

ggsave(file="graphs/recombination.svg")

### Mutation 4f ###

mutation_factor = cut(k4f,c(0,0.12,0.20,0.75))

Mutation2=na.omit(mutation_factor)

df <- data.frame(mutation_factor,k0f)
df <- na.omit(df)

mutation <- df$mutation_factor
k0f <- df$k0f

qplot(mutation, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=Mutation2
) +
  xlab("K4") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) +  
  coord_cartesian(ylim=c(0.0, 0.16)
  ) 

ggsave(file="graphs/mutationk4.svg")

### Mutation ins ###

k0f <- (data[["div_0f"]]+1)/(data[["mdyak_0f"]])
mutation_factor_ins = cut(kins,c(0,0.18,0.25,0.7))
mutation_ins=na.omit(mutation_factor_ins)

df <- data.frame(mutation_factor_ins,k0f)
df <- na.omit(df)

mutation_ins_f <- df$mutation_factor_ins
k0f <- df$k0f

qplot(mutation_ins_f, k0f,
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=mutation_ins
) +
  xlab("Ki") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) +  
  coord_cartesian(ylim=c(0.0, 0.12)
  ) 

ggsave(file="graphs/mutation_ins.svg")

### Number of exons ###

numexons_factor = cut(data[["num_exons"]],c(0,5,9,114))

qplot(numexons_factor, k0f,
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=numexons_factor
) +
  xlab("Number of exons") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) +  
  coord_cartesian(ylim=c(0.0, 0.14)
  ) 

ggsave(file="graphs/exons.svg")

### Distance ### 

density = cut(data[["distance"]],c(0,128,703,200000))

qplot(density, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=density
) +
  xlab("Density") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) +  
  coord_cartesian(ylim=c(0.0, 0.11)
  ) 

ggsave(file="graphs/distance.svg")




