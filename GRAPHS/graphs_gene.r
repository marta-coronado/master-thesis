rm(list = ls())

data<-read.table(file="GENE",header=TRUE,sep="\t")  
data<-read.table(file="GENE.SHORTINTRON.65_8-30",header=TRUE,sep="\t")

data <- na.omit(data)
data <- subset(data,comeron_100kb > 0)
nrow(data) 

m <- data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]]
k4f <- (data[["div_4f"]])/(data[["mdyak_4f"]])
k0f <- (data[["div_0f"]]+1)/(data[["mdyak_0f"]])
kins <- (data[["div_ins"]])/(data[["mdyak_ins"]])

### Distance ### 

density = cut(data[["distance"]],c(0,554,1979,140000))

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
  coord_cartesian(ylim=c(0.0, 0.13)
) 

ggsave(file="graphs/distance.svg")

### Size ###

m_factor = cut(m,c(0,888,1743,55353))

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
  coord_cartesian(ylim=c(0.0, 0.16)
) 

ggsave(file="graphs/size.svg")

### Transcripts ###

num_transcripts_factor = cut(data[["num_transcripts"]],c(0,1,2,75))

qplot(num_transcripts_factor, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=num_transcripts_factor
) +
  xlab("Number of transcripts per gene") + 
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

ggsave(file="graphs/transcripts.svg")

### Exons ###

num_exons_factor = cut(data[["num_exons"]],c(0,3,8,114))

qplot(num_exons_factor, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=num_exons_factor
) +
  xlab("Number of exons per gene") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(0.0, 0.135)
) 

ggsave(file="graphs/exons.svg")

### Messenger Complexity ###

mess_com = cut(data[["num_transcripts"]]/data[["num_exons"]],c(0,0.33,0.66,8))

qplot(mess_com, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=mess_com
) +
  xlab("Messenger complexity") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(0.0, 0.135)
) 

ggsave(file="graphs/mcomp.svg")

### Expression ###

expression_factor = cut(data[["max_dev"]],c(0,1.477,1.87,4.4))

qplot(expression_factor, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=expression_factor
) +
  xlab("Expression") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(0.0, 0.15)
) 

ggsave(file="graphs/expression.svg")

### Breadth ###

breadth_factor = cut(data[["bias_dev"]],c(0,0.27,0.54,1))

qplot(breadth_factor, k0f, 
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=breadth_factor
) +
  xlab("Expression bias") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(angle=360, vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(angle=360, vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) + 
  coord_cartesian(ylim=c(0.0, 0.19)
) 

ggsave(file="graphs/breadth.svg")

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
  coord_cartesian(ylim=c(0.0, 0.13)
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
  coord_cartesian(ylim=c(0.0, 0.22)
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
  coord_cartesian(ylim=c(0.0, 0.17)
) 

ggsave(file="graphs/domain.svg")

### Recombination ###

recombination_factor = cut(data[["comeron_100kb"]],c(0,1.44,2.89,15))

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

k0f <- (data[["div_0f"]]+1)/(data[["mdyak_0f"]])
mutation_factor = cut(k4f,c(0,0.16,0.20,0.75)) #ojo le he quitado na omit a k4f...
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
  coord_cartesian(ylim=c(0.0, 0.18)
) 

ggsave(file="graphs/mutationk4.svg")

### mutation ins ###
# leer el otro
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

### Theta ###

var = 1
cal = 0
dom = 0
lines = 128
while (var < lines)
{
  cal = 1/var
  dom = dom + cal
  var = var + 1
}
dom

K_watterson <- (data[["seg_ins"]])/(data[["mdmel_ins"]])
k0f <- (data[["div_0f"]]+1)/(data[["mdyak_0f"]])
watterson = K_watterson/dom

watterson_factor = cut(watterson,c(0,0.005,0.013,0.087))

df <- data.frame(watterson_factor,k0f)
df <- na.omit(df)

watterson_factor <- df$watterson_factor
k0f <- df$k0f

qplot(watterson_factor, k0f,
      geom=("boxplot"), 
      outlier.color = NA, outlier.shape = NA, 
      fill=watterson_factor
) +
  xlab("Watterson estimator") + 
  ylab("Kz") + 
  scale_fill_brewer(type="seq", palette=3) + 
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold", color="grey20", size=12), 
        axis.title.y=element_text(face="bold", color="grey20", size=12), 
        axis.text.x=element_text(vjust=0.5, size=9, color="grey50"),
        axis.text.y=element_text(vjust=0.5, size=9, color="grey50")
  )  + 
  scale_x_discrete(labels=c("Low","Medium","High")) +  
  coord_cartesian(ylim=c(0.0, 0.10)
  ) 

ggsave(file="graphs/theta.svg")
