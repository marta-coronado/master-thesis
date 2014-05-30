library(car)
library(lattice)
library(ppcor)
library(pls) 

rm(list = ls())


#data<-read.table(file="GENE",header=TRUE,sep="\t")		
#data<-read.table(file="GENE.SHORTINTRON.65",header=TRUE,sep="\t")
data<-read.table(file="GENE.SHORTINTRON.65_8-30",header=TRUE,sep="\t")
data <- na.omit(data)

data <- subset(data,comeron_100kb > 0)
nrow(data) 

data <- subset(data,chr == "2L")
data <- subset(data,chr == "2R")
data <- subset(data,chr == "3L")
data <- subset(data,chr == "3R")
data <- subset(data,chr == "X")
data <- subset(data,chr != "X")

m <- data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]]
k4f <- (data[["div_4f"]])/(data[["mdyak_4f"]])
k0f <- (data[["div_0f"]])/(data[["mdyak_0f"]])
k0f <- (data[["div_0f"]]+1)/(data[["mdyak_0f"]])
kins <- (data[["div_ins"]])/(data[["mdyak_ins"]])


### ASSESSIN COEVOLUTION ###

gene <- data.frame(
                                Transcripts = data[["num_transcripts"]],        
                                Exons = data[["num_exons"]],
                                m = m,
                                Distance = data[["distance"]],
                                Breadth = data[["bias_dev"]],
                                Expression = data[["max_dev"]],
                                Ki = kins,
                                AAsubstitutions = k0f,
                                Mutation = k4f,   
                                #Chr_state = data[["state"]],                                
                                #Chr_domain = data[["domain"]],
                                #Chromosome = data[["chr"]],                                                                                          
                                Recombination=(data[["comeron_100kb"]]),
                                mess_com = (data[["num_transcripts"]]/data[["num_exons"]])
                          )
gene <- na.omit(gene)
nrow(gene)
spcor(gene,method=c("spearman"))

Transcripts = data[["num_transcripts"]]
Exons = data[["num_exons"]]
Distance = data[["distance"]]
Breadth = data[["bias_dev"]]
Expression = data[["max_dev"]]
AAsubstitutions = k0f
Mutation = k4f
Chr_state = data[["state"]]                               
Chr_domain = data[["domain"]]
Chromosome = data[["chr"]]                                                                                        
Recombination=(data[["comeron_100kb"]]) 
mess_com = (data[["num_transcripts"]]/data[["num_exons"]])

	### Distance ### 
	summary(data[["distance"]])
	quantile(data[["distance"]],(0:3)/3)	
	density = cut(data[["distance"]],c(0,554,1979,140000))

	with(data,tapply(m,density,sum))
	summary(density)
		
		boxplot(k0f~density,outline=F,xlab="Distance",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,data[["distance"]],method="spearman") 
		cor.test(data[["distance"]],m,method="spearman") 
		cor.test(k0f,m,method="spearman") 
		cor.test(k0f,data[["distance"]]/(m/data[["num_exons"]]),method="spearman") 

		kruskal.test(k0f~density)	

	### Size ###
	summary(m)
	quantile(m,(0:3)/3)	
		m_factor = cut(m,c(0,888,1743,55353))
		summary(m_factor)
		with(data,tapply(m,m_factor,sum))
		
		boxplot(k0f~m_factor,outline=F,xlab="Gene Size",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,m,method="spearman")
		kruskal.test(k0f~m_factor)

	### Transcripts ###
quantile(data$num_transcripts,(0:3)/3)	

hist(data$num_transcripts,breaks=100,xlim=c(0,20))
	summary(data[["num_transcripts"]])
	with(data,tapply(m,data[["num_transcripts"]],sum))
		
		num_transcripts_factor = cut(data[["num_transcripts"]],c(0,1,2,75))
		summary(num_transcripts_factor)
		with(data,tapply(m,num_transcripts_factor,sum))
		
		boxplot(k0f~num_transcripts_factor,outline=F,xlab="Number of Transcripts/Gene",ylab="Ka")
		boxplot(k0f~data$num_transcripts,outline=F,xlab="Number of Transcripts/Gene",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,data$num_transcripts,method="spearman")
		kruskal.test(k0f~num_transcripts_factor)

	### Exons ###
quantile(data$num_exons,(0:3)/3)  
hist(data$num_exons,breaks=100,xlim=c(0,20))

	summary(data[["num_exons"]])	
	with(data,tapply(m,data[["num_exons"]],sum))
	
		num_exons_factor = cut(data[["num_exons"]],c(0,3,8,114))
		summary(num_exons_factor)
		with(data,tapply(m,num_exons_factor,sum))
		
		boxplot(k0f~num_exons_factor,outline=F,xlab="Number of Exons/Gene",ylab="Ka")
		boxplot(k0f~data$num_exons,outline=F,xlab="Number of Exons/Gene",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,data$num_exons,method="spearman")
		kruskal.test(k0f~num_exons_factor)

	### Messenger Complexity ###
	plot(data[["num_transcripts"]],data[["num_exons"]])
	quantile(data[["num_transcripts"]]/data[["num_exons"]],(0:3)/3)
	summary(data[["num_transcripts"]]/data[["num_exons"]])
	hist(data[["num_transcripts"]]/data[["num_exons"]])
	
	mess_com = cut(data[["num_transcripts"]]/data[["num_exons"]],c(0,0.33,0.66,8))
	tapply(m,mess_com,sum)
	summary(mess_com) 
		
		boxplot(k0f~mess_com,outline=F,xlab="Messenger Complexity",ylab="Ka")
		abline(h=median(k0f),col="black")
		boxplot(data[["num_exons"]]~mess_com,outline=F,xlab="Messenger Complexity",ylab="#Exons")
		boxplot(data[["num_transcripts"]]~mess_com,outline=F,xlab="Messenger Complexity",ylab="#mRNAs")

		cor.test(k0f,data[["num_transcripts"]]/data[["num_exons"]],method="spearman")
		cor.test(data[["num_transcripts"]],data[["num_exons"]],method="spearman")

		kruskal.test(k0f~mess_com)	
		
	### Expression ###
	plot(data[["bias_dev"]],data[["max_dev"]])
	summary(data[["max_dev"]])
	quantile(data[["max_dev"]],(0:3)/3)
	
		expression_factor = cut(data[["max_dev"]],c(0,1.477,1.87,4.4))
		summary(expression_factor)
		tapply(m,expression_factor,sum)
		
		boxplot(k0f~expression_factor,outline=F,xlab="Expression",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,data[["max_dev"]],method="spearman")
		kruskal.test(k0f~expression_factor)

	### Breadth ###
	summary(data[["bias_dev"]])
	quantile(data[["bias_dev"]],(0:3)/3)
	
		breadth_factor = cut(data[["bias_dev"]],c(0,0.27,0.54,1))
		summary(breadth_factor)
		tapply(m,breadth_factor,sum)
		
		boxplot(k0f~breadth_factor,outline=F,xlab="Expression Breadth",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,data[["bias_dev"]],method="spearman")
		kruskal.test(k0f~breadth_factor)

	### Chromosome ###
	
	tapply(m,data[["chr"]],sum)
	summary(data[["chr"]]) 
	
	boxplot(k0f~data[["chr"]],outline=F,xlab="Chromosome",ylab="Ka")
		
	### Chromatin State and Chromatin Domain ###
	
	tapply(m,data[["state"]],sum)
	tapply(m,data[["domain"]],sum)
	summary(data[["state"]]) 
	summary(data[["domain"]]) 
	
	boxplot(k0f~data[["state"]],outline=F,xlab="Chromatin State",ylab="Ka")
	boxplot(k0f~data[["domain"]],outline=F,xlab="Chromatin Domain",ylab="Ka")
		
	### Recombination ###
	summary(data[["comeron_100kb"]])
	quantile(data[["comeron_100kb"]],(0:3)/3)
	
		recombination_factor = cut(data[["comeron_100kb"]],c(0,1.44,2.89,15))
		summary(recombination_factor)
		tapply(m,recombination_factor,sum)
		
		boxplot(k0f~recombination_factor,outline=F,xlab="Recombination Rate",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,data[["comeron_100kb"]],method="spearman")
		kruskal.test(k0f~recombination_factor)

	### mutation 4f ###
	summary(k4f)
	quantile(k4f,(0:3)/3)
	
		mutation_factor = cut(k4f,c(0,0.16,0.20,0.75))
		summary(mutation_factor)
		tapply(m,mutation_factor,sum)
		
		boxplot(k0f~mutation_factor,outline=F,xlab="Mutation Rate",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,k4f,method="spearman")
		kruskal.test(k0f~mutation_factor)	

	### mutation ins ###
	summary(kins)
	quantile(kins,(0:3)/3)
	
		mutation_factor = cut(kins,c(0,0.18,0.25,0.7))
		summary(mutation_factor)
		tapply(m,mutation_factor,sum)
		
		boxplot(k0f~mutation_factor,outline=F,xlab="Mutation Rate",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,kins,method="spearman")
		kruskal.test(k0f~mutation_factor)			
	
##### THETA #####
 
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
dom #5.4, media harmónica

K_watterson <- (data[["seg_ins"]])/(data[["mdmel_ins"]])

watterson = K_watterson/dom
head(watterson)

summary(watterson)
quantile(watterson ,(0:3)/3)

watterson_factor = cut(watterson,c(0,0.005,0.013,0.087))
summary(watterson_factor)
tapply(m,watterson_factor,sum)

boxplot(k0f~watterson_factor,outline=F,xlab="Watterson estimator",ylab="Ka")
abline(h=median(k0f),col="black")

cor.test(k0f,watterson,method="spearman")
kruskal.test(k0f~watterson_factor)

qplot(Chr_domain,Recombination, geom=("boxplot"),fill=Chromosome, binwidth=3) + xlab("Chromatin domain") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right") + ggtitle("R")
