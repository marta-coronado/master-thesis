library(car)
library(lattice)
library(ppcor)
library(pls) 

rm(list = ls())

data<-read.table(file="EXON",header=TRUE,sep="\t")
#data<-read.table(file="EXON.SHORTINTRON.65",header=TRUE,sep="\t")
#data<-read.table(file="EXON.SHORTINTRON.65_8-30",header=TRUE,sep="\t")

data<-na.omit(data) 
data <- subset(data,comeron_100kb > 0) # Con recombincion
summary(data) 
nrow(data) 

data <- subset(data,chr == "2L")
data <- subset(data,chr == "2R")
data <- subset(data,chr == "3L")
data <- subset(data,chr == "3R")
data <- subset(data,chr == "X") # X
data <- subset(data,chr != "X") # Autosomas

m <- data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]]
k4f <- (data[["div_4f"]])/(data[["mdyak_4f"]])
k0f <- (data[["div_0f"]]+1)/(data[["mdyak_0f"]])
k0f <- (data[["div_0f"]])/(data[["mdyak_0f"]])
kins <- (data[["div_ins"]])/(data[["mdyak_ins"]])

exon <- data.frame(
  Transcripts = data[["num_transcripts_exon"]] ,        
  Erep = (data[["num_transcripts_exon"]]/data[["num_transcripts_gene"]]),
  m = m,
  Distance = data[["distance"]],
  Breadth = data[["bias_dev"]],
  Expression = data[["max_dev"]],
  Ki = kins,
   AAsubstitutions = k0f,
  Mutation = k4f,   
  
Recombination=(data[["comeron_100kb"]]),
  Order = data[["order"]]
  #Chr_state = data[["state"]],                                
  #Chr_domain = data[["domain"]],
  #Chromosome = data[["chr"]],                                                                                          
 
  )                             
)
exon <- na.omit(exon)
nrow(exon)
spcor(exon,method=c("spearman"))



Transcripts = data[["num_transcripts_exon"]]     
Erep = (data[["num_transcripts_exon"]]/data[["num_transcripts_gene"]])
Distance = data[["distance"]]
AAsubstitutions = k0f
Mutation = k4f
Order = data[["order"]]
Chr_state = data[["state"]]                              
Chr_domain = data[["domain"]]
Chromosome = data[["chr"]]                                                                                        
Recombination=(data[["comeron_100kb"]])  
Num_Exons = data[["num_exons"]]
Breadth = data[["bias_dev"]]
Expression = data[["max_dev"]]

	### Size ###
	quantile(m,(0:3)/3)	
		
		m_factor = cut(m,c(0,149,328,14000))
		summary(m_factor)
		with(data,tapply(m,m_factor,sum))
		
		boxplot(k0f~m_factor,outline=F,xlab="Exon Size (bp)",ylab="Ka (All Exons)")
		boxplot(k0f~m_factor,outline=F,xlab="Exon Size (bp)",ylab="Ka (All Exons dN+1)")
		boxplot(k0f~m_factor,outline=F,xlab="Exon Size (bp)",ylab="Ka (Exons dN > 0)")
		abline(h=median(k0f),col="black")
		
		cor.test(log10(k0f),log10(m),method="spearman")
		plot(log10(m),log10(k0f))
		kruskal.test(k0f~m_factor)


	### Transcripts ###
	
	with(data,tapply(m,num_transcripts_exon,sum))
		
		num_transcripts_exon_factor = cut(data[["num_transcripts_exon"]],c(0,1,2,75))
		summary(num_transcripts_exon_factor)
		with(data,tapply(m,num_transcripts_exon_factor,sum))
		
		boxplot(k0f~num_transcripts_exon_factor,outline=F,xlab="Number of Transcripts/Exon",ylab="Ka")
		boxplot(k0f~data$num_transcripts_exon,outline=F,xlab="Number of Transcripts/Exon",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,data$num_transcripts_exon,method="spearman")
		kruskal.test(k0f~num_transcripts_exon_factor)


	### Inclusion Levels ### 
	
	Erep <-(data[["num_transcripts_exon"]]/data[["num_transcripts_gene"]])
	alternative = cut(Erep,c(0,0.99,1))
	with(data,tapply(m,alternative,sum))
	summary(alternative)
		
		boxplot(k0f~alternative,outline=F,xlab="Inclusion Levels",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,Erep,method="spearman") 
		kruskal.test(k0f~alternative)

	### Order ####
#	order = cut(data[["order"]],c(0,1,2,4,112))
order = cut(data[["order"]],c(0,1,112))
	tapply(m,order,sum)
	summary(order) 

		boxplot(k0f~order,outline=F,xlab="Order",ylab="Ka")
		abline(h=median(k0f),col="black")
		cor.test(k0f,data[["order"]],method="spearman") 
		kruskal.test(k0f~order)	
				
	### Chromatin State and Chromatin Domain ###
	
	boxplot(k0f~data[["state"]],outline=F,xlab="Chromatin State",ylab="Ka")
	abline(h=0.01235,col="black")
	boxplot(k0f~data[["domain"]],outline=F,xlab="Chromatin Domain",ylab="Ka")
	abline(h=0.01235,col="black")
	
	### Density ### 
	summary(data[["distance"]]/m)
	quantile(data[["distance"]]/m,(0:3)/3)	
	density = cut(data[["distance"]]/m,quantile(data[["distance"]]/m,(0:3)/3))
	
	summary(data[["distance"]])
	quantile(data[["distance"]],(0:3)/3)	
	#density = cut(data[["distance"]],c(0,100,1000,1000000))
	density = cut(data[["distance"]],c(0,128,703,200000))
	#density = cut(log10(data[["distance"]]/m),c(-2,-0.5,1,4))

	with(data,tapply(m,density,sum))
	summary(density)
		
		boxplot(k0f~density,outline=F,xlab="Density",ylab="Ka")
		abline(h=median(k0f),col="black")
		
		cor.test(k0f,data[["distance"]],method="spearman") 
		cor.test(data[["distance"]],m,method="spearman") 
				cor.test(k0f,m,method="spearman") 

		kruskal.test(k0f~density)	


### Recombination ###
summary(data[["comeron_100kb"]])
quantile(data[["comeron_100kb"]],(0:3)/3)

recombination_factor = cut(data[["comeron_100kb"]],c(0,1.42,2.9,15))
summary(recombination_factor)
tapply(m,recombination_factor,sum)

boxplot(k0f~recombination_factor,outline=F,xlab="Recombination Rate",ylab="Ka")
abline(h=median(k0f),col="black")

cor.test(k0f,data[["comeron_100kb"]],method="spearman")
kruskal.test(k0f~recombination_factor)
		
### mutation 4f ###
summary(k4f)
quantile(k4f,(0:3)/3)

mutation_factor = cut(k4f,c(0,0.12,0.20,0.75))
summary(mutation_factor)
tapply(m,mutation_factor,sum)

boxplot(k0f~mutation_factor,outline=F,xlab="Mutation Rate",ylab="Ka")
abline(h=median(k0f),col="black")

cor.test(k0f,k4f,method="spearman")
kruskal.test(k0f~mutation_factor)	

### mutation ins ###
summary(kins)
quantile(kins,(0:3)/3)

mutation_factor = cut(kins,c(0,0.19,0.25,0.7))
summary(mutation_factor)
tapply(m,mutation_factor,sum)

boxplot(k0f~mutation_factor,outline=F,xlab="Mutation Rate",ylab="Ka")
abline(h=median(k0f),col="black")

cor.test(k0f,kins,method="spearman")
kruskal.test(k0f~mutation_factor)		



### Chromosome ###

tapply(m,data[["chr"]],sum)
summary(data[["chr"]]) 

boxplot(k0f~data[["chr"]],outline=F,xlab="Chromosome",ylab="Ka")

qplot(Chr_domain,Recombination, geom=("boxplot"),fill=Chromosome, binwidth=3) + xlab("Chromatin domain") + ylab("Recombination") + scale_fill_brewer(type="seq", palette=3) + theme(axis.title=element_text(face="bold",size="12", color="black"), legend.position="right") + ggtitle("R")
kruskal.test(Chr_domain~Recombination)
ggsave(file="rec_dom.png")

### number of exons ### 
summary(data[["num_exons"]])
quantile(data[["num_exons"]],(0:3)/3)

numexons_factor = cut(data[["num_exons"]],c(0,5,9,114))
summary(numexons_factor)
tapply(m,numexons_factor,sum)

boxplot(k0f~numexons_factor,outline=F,xlab="Number of exons",ylab="Ka")
abline(h=median(k0f),col="black")

cor.test(k0f,data[["num_exons"]],method="spearman")
kruskal.test(k0f~numexons_factor)

### Expression ###
plot(data[["bias_dev"]],data[["max_dev"]])
summary(data[["max_dev"]])
quantile(data[["max_dev"]],(0:3)/3)

expression_factor = cut(data[["max_dev"]],c(0,1.45,1.81,4.4))
summary(expression_factor)
tapply(m,expression_factor,sum)

boxplot(k0f~expression_factor,outline=F,xlab="Expression",ylab="Ka")
abline(h=median(k0f),col="black")

cor.test(k0f,data[["max_dev"]],method="spearman")
kruskal.test(k0f~expression_factor)

### Breadth ###
summary(data[["bias_dev"]])
quantile(data[["bias_dev"]],(0:3)/3)

breadth_factor = cut(data[["bias_dev"]],c(0,0.25,0.47,1))
summary(breadth_factor)
tapply(m,breadth_factor,sum)

boxplot(k0f~breadth_factor,outline=F,xlab="Expression Breadth",ylab="Ka")
abline(h=median(k0f),col="black")

cor.test(k0f,data[["bias_dev"]],method="spearman")
kruskal.test(k0f~breadth_factor)
