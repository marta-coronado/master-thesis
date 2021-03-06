library(car)
library(lattice)
library(ppcor)

rm(list = ls())
# gene features
data<-read.table(file="gene_M",header=TRUE,sep="\t")  # breadth, expression
data<-read.table(file="theta",header=TRUE,sep="\t")   # theta 
data<-read.table(file="gene_A",header=TRUE,sep="\t")  # chr  recombination	mutation	 
data<-read.table(file="gene_E",header=TRUE,sep="\t")  #  
data<-read.table(file="gene_F",header=TRUE,sep="\t")  # chr  recombination	mutation	state
data<-read.table(file="gene_G",header=TRUE,sep="\t")  #  chr  recombination	mutation	domain


####
# exon features
data<-read.table(file="order1",header=TRUE,sep="\t") 
data<-read.table(file="size1",header=TRUE,sep="\t") 
data<-read.table(file="transcript1",header=TRUE,sep="\t") 
data<-read.table(file="erep1",header=TRUE,sep="\t") 
data<-read.table(file="domain1",header=TRUE,sep="\t") 

# div0 fold > 0
data<-read.table(file="order2solo",header=TRUE,sep="\t") 
data<-read.table(file="size2",header=TRUE,sep="\t") 
data<-read.table(file="transcript2",header=TRUE,sep="\t") 
data<-read.table(file="erep2",header=TRUE,sep="\t") 
data<-read.table(file="domain2",header=TRUE,sep="\t")
data<-read.table(file="distance",header=TRUE,sep="\t") 
data<-read.table(file="expr_ex",header=TRUE,sep="\t") 
data<-read.table(file="breadth_ex",header=TRUE,sep="\t") 

data <- subset(data,chr == "X") # X
data <- subset(data,chr != "X") # Autosomas

data[["mutation"]] = factor(data[["mutation"]], levels=c("low","medium","high"))
data[["recombination"]] = factor(data[["recombination"]], levels=c("low","medium","high"))
data[["order"]] = factor(data[["order"]], levels=c("low","medium","high","sup"))
data[["order"]] = factor(data[["order2"]], levels=c("low","high"))
data[["num_exon"]] = factor(data[["exons"]], levels=c("low","medium","high"))
data[["size"]] = factor(data[["size"]], levels=c("low","medium","high"))
data[["transcript"]] = factor(data[["transcript"]], levels=c("low","medium","high"))
data[["erep"]] = factor(data[["erep"]], levels=c("low","high"))
data[["domain"]] = factor(data[["domain"]], levels=c("active","both","inactive"))
data[["distance"]] = factor(data[["distance"]], levels=c("low","medium","high"))
data[["expression"]] = factor(data[["expression"]], levels=c("low","medium","high"))
data[["breadth"]] = factor(data[["breadth"]], levels=c("low","medium","high"))


# sólo autosomas
data <- subset(data,chr != "X")

##### features
## order
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["order"]],main="Order",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["order"]],main="Order",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["order"]],main="Order",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["order"]],main="Order",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["order"]],main="Order",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["order"]],main="Order",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

# size
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["size"]],main="Size",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["size"]],main="Size",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["size"]],main="Size",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["size"]],main="Size",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["size"]],main="Size",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["size"]],main="Size",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

# transcripts
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcript"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcript"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcript"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcript"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcript"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcript"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["transcript"]],main="Transcript",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["transcript"]],main="Transcript",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["transcript"]],main="Transcript",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["transcript"]],main="Transcript",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["transcript"]],main="Transcript",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["transcript"]],main="Transcript",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

# inclusion levels
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["erep"]],main="Inclusion levels",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["erep"]],main="Inclusion levels",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["erep"]],main="Inclusion levels",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["erep"]],main="Inclusion levels",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["erep"]],main="Inclusion levels",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["erep"]],main="Inclusion levels",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

# domain
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["domain"]],main="Domain",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["domain"]],main="Domain",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["domain"]],main="Domain",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["domain"]],main="Domain",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["domain"]],main="Domain",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["domain"]],main="Domain",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

# order and num_exons
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["order"]],main="Order",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["order"]],main="Order",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["order"]],main="Order",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["order"]],main="Order",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["order"]],main="Order",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["order"]],main="Order",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

boxplot(data[["Ka"]]~data[["num_exon"]],main="Number of exons",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["num_exon"]],main="Number of exons",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["num_exon"]],main="Number of exons",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["num_exon"]],main="Number of exons",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["num_exon"]],main="Number of exons",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["num_exon"]],main="Number of exons",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")


#num exons

mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["num_exon"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]]  + data[["num_exon"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]]  + data[["num_exon"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["num_exon"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["num_exon"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["num_exon"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["num_exon"]],main="Number of exons",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["num_exon"]],main="Number of exons",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["num_exon"]],main="Number of exons",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["num_exon"]],main="Number of exons",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["num_exon"]],main="Number of exons",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["num_exon"]],main="Number of exons",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

#distance

mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]]  + data[["distance"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]]  + data[["distance"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["distance"]],main="Distance",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["distance"]],main="Distance",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["distance"]],main="Distance",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["distance"]],main="Distance",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["distance"]],main="Distance",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["distance"]],main="Distance",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")


# order only 2 categories
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["num_exon"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["order"]],main="Order",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["order"]],main="Order",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["order"]],main="Order",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["order"]],main="Order",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["order"]],main="Order",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["order"]],main="Order",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

boxplot(data[["Ka"]]~data[["num_exon"]],main="Number of exons",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["num_exon"]],main="Number of exons",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["num_exon"]],main="Number of exons",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["num_exon"]],main="Number of exons",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["num_exon"]],main="Number of exons",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["num_exon"]],main="Number of exons",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

# expression
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["expression"]],main="Expression",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["expression"]],main="Expression",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["expression"]],main="Expression",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["expression"]],main="Expression",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["expression"]],main="Expression",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["expression"]],main="Expression",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

# breadth
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["breadth"]],main="Expression bias",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["breadth"]],main="Expression bias",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["breadth"]],main="Expression bias",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["breadth"]],main="Expression bias",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["breadth"]],main="Expression bias",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["breadth"]],main="Expression bias",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

# order and size
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["size"]])
mod2=lm(data[["Ks"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["size"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["size"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["size"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["size"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["order"]] + data[["size"]])
summary(mod2) 
Anova(mod2,type="2") 
# boxplots
boxplot(data[["Ka"]]~data[["order"]],main="Order",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["order"]],main="Order",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["order"]],main="Order",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["order"]],main="Order",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["order"]],main="Order",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["order"]],main="Order",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")

boxplot(data[["Ka"]]~data[["size"]],main="Size",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["size"]],main="Size",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["size"]],main="Size",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["alpha"]]~data[["size"]],main="Size",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["omega_A"]]~data[["size"]],main="Size",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["Ks"]]~data[["size"]],main="Size",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")




### Polymorphism and m ###

theta <- data[["segins"]]/data[["mdmelins"]]/5.42533

boxplot(theta~data[["recombination"]],main="Recombination",ylab="Theta",outline=F)
boxplot(theta~data[["mutation"]],main="Mutation",ylab="Theta",outline=F)
boxplot(theta~data[["chr"]],main="Chromosomes",ylab="Theta",outline=F)

boxplot(data[["mdmel0f"]]~data[["recombination"]],main="Recombination",ylab="Supergene size",outline=F)
boxplot(data[["mdmel0f"]]~data[["mutation"]],main="Mutation",ylab="Supergene size",outline=F)
boxplot(data[["mdmel0f"]]~data[["chr"]],main="Chromosome",ylab="Supergene size",outline=F)


mod2=lm(Ka_pos ~recombination + mutation + chr + theta + mdmel0f,data)
mod2=lm(Ka_pos ~recombination + mutation + chr,data)

summary(mod2) 
Anova(mod2,type="2") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))
	
	
### Alpha ###
summary(data[["alpha"]])
mod=lm(alpha ~1,data)
mod2=step(mod, scope=alpha ~ recombination * mutation + chr) 
mod2=lm(alpha ~recombination * mutation * chr,data)

summary(mod2) 
Anova(mod2,type="3") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))

boxplot(1-data[["alpha"]]~data[["recombination"]],main="Recombination",ylab="1-Alpha",outline=F)
	
boxplot(data[["alpha"]]~data[["recombination"]],main="Recombination",ylab="Alpha",outline=F)
boxplot(data[["alpha"]]~data[["mutation"]],main="Mutation",ylab="Alpha",outline=F)
boxplot(data[["alpha"]]~data[["chr"]],main="Chromosomes",ylab="Alpha",outline=F)
kruskal.test(data[["alpha"]]~data[["recombination"]])

### Omega_A ###

mod=lm(omega_A ~1,data)
mod2=step(mod, scope=omega_A ~ recombination * mutation + chr) 
mod2=lm(omega_A ~recombination + mutation + chr,data)
mod2=aov(omega_A ~recombination + mutation + chr,data)

TukeyHSD(mod2)
summary(mod2) 
Anova(mod2,type="2") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))

boxplot(data[["Ka_neg"]]/data[["Ks"]]~data[["recombination"]],main="Recombination",ylab="omega_D",outline=F)

boxplot(data[["omega_A"]]~data[["recombination"]],main="Recombination",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["mutation"]],main="Mutation",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["chr"]],main="Chromosomes",ylab="omega_A",outline=F)
kruskal.test(data[["omega_A"]]~data[["recombination"]])
cor.test(data[["omega_A"]],data[["recombination"]])

### Ka+ ###

mod=lm(Ka_pos~1,data)
mod2=step(mod, scope=Ka_pos ~ recombination * mutation * chr) 
mod2=lm(Ka_pos~recombination + mutation + chr,data)

summary(mod2) 
Anova(mod2,type="3") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))
	
boxplot(data[["Ka_pos"]]~data[["recombination"]],main="Recombination",ylab="Ka+")
boxplot(data[["Ka_pos"]]~data[["mutation"]],main="Mutation",ylab="Ka+")
boxplot(data[["Ka_pos"]]~data[["chr"]],main="Chromosomes",ylab="Ka+")
kruskal.test(data[["Ka_pos"]]~data[["recombination"]])

### Ka ###

mod=lm(Ka~1,data)
mod2=step(mod, scope=Ka ~ recombination * mutation * chr) 
mod2=lm(Ka~recombination + mutation + chr,data)

summary(mod2) 
Anova(mod2,type="3") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))
	
boxplot(data[["Ka"]]~data[["recombination"]],main="Recombination",ylab="Ka")
boxplot(data[["Ka"]]~data[["mutation"]],main="Mutation",ylab="Ka")
boxplot(data[["Ka"]]~data[["chr"]],main="Chromosomes",ylab="Ka")

### Ka- ###

Ka_neg <- data[["Ka"]]-data[["Ka_pos"]]

mod=lm(Ka_neg~1,data)
mod2=step(mod, scope=Ka_neg ~ data$recombination * data$mutation + data$chr) 
mod2=lm(Ka_neg~ data$recombination + data$mutation + data$chr)
mod2=aov(Ka_neg/data[["Ks"]]~ data$recombination + data$mutation + data$chr)
mod2=aov(Ka_neg~ data$recombination + data$mutation + data$chr)

TukeyHSD(mod2)

summary(mod2) 
Anova(mod2,type="2") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))
	
boxplot(Ka_neg~data[["recombination"]],main="Recombination",ylab="Ka-/0")
boxplot(Ka_neg~data[["mutation"]],main="Mutation",ylab="Ka-/0")
boxplot(Ka_neg~data[["chr"]],main="Chromosomes",ylab="Ka-/0")
kruskal.test(data[[Ka_neg]]~data[["recombination"]])
boxplot(Ka_neg/data$Ks~data[["mutation"]],main="Mutation",ylab="omega_D")


### Ks ###

mod=lm(data[["Ks"]]~1,data)
mod2=step(mod, scope=data[["Ks"]] ~ data$recombination * data$mutation + data$chr) 
mod2=lm(data[["Ks"]]~ data$recombination + data$chr)

summary(mod2) 
Anova(mod2,type="3") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))
	
boxplot(data[["Ks"]]~data[["recombination"]],main="Recombination",ylab="Ks")
boxplot(data[["Ks"]]~data[["chr"]],main="Chromosomes",ylab="Ks")

####

data<-read.table(file="gene_E",header=TRUE,sep="\t") #gene features
data<-read.table(file="gene_G",header=TRUE,sep="\t") #gene chromatin domain
data<-read.table(file="gene_F",header=TRUE,sep="\t") #gene chromatin state
data<-read.table(file="gene_M",header=TRUE,sep="\t") #gene MARTA

nrow(data)
summary(data)
data <- na.omit(data) 

### FEATURES ###
data <- subset(data, mdmelins > 1379)
data[["mutation"]] = factor(data[["mutation"]], levels=c("low","medium","high"))
data[["recombination"]] = factor(data[["recombination"]], levels=c("low","medium","high"))
data[["transcripts"]] = factor(data[["transcripts"]], levels=c("low","medium","high"))
data[["size"]] = factor(data[["size"]], levels=c("low","medium","high"))
data[["distance"]] = factor(data[["distance"]], levels=c("low","medium","high"))
data[["exons"]] = factor(data[["exons"]], levels=c("low","medium","high"))
data[["breadth"]] = factor(data[["breadth"]], levels=c("low","medium","high"))
data[["expression"]] = factor(data[["expression"]], levels=c("low","medium","high"))
data[["mcomplexity"]] = factor(data[["mcomplexity"]], levels=c("low","medium","high"))
theta <- data[["segins"]]/data[["mdmelins"]]/5.42533

### DOMAIN ###
data <- subset(data, mdmelins > 1064)
data[["mutation"]] = factor(data[["mutation"]], levels=c("low","medium","high"))
data[["recombination"]] = factor(data[["recombination"]], levels=c("low","medium","high"))
data[["domain"]] = factor(data[["domain"]], levels=c("active","both","inactive"))
theta <- data[["segins"]]/data[["mdmelins"]]/5.42533

### STATE ###
data[["mutation"]] = factor(data[["mutation"]], levels=c("low","medium","high"))
data[["recombination"]] = factor(data[["recombination"]], levels=c("low","medium","high"))
data[["state"]] = factor(data[["state"]], levels=c("A","B","C","D","E","F","G","H","I"))


### FEATURES ANALYSIS ###

mod=lm(Ka_pos~1,data)
mod2=step(mod, scope= Ka_pos ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["breadth"]]) 

mod2=lm(data[["Ka_pos"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["transcripts"]] + data[["Ka_neg"]])
mod2=lm(data[["Ka_pos"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["exons"]] + data[["Ka_neg"]])
mod2=lm(data[["Ka_pos"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["size"]] + data[["Ka_neg"]])
mod2=lm(data[["Ka_pos"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["distance"]] + data[["Ka_neg"]])

mod2=lm(data[["Ka_pos"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["breadth"]] + data[["Ka_neg"]])
mod2=lm(data[["Ka_neg"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["breadth"]])

mod2=lm(data[["Ka_pos"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["expression"]] + data[["Ka_neg"]])
mod2=lm(data[["Ka_pos"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["mcomplexity"]] + data[["Ka_neg"]])


summary(mod2) 
Anova(mod2,type="2") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))
	
cor.test(data[["Ka_pos"]],data[["Ka_neg"]],method="spearman")	#-0.26***
cor.test(data[["Ka_pos"]],data[["Ka"]],method="spearman") #+0.45
cor.test(data[["Ka_neg"]],data[["Ka"]],method="spearman") #+0.67

gene <- data.frame(
                                Ks = data[["Ks"]],
                                Ka = data[["Ka"]],                                
                                Ka_pos = data[["Ka_pos"]],
                                Ka_neg = data[["Ka_neg"]]                                                                                       
                          )
gene <- na.omit(gene)
nrow(gene)
spcor(gene,method=c("spearman"))


gene <- data.frame(
                                Breadth = data[["breadth"]],
                                Ka_pos = data[["Ka_pos"]],
                                Ka_neg = data[["Ka_neg"]]                                                                                       
                          )
gene <- na.omit(gene)
nrow(gene)
spcor(gene,method=c("spearman"))
cor.test(gene[["Ka_pos"]],gene[["Ka_neg"]],method="spearman") #ns

boxplot(gene[["Ka_pos"]]~gene[["Breadth"]])
boxplot(gene[["Ka_neg"]]~gene[["Breadth"]])
quantile(gene[["Ka_neg"]],(0:3)/3)
Ka_neg = cut(gene[["Ka_neg"]],quantile(gene[["Ka_neg"]],(0:3)/3))
boxplot(gene[["Ka_pos"]]~Ka_neg)

### DOMAIN ANALYSIS ###

mod=lm(Ka_pos~1,data)
mod2=step(mod, scope= Ka_pos ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["domain"]]) 
mod2=lm(data[["Ka_pos"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["recombination"]] + data[["mutation"]] + data[["chr"]] + data[["domain"]])

### MARTA ANALYSIS ###
# efecte de theta:
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["theta"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["theta"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["theta"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["theta"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["theta"]])
summary(mod2) 
Anova(mod2,type="2") 
# gràfics per theta
boxplot(data[["Ka_pos"]]~data[["theta"]],main="Theta",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["theta"]],main="Theta",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["theta"]],main="Theta",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["theta"]],main="Theta",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["theta"]],main="Theta",ylab="alpha",outline=F)


# efecte de breadth
mod2=lm(data[["Ka"]] ~ data[["breadth"]] + data[["expression"]] + data[["chr"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["breadth"]] + data[["expression"]] + data[["chr"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["breadth"]] + data[["expression"]] + data[["chr"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["breadth"]] + data[["expression"]] + data[["chr"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["breadth"]] + data[["expression"]] + data[["chr"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

boxplot(data[["Ka"]]~data[["breadth"]],main="breadth",ylab="Ka",outline=F)
boxplot(data[["Ka_pos"]]~data[["breadth"]],main="breadth",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_neg"]]~data[["breadth"]],main="breadth",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["breadth"]],main="breadth",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["breadth"]],main="breadth",ylab="alpha",outline=F)



summary(mod2) 
Anova(mod2,type="2") 
	par(mfrow=c(2,2))
	plot(mod2)
	par(mfrow=c(1,1))	
	
### BOXPLOTS ###

boxplot(data[["omega_A"]]~data[["recombination"]],main="Recombination",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["mutation"]],main="Mutation",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["chr"]],main="Chromosomes",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")
boxplot(data[["omega_A"]]~data[["transcripts"]],main="Transcripts",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["size"]],main="Protein Lenght",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["distance"]],main="Distance",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["exons"]],main="Exons",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["breadth"]],main="Breadht",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["expression"]],main="Expression",ylab="omega_A",outline=F)
boxplot(data[["omega_A"]]~data[["mcomplexity"]],main="#Transcripts/#Exons",ylab="omega_A",outline=F)

boxplot(data[["omega_A"]]~data[["domain"]],main="Chromatin Domain",ylab="omega_A",outline=F)

boxplot(data[["omega_A"]]~data[["state"]],main="Chromatin State",ylab="omega_A",outline=F)
abline(h=median(data[["omega_A"]]),col="black")

boxplot(data[["alpha"]]~data[["recombination"]],main="Recombination",ylab="alpha",outline=F)
boxplot(data[["alpha"]]~data[["mutation"]],main="Mutation",ylab="alpha",outline=F)
boxplot(data[["alpha"]]~data[["chr"]],main="Chromosomes",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")
boxplot(data[["alpha"]]~data[["transcripts"]],main="Transcripts",ylab="alpha",outline=F)
boxplot(data[["alpha"]]~data[["size"]],main="Protein Lenght",ylab="alpha",outline=F)
boxplot(data[["alpha"]]~data[["distance"]],main="Distance",ylab="alpha",outline=F)
boxplot(data[["alpha"]]~data[["exons"]],main="Exons",ylab="alpha",outline=F)
boxplot(data[["alpha"]]~data[["breadth"]],main="Breadht",ylab="alpha",outline=F)
boxplot(data[["alpha"]]~data[["expression"]],main="Expression",ylab="alpha",outline=F)
boxplot(data[["alpha"]]~data[["mcomplexity"]],main="#Transcripts/#Exons",ylab="alpha",outline=F)

boxplot(data[["alpha"]]~data[["domain"]],main="Chromatin Domain",ylab="alpha",outline=F)

boxplot(data[["alpha"]]~data[["state"]],main="Chromatin State",ylab="alpha",outline=F)
abline(h=median(data[["alpha"]]),col="black")

boxplot(data[["Ka_pos"]]~data[["recombination"]],main="Recombination",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_pos"]]~data[["mutation"]],main="Mutation",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_pos"]]~data[["chr"]],main="Chromosomes",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")
boxplot(data[["Ka_pos"]]~data[["transcripts"]],main="Transcripts",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_pos"]]~data[["size"]],main="Protein Lenght",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_pos"]]~data[["distance"]],main="Distance",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_pos"]]~data[["exons"]],main="Exons",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_pos"]]~data[["breadth"]],main="Breadht",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_pos"]]~data[["expression"]],main="Expression",ylab="Ka_pos",outline=F)
boxplot(data[["Ka_pos"]]~data[["mcomplexity"]],main="#Transcripts/#Exons",ylab="Ka_pos",outline=F)

boxplot(data[["Ka_pos"]]~data[["domain"]],main="Chromatin Domain",ylab="Ka_pos",outline=F)

boxplot(data[["Ka_pos"]]~data[["state"]],main="Chromatin State",ylab="Ka_pos",outline=F)
abline(h=median(data[["Ka_pos"]]),col="black")

boxplot(data[["Ka_neg"]]~data[["recombination"]],main="Recombination",ylab="Ka_neg",outline=F)
boxplot(data[["Ka_neg"]]~data[["mutation"]],main="Mutation",ylab="Ka_neg",outline=F)
boxplot(data[["Ka_neg"]]~data[["chr"]],main="Chromosomes",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")
boxplot(data[["Ka_neg"]]/data[["Ks"]]~data[["chr"]],main="Chromosomes",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]/data[["Ks"]]),col="black")
boxplot(data[["Ka_neg"]]~data[["transcripts"]],main="Transcripts",ylab="Ka_neg",outline=F)
boxplot(data[["Ka_neg"]]~data[["size"]],main="Protein Lenght",ylab="Ka_neg",outline=F)
boxplot(data[["Ka_neg"]]~data[["distance"]],main="Distance",ylab="Ka_neg",outline=F)
boxplot(data[["Ka_neg"]]~data[["exons"]],main="Exons",ylab="Ka_neg",outline=F)
boxplot(data[["Ka_neg"]]~data[["breadth"]],main="Breadht",ylab="Ka_neg",outline=F)
boxplot(data[["Ka_neg"]]~data[["expression"]],main="Expression",ylab="Ka_neg",outline=F)
boxplot(data[["Ka_neg"]]~data[["mcomplexity"]],main="#Transcripts/#Exons",ylab="Ka_neg",outline=F)

boxplot(data[["Ka_neg"]]~data[["domain"]],main="Chromatin Domain",ylab="Ka_neg",outline=F)

boxplot(data[["Ka_neg"]]~data[["state"]],main="Chromatin State",ylab="Ka_neg",outline=F)
abline(h=median(data[["Ka_neg"]]),col="black")

boxplot(data[["Ka"]]~data[["recombination"]],main="Recombination",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["mutation"]],main="Mutation",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["chr"]],main="Chromosomes",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["transcripts"]],main="Transcripts",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["size"]],main="Protein Lenght",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["distance"]],main="Distance",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["exons"]],main="Exons",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["breadth"]],main="Breadht",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["expression"]],main="Expression",ylab="Ka",outline=F)
boxplot(data[["Ka"]]~data[["mcomplexity"]],main="#Transcripts/#Exons",ylab="Ka",outline=F)

boxplot(data[["Ka"]]~data[["domain"]],main="Chromatin Domain",ylab="Ka",outline=F)

boxplot(data[["Ka"]]~data[["state"]],main="Chromatin State",ylab="Ka",outline=F)
abline(h=median(data[["Ka"]]),col="black")

boxplot(data[["Ks"]]~data[["recombination"]],main="Recombination",ylab="Ks",outline=F)
boxplot(data[["Ks"]]~data[["mutation"]],main="Mutation",ylab="Ks",outline=F)
boxplot(data[["Ks"]]~data[["chr"]],main="Chromosomes",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")
boxplot(data[["Ks"]]~data[["transcripts"]],main="Transcripts",ylab="Ks",outline=F)
boxplot(data[["Ks"]]~data[["size"]],main="Protein Lenght",ylab="Ks",outline=F)
boxplot(data[["Ks"]]~data[["distance"]],main="Distance",ylab="Ks",outline=F)
boxplot(data[["Ks"]]~data[["exons"]],main="Exons",ylab="Ks",outline=F)
boxplot(data[["Ks"]]~data[["breadth"]],main="Breadht",ylab="Ks",outline=F)
boxplot(data[["Ks"]]~data[["expression"]],main="Expression",ylab="Ks",outline=F)
boxplot(data[["Ks"]]~data[["mcomplexity"]],main="#Transcripts/#Exons",ylab="Ks",outline=F)

boxplot(data[["Ks"]]~data[["domain"]],main="Chromatin Domain",ylab="Ks",outline=F)

boxplot(data[["Ks"]]~data[["state"]],main="Chromatin State",ylab="Ks",outline=F)
abline(h=median(data[["Ks"]]),col="black")
