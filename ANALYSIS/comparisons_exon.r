library(car)
library(lattice)
library(ppcor)
library(ggplot2)

rm(list = ls())

### FILES EXON ###
data<-read.table(file="order2size",header=TRUE,sep="\t")
data<-read.table(file="order2tran",header=TRUE,sep="\t")
data<-read.table(file="order2erep",header=TRUE,sep="\t")
data<-read.table(file="order2domain",header=TRUE,sep="\t")
data<-read.table(file="order2dist",header=TRUE,sep="\t")
data<-read.table(file="order2expr",header=TRUE,sep="\t")
data<-read.table(file="order2breadth",header=TRUE,sep="\t")

data<-read.table(file="order2solo",header=TRUE,sep="\t")

data<-read.table(file="exprtran",header=TRUE,sep="\t")
data<-read.table(file="exprsize",header=TRUE,sep="\t")
data<-read.table(file="exprerep",header=TRUE,sep="\t")
data<-read.table(file="exprbreadth",header=TRUE,sep="\t")
data<-read.table(file="exprdomain",header=TRUE,sep="\t")
data<-read.table(file="exprexon",header=TRUE,sep="\t")
data<-read.table(file="exprdist",header=TRUE,sep="\t")

data<-read.table(file="breadthsize",header=TRUE,sep="\t")
data<-read.table(file="breadthtran",header=TRUE,sep="\t")
data<-read.table(file="breadtherep",header=TRUE,sep="\t")
data<-read.table(file="breadthdomain",header=TRUE,sep="\t")
data<-read.table(file="breadthexon",header=TRUE,sep="\t")
data<-read.table(file="breadthdist",header=TRUE,sep="\t")

data<-read.table(file="sizetran",header=TRUE,sep="\t")
data<-read.table(file="sizeorder_BIEN",header=TRUE,sep="\t")
data<-read.table(file="sizerep",header=TRUE,sep="\t")
data<-read.table(file="sizedom_BIEN",header=TRUE,sep="\t")
data<-read.table(file="sizedist",header=TRUE,sep="\t")
data<-read.table(file="sizeexon",header=TRUE,sep="\t")

data<-read.table(file="tranorder_BIEN",header=TRUE,sep="\t")
data<-read.table(file="tranerep",header=TRUE,sep="\t")
data<-read.table(file="trandom",header=TRUE,sep="\t")
data<-read.table(file="trandist",header=TRUE,sep="\t")
data<-read.table(file="tranexon",header=TRUE,sep="\t")

data<-read.table(file="ordererep_BIEN",header=TRUE,sep="\t")
data<-read.table(file="orderdom_BIEN",header=TRUE,sep="\t")
data<-read.table(file="orderdist_BIEN",header=TRUE,sep="\t")
data<-read.table(file="orderexon_BIEN",header=TRUE,sep="\t")

data<-read.table(file="erepdom",header=TRUE,sep="\t")
data<-read.table(file="erepdist",header=TRUE,sep="\t")
data<-read.table(file="erepexon",header=TRUE,sep="\t")

data<-read.table(file="domdist",header=TRUE,sep="\t")
data<-read.table(file="domexon",header=TRUE,sep="\t")

data<-read.table(file="distexon",header=TRUE,sep="\t")

data_subset <- subset(data,breadth=="low")

mod2=lm(data_subset[["Ka_pos"]] ~ data_subset[["chr"]] + data_subset[["recombination"]] + data_subset[["mutation"]] + data_subset[["erep"]])
mod2=lm(data_subset[["Ka_neg"]] ~ data_subset[["chr"]] + data_subset[["recombination"]] + data_subset[["mutation"]] + data_subset[["erep"]])
summary(mod2) 
Anova(mod2,type="2") 

### FEATURES ###

data[["mutation"]] = factor(data[["mutation"]], levels=c("low","medium","high"))
data[["recombination"]] = factor(data[["recombination"]], levels=c("low","medium","high"))
data[["order"]] = factor(data[["order2"]], levels=c("low","high"))
data[["size"]] = factor(data[["size"]], levels=c("low","medium","high"))
data[["transcripts"]] = factor(data[["transcripts"]], levels=c("low","medium","high"))
data[["distance"]] = factor(data[["distance"]], levels=c("low","medium","high"))
data[["exons"]] = factor(data[["exons"]], levels=c("low","medium","high"))
data[["expression"]] = factor(data[["expression"]], levels=c("low","medium","high"))
data[["breadth"]] = factor(data[["breadth"]], levels=c("low","medium","high"))
data[["erep"]] = factor(data[["erep"]], levels=c("low","high"))
data[["domain"]] = factor(data[["domain"]], levels=c("active","both","inactive"))
data[["order"]] = factor(data[["order"]], levels=c("low","medium","high","sup"))

## 1. SIZE TRANSCRIPTS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
summary(mod2) 
Anova(mod2,type="2") 

## 2. SIZE ORDER
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2") 

## 3. SIZE EREP
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["erep"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["erep"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["erep"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["erep"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["erep"]])
summary(mod2) 
Anova(mod2,type="2") 

## 4. SIZE DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 5. SIZE DISTANCE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
summary(mod2) 
Anova(mod2,type="2") 

## 6. SIZE EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
summary(mod2) 
Anova(mod2,type="2") 

## 7. TRANSCRIPTS ORDER
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2") 

## 8. TRANSCRIPTS INCLUSION LEVELS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["erep"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["erep"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["erep"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["erep"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["erep"]])
summary(mod2) 
Anova(mod2,type="2")

## 9. TRANSCRIPTS DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 10. TRANSCRIPTS DISTANCE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
summary(mod2) 
Anova(mod2,type="2") 

## 11. TRANSCRIPTS EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
summary(mod2) 
Anova(mod2,type="2") 

## 12. ORDER INCLUSION LEVELS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2")

## 13. ORDER DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2")

## 14. ORDER DISTANCE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2")

## 15. ORDER EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2") 

## 15. ORDER EXPRESSION
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2") 

## 15. ORDER BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["order"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["order"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["order"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["order"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["order"]])
summary(mod2) 
Anova(mod2,type="2") 

## 16. INCLUSION LEVELS DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["erep"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 17. DISTANCE INCLUSION LEVELS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["erep"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["erep"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["erep"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["erep"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["erep"]])
summary(mod2) 
Anova(mod2,type="2") 

## 18. EXONS INCLUSION LEVELS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["erep"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["erep"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["erep"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["erep"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["erep"]])
summary(mod2) 
Anova(mod2,type="2") 

## 19. DISTANCE DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 20. DISTANCE EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
summary(mod2) 
Anova(mod2,type="2") 

## 21. EXONS DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 1. EXPRESSION TRANSCRIPTS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["transcripts"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["transcripts"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["transcripts"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["transcripts"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["transcripts"]])
summary(mod2) 
Anova(mod2,type="2") 

## 2. EXPRESSION SIZE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["size"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["size"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["size"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["size"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["size"]])
summary(mod2) 
Anova(mod2,type="2") 

## 3. EXPRESSION EREP
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["erep"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["erep"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["erep"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["erep"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["erep"]])
summary(mod2) 
Anova(mod2,type="2") 


## 4. EXPRESSION BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 


## 4. EXPRESSION DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 6. EXPRESSION EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["exons"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["exons"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["exons"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["exons"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["exons"]])
summary(mod2) 
Anova(mod2,type="2") 

## 5. EXPRESSION DISTANCE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["distance"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["distance"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["distance"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["distance"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["distance"]])
summary(mod2) 
Anova(mod2,type="2") 

## 2. BREADTH SIZE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["size"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["size"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["size"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["size"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["size"]])
summary(mod2) 
Anova(mod2,type="2") 

## 1. BREADTH TRANSCRIPTS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["transcripts"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["transcripts"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["transcripts"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["transcripts"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["transcripts"]])
summary(mod2) 
Anova(mod2,type="2") 

## 3. BREADTH EREP
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["erep"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["erep"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["erep"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["erep"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["erep"]])
summary(mod2) 
Anova(mod2,type="2") 

## 4. BREADTH DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 6. BREADTH EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["exons"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["exons"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["exons"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["exons"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["exons"]])
summary(mod2) 
Anova(mod2,type="2") 

## 5. BREADTH DISTANCE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["distance"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["distance"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["distance"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["distance"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["breadth"]] + data[["distance"]])
summary(mod2) 
Anova(mod2,type="2") 






# boxplots
boxplot(data[["Ka_pos"]]~data[["transcripts"]],main="Transcripts",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["transcripts"]],main="Transcripts",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["transcripts"]],main="Transcripts",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["transcripts"]],main="Transcripts",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["transcripts"]],main="Transcripts",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["transcripts"]],main="Transcripts",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["exons"]],main="Number of exons",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["exons"]],main="Number of exons",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["exons"]],main="Number of exons",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["exons"]],main="Number of exons",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["exons"]],main="Number of exons",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["exons"]],main="Number of exons",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["size"]],main="Size",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["size"]],main="Size",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["size"]],main="Size",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["size"]],main="Size",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["size"]],main="Size",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["size"]],main="Size",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["breadth"]],main="Breadth",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["breadth"]],main="Breadth",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["breadth"]],main="Breadth",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["breadth"]],main="Breadth",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["breadth"]],main="Breadth",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["breadth"]],main="Breadth",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["distance"]],main="Distance",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["distance"]],main="Distance",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["distance"]],main="Distance",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["distance"]],main="Distance",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["distance"]],main="Distance",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["distance"]],main="Distance",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["order"]],main="Order",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["order"]],main="Order",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["order"]],main="Order",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["order"]],main="Order",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["order"]],main="Order",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["order"]],main="Order",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["erep"]],main="Inclusion levels",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["erep"]],main="Inclusion levels",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["erep"]],main="Inclusion levels",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["erep"]],main="Inclusion levels",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["erep"]],main="Inclusion levels",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["erep"]],main="Inclusion levels",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["domain"]],main="Domain",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["domain"]],main="Domain",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["domain"]],main="Domain",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["domain"]],main="Domain",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["domain"]],main="Domain",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["domain"]],main="Domain",ylab="Ks",outline=F)
