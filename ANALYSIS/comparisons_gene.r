library(car)
library(lattice)
library(ppcor)

rm(list = ls())

### FILES GENE ###

data<-read.table(file="sizebr",header=TRUE,sep="\t")
data<-read.table(file="tranbr",header=TRUE,sep="\t")
data<-read.table(file="distbr",header=TRUE,sep="\t")
data<-read.table(file="exonbr",header=TRUE,sep="\t")
data<-read.table(file="exprbr",header=TRUE,sep="\t")
data<-read.table(file="dombr",header=TRUE,sep="\t")
data<-read.table(file="mcombr",header=TRUE,sep="\t")

data<-read.table(file="sizetran",header=TRUE,sep="\t")
data<-read.table(file="sizedist",header=TRUE,sep="\t")
data<-read.table(file="sizeexon",header=TRUE,sep="\t")
data<-read.table(file="sizeexpr",header=TRUE,sep="\t")
data<-read.table(file="sizedom",header=TRUE,sep="\t")
data<-read.table(file="sizemcomp",header=TRUE,sep="\t")

data<-read.table(file="trandist",header=TRUE,sep="\t")
data<-read.table(file="tranexon",header=TRUE,sep="\t")
data<-read.table(file="tranexpr",header=TRUE,sep="\t")
data<-read.table(file="trandom",header=TRUE,sep="\t")
data<-read.table(file="tranmcomp",header=TRUE,sep="\t")

data<-read.table(file="distexon",header=TRUE,sep="\t")
data<-read.table(file="distexpr",header=TRUE,sep="\t")
data<-read.table(file="distdom",header=TRUE,sep="\t")
data<-read.table(file="distmcomp",header=TRUE,sep="\t")

data<-read.table(file="exonexpr",header=TRUE,sep="\t")
data<-read.table(file="exondom",header=TRUE,sep="\t")
data<-read.table(file="exonmcomp",header=TRUE,sep="\t")

data<-read.table(file="exprdom",header=TRUE,sep="\t")
data<-read.table(file="exprmcomp",header=TRUE,sep="\t")

data<-read.table(file="dommcomp",header=TRUE,sep="\t")

### FEATURES ###

data[["mutation"]] = factor(data[["mutation"]], levels=c("low","medium","high"))
data[["recombination"]] = factor(data[["recombination"]], levels=c("low","medium","high"))
data[["breadth"]] = factor(data[["breadth"]], levels=c("low","medium","high"))
data[["size"]] = factor(data[["size"]], levels=c("low","medium","high"))
data[["transcripts"]] = factor(data[["transcripts"]], levels=c("low","medium","high"))
data[["distance"]] = factor(data[["distance"]], levels=c("low","medium","high"))
data[["exons"]] = factor(data[["exons"]], levels=c("low","medium","high"))
data[["expression"]] = factor(data[["expression"]], levels=c("low","medium","high"))
data[["mcomplexity"]] = factor(data[["mcomplexity"]], levels=c("low","medium","high"))
data[["domain"]] = factor(data[["domain"]], levels=c("active","both","inactive"))

## 1. SIZE BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 

## 2. TRANSCRIPTS BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 

## 3. DISTANCE BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 

## 4. EXONS BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 

## 5. EXPRESSION BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 

## 6. DOMAIN BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 

## 7. MESSENGER COMPLEXITY BREADTH
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["mcomplexity"]] + data[["breadth"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["mcomplexity"]] + data[["breadth"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["mcomplexity"]] + data[["breadth"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["mcomplexity"]] + data[["breadth"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["mcomplexity"]] + data[["breadth"]])
summary(mod2) 
Anova(mod2,type="2") 

## 8. SIZE TRANSCRIPTS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["transcripts"]])
summary(mod2) 
Anova(mod2,type="2") 

## 9. SIZE DISTANCE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["distance"]])
summary(mod2) 
Anova(mod2,type="2") 

## 10. SIZE EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["exons"]])
summary(mod2) 
Anova(mod2,type="2") 

## 11. SIZE EXPRESION
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["expression"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["expression"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["expression"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["expression"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["expression"]])
summary(mod2) 
Anova(mod2,type="2") 


## 12. SIZE DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 13. SIZE MESSENGER-COMPLEXITY
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["mcomplexity"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["mcomplexity"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["size"]] + data[["mcomplexity"]])
summary(mod2) 
Anova(mod2,type="2") 

## 14. TRANSCRIPTS DISTANCE
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["distance"]])
summary(mod2) 
Anova(mod2,type="2") 

## 15. TRANSCRIPTS EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["exons"]])
summary(mod2) 
Anova(mod2,type="2") 

## 16. TRANSCRIPTS EXPRESSION
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["expression"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["expression"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["expression"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["expression"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["expression"]])
summary(mod2) 
Anova(mod2,type="2") 

## 17. TRANSCRIPTS DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 18. TRANSCRIPTS MESSENGER COMPLEXITY
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["mcomplexity"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["mcomplexity"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["transcripts"]] + data[["mcomplexity"]])
summary(mod2) 
Anova(mod2,type="2") 

## 19. DISTANCE EXONS
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["exons"]])
summary(mod2) 
Anova(mod2,type="2") 


## 20. DISTANCE EXPRESSION
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["expression"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["expression"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["expression"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["expression"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["expression"]])
summary(mod2) 
Anova(mod2,type="2") 

## 21. DISTANCE DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 22. DISTANCE MESSENGER COMPLEXITY
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["mcomplexity"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["mcomplexity"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["distance"]] + data[["mcomplexity"]])
summary(mod2) 
Anova(mod2,type="2") 

## 23. EXONS EXPRESSION
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["expression"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["expression"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["expression"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["expression"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["expression"]])
summary(mod2) 
Anova(mod2,type="2") 

## 24. EXONS DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 25. EXONS MESSENGER COMPLEXITY
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["mcomplexity"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["mcomplexity"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["exons"]] + data[["mcomplexity"]])
summary(mod2) 
Anova(mod2,type="2") 

## 26. EXPRESSION DOMAIN
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["domain"]])
summary(mod2) 
Anova(mod2,type="2") 

## 27. EXPRESSION MESSENGER COMPLEXITY
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["mcomplexity"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["mcomplexity"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["expression"]] + data[["mcomplexity"]])
summary(mod2) 
Anova(mod2,type="2")

## 28. DOMAIN MESSENGER COMPLEXITY
# analysis
mod2=lm(data[["Ka"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_pos"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["mcomplexity"]])
mod2=lm(data[["Ka_neg"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["mcomplexity"]])
mod2=lm(data[["alpha"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["mcomplexity"]])
mod2=lm(data[["omega_A"]] ~ data[["chr"]] + data[["recombination"]] + data[["mutation"]] + data[["domain"]] + data[["mcomplexity"]])
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

boxplot(data[["Ka_pos"]]~data[["mcomplexity"]],main="Messenger complexity",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["mcomplexity"]],main="Messenger complexity",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["mcomplexity"]],main="Messenger complexity",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["mcomplexity"]],main="Messenger complexity",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["mcomplexity"]],main="Messenger complexity",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["mcomplexity"]],main="Messenger complexity",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["expression"]],main="Expression",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["expression"]],main="Expression",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["expression"]],main="Expression",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["expression"]],main="Expression",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["expression"]],main="Expression",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["expression"]],main="Expression",ylab="Ks",outline=F)

boxplot(data[["Ka_pos"]]~data[["domain"]],main="Domain",ylab="Ka_pos",outline=F)
boxplot(data[["Ka"]]~data[["domain"]],main="Domain",ylab="Ka",outline=F)
boxplot(data[["Ka_neg"]]~data[["domain"]],main="Domain",ylab="Ka_neg",outline=F)
boxplot(data[["omega_A"]]~data[["domain"]],main="Domain",ylab="omega_A",outline=F)
boxplot(data[["alpha"]]~data[["domain"]],main="Domain",ylab="alpha",outline=F)
boxplot(data[["Ks"]]~data[["domain"]],main="Domain",ylab="Ks",outline=F)