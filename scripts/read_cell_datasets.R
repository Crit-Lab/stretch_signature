
library(limma)
library(MetaIntegrator)
library(COCONUT)
library(oligo)
library(tidyverse)
library(ggpubr)

source("scripts/COCONUT_tester.R")

# Read datasets with gene expression in stretched cells #############

## GSE59128 ----

x<-read.ilmn("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE59128_RAW/GSE59128_non-normalized.txt", probeid = "PROBE_ID") 
x<-x[,c(1:60)] 
y<-neqc(x) 
expr<-as.matrix(y)

class<-c(rep(0,12), rep(1,28), rep(0,12), rep(1,8)) #¿De dónde se obtienen estos valores??
cyclic<-c(rep(0,24), rep(1,16), rep(0,20))
second_hit<-c(rep(0,12), rep(1, 12), rep(0,8), rep(1,8), rep(0,12), rep(1,8))
names(class)<-colnames(expr)

annotation<-read.csv2("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE59128_RAW/annotation_GPL10558.csv", skip=8, stringsAsFactors = FALSE)

keys<-annotation$Symbol[match(rownames(expr), annotation$Probe_Id)]

names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=cyclic, second_hit=second_hit)
rownames(pheno)<-colnames(expr)

formattedName<-"GSE59128"
pheno$dataset <- formattedName

dataObj1<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj1, "Dataset")
boxplot(dataObj1$expr)


## GSE16650 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE16650_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(0,rep(1,5),0,rep(1,5))
cyclic<-c(0,1,0,1,0,1,0,1,0,1,0,1)
second_hit<-c(0,0,1,1,1,1,0,0,1,1,1,1)
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

annotations<-read.csv("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE34789_RAW/HG-U133_Plus_2.na36.annot.csv", skip=25, stringsAsFactors = FALSE)

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=cyclic, second_hit=second_hit)
rownames(pheno)<-files

formattedName<-"GSE16650"
pheno$dataset <- formattedName

dataObj2<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj2, "Dataset") 
boxplot(dataObj2$expr)


## GSE34789 and GSE27128 (Both Calu-3, same group) ----

files1<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE34789_RAW",
                      listGzipped = TRUE,
                      full.names=TRUE)
files2<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE27128_RAW",
                      listGzipped = TRUE,
                      full.names=TRUE)
files<-c(files1, files2)
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(0,1,1,1,1,0,0,0)
cyclic<-class
second_hit<-rep(0,8)
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=cyclic, second_hit=second_hit)
rownames(pheno)<-files

formattedName<-"GSE27128"
pheno$dataset <- formattedName

dataObj3<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj3, "Dataset")
boxplot(dataObj3$expr)


## GSE1541 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE1541_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(rep(0,4), rep(1,16))
cyclic<-c(rep(0,8), rep(1, 4), rep(0,4), rep(1,4))
second_hit<-c(rep(0,4), rep(1,4), rep(0,4), rep(1,8))
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

annotations<-read.csv("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE1541_RAW/HG-Focus.na36.annot.csv", skip=25, stringsAsFactors = FALSE)

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=cyclic, second_hit=second_hit)
rownames(pheno)<-files

formattedName<-"GSE1541"
pheno$dataset <- formattedName

dataObj4<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj4, "Dataset")
boxplot(dataObj4$expr)


## GSE3541 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE3541_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(0,0,0,1,1,1)
cyclic<-class
second_hit<-rep(0,6)
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

annotations<-read.csv("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE3541_RAW/Rat230_2.na36.annot.csv", skip=21, stringsAsFactors = FALSE)

keys<-toupper(annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)])
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=cyclic, second_hit=second_hit)
rownames(pheno)<-files

formattedName<-"GSE3541"
pheno$dataset <- formattedName

dataObj5<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj5, "Dataset")
boxplot(dataObj5$expr)


# Metaintegrator & COCONUT ################################

## Metaintegrator ----

discovery_datasets<-list(dataObj1, dataObj2, dataObj3, dataObj4, dataObj5)
names(discovery_datasets) = c(dataObj1$formattedName, 
                              dataObj2$formattedName, 
                              dataObj3$formattedName, 
                              dataObj4$formattedName, 
                              dataObj5$formattedName)
MetaObj_genes_cells=list() 
MetaObj_genes_cells$originalData <- discovery_datasets

checkDataObject(MetaObj_genes_cells, "Meta", "Pre-Analysis")


## COCONUT ----

coconut_genes_cells<-coconutMetaIntegrator(MetaObj_genes_cells)
coco_out_genes_cells<-combineCOCOoutput(coconut_genes_cells)

testplots <- COCONUT_tester(coconut_genes_cells, MetaObj_genes_cells)
pdf("plots/SUP13_COCONUT_genes_cells.pdf",
    paper = "a4r")
ggarrange(testplots$raw_controls, testplots$raw_cases, testplots$COCONUT_controls, testplots$COCONUT_cases,
          labels = c("Raw controls", "Raw cases", "Normalized controls", "Normalized cases"),
          hjust = -1,
          nrow = 2,
          ncol = 2)
dev.off()


save(coconut_genes_cells, file = "coconut_genes_cells.RData")
save(coco_out_genes_cells, file = "coco_out_genes_cells.RData")


rm(affyRaw, annotation, annotations, discovery_datasets,eset,expr, MetaObj_genes_cells, pheno, x, y, class, cyclic, files, files1, files2, formattedName, keys, second_hit)
rm(list=ls(pattern = "dataObj"))

