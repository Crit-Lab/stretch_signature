
library(oligo)
library(mogene20sttranscriptcluster.db)
library(mogene10sttranscriptcluster.db)
library(mouse430a2.db)
library(GEOquery)
library(m20kcod.db)
library(mgu74av2.db)
library(rgu34a.db)
library(MetaIntegrator)
library(COCONUT)

# Load datasets from animal experiments ####

## GSE121550 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE121550_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(1,1,1,0,0,0,0,0,0,1,1,1)
stretch<-class
second_hit<-rep(0,12)
vt<-rep(NA, 12)
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

annotations<-data.frame(SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=","), 
                        DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=","),
                        ENTREZID=sapply(contents(mogene20sttranscriptclusterENTREZID), paste, collapse=","),
                        ENSEMBLID=sapply(contents(mogene20sttranscriptclusterENSEMBL), paste, collapse=","))

keys<-as.character(annotations$SYMBOL[match(rownames(expr), rownames(annotations))])
names(keys)<-rownames(expr)


pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE121550"
pheno$dataset <- formattedName

dataObj1<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj1, "Dataset")
boxplot(dataObj1$expr)


## GSE85269 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE85269_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
files<-files[c(1:3, 7:9)]
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(0,0,0,1,1,1)
stretch<-class
second_hit<-rep(0,6)
vt<-rep(NA,6)
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

keys<-as.character(annotations$SYMBOL[match(rownames(expr), rownames(annotations))])
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE85269"
pheno$dataset <- formattedName

dataObj2<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj2, "Dataset")
boxplot(dataObj2$expr)


## GSE18341 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE18341_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(0,0,0,0, 1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1) #Only pure controls
stretch<-c(rep(0,4), rep(1, 3), rep(0,4), rep(1,3), rep(0,4), rep(1,4), rep(0,4), rep(1, 4))
second_hit<-c(rep(0,7), rep(1,7), rep(0,8), rep(1, 8))
vt<-c(rep(0,4), rep(15, 3), rep(0,4), rep(15,3), rep(0,4), rep(15,4), rep(0,4), rep(15, 4))
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

annotations<-read.csv("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE18341_RAW/Mouse430_2.na36.annot.csv", skip=22, stringsAsFactors = FALSE)

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE18341"
pheno$dataset <- formattedName

dataObj3<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj3, "Dataset")
boxplot(dataObj3$expr)


## GSE11434 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE11434_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(1,1, 1,1,1, 0,0,0,0,0)
stretch<-class
second_hit<-rep(0,10)
vt<-rep(NA, 10)
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE11434"
pheno$dataset <- formattedName

dataObj4<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj4, "Dataset")
boxplot(dataObj4$expr)


## GSE9368 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE9368_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
files<-files[1:6]
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(0,0,0,1,1,1)
stretch<-class
second_hit<-rep(0,6)
vt<-c(0,0,0,30,30,30)
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE9368"
pheno$dataset <- formattedName

dataObj5<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj5, "Dataset")
boxplot(dataObj5$expr)


## GSE9314 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE9314_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
files<-files[1:8]
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(0,0,0,0,1,1,1,1)
stretch<-class
second_hit<-rep(0,8)
vt<-c(0,0,0,0,30,30,30,30)
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE9314"
pheno$dataset <- formattedName

dataObj6<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj6, "Dataset")
boxplot(dataObj6$expr)


## GSE86229 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE86229_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
files<-files[c(1:5,16:25)]
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(rep(0,5), rep(1, 10))
stretch<-class
second_hit<-rep(0,15)
vt<-c(rep(0,5), rep(35, 10))
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

annotations<-data.frame(SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=","), 
                        DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=","),
                        ENTREZID=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=","),
                        ENSEMBLID=sapply(contents(mogene10sttranscriptclusterENSEMBL), paste, collapse=","))

keys<-as.character(annotations$SYMBOL[match(rownames(expr), rownames(annotations))])
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE86229"
pheno$dataset <- formattedName

dataObj7<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj7, "Dataset")
boxplot(dataObj7$expr)

## GSE31678 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE31678_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
files<-files[c(1,2,6:8)]

affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(rep(0,2), rep(1, 3))
stretch<-class
second_hit<-rep(0, 5)
vt<-c(rep(0,2), rep(18,3))
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

annotations<-read.csv("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/coding/GSE3541_RAW/Rat230_2.na36.annot.csv", skip=21, stringsAsFactors = FALSE)

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE31678"
pheno$dataset <- formattedName

dataObj8<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj8, "Dataset")
boxplot(dataObj8$expr)


## GSE7041 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE7041_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
files<-files[c(4:6,10:12)]

affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(rep(0,3), rep(1, 3))
stretch<-class
second_hit<-rep(0,6)
vt<-c(rep(0,3), rep(20, 3))
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

keys<-annotations$Gene.Symbol[match(rownames(expr), annotations$Probe.Set.ID)]
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE7041"
pheno$dataset <- formattedName

dataObj9<-list(class=class,
               expr=expr,
               keys=keys,
               pheno=pheno,
               formattedName=formattedName)

checkDataObject(dataObj8, "Dataset")
boxplot(dataObj9$expr)


## GSE9208 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE9208_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)
files<-files[c(1:5,7)]
affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)

class<-c(rep(0,3), rep(1, 3))
stretch<-class
second_hit<-rep(0, 6)
vt<-c(rep(0,3), rep(30, 3))
names(class)<-files

expr<-exprs(eset)
colnames(expr)<-files

#BiocManager::install("mouse430a2.db")

annotations<-data.frame(SYMBOL=sapply(contents(mouse430a2SYMBOL), paste, collapse=","), 
                        DESC=sapply(contents(mouse430a2GENENAME), paste, collapse=","),
                        ENTREZID=sapply(contents(mouse430a2ENTREZID), paste, collapse=","),
                        ENSEMBLID=sapply(contents(mouse430a2ENSEMBL), paste, collapse=","))

keys<-as.character(annotations$SYMBOL[match(rownames(expr), rownames(annotations))])
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files

formattedName<-"GSE9208"
pheno$dataset <- formattedName

dataObj10<-list(class=class,
                expr=expr,
                keys=keys,
                pheno=pheno,
                formattedName=formattedName)

checkDataObject(dataObj10, "Dataset")
boxplot(dataObj10$expr)


## GSE7742 ----

codset<-getGEO("GSE7742", destdir = "C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE7742_RAW")
eset<-codset[[1]]
eset<-eset[,1:7]

class<-c(rep(0,3), rep(1, 4))
stretch<-class
second_hit<-rep(0, 7)
vt<-c(rep(0,3), rep(10, 4))
names(class)<-eset@phenoData@data$geo_accession

expr<-exprs(eset)
rownames(expr)<-eset@featureData@data$ENTREZ_ID
colnames(expr)<-eset@phenoData@data$geo_accession

expr<-expr[!is.na(rownames(expr)),]

annotations<-data.frame(SYMBOL=sapply(contents(m20kcodSYMBOL), paste, collapse=","), 
                        DESC=sapply(contents(m20kcodGENENAME), paste, collapse=","),
                        ENTREZID=sapply(contents(m20kcodENTREZID), paste, collapse=","),
                        ENSEMBLID=sapply(contents(m20kcodENSEMBL), paste, collapse=","))

keys<-as.character(annotations$SYMBOL[match(rownames(expr), annotations$ENTREZID)])
names(keys)<-rownames(expr)


pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-eset@phenoData@data$geo_accession

formattedName<-"GSE7742"
pheno$dataset <- formattedName

dataObj11<-list(class=class,
                expr=expr,
                keys=keys,
                pheno=pheno,
                formattedName=formattedName)

checkDataObject(dataObj11, "Dataset")
boxplot(dataObj11$expr)


## GSE2411 ----

affyset<-getGEO("GSE2411", destdir = "C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE2411_RAW")
eset<-affyset[[1]]

class<-c(rep(0,6), rep(1, 18))
stretch<-c(rep(0,6), rep(1, 6), rep(0,6), rep(1, 6))
second_hit<-c(rep(0,12), rep(1,12))
vt<-c(rep(0,6), rep(10, 6), rep(0,6), rep(10, 6))
names(class)<-rownames(eset@phenoData@data)

expr<-log2(exprs(eset))
colnames(expr)<-rownames(eset@phenoData@data)

keys<-eset@featureData@data$`Gene Symbol`
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-rownames(eset@phenoData@data)

formattedName<-"GSE2411"
pheno$dataset <- formattedName

dataObj12<-list(class=class,
                expr=expr,
                keys=keys,
                pheno=pheno,
                formattedName=formattedName)

checkDataObject(dataObj12, "Dataset")
boxplot(dataObj12$expr)


## GSE2368 ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/animals/coding/GSE2368_RAW",
                     listGzipped = TRUE,
                     full.names=TRUE)

#Mouse samples
affyRaw<-read.celfiles(files[1:4])
eset<-rma(affyRaw)

class<-c(0,0,1,1)
stretch<-class
second_hit<-rep(0,4)
vt<-c(0,0,35,35)
names(class)<-files[1:4]

expr<-exprs(eset)
colnames(expr)<-files[1:4]

annotations<-data.frame(SYMBOL=sapply(contents(mgu74av2SYMBOL), paste, collapse=","), 
                        DESC=sapply(contents(mgu74av2GENENAME), paste, collapse=","),
                        ENTREZID=sapply(contents(mgu74av2ENTREZID), paste, collapse=","),
                        ENSEMBLID=sapply(contents(mgu74av2ENSEMBL), paste, collapse=","))

keys<-as.character(annotations$SYMBOL[match(rownames(expr), rownames(annotations))])
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files[1:4]

formattedName<-"GSE2368_mouse"
pheno$dataset <- formattedName

dataObj13<-list(class=class,
                expr=expr,
                keys=keys,
                pheno=pheno,
                formattedName=formattedName)

checkDataObject(dataObj13, "Dataset")
boxplot(dataObj13$expr)

#rat

affyRaw<-read.celfiles(files[5:8])
eset<-rma(affyRaw)

class<-c(0,0,1,1)
stretch<-class
second_hit<-rep(0,4)
vt<-c(0,0,12,12)
names(class)<-files[5:8]

expr<-exprs(eset)
colnames(expr)<-files[5:8]

annotations<-data.frame(SYMBOL=sapply(contents(rgu34aSYMBOL), paste, collapse=","), 
                        DESC=sapply(contents(rgu34aGENENAME), paste, collapse=","),
                        ENTREZID=sapply(contents(rgu34aENTREZID), paste, collapse=","),
                        ENSEMBLID=sapply(contents(rgu34aENSEMBL), paste, collapse=","))

keys<-as.character(annotations$SYMBOL[match(rownames(expr), rownames(annotations))])
names(keys)<-rownames(expr)

pheno<-data.frame(disease=class, stretch=stretch, second_hit=second_hit, vt=vt)
rownames(pheno)<-files[5:8]

formattedName<-"GSE2368_rat"
pheno$dataset <- formattedName

dataObj14<-list(class=class,
                expr=expr,
                keys=keys,
                pheno=pheno,
                formattedName=formattedName)

checkDataObject(dataObj14, "Dataset")
boxplot(dataObj14$expr)


# Metaintegrator & COCONUT ################################

## Metaintegrator ----

validation_datasets_animal<-mget(ls(pattern="dataObj"))#Get multiple documents based on an index, type (optional) and ids

names(validation_datasets_animal) = sapply(validation_datasets_animal, FUN=function(x) x$formattedName)
MetaObj_validation_animal=list() 
MetaObj_validation_animal$originalData <- validation_datasets_animal

checkDataObject(MetaObj_validation_animal, "Meta", "Pre-Analysis")


## COCONUT ----

coconutRes_validation_animal <- coconutMetaIntegrator(MetaObj_validation_animal)
coco_out_validation_animal <- combineCOCOoutput(coconutRes_validation_animal)


save(coconutRes_validation_animal, file = "coconutRes_validation_animal.RData")
save(coco_out_validation_animal, file = "coco_out_validation_animal.RData")


rm(affyRaw, affyset, annotations, codset, eset, expr, MetaObj_validation_animal, pheno, validation_datasets_animal, class, files, formattedName, keys,
   second_hit, stretch, vt)
rm(list=ls(pattern = "dataObj"))


