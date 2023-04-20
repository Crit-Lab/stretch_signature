
library(limma)
library(stringr)
library(MetaIntegrator)
library(preprocessCore)
library(oligo)
library(MetaIntegrator)
library(COCONUT)

# Read miRNA datasets ###############################

## GSE36256: rat cell line ----

targets<-readTargets("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/code/targets_GSE36256.txt") 
f<-function(x) as.numeric(x$Flags>-99)

RG<-read.maimages(targets, source="genepix", wt.fun=f) 
RG<-limma::backgroundCorrect(RG, method="minimum") 
MA<-normalizeWithinArrays(RG, method="median") 
MA<-normalizeBetweenArrays(MA, method="Rquantile") 

expr<-exprs.MA(MA)[,c(2,4,6,8,10,12,14,16,18,20,22,24)] 

annotations<-read.csv2("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/code/annotation_GPL7722.csv", skip=18, stringsAsFactors = FALSE)

rownames(expr)<-MA$genes$ID
colnames(expr)<-targets$FileName

t<-MA$genes$Name
keep <- str_detect(t, "rno")

t[grep("-",t)]<-str_sub(str_extract(t[grep("-",t)], "-.*"),2)
t[grep("/",t)]<-str_extract(t[grep("/",t)], ".*(?=/)")
t<-sub("/.*", "", t)

#Further trimming
t<-sub("\\*", "", t)
t<-sub("\\-5p", "", t)
t<-sub("\\-3p", "", t)

keys<-t[keep]
names(keys)<-MA$genes$ID[keep]

expr <- expr[keep,]

pheno<-data.frame(stretch=c(1,1,0,1,1,0,1,1,0,1,1,0))
rownames(pheno)<-targets$FileName

class<-pheno$stretch
names(class)<-targets$FileName

formattedName<-"GSE36256"
pheno$dataset <- formattedName

mirDataObj1 <- list(class=class,
                    expr=expr,
                    keys=keys,
                    pheno=pheno,
                    formattedName=formattedName)

checkDataObject(mirDataObj1, "Dataset")


## GSE75100: human cell line ----

targets<-readTargets("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/code/targets_GSE75100.txt")
RG<-read.maimages(targets, source="genepix", wt.fun=f, green.only = TRUE)
RG<-limma::backgroundCorrect(RG, method="normexp")
MA<-normalizeBetweenArrays(RG, method="quantile")

class<-c(0,1,0,1,0,1)
names(class)<-targets$FileName

expr<-matrix(log2(MA$E), ncol=nrow(targets))
expr<-expr+1

rownames(expr)<-MA$genes$ID
colnames(expr)<-targets$FileName

t<-MA$genes$Name
keep <- str_detect(t, "hsa")

t[grep("-",t)]<-str_sub(str_extract(t[grep("-",t)], "-.*"),2)
t[grep("/",t)]<-str_extract(t[grep("/",t)], ".*(?=/)")
t<-sub("/.*", "", t)
t<-sub(";.*", "", t)

#Further trimming
t<-sub("\\*", "", t)
t<-sub("\\-5p", "", t)
t<-sub("\\-3p", "", t)

keys<-t[keep]
names(keys)<-MA$genes$ID[keep]

expr <- expr[keep,]

pheno<-data.frame(stretch=c(0,1,0,1,0,1))
rownames(pheno)<-targets$FileName

formattedName<-"GSE75100"
pheno$dataset <- formattedName
mirDataObj2<-list(class=class,
                  expr=expr,
                  keys=keys,
                  pheno=pheno,
                  formattedName=formattedName)

checkDataObject(mirDataObj2, "Dataset")


## GSE131645: mouse cell line ----

files<-list.celfiles(path="C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/raw data/cells/non_coding/GSE131645_RAW", listGzipped = TRUE, full.names=TRUE)
files<-files[grep("Cell_",files)]
files<-files[1:6]

affyRaw<-read.celfiles(files)
eset<-rma(affyRaw)                            
class<-c(0,0,1,1,0,0)
names(class)<-files

expr<-exprs(eset)
annotations<-read.csv("C:/Users/cecil/Desktop/Crit lab/proyectos/2019 firma stretch/Stretch signature/code/miRNA-4_0-st-v1.annotations.20160922.csv", skip=4, stringsAsFactors = FALSE)
t<-annotations$Transcript.ID.Array.Design. #Diseño del array
keep <- str_detect(t, "mmu")


t[grep("-",t)]<-str_sub(str_extract(t[grep("-",t)], "-.*"),2)

#Further trimming
t<-sub("\\-5p", "", t)
t<-sub("\\-3p", "", t)


annotations$Transcript.ID.Array.Design.<-t

rownames(expr)<-annotations$Probe.Set.ID
colnames(expr)<-files
expr <- expr[keep,]

keys<-annotations$Transcript.ID.Array.Design.
names(keys)<-annotations$Probe.Set.ID
keys <- keys[keep]

pheno<-data.frame(stretch=c(0,0,1,1,0,0))
rownames(pheno)<-files

formattedName<-"GSE131645"
pheno$dataset <- formattedName

mirDataObj3<-list(class=class,
                  expr=expr,
                  keys=keys,
                  pheno=pheno,
                  formattedName=formattedName)

checkDataObject(mirDataObj3, "Dataset")

# Metaintegrator & COCONUT ###############################


## Metaintegrator ----

datasets <- list(mirDataObj1, mirDataObj2, mirDataObj3)
names(datasets) <- c(mirDataObj1$formattedName, 
                     mirDataObj2$formattedName, 
                     mirDataObj3$formattedName)
mirMetaObj <- list() 
mirMetaObj$originalData <- datasets

checkDataObject(mirMetaObj, "Meta", "Pre-Analysis")

## COCONUT ----

coconut_miRNA <- coconutMetaIntegrator(mirMetaObj) 
coco_out_miRNA <- combineCOCOoutput(coconut_miRNA) 

save(coconut_miRNA, file = "coconut_miRNA.RData")
save(coco_out_miRNA, file = "coco_out_miRNA.RData")

rm(affyRaw, annotations,datasets, eset, expr, MA, mirMetaObj, pheno, RG, targets, class, files, formattedName, keys, keep, t, f)
rm(list=ls(pattern = "mirDataObj"))
