########################## FINAL ANALYSIS V2 ###################################

library(tidyverse)
library(limma)
library(stringr)
library(MetaIntegrator)
library(COCONUT)
library(preprocessCore)
library(oligo)
library(multiMiR)
library(pROC)
library(FSelectorRcpp)
library(pheatmap)
library(viridis)
library(eulerr)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library("org.Hs.eg.db", character.only = TRUE)
library(DOSE)
library(ungeviz)
library(DESeq2)
library(pspearman)

source("scripts/gm_mean_v2.R")


# miRNA signature ##############################################################

## 1. Read miRNA arrays ####
 
# Output from:
# source("read_mir_data.R")

load("data objects/coconut_miRNA.RData")
load("data objects/coco_out_miRNA.RData")

## 2. limma ####

miRNAs <- as.matrix(coco_out_miRNA$genes)
dim(miRNAs)
design <- model.matrix(~factor(coco_out_miRNA$class.cntl0.dis1))
fit <- lmFit(miRNAs, design)
fit <- eBayes(fit)

complete_miRNAs <- topTable(fit, number=500000, adjust="fdr")

tt_miRNAs<-topTable(fit, number=50, adjust="fdr", p.value=0.1) 
idx_miRNAs <- ifelse(tt_miRNAs$logFC > 0, "UP", "DOWN")
names(idx_miRNAs) <- rownames(tt_miRNAs)

idx_miRNAs
# miR-383 miR-146b miR-181b  miR-26b  miR-877 miR-130b 
# "UP"   "DOWN"   "DOWN"   "DOWN"     "UP"     "UP" 

## 3. Plots ####

### B. Heatmap

legend <- data.frame(file=colnames(miRNAs), stretch=coco_out_miRNA$class.cntl0.dis1)
annotation <- data.frame(stretch=as.factor(legend$stretch))
rownames(annotation) <- colnames(miRNAs)

ann_colors <- list(stretch=c("0"="#ECC192", "1"="#507779"))

pheatmap(miRNAs[names(idx_miRNAs),order(legend$stretch)],
         color=viridis(n=25, option="E"),
         cellwidth = 15,
         cellheight = 15,
         # annotation_col=annotation,   Not pretty
         scale="row",
         cluster_cols = FALSE,
         #cutree_cols = 2,  #could be
         clustering_method = "complete",
         labels_col = legend$stretch[order(legend$stretch)],
         filename="plots/Fig2_heatmap_miRNAs.pdf")
dev.off()

### C. Metascore of miRNAs

df <- data.frame(score = apply(data.frame(miRNAs[rownames(miRNAs) %in% names(idx_miRNAs)[idx_miRNAs == "UP"],]),
                               MARGIN = 2,
                               FUN = gm_mean) -
                   apply(data.frame(miRNAs[rownames(miRNAs) %in% names(idx_miRNAs)[idx_miRNAs == "DW"],]),
                         MARGIN = 2,
                         FUN = gm_mean),
                 condition = as.factor(coco_out_miRNA$class.cntl0.dis1))


colors_two<-c("0"="#ECC192",
              "1"="#507779",
              "2"="#D9843A",
              "3"="#32474A")


ggplot(df, aes(x = condition, y = score))+
  geom_violin(aes(fill = condition, alpha=0.5), col=NA)+
  geom_jitter(aes(col = condition), width=0.15)+
  theme_bw(base_size = 24)+
  theme(legend.position = "none", aspect.ratio = 1.618)+
  scale_fill_manual(values=colors_two)+
  scale_color_manual(values=colors_two)+
  labs(y="Meta-score", x=NULL)
ggsave("plots/Fig2_metascore_miRNAs.pdf", dpi=300, useDingbats=FALSE)  
dev.off()

shapiro.test(df$score)
t.test(score~condition, data = df)

### C. AUC


stretch.glm <- glm(condition ~ score, family=binomial, data = df)
probs<-predict(stretch.glm, df, type="response")
roc.obj<-roc(df$condition, probs, ci=TRUE)
roc.obj


### Sup: miRNA expression per dataset

pheno <- coco_out_miRNA$pheno
pheno <- bind_cols(pheno, t(coco_out_miRNA$genes)[, names(idx_miRNAs)])

ann_colors<-c("0"="#ECC192", "1"="#507779")

pheno %>% pivot_longer(cols = 4:9, names_to = "gene", values_to = "value") %>% 
  ggplot( aes(x = dataset, y = value, fill = as.factor(stretch), col = as.factor(stretch)))+
  geom_violin(aes(fill=as.factor(stretch), alpha=0.5), col=NA)+
  geom_point(aes(col=as.factor(stretch)), position = position_jitterdodge())+
  theme_bw(base_size = 10)+
  theme(legend.position = "none", aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = NULL, y = NULL)+
  scale_fill_manual(values=ann_colors)+
  scale_color_manual(values=ann_colors)+
  facet_wrap(~gene, scales = "free", nrow = 2)
ggsave("plots/SUP_miRNAs_per_dataset.pdf", dpi=300, useDingbats=FALSE)  
dev.off()


# miRNA targets ################################################################

targets_miRNAs <- list_multimir("mirna")
target_ids <- targets_miRNAs[grep(paste(names(idx_miRNAs), "($|\\-3p$|\\-5p$)",sep="", collapse="|"), targets_miRNAs$mature_mirna_id),] 

gene_targets_human <- get_multimir(org="hsa",
                                   mirna=target_ids$mature_mirna_id[target_ids$org=="hsa"],
                                   table="validated", 
                                   summary=TRUE)

gene_targets_mouse <- get_multimir(org="mmu",
                                   mirna=target_ids$mature_mirna_id[target_ids$org=="mmu"],
                                   table="validated", 
                                   summary=TRUE)

gene_targets_rat <- get_multimir(org="rno",
                                 mirna=target_ids$mature_mirna_id[target_ids$org=="rno"],
                                 table="validated", 
                                 summary=TRUE)

gene_targets <- rbind(gene_targets_human@data, gene_targets_mouse@data, gene_targets_rat@data)

t <- gene_targets$mature_mirna_id
t[grep("-",t)] <- str_sub(str_extract(t[grep("-",t)], "-.*"),2)
t[grep("/",t)] <- str_extract(t[grep("/",t)], ".*(?=/)")
t<-sub("/.*", "", t)

#Further trimming
t <- sub("\\*", "", t)
t <- sub("\\-5p", "", t)
t <- sub("\\-3p", "", t)

gene_targets$trimmed_names <- t
gene_targets$species <- substr(gene_targets$mature_mirna_id,1,3)
gene_targets$target_symbol <- toupper(gene_targets$target_symbol)

length(unique(gene_targets$target_symbol))
# 7980

rm(gene_targets_human, gene_targets_mouse, gene_targets_rat, target_ids, targets_miRNAs)


# Gene signature in cells ######################################################

## 1. Read cell datasets ####

# Output from:
# source("read_cell_datasets.R")

load("data objects/coconut_genes_cells.RData")
load("data objects/coco_out_genes_cells.RData")

## 2. Intersect genes with miRNA targets ####

length(rownames(coco_out_genes_cells$genes)) # gene universe
#5396

length(intersect(rownames(coco_out_genes_cells$genes), unique(gene_targets$target_symbol)))
# 2710 genes regulated by themiRNA signature can be found in the cell datasets


## 3. Limma for the target genes ####

genes<-coco_out_genes_cells$genes
genes_sel<-genes[rownames(genes) %in% intersect(rownames(coco_out_genes_cells$genes), unique(gene_targets$target_symbol)),]

stretch<-coco_out_genes_cells$pheno$stretch
second_hit<-coco_out_genes_cells$pheno$second_hit
group<-stretch+2*second_hit

stretch<-as.factor(stretch)
second_hit<-as.factor(second_hit)
group<-factor(group, labels=c("Control", "Stretch", "Other", "Two_hits"))
table(group)

design<-model.matrix(~stretch+second_hit)
colnames(design)
fit<-lmFit(genes_sel, design)
fit2<-eBayes(fit)
complete_genes<-topTable(fit2, coef=2, number=5000000)
tt_genes_cells<-topTable(fit2, coef=2, number=5000000, adjust="fdr", p.value=0.01)
dim(tt_genes_cells)
# 451 target genes differentially expressed

tt_genes_second_hit<-topTable(fit2, coef=3, number=5000000, adjust="fdr", p.value=0.01)
dim(tt_genes_second_hit)
# 65 for the second hit

idx_genes_cells<-rownames(tt_genes_cells)

## 4. Plots ####

### A. Venn diagram

input_venn <- c(Stretch=nrow(tt_genes_cells), 
                Second_hit=nrow(tt_genes_second_hit),
                "Stretch&Second_hit"=sum(rownames(tt_genes_cells) %in% rownames(tt_genes_second_hit)))
pdf("plots/Fig3_venn_diagram.pdf")
plot(euler(input_venn),
     fills=c("#507779", "#D9843A"),
     quantities=TRUE)
dev.off()

### B. volcano plot of the gene universe

fit<-lmFit(data.frame(genes), design)
fit2<-eBayes(fit)
gene_universe<-topTable(fit2, coef=2, number=5000000)
gene_universe

ggplot(gene_universe, aes(x = logFC, y = -log10(adj.P.Val)))+
  theme_bw()+
  geom_point(col = "#ECC192")+
  geom_point(data = gene_universe[rownames(gene_universe) %in% rownames(tt_genes_cells),], 
             aes(x = logFC, y = -log10(adj.P.Val)),
             col = "#507779")+
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed")
ggsave("plots/Fig3_volcanoplot_cell_universe.pdf", dpi=300, useDingbats=FALSE)  
dev.off()

### C. Overrepresentation analysis

gene_list <- tt_genes_cells$t
names(gene_list) <- rownames(tt_genes_cells)
gene_list <- gene_list[order(gene_list, decreasing = TRUE)]

res <- enrichGO(gene = names(gene_list),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                minGSSize = 10, 
                maxGSSize = 1000)

int.GOs <- res@result %>% filter(p.adjust<0.001)
int.GOs <- int.GOs$Description
length(int.GOs)

res2 <- pairwise_termsim(res, showCategory = 1000)
pdf("plots/Fig3_treeplot.pdf",
    paper = "a4r")
treeplot(res2, 
         showCategory = int.GOs,
         nCluster = 7,
         color = NULL,
         cex_category = 0.5,
         offset = 30,
         offset_tiplab = 0.5,
         group_color = viridis(8)[-8],
         fontsize = 3.5)+
  theme(aspect.ratio = 1/1.62)
dev.off()

# Gene signature in animals ####################################################

## 1. Read animal datsets ####

# Output from:
# source("read_animal_datasets.R")

load("data objects/coconutRes_validation_animal.RData")
load("data objects/coco_out_validation_animal.RData")

## 2. Score and ROC curve with all DF genes in cells ####

positive_genes_animal<-rownames(tt_genes_cells)[tt_genes_cells$logFC>0]
negative_genes_animal<-rownames(tt_genes_cells)[tt_genes_cells$logFC<0]

genes<-as.matrix(coco_out_validation_animal$genes)
rownames(genes)<-toupper(rownames(genes))

dim(genes)
dim(genes[rownames(genes) %in% rownames(tt_genes_cells),])
# 144 DF genes in cells can be found in the animal metaanalysis

sel.genes <- genes[rownames(genes) %in% rownames(tt_genes_cells),]


## 3. Plots before greedy ####

### A. metascore

df <- data.frame(score = apply(data.frame(sel.genes[rownames(sel.genes) %in% positive_genes_animal, ]),
                               MARGIN = 2,
                               FUN = gm_mean) -
                   apply(data.frame(sel.genes[rownames(sel.genes) %in% negative_genes_animal, ]),
                         MARGIN = 2,
                         FUN = gm_mean),
                 condition = coco_out_validation_animal$pheno$stretch,
                 second_hit = coco_out_validation_animal$pheno$second_hit,
                 group = as.factor(coco_out_validation_animal$pheno$stretch + 2*coco_out_validation_animal$pheno$second_hit),
                 vt = coco_out_validation_animal$pheno$vt)


colors_two<-c("0"="#ECC192",
              "1"="#507779",
              "2"="#D9843A",
              "3"="#32474A")


ggplot(df, aes(x = group, y = score))+
  geom_violin(aes(fill = group, alpha=0.5), col=NA)+
  geom_jitter(aes(col = group), width=0.15)+
  theme_bw(base_size = 24)+
  theme(legend.position = "none", aspect.ratio = 1)+
  scale_fill_manual(values=colors_two)+
  scale_color_manual(values=colors_two)+
  labs(y="Meta-score", x=NULL)
ggsave("plots/Fig4_metascore_animal_allgenes.pdf", dpi=300, useDingbats=FALSE)  
dev.off()

aggregate(score~group, data=df, FUN=summary)
model.aov.animal<-aov(score~as.factor(condition):as.factor(second_hit), data=df)
summary(model.aov.animal)
TukeyHSD(model.aov.animal)

### B. ROC curve

model_animal<-glm(as.factor(condition)~score, family=binomial(link=logit), data=df)
probs<-predict(model_animal, df, type="response")
roc.obj<-roc(df$condition, probs, ci=TRUE)
roc.obj
#AUC=0.878

pdf("plots/Fig4_ROC_animal_allgenes.pdf")
ggroc(roc.obj, col="#507779", size=2, legacy.axes = TRUE)+
  theme_bw(base_size = 24)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "darkgrey", 
               linetype = "dashed")
dev.off()


### C. Quantitative analysis

breaks <- c(-1, 9, 21, 39)
df$vt[df$vt==0] <- 6
df$vt_group <- cut(df$vt, 
                   breaks = breaks,
                   labels=c("Spont", "10-20", ">20"),
                   ordered_result=TRUE)

ggplot(df[!is.na(df$vt_group),], aes(x = vt_group, y = score))+
  theme_bw(base_size = 24)+
  labs(x="Tidal volume (ml/Kg)", y="Meta-Score")+
  geom_jitter(aes(col = vt_group), width = 0.05)+
  scale_color_manual(values=c("#F5DCC1", "#92ABAD", "#5B8283"))+
  stat_boxplot(geom = "errorbar", width = 0.3)+
  stat_summary(geom = "hpline", fun = "median")+
  theme(legend.position = "none", aspect.ratio = 1)
ggsave("plots/Fig4_Metascore_allgenes_vt_group.pdf", useDingbats=FALSE)
dev.off()

cor.sp<-spearman.test(df$score, df$vt_group)
cor.sp

## 4. greedy selection of genes ####

greedy_set <- data.frame(t(sel.genes))
stretch <-  coco_out_validation_animal$pheno$stretch

evaluator_AUC <- function(attributes, data, condition = stretch) {
  
  df <- data.frame(score = apply(data.frame(data[,colnames(data) %in% attributes[attributes %in% positive_genes_animal]]),
                                 MARGIN = 1,
                                 FUN = gm_mean) -
                     apply(data.frame(data[,colnames(data) %in% attributes[attributes %in% negative_genes_animal]]),
                           MARGIN = 1,
                           FUN = gm_mean),
                   condition = condition)
  
  stretch.glm <- glm(condition ~ score, family=binomial, data = df)
  probs <- predict(stretch.glm, df, type="response")
  roc.obj <- roc(df$condition, probs, ci=FALSE)
  return(as.numeric(roc.obj$auc))
}

genenames <- colnames(greedy_set)

t<-feature_search(attributes= genenames,
                  fun = evaluator_AUC,
                  data = greedy_set,
                  mode = "greedy",
                  type = "forward")

t2<-as.matrix(t(t$best))
t2<-t2[t2[,1]>0,]  
greedy_genes<-names(t2) 
greedy_genes


## 5. plots after greedy ####

### D. Heatmap of greedy genes

group <- coco_out_validation_animal$pheno$stretch + 2*coco_out_validation_animal$pheno$second_hit

mat <- greedy_set[,colnames(greedy_set) %in% greedy_genes]
mat <- mat[order(group),]

pheatmap(t(mat),
         color=viridis(n=25, option="E"),
         cellwidth = 1.5,
         cellheight = 15,
         scale="row",
         cluster_cols = FALSE,
         clustering_method = "complete",
         labels_col = group[order(group)],
         #labels_row = rownames(genes_long),
         fontsize_col=2,
         filename="plots/Fig4_heatmap_greedy_genes.pdf")


### E. Metascore

df <- data.frame(score = apply(data.frame(greedy_set[,colnames(greedy_set) %in% greedy_genes[greedy_genes %in% positive_genes_animal]]),
                               MARGIN = 1,
                               FUN = gm_mean) -
                   apply(data.frame(greedy_set[,colnames(greedy_set) %in% greedy_genes[greedy_genes %in% negative_genes_animal]]),
                         MARGIN = 1,
                         FUN = gm_mean),
                 condition = coco_out_validation_animal$pheno$stretch,
                 second_hit = coco_out_validation_animal$pheno$second_hit,
                 group = as.factor(coco_out_validation_animal$pheno$stretch + 2*coco_out_validation_animal$pheno$second_hit),
                 vt = coco_out_validation_animal$pheno$vt)


colors_two<-c("0"="#ECC192",
              "1"="#507779",
              "2"="#D9843A",
              "3"="#32474A")


ggplot(df, aes(x = group, y = score))+
  geom_violin(aes(fill = group, alpha=0.5), col=NA)+
  geom_jitter(aes(col = group), width=0.15)+
  theme_bw(base_size = 24)+
  theme(legend.position = "none", aspect.ratio = 1)+
  scale_fill_manual(values=colors_two)+
  scale_color_manual(values=colors_two)+
  labs(y="Meta-score", x=NULL)
ggsave("plots/Fig4_metascore_greedy_genes.pdf", dpi=300, useDingbats=FALSE)  

aggregate(score~group, data=df, FUN=summary)
model.aov.animal<-aov(score~as.factor(condition):as.factor(second_hit), data=df)
summary(model.aov.animal)
TukeyHSD(model.aov.animal)

### F. ROC curve

stretch.glm <- glm(condition ~ score, family=binomial, data = df)
probs<-predict(stretch.glm, df, type="response")
roc.obj<-roc(df$condition, probs, ci=TRUE)
roc.obj

pdf("plots/Fig4_ROC_greedy_genes.pdf")
ggroc(roc.obj, col="#507779", size=2, legacy.axes = TRUE)+
  theme_bw(base_size = 24)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "darkgrey", 
               linetype = "dashed")
dev.off()


### G. Quantitative analysis

breaks <- c(-1, 9, 21, 39)
df$vt[df$vt == 0] <- 6
df$vt_group <- cut(df$vt, 
                   breaks = breaks,
                   labels=c("Spont", "10-20", ">20"),
                   ordered_result=TRUE)

ggplot(df[!is.na(df$vt_group),], aes(x = vt_group, y = score))+
  theme_bw(base_size = 24)+
  labs(x="Tidal volume (ml/Kg)", y="Meta-Score")+
  geom_jitter(aes(col = vt_group), width = 0.05)+
  scale_color_manual(values=c("#F5DCC1", "#92ABAD", "#5B8283"))+
  stat_boxplot(geom = "errorbar", width = 0.3)+
  stat_summary(geom = "hpline", fun = "median")+
  theme(legend.position = "none", aspect.ratio = 1)
ggsave("plots/Fig4_Metascore_greedygenes_vt_group.pdf", useDingbats=FALSE)
dev.off()

cor.sp<-spearman.test(df$score, df$vt_group)
cor.sp

### Sup: gene expression per dataset

pheno <- coco_out_validation_animal$pheno
pheno <- bind_cols(pheno, t(coco_out_validation_animal$genes)[, toupper(rownames(coco_out_validation_animal$genes)) %in% greedy_genes])

ann_colors<-c("0"="#ECC192", "1"="#507779")

pheno %>% pivot_longer(cols = 7:12, names_to = "gene", values_to = "value") %>% 
  ggplot( aes(x = dataset, y = value, fill = as.factor(stretch), col = as.factor(stretch)))+
  geom_violin(aes(fill=as.factor(stretch), alpha=0.5), col=NA)+
  geom_point(aes(col=as.factor(stretch)), position = position_jitterdodge())+
  theme_bw(base_size = 10)+
  theme(legend.position = "none", aspect.ratio = 1/1.618,
        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = NULL, y = NULL)+
  scale_fill_manual(values=ann_colors)+
  scale_color_manual(values=ann_colors)+
  facet_wrap(~gene, scales = "free", nrow = 2)
ggsave("plots/SUP_animal_genes_per_dataset.pdf", useDingbats=FALSE)
dev.off()

### Same thing for cells

pheno <- coco_out_genes_cells$pheno
pheno <- bind_cols(pheno, t(coco_out_genes_cells$genes)[, toupper(rownames(coco_out_genes_cells$genes)) %in% greedy_genes])

ann_colors<-c("0"="#ECC192", "1"="#507779")

pheno %>% pivot_longer(cols = 6:11, names_to = "gene", values_to = "value") %>% 
  ggplot( aes(x = dataset, y = value, fill = as.factor(stretch), col = as.factor(stretch)))+
  geom_violin(aes(fill=as.factor(stretch), alpha=0.5), col=NA)+
  geom_point(aes(col=as.factor(stretch)), position = position_jitterdodge())+
  theme_bw(base_size = 10)+
  theme(legend.position = "none", aspect.ratio = 1/1.618,
        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = NULL, y = NULL)+
  scale_fill_manual(values=ann_colors)+
  scale_color_manual(values=ann_colors)+
  facet_wrap(~gene, scales = "free", nrow = 2)
ggsave("plots/SUP_cells_genes_per_dataset.pdf", useDingbats=FALSE)
dev.off()

# External validations ##################################################

## GSE114132 ----

load("data objects/GSE114132_rawcounts.RData")

# filtering genes with low counts
dim(data)
data <- data[rowSums(data) >= ncol(data),]
dim(data)

samples <- data.frame(samples = colnames(data), condition = as.factor(rep(c("control", "VILI"), each = 5)))
rownames(samples) <- samples$samples

# Normalizing counts with Deseq2
dds <- DESeqDataSetFromMatrix(data, samples, design = ~condition)
dds <- DESeq(dds)

int.genes <- counts(dds, normalized = TRUE)[toupper(rownames(counts(dds))) %in% greedy_genes,]

samples$score <- apply(int.genes, MARGIN = 2, FUN = gm_mean)


ggplot(samples, aes(x = condition, y = score))+
  geom_violin(aes(fill = condition, alpha=0.5), col=NA)+
  geom_jitter(aes(col = condition), width=0.15, size = 3)+
  theme_bw(base_size = 24)+
  theme(legend.position = "none", aspect.ratio = 1.618)+
  scale_fill_manual(values = c("#ECC192","#507779" ))+
  scale_color_manual(values = c("#ECC192","#507779" ))+
  labs(y="RNA score", x=NULL)

ggsave("plots/Fig5_GSE114132_score.pdf", useDingbats=FALSE)
dev.off()

shapiro.test(samples$score)
t.test(score~condition, samples)



stretch.glm <- glm(condition ~ score, family=binomial, data = samples)
probs<-predict(stretch.glm, samples, type="response")
roc.obj<-roc(samples$condition, probs, ci=FALSE)
as.numeric(roc.obj$auc)

pdf("plots/Fig5_GSE114132_ROC.pdf")
ggroc(roc.obj, col="#507779", size=2, legacy.axes = TRUE)+
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "darkgrey", 
               linetype = "dashed")

dev.off()


## Ex vivo lungs RNA and miRNA ##################################################

### miRNA ####

#### 1. Read mirDeep data #####

load("data objects/exvivo_lungs_miRNA_counts.RData")

numbers <- str_replace(names(idx_miRNAs), "miR", "")

int.genes <- data[str_detect(data$precursor, paste(numbers, collapse = "|")),]
rownames(int.genes) <- int.genes$precursor
int.genes <- int.genes[, 11:16]

#managing the two locus of mir-181b by adding them up
int.genes["hsa-mir-181b-1",] <- int.genes["hsa-mir-181b-1",]+int.genes["hsa-mir-181b-2",]
int.genes <- int.genes[-4,]

load("data objects/exvivo_lung_sample_data.RData")

samples$score <- apply(int.genes[rownames(int.genes) %in% paste0("hsa-", tolower(names(idx_miRNAs)[idx_miRNAs == "UP"])),], 
                       MARGIN = 2, 
                       FUN = gm_mean)-
  apply(int.genes[rownames(int.genes) %in% paste0("hsa-", tolower(names(idx_miRNAs)[idx_miRNAs == "DOWN"])),], 
        MARGIN = 2, 
        FUN = gm_mean)

#### 2. Plots #####

## A. score

ggplot(samples, aes(x = condition, y = score))+
  geom_violin(aes(fill = condition, alpha=0.5), col=NA)+
  geom_jitter(aes(col = condition), width=0.15, size = 3)+
  theme_bw(base_size = 24)+
  theme(legend.position = "none", aspect.ratio = 1.618)+
  scale_fill_manual(values = c("#ECC192","#507779" ))+
  scale_color_manual(values = c("#ECC192","#507779" ))+
  labs(y="miRNA score", x=NULL)
ggsave("plots/Fig5_exvivo_miRNA_score.pdf", useDingbats=FALSE)
dev.off()


## B. ROC curve

samples$comp <- c(0,0,1,1,0,1)

stretch.glm <- glm(comp ~ score, family=binomial, data = samples)
probs<-predict(stretch.glm, samples, type="response")
roc.obj<-roc(samples$comp, probs, ci=TRUE)
roc.obj

pdf("plots/Fig5_exvivo_miRNA_ROC.pdf")
ggroc(roc.obj, col="#507779", size=2, legacy.axes = TRUE)+
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "darkgrey", 
               linetype = "dashed")
dev.off()

### RNA ####

load("data objects/exvivo_lungs_mRNA_normcounts.RData")

load("data objects/exvivo_lung_sample_data.RData")
rownames(samples) <- colnames(norm.counts)

int.genes <- as.data.frame(t(norm.counts[rownames(norm.counts) %in% greedy_genes, ]))
samples <- cbind(samples, int.genes)
samples$score <- apply(samples[,colnames(samples) %in% greedy_genes], MARGIN = 1, FUN = gm_mean)


#### 2. Plots #####

## A. score

ggplot(samples, aes(x = condition, y = score))+
  geom_violin(aes(fill = condition, alpha=0.5), col=NA)+
  geom_jitter(aes(col = condition), width=0.15, size = 3)+
  theme_bw(base_size = 24)+
  theme(legend.position = "none", aspect.ratio = 1.618)+
  scale_fill_manual(values = c("#ECC192","#507779" ))+
  scale_color_manual(values = c("#ECC192","#507779" ))+
  labs(y="RNA score", x=NULL)
ggsave("plots/Fig5_exvivo_RNA_score.pdf", useDingbats=FALSE)
dev.off()


## B. ROC curve

samples$comp <- c(0,0,1,1,0,1)

stretch.glm <- glm(comp ~ score, family=binomial, data = samples)
probs<-predict(stretch.glm, samples, type="response")
roc.obj<-roc(samples$comp, probs, ci=TRUE)
roc.obj

pdf("plots/Fig5_exvivo_RNA_ROC.pdf")
ggroc(roc.obj, col="#507779", size=2, legacy.axes = TRUE)+
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "darkgrey", 
               linetype = "dashed")
dev.off()


## miRNA signature in BALF ######################################################

load("data objects/BALF_miRNAs.RData")

mirnas_abundance <- mirnas_balf
mirnas_abundance[,3:14] <- mirnas_abundance[,3:14]/mirnas_abundance$`Total Counts`

mirnas_abundance <- mirnas_abundance %>% 
  mutate(`miR-146b` = `miR-146b-3p` + `miR-146b-5p`,
         `miR-181b` = `miR-181b-3p` + `miR-181b-5p`,
         `miR-26b` = `miR-26b-3p` + `miR-26b-5p`,
         `miR-383` = `miR-383-3p` + `miR-383-5p`,
         `miR-877` = `miR-877-3p` + `miR-877-5p`,
         `miR-130b` = `miR-130b-3p` + `miR-130b-5p`) %>% 
  dplyr::select(- c(`miR-146b-3p`, `miR-146b-5p`, `miR-181b-3p`, `miR-181b-5p`,
             `miR-26b-3p`, `miR-26b-5p`, `miR-383-3p`, `miR-383-5p`,
             `miR-877-3p`, `miR-877-5p`, `miR-877-3p`, `miR-877-5p`,
             `miR-130b-3p`, `miR-130b-5p`))

mirnas_abundance$score <- apply(mirnas_abundance[,str_which(colnames(mirnas_abundance), pattern = paste(names(idx_miRNAs)[idx_miRNAs == "UP"], collapse = "|"))], 
                           MARGIN = 1, 
                           FUN = gm_mean)-
  apply(mirnas_abundance[,str_which(colnames(mirnas_abundance), pattern = paste(names(idx_miRNAs)[idx_miRNAs == "DOWN"], collapse = "|"))], 
        MARGIN = 1, 
        FUN = gm_mean)

colors_two<-c("3"="#ECC192",
              "6"="#507779")

ggplot(mirnas_abundance[!is.na(mirnas_abundance$strain),], aes(x=as.factor(vt), y=score))+
  geom_line(aes(group=Patient), linetype=2, col="gray")+
  geom_point(aes(col=as.factor(vt)), size=3)+
  theme_bw(base_size = 24)+
  theme(legend.position = "none", aspect.ratio = 1.618, strip.background = element_blank())+
  scale_fill_manual(values=colors_two)+
  scale_color_manual(values=colors_two)+
  scale_x_discrete(labels=c("3 ml/Kg", "6 ml/Kg"))+
  labs(y="miRNA score", x=NULL)+
  facet_grid(~(delta_strain<0))
ggsave("plots/Fig5_BALF_metascores.pdf", useDingbats=FALSE)
dev.off()

impares <- seq(1,24,2)
delta_score <- mirnas_abundance$score[impares] - mirnas_abundance$score[impares+1]
roc.curve <- roc(as.factor(mirnas_abundance$delta_strain[impares]<0)~delta_score, ci=TRUE)
as.numeric(roc.curve$auc)
roc.curve


shapiro.test(mirnas_abundance$score[!is.na(mirnas_abundance$strain) & mirnas_abundance$delta_strain>=0])
shapiro.test(mirnas_abundance$score[!is.na(mirnas_abundance$strain) & mirnas_abundance$delta_strain<0])

with(mirnas_abundance[!is.na(mirnas_abundance$strain) & mirnas_abundance$delta_strain>=0,],
     t.test(score[impares], score[impares+1], paired=TRUE))
with(mirnas_abundance[!is.na(mirnas_abundance$strain) & mirnas_abundance$delta_strain<0,],
     t.test(score[impares], score[impares+1], paired=TRUE))


pdf("Fig5_BALF_ROC.pdf")
ggroc(roc.curve, col="#507779", size=2, legacy.axes = TRUE)+
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "darkgrey", 
               linetype = "dashed")
dev.off()


## miRNA signature in COVID19 patients ########################################

# GSE197258

covid_data <- read.table("data objects/GSE197258_raw_norm_counts.csv", sep = ",")
colnames(covid_data) <- covid_data[1,]
covid_data <- covid_data[-1,]
covid_norm <- covid_data[,52:97]
covid_norm <- data.frame(apply(covid_norm, MARGIN = 2, FUN = as.numeric))
covid_norm$precursor <- covid_data$precursor

load("data objects/covid_miRNA_sample_data.RData")
samples$clasif <- ifelse((samples$PEEP_2 - samples$PEEP_1)>= 0, "2", "1")
samples$delta_pCO2 <- samples$pCO2_2-samples$pCO2_1
samples$pCO2_by_PEEP <- ifelse(samples$clasif == 1, samples$delta_pCO2*-1, samples$delta_pCO2)


numbers <- str_replace(names(idx_miRNAs), "miR", "")

int.data <- covid_norm[str_detect(covid_norm$precursor, paste(numbers, collapse ="|")),]
rownames(int.data) <- int.data$precursor
int.data <- int.data[,-47]

# Manage the two locus of 181b
int.data["hsa-mir-181b-1",] <- int.data["hsa-mir-181b-1",] + int.data["hsa-mir-181b-2",]
rownames(int.data)[rownames(int.data) == "hsa-mir-181b-1"] <- "hsa-mir-181b"
int.data <- int.data[-4,]

samples$score <- apply(int.data[rownames(int.data) %in% c("hsa-mir-383","hsa-mir-877", "hsa-mir-130b"),],
               MARGIN = 2,
               FUN = gm_mean) -
  apply(int.data[rownames(int.data) %in% c("hsa-mir-146b","hsa-mir-181b", "hsa-mir-26b"),],
        MARGIN = 2,
        FUN = gm_mean) 


#### 2. Plots #####

samples <-  samples[!is.na(samples$pCO2_by_PEEP),]

## A. score

ggplot(samples, aes(x = as.factor(pCO2_by_PEEP > 0), y = score))+
  geom_violin(aes(fill = as.factor(pCO2_by_PEEP > 0), alpha=0.5), col=NA)+
  geom_jitter(aes(col = as.factor(pCO2_by_PEEP > 0)), width=0.15, size = 3)+
  theme_bw(base_size = 24)+
  theme(legend.position = "none", aspect.ratio = 1.618)+
  scale_fill_manual(values = c("#ECC192","#507779" ))+
  scale_color_manual(values = c("#ECC192","#507779" ))+
  labs(y="Meta-score", x=NULL)


ggsave("plots/Fig5_covid19_score.pdf", useDingbats=FALSE)
dev.off()

t.test(score~as.factor(pCO2_by_PEEP > 0), samples)


## B. ROC curve

samples$comp <- ifelse(samples$pCO2_by_PEEP > 0, 1, 0)

stretch.glm <- glm(comp ~ score, family=binomial, data = samples)
probs<-predict(stretch.glm, samples, type="response")
roc.obj<-roc(samples$comp, probs, ci=TRUE)
roc.obj

pdf("plots/Fig5_covid19_ROC.pdf")
ggroc(roc.obj, col="#507779", size=2, legacy.axes = TRUE)+
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color = "darkgrey", 
               linetype = "dashed")
dev.off()


