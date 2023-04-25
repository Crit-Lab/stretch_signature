################ AUC distribution of random signatures #########################

random_AUCs<-vector(mode="numeric", length=10000)

for(i in 1:10000) {
  all_animal_genes <- c(positive_genes_animal, negative_genes_animal)
  random.names <- sample(all_animal_genes, 6)
  random.genes <- ifelse(random.names %in% positive_genes_animal, "UP", "DOWN")
  names(random.genes) <- random.names
  
  df <- data.frame(score = apply(data.frame(sel.genes[rownames(sel.genes) %in% names(random.genes)[random.genes == "UP"], ]),
                                 MARGIN = 2,
                                 FUN = gm_mean) -
                     apply(data.frame(sel.genes[rownames(sel.genes) %in% names(random.genes)[random.genes == "DOWN"], ]),
                           MARGIN = 2,
                           FUN = gm_mean),
                   condition = coco_out_validation_animal$pheno$stretch,
                   second_hit = coco_out_validation_animal$pheno$second_hit,
                   group = as.factor(coco_out_validation_animal$pheno$stretch + 2*coco_out_validation_animal$pheno$second_hit),
                   vt = coco_out_validation_animal$pheno$vt)
  
  
  
  model_animal<-glm(as.factor(condition)~score, family=binomial(link=logit), data=df)
  probs<-predict(model_animal, df, type="response")
  roc.obj<-roc(df$condition, probs, ci=TRUE)
  random_AUCs[i]<-as.numeric(roc.obj$auc)
}


ggplot(NULL, aes(x=random_AUCs))+
  geom_density(fill="darkblue", alpha=0.5)+
  xlim(0.5,1)+
  theme_bw(base_size = 24)+
  theme(aspect.ratio=1/1.618)+
  labs(x="AUC w/ 6 random genes", y="Density")
ggsave("plots/SUP_random_greedy_signature.pdf", useDingbats=FALSE)
sum(random_AUCs>0.97)/length(random_AUCs)
# 0.0011
sum(random_AUCs>0.95)/length(random_AUCs)
# 0.0055
