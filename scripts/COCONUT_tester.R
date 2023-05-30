COCONUT_tester<-function(coco_object, original_object) {
  
  ggDensities<-function(x) {
    ggplot(x, aes(x=value, col=key))+
      geom_density()+
      theme_bw(base_size=24)+
      labs(x="Gene expression", y="Density")+
      xlim(0,16)+
      theme(legend.position = "none", aspect.ratio=1/1.618)
  }
  
  COCONUT_controls<-Reduce(rbind, 
                           lapply(
                             lapply(coco_object$controlList$GSEs, 
                                    function(x) as.data.frame(x$genes)), gather))
  COCONUT_cases<-Reduce(rbind,
                        lapply(
                          lapply(coco_object$COCONUTList, 
                                 function(x) as.data.frame(x$genes)), gather))
  raw_controls<-Reduce(rbind, 
                       lapply(
                         lapply(original_object$originalData, 
                                function(x) as.data.frame(x$expr[x$keys %in% rownames(coco_object$COCONUTList[[1]]$genes),x$class==0])), gather))
  raw_cases<-Reduce(rbind, 
                    lapply(
                      lapply(original_object$originalData, 
                             function(x) as.data.frame(x$expr[x$keys %in% rownames(coco_object$COCONUTList[[1]]$genes),x$class==1])), gather))
  
  plotlist<-list(raw_controls=raw_controls, 
                 raw_cases=raw_cases, 
                 COCONUT_controls=COCONUT_controls, 
                 COCONUT_cases=COCONUT_cases)
  plots<-lapply(plotlist, ggDensities)
  return(plots)
}

