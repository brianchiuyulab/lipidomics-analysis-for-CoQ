# Exploratory PCA views; PC1-vs-PC5 was chosen after inspecting PCs1-10.
# These views do not replace the ordinary PC1-vs-PC2 overview.
options(stringsAsFactors=FALSE,scipen=999)
needed <- c("analysis_pareto", "sample_metadata", "analysis_root")
if (!all(vapply(needed, exists, logical(1), inherits=TRUE))) {
  candidates <- c("01_load_data_and_PCA.R", file.path("RStudio","01_load_data_and_PCA.R"),
                  file.path("Code availability","01_load_data_and_PCA.R"))
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f)) stop("Run Code 01 from the analysis project")
  source(f)
}
out <- file.path(analysis_root,"Exploratory","PCA")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
model <- prcomp(t(analysis_pareto),center=FALSE,scale.=FALSE)
variance <- 100*model$sdev^2/sum(model$sdev^2)
description <- rbindlist(lapply(1:10,function(k) {
  fit <- lm(model$x[,k] ~ day + passage, data=sample_metadata)
  sums <- anova(fit)[["Sum Sq"]]
  data.table(component=k,variance_percent=variance[k],
    day_SS_percent=100*sums[1]/sum(sums),passage_SS_percent=100*sums[2]/sum(sums),
    residual_SS_percent=100*sums[3]/sum(sums))
}))
fwrite(description,file.path(out,"PCA_component_description.csv"),bom=TRUE)
d <- as.data.table(model$x[,1:10,drop=FALSE],keep.rownames="biological_sample")
d <- sample_metadata[d,on="biological_sample"]
fwrite(d,file.path(out,"PCA_scores.csv"),bom=TRUE)
colors <- c("Early Passage"="#3B6FB6","Middle Passage"="#E69F00","Late Passage"="#B24C63")
shapes <- c(D0=21,D1=22,D3=24,D6=25,D9=8)
for (y in c(2L,5L)) {
  p <- ggplot(d,aes(PC1,.data[[paste0("PC",y)]],color=passage,shape=day))+
    stat_ellipse(aes(x=PC1,y=.data[[paste0("PC",y)]],color=passage,group=passage),
      inherit.aes=FALSE,type="t",level=.95,linewidth=.6,show.legend=FALSE)+
    geom_point(size=3.2,stroke=.9,fill="white")+
    scale_color_manual(values=colors)+scale_shape_manual(values=shapes)+
    labs(title=sprintf("PCA: PC1 vs PC%d",y),
      x=sprintf("PC1 (%.1f%% variance)",variance[1]),
      y=sprintf("PC%d (%.1f%% variance)",y,variance[y]),color="Passage",shape="Day")+
    theme_bw(base_size=13)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"),axis.title=element_text(face="bold"),
      legend.title=element_text(face="bold"))
  name <- sprintf("PCA_PC1_PC%d_passage_ellipses",y)
  ggsave(file.path(out,paste0(name,".png")),p,width=7.5,height=5.5,dpi=300,bg="white")
  ggsave(file.path(out,paste0(name,".svg")),p,width=7.5,height=5.5)
}
print(description)
