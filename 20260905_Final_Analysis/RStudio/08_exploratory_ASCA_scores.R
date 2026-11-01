# Exploratory ASCA-guided factor-focused scores; not global PCA or standard PLS-DA.
options(stringsAsFactors=FALSE,scipen=999)
needed <- c("analysis_pareto","sample_metadata","analysis_root","input_file")
if (!all(vapply(needed,exists,logical(1),inherits=TRUE))) {
  candidates <- c("01_load_data_and_PCA.R",file.path("RStudio","01_load_data_and_PCA.R"),
                  file.path("Code availability","01_load_data_and_PCA.R"))
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f)) stop("Run Code 01 from the analysis project")
  source(f)
}
out <- file.path(analysis_root,"Exploratory","ASCA")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
X <- t(analysis_pareto)
meta <- copy(sample_metadata)
stopifnot(nrow(X)==45L,ncol(X)==485L,
 identical(rownames(X),meta$biological_sample),
 all(meta[,.N,by=.(day,passage)]$N==3L),
 unname(tools::md5sum(input_file))=="0c1aa6b5304897d2178794ee9bc62a31")
grand <- colMeans(X)
centered <- sweep(X,2,grand,"-")
mean_rows <- function(groups) {
 answer <- X
 for(g in unique(groups)) answer[groups==g,] <- matrix(colMeans(X[groups==g,,drop=FALSE]),
 sum(groups==g),ncol(X),byrow=TRUE)
 answer
}
D <- sweep(mean_rows(meta$day),2,grand,"-")
P <- sweep(mean_rows(meta$passage),2,grand,"-")
cell_mean <- mean_rows(interaction(meta$day,meta$passage,drop=TRUE))
I <- sweep(cell_mean,2,grand,"-")-D-P
E <- X-cell_mean
stopifnot(max(abs(centered-D-P-I-E))<1e-10)
totalSS <- sum(centered^2)
effects <- list(Day=D,Passage=P,Interaction=I,Within_condition=E)
factor_summary <- data.table(effect=names(effects),
 SS=vapply(effects,function(z)sum(z^2),numeric(1)))
factor_summary[,percent_total_X_SS:=100*SS/totalSS]
stopifnot(abs(sum(factor_summary$SS)-totalSS)<1e-7)
# Independent model-matrix check of main and interaction effect matrices.
mm <- model.matrix(~day*passage,data=meta,
 contrasts.arg=list(day=contr.sum(5),passage=contr.sum(3)))
fit <- lm.fit(mm,X)
coef <- fit$coefficients
assignment <- attr(mm,"assign")
fitted_D <- mm[,assignment==1,drop=FALSE]%*%coef[assignment==1,,drop=FALSE]
fitted_P <- mm[,assignment==2,drop=FALSE]%*%coef[assignment==2,,drop=FALSE]
fitted_I <- mm[,assignment==3,drop=FALSE]%*%coef[assignment==3,,drop=FALSE]
stopifnot(max(abs(D-fitted_D))<1e-9,max(abs(P-fitted_P))<1e-9,
 max(abs(I-fitted_I))<1e-9,max(abs(E-fit$residuals))<1e-9)
day_fit <- prcomp(D,center=FALSE,scale.=FALSE)
passage_fit <- prcomp(P,center=FALSE,scale.=FALSE)
ld <- day_fit$rotation[,1]
lp <- passage_fit$rotation[,1]
# Signs only: later days rightward; Late higher than Early.
if(cor(as.numeric(meta$day_numeric),as.vector(D%*%ld))<0)ld <- -ld
passage_order <- match(as.character(meta$passage),c("Early Passage","Middle Passage","Late Passage"))
if(cor(passage_order,as.vector(P%*%lp))<0)lp <- -lp
day_within <- 100*sum((D%*%ld)^2)/sum(D^2)
passage_within <- 100*sum((P%*%lp)^2)/sum(P^2)
scores <- copy(meta)
scores[,day_effect_score:=as.vector(D%*%ld)]
scores[,passage_effect_score:=as.vector(P%*%lp)]
# Separate covariate adjustments: retain the interaction and all within-cell deviations.
# Each coordinate uses a different adjusted matrix. This is not a 2-PC plot.
scores[,day_score:=as.vector((centered-P)%*%ld)]
scores[,passage_score:=as.vector((centered-D)%*%lp)]
scores[,day_unadjusted:=as.vector(centered%*%ld)]
scores[,passage_unadjusted:=as.vector(centered%*%lp)]
stopifnot(max(abs(scores$day_score-as.vector((D+I+E)%*%ld)))<1e-9,
 max(abs(scores$passage_score-as.vector((P+I+E)%*%lp)))<1e-9)
# All residual and interaction terms remain; no manually moved points or jitter.
fwrite(scores,file.path(out,"factor_focused_sample_scores.csv"),bom=TRUE)
fwrite(factor_summary,file.path(out,"ASCA_factor_variance.csv"),bom=TRUE)
component_summary <- data.table(component=c("Day component 1","Passage component 1"),
 percent_within_its_effect=c(day_within,passage_within),
 percent_of_total_X_from_fitted_effect=c(100*sum((D%*%ld)^2)/totalSS,
 100*sum((P%*%lp)^2)/totalSS))
fwrite(component_summary,file.path(out,"ASCA_first_components.csv"),bom=TRUE)
fwrite(data.table(met_id=rownames(day_fit$rotation),day_loading=ld,passage_loading=lp),
 file.path(out,"ASCA_first_component_loadings.csv"),bom=TRUE)
colors <- c("Early Passage"="#3B6FB6","Middle Passage"="#E69F00","Late Passage"="#B24C63")
shapes <- c(D0=21,D1=22,D3=24,D6=25,D9=8)
p <- ggplot(scores,aes(day_score,passage_score,color=passage,shape=day))+
 stat_ellipse(aes(x=day_score,y=passage_score,color=passage,group=passage),
 inherit.aes=FALSE,type="t",level=.95,linewidth=.6,show.legend=FALSE)+
 geom_point(size=3.3,stroke=.9,fill="white")+
 scale_color_manual(values=colors)+scale_shape_manual(values=shapes)+
 labs(title="Factor-focused scores: Day and Passage",
 x="Day component 1",y="Passage component 1",color="Passage",shape="Day",
 caption=paste("Axes learned separately from Day and Passage effects (ASCA).",
 "Opposite main effect adjusted; interaction and within-condition variation retained.",sep="\n"))+
 theme_bw(base_size=13)+theme(panel.grid.minor=element_blank(),
 plot.title=element_text(face="bold",size=15),axis.title=element_text(face="bold",size=14),
 legend.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=9))
ggsave(file.path(out,"Day_horizontal_Passage_vertical_ASCA.png"),p,width=8,height=6,dpi=300,bg="white")
ggsave(file.path(out,"Day_horizontal_Passage_vertical_ASCA.svg"),p,width=8,height=6)
# Standard effect score plots, explicitly showing effect + within-condition residuals.
standard <- rbindlist(list(
 data.table(meta,axis1=as.vector((D+E)%*%day_fit$rotation[,1]),
 axis2=as.vector((D+E)%*%day_fit$rotation[,2]),effect="Day effect + residual"),
 data.table(meta,axis1=as.vector((P+E)%*%passage_fit$rotation[,1]),
 axis2=as.vector((P+E)%*%passage_fit$rotation[,2]),effect="Passage effect + residual")))
sp <- ggplot(standard,aes(axis1,axis2,color=passage,shape=day))+
 geom_point(size=2.5,stroke=.8,fill="white")+facet_wrap(~effect,scales="free")+
 scale_color_manual(values=colors)+scale_shape_manual(values=shapes)+
 labs(x="Effect component 1",y="Effect component 2",color="Passage",shape="Day")+
 theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),strip.text=element_text(face="bold"))
ggsave(file.path(out,"ASCA_separate_effect_diagnostics.png"),sp,width=11,height=4.5,dpi=240,bg="white")
audit <- data.table(metric=c("raw_input_md5","features","samples","labels_used_to_learn_axes",
 "largest_ANOVA_reconstruction_error","Day_loading_dot_Passage_loading",
 "Day_first_component_within_Day_percent","Passage_first_component_within_Passage_percent"),
 value=as.character(c(unname(tools::md5sum(input_file)),485,45,TRUE,
 max(abs(centered-D-P-I-E)),sum(ld*lp),day_within,passage_within)))
fwrite(audit,file.path(out,"input_and_method_audit.csv"),bom=TRUE)
print(factor_summary);print(component_summary)
print(scores[,.(day_score=mean(day_score),passage_score=mean(passage_score)),by=.(day,passage)])
cat("Input and two-factor decomposition verified; exploratory plots written to:",out,"\n")
