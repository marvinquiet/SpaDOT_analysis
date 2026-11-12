# Fig. S45. Usage: Rscript render_s45.R <prepared-input-dir> <new-output-dir>
# Plot structure follows the user's plot_OT_transition_metrics.R:
# dodged identity bars, log2FC over SpaDOT, free-x facets, theme_bw.
# Only simulation faceting, seed aggregation, palette, and missingness are added.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
input <- args[[1]]; output <- args[[2]]
if(dir.exists(output)) stop('Output collision')
dir.create(output,recursive=TRUE)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
x <- read.csv(file.path(input,'transition_accuracy_by_scenario_method.csv'))
fc <- read.csv(file.path(input,'SpaDOT_fold_changes.csv'))
x <- x[x$mode=='single' | x$method=='SpaDOT',]
colors <- c(SpaDOT='#c1121f',PRECAST='#68a756',SEDR='#f09232',spaGCN='#4c2623',GraphST='#165451',spatialPCA='#d3b300',spatiAlign='#ded6c1',spaVAE='#0082c0')
metrics <- c('annotation_naive_accuracy')
pdf(file.path(output,'S45_Simulation_transition_accuracy.pdf'),width=18,height=15,onefile=TRUE)
for(metric in metrics){
  df <- x %>% group_by(deg,proportion,transition,method) %>%
    summarise(naive_accuracy=if(all(is.na(.data[[metric]]))) NA_real_ else mean(.data[[metric]],na.rm=TRUE),
              n=sum(!is.na(.data[[metric]])),.groups='drop')
  # Compare only conditions with all three independent simulation seeds.
  df$naive_accuracy[df$n!=3] <- NA_real_
  gap <- fc[fc$metric==metric,c('deg','proportion','transition','fold')]
  gap$log2FC <- log2(gap$fold)
  df <- left_join(df,gap,by=c('deg','proportion','transition'))
  ranking <- df %>% group_by(method) %>% summarise(avg=mean(naive_accuracy,na.rm=TRUE),.groups='drop') %>% arrange(desc(avg))
  df$method <- factor(df$method,levels=c('SpaDOT',setdiff(ranking$method,'SpaDOT')))
  df$Data <- factor(df$deg,levels=c(100,200,500),labels=c('DEG 100','DEG 200','DEG 500'))
  df$Proportion <- factor(df$proportion,levels=c(.5,.6,.7,.8,.9),labels=paste0('p = ',c(.5,.6,.7,.8,.9)))
  df$label <- factor(df$transition,levels=c('t0_t1','t1_t2'),labels=c('t0 -> t1','t1 -> t2'))
  g <- ggplot(df,aes(x=label,y=naive_accuracy,fill=method)) +
    geom_bar(stat='identity',position=position_dodge(),na.rm=TRUE) +
    geom_text(data=subset(df,method=='SpaDOT' & is.finite(log2FC)),
              aes(label=paste0('log2FC=',sprintf('%.2f',log2FC)),color=log2FC>0),
              vjust=-.6,size=4,fontface='bold') +
    scale_color_manual(values=c('TRUE'='green','FALSE'='red'),guide='none') +
    scale_fill_manual(values=colors,drop=FALSE) +
    facet_grid(Proportion ~ Data,scales='free_x',space='free_x') +
    scale_y_continuous(limits=c(0,1.15),breaks=seq(0,1,.25),expand=expansion(mult=c(0,.1))) +
    theme_bw() +
    theme(strip.text=element_text(size=15,face='bold'),strip.background=element_blank(),
          axis.title=element_text(size=16,face='bold'),axis.text.x=element_text(size=14,angle=45,hjust=1),
          axis.text.y=element_text(size=13),legend.title=element_text(size=15,face='bold'),
          legend.text=element_text(size=14),legend.position='right') +
    guides(fill=guide_legend(ncol=1)) +
    labs(x='Transitions',y='Transition Accuracy',fill='Methods')
  text_layer <- ggplot_build(g)$data[[2]]
  stopifnot(all(text_layer$colour == 'green'))
  stopifnot(identical(unname(col2rgb('green')[,1]), c(0L,255L,0L)))
  write.csv(text_layer[,c('label','colour')],file.path(output,'fold_label_colors.csv'),row.names=FALSE)
  print(g)
  ggsave(file.path(output,paste0(metric,'.png')),g,width=18,height=15,dpi=140)
  write.csv(df,file.path(output,paste0(metric,'_plot_data.csv')),row.names=FALSE)
}
dev.off()
