#line plots
#data
#trait           model           gen             corr            regr            pi             
#eah             BayesCPi        4                0.3286         0.8895         0.944998      
#eah             BayesCPi        5                0.3456         0.9439         0.960500      
#eah             BayesCPi        6                0.3565         1.0356         0.918984       
#eah             BayesB95        4                0.3316         0.6764         0.950000      
#eah             BayesB95        5                0.3323         0.6899         0.950000      
#eah             BayesB95        6                0.3455         0.7769         0.950000       
#

library(ggplot2)
pdf("eah1.pdf",width=8,height=6)
eah = read.table("eah.txt", header=TRUE)
plot = ggplot(subset(eah, model %in% c("QTLmodel","BayesNPi_allin","BayesCPi","BayesNPi","ssGBLUP")), aes(gen, corr, col=model))
plot = plot + stat_summary(fun.y=mean, geom="point",size=3)
plot = plot + stat_summary(fun.y=mean, geom="line",size=1)
#plot = plot + scale_colour_brewer(palette="Set1")
plot = plot + ylab("Correlation between GEBV and Phenotypes\n") + xlab("Number of Pedigree Generations")
plot = plot + theme_bw() + theme(legend.title=element_blank()) + theme(legend.key = element_blank())
plot = plot + scale_color_manual(values=c("#E69F00", "#009E73", "#0072B2", "red", "#CC79A7"))
#plot = plot + theme_set(theme_gray(base_size=2))
plot
dev.off()
