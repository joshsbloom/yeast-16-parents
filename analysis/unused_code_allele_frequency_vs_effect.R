fisher.test(rbind(table((derived.eff)<(-.5)),
                  c(table(ancestral.eff<(-.5)),0)))


fisher.test(rbind(table((derived.eff)>(.5)),
                  c(table(ancestral.eff>(.5)))))





t(table(abs(r2MeanEffect$MeanAbsBeta)>.1, r2MeanEffect$maf1012>.01))


split(abs(r2MeanEffect$MeanAbsBeta), r2MeanEffect$maf2>.01)
wilcox.test(abs(r2MeanEffect$MeanAbsBeta)[r2MeanEffect$maf2>.95], abs(r2MeanEffect$MeanAbsBeta)[r2MeanEffect$maf2<.05])
median(abs(r2MeanEffect$MeanAbsBeta)[r2MeanEffect$maf2<.05], na.rm=T)
median(abs(r2MeanEffect$MeanAbsBeta)[r2MeanEffect$maf2>.95], na.rm=T)






fisher.test(rbind(table(abs(derived.eff)>.5),
      table(abs(ancestral.eff)>.5)))

fisher.test(rbind(table((derived.eff)>.5),
                  c(table((ancestral.eff)>.5),0) ))
fisher.test(rbind(table((derived.eff)<(-.5)),
                  c(table((ancestral.eff)<(-.5)),0) ))
#



#r2MeanEffectBig=r2MeanEffect[abs(r2MeanEffect$MeanAbsBeta)>.1,]
#fisher.test(sapply(split(sign(r2MeanEffectBig$MeanAbsBeta), r2MeanEffectBig$maf2>.01), table))
ggplot(r2MeanEffect, aes(x=maf1012>.01,y=abs(MeanAbsBeta)))+
    geom_violin(width=1.3)+
    geom_quasirandom(alpha = .5, width = 0.5)

    #geom_jitter(position='jitter',alpha=.3, size=1)+

#x11()
#plot(r2MeanEffect$maf2, r2MeanEffect$MeanAbsBeta, xlab='ancestral allele frequency', ylab='beta', col='#00000066')

# plots 
# unfolded allele frequency spectrum
Fig3C=ggplot(r2MeanEffect)+geom_point(alpha=.4, size=.75, aes(x=maf2,y=MeanAbsBeta))+
     #scale_color_viridis()+
    scale_x_continuous(name='Derived allele frequency', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='Average effect, SD units', limits=c(-1,1), breaks=seq(-1,1,.1), expand=c(0,0))+theme_bw()
#ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/unfolded_spectrum.png', width=11,height=8)
#ggsave('/data/rrv2/Figures_and_Tables/Figs/Figure4C.png', width=7,height=5)

library(Hmisc)
Fig3B=ggplot(r2MeanEffect)+geom_point(alpha=.4, size=.75, aes(x=maf1012,y=abs(MeanAbsBeta)))+ #,color=density))+
    #scale_color_viridis(option = "inferno", direction=1,end=1)+
    scale_x_continuous(name='Minor allele frequency', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='Average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()+
    theme(panel.grid.major = element_line(colour = "#80808022"))+
    theme(panel.grid.minor = element_line(colour = "#80808022"))+
    theme(axis.text.x=element_text(size=rel(1.1),color='black'),
          axis.text.y=element_text(size=rel(1.1),color='black'))+
    stat_summary_bin(aes(x=maf1012,y=abs(MeanAbsBeta)), breaks=cut2(r2MeanEffect$maf1012,g=42, onlycuts=T),  col='red')
Fig3Blank=ggplot(r2MeanEffect)+geom_blank()+theme_bw()
ggarrange(Fig3A,Fig3Blank,Fig3B, Fig3C, ncol=2, nrow=2, labels=c('A','','B', 'C'))
#ggsave('/home/jbloom/Dropbox/RR/Figures\ and\ Tables/Figure3A.png', width=11, height=11)


#+
#    geom_vline(xintercept=.01, col='red')
#ggsave('/data/rrv2/Figures_and_Tables/Figs/Figure4B.png', width=10,height=7)

ggarrange(Fig4A, labels='A')
ggsave('/home/jbloom/Dropbox/RR/Figures\ and\ Tables/Figure3A.png', width=12,height=5)

ggarrange(Fig4B,Fig4C, ncol=2, labels=c('B', 'C'))
ggsave('/home/jbloom/Dropbox/RR/Figures\ and\ Tables/Figure3BC.png', width=12,height=5)

label_trait=utraits
names(label_trait)=utraits.orig
ggplot(r2MeanEffect)+geom_point(alpha=.4, size=.75, aes(x=maf1012,y=abs(MeanAbsBeta)))+ #,color=density))+
    #scale_color_viridis(option = "inferno", direction=1,end=1)+
    scale_x_continuous(name='Minor allele frequency', breaks=seq(0,1,.1)) +
    scale_y_continuous(name='Average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.2), expand=c(0,0))+theme_bw()+
    facet_wrap(~trait, ncol=4, labeller=as_labeller(label_trait))+
    theme(panel.grid.major = element_line(colour = "#80808022"))+
    theme(panel.grid.minor = element_line(colour = "#80808022"))
ggsave('/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/SupplementaryFigure5.pdf', width=8,height=11)



ggsave('/home/jbloom/Dropbox/RR/Figures\ and\ Tables/Figure4.png')



ggplot(r2MeanEffect)+geom_point(alpha=.4, size=.8, aes(x=log10(maf1012),y=(abs(MeanAbsBeta))))+ 
    theme_bw()+
    scale_y_continuous(name='average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1), expand=c(0,0))+
    scale_x_continuous(name='log10(minor allele frequency)', breaks=log10(c(.001,.01, .05,.1,.2,.5)), labels=c(c(.001,.01, .05,.1,.2,.5)))+
    theme(panel.grid.major = element_line(colour = "#80808022"))+
    theme(panel.grid.minor = element_line(colour = "#80808022"))
    #,color=density))+
    #scale_color_viridis(option = "inferno", direction=1,end=1)+
    #scale_x_continuous(name='folded allele frequency spectrum (minor allele frequency)', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    #scale_y_continuous(name='average absolute effect, SD units') + theme_bw()+
    #, limits=c(0,1), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()+
    #+
    #    geom_vline(xintercept=.01, col='red')
ggsave('/data/rrv2/Figures_and_Tables/Figs/Figure4BLog.png', width=10,height=7)

#rre=sapply(seq(0.1,.9,.1), function(x)
#      table(r2MeanEffect$maf1012>.01, abs(r2MeanEffect$MeanAbsBeta)>x))
#rownames(rre)=c('maf<.01_eff<T', 'maf<.01_eff>T', 'maf>.01_eff<T', 'maf<.01_eff>T')
#barplot(rre, beside=T, legend=rownames(rre), names.arg=seq(0.1,.9,.1), xlab='T')
#ggarrange(af,uaf,ncol=1, nrow=2, labels=c('a','b') )

fisher.test(table(r2MeanEffect$maf1012<.01,abs(r2MeanEffect$MeanAbsBeta)>.3))


t(table(abs(r2MeanEffect$MeanAbsBeta)>.1, r2MeanEffect$maf1012>.01))
t(table(abs(r2MeanEffect$MeanAbsBeta)>.2, r2MeanEffect$maf1012>.01))
t(table(abs(r2MeanEffect$MeanAbsBeta)>.3, r2MeanEffect$maf1012>.01))
t(table(abs(r2MeanEffect$MeanAbsBeta)>.4, r2MeanEffect$maf1012>.01))
t(table(abs(r2MeanEffect$MeanAbsBeta)>.5, r2MeanEffect$maf1012>.01))


table(abs(r2MeanEffect$MeanAbsBeta)>.3163, r2MeanEffect$maf1012>.05)

table(abs(r2MeanEffect$MeanAbsBeta)>.1, r2MeanEffect$maf1012>.01)

r2MeanEffect$absBeta=abs(r2MeanEffect$MeanAbsBeta)
rre=(sapply(seq(0.1,.4,.1), function(x)
      table(r2MeanEffect$maf1012>.01, r2MeanEffect$absBeta>x)))
rownames(rre)=c('rare_small', 'common_small', 'rare_large', 'common_large')
#barplot(rre, beside=T, legend=rownames(rre), names.arg=seq(0.1,.4,.1), xlab='T')
rre2=rre[c(1,3,2,4),]
par(mfrow=c(2,1))
barplot(rre2, beside=T, legend=rownames(rre2), names.arg=paste('>' ,seq(0.1,.4,.1)), xlab='effect size threshold', ylab='count', main='rare is maf<.01')

nf1=colSums(rre2[1:2,])
nf2=colSums(rre2[3:4,])
rre3=rbind(t(t(rre2[1:2,])/nf1), t(t(rre2[3:4,])/nf2))
barplot(rre3, beside=T, legend=rownames(rre3), names.arg=paste('>' ,seq(0.1,.4,.1)), xlab='effect size threshold', ylab='proportion')

# normalize bars total count for rare and for common ...and then two bars ... 
library(Hmisc)
test=cut2(r2MeanEffect$maf1012, g=7)
test=cut(r2MeanEffect$maf1012, b=20)

par(mfrow=c(4,1))
barplot(table(r2MeanEffect$absBeta>.1, test),main=' effect > .05',legend=c('less than', 'greater than'), xlab='maf bin')
barplot(table(r2MeanEffect$absBeta>.1, test),main='effect > .1',legend=c('less than', 'greater than'), xlab='maf bin')
barplot(table(r2MeanEffect$absBeta>.2, test),main='effect > .2',legend=c('less than', 'greater than'), xlab='maf bin')
barplot(table(r2MeanEffect$absBeta>.3, test),main='effect > .3',legend=c('less than', 'greater than'), xlab='maf bin')

R>  table(r2MeanEffect$absBeta>.1, r2MeanEffect$maf1012>.01)
 ###                effect
#                <.1 >.1
#              FALSE TRUE
#AF<.01  FALSE   388 2523
#AF>.01  TRUE    596 1045
sum(r2MeanEffect$maf1012<.01 & r2MeanEffect$absBeta<.1)
#[1] 388
sum(r2MeanEffect$maf1012<.01 & r2MeanEffect$absBeta>.1)
#[1] 596
sum(r2MeanEffect$maf1012>.01 & r2MeanEffect$absBeta>.1)
#[1] 1045
sum(r2MeanEffect$maf1012>.01 & r2MeanEffect$absBeta<.1)
#[1] 2523

sum(r2MeanEffect$maf1012<.05 & r2MeanEffect$absBeta<.1)
#[1] 388
sum(r2MeanEffect$maf1012<.05 & r2MeanEffect$absBeta>.1)
#[1] 596
sum(r2MeanEffect$maf1012>.05 & r2MeanEffect$absBeta>.1)
#[1] 1045
sum(r2MeanEffect$maf1012>.05 & r2MeanEffect$absBeta<.1)


ggplot(r2MeanEffect, aes(x=(maf2), y=MeanAbsBeta))+geom_point(alpha=.3, size=.75)+facet_wrap(~trait)+
       scale_x_continuous(name='unfolded allele frequency spectrum (derived allele frequency)', breaks=seq(0,1,.05), expand=c(0.01,0)) +
       scale_y_continuous(name='average absolute effect, SD units', limits=c(-1,1), breaks=seq(-1,1,.1), expand=c(0,0))

    scale_x_continuous(name='peak marker minor allele frequency from 1012 yeast panel', breaks=seq(0,.5,.1)) +
    scale_y_continuous(name='average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1))+
    geom_smooth(method='lm', formula=y~x, size=.25)+theme_bw()
test=cbind(ifelse(r2MeanEffect$maf2>.01, 'daf >.01', 'daf<.01'), sign(r2MeanEffect$MeanAbsBeta))
table(data.frame(test))

x=(sapply(split(sign(r2MeanEffect$MeanAbsBeta), cut(r2MeanEffect$maf2, 100)), table))
#test=cbind(r2MeanEffect$maf2<.05, r2MeanEffect$MeanAbsBeta<0)
#table(data.frame(test))
#              -   +
#   af<.05   1225 1130
#   af>.05    712  928
library(viridis)
ggplot(r2MeanEffect)+geom_jitter(alpha=.6, size=.75, aes(x=crossCount,y=MeanAbsBeta,color=density))+scale_color_viridis()+
    scale_x_continuous(name='# of crosses a peak marker segregates in', breaks=seq(0,14,2)) +
    scale_y_continuous(name='average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1))

ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/04_jointEffectDist_perCross.png', width=11,height=8)

# all effects per QTL
nlsfit= nlsfit(data.frame(r2$maf1012,abs(r2$betas)), model=6)
a <- nlsfit$Parameters[row.names(nlsfit$Parameters) == 'coefficient a',]
b <- nlsfit$Parameters[row.names(nlsfit$Parameters) == 'coefficient b',]
ggplot(r2)+geom_point(alpha=.5, size=.75, aes(x=maf1012,y=abs(betas),color=density))+scale_color_viridis()+
    scale_x_continuous(name='peak marker minor allele frequency from 1012 yeast panel', breaks=seq(0,.5,.05)) +
    scale_y_continuous(name='average absolute effect, SD units', breaks=seq(0,1,.1))+
    stat_function(fun=function(x) a*exp(b*x), colour = "red")+
     stat_smooth(method='lm', aes(x=maf1012,y=abs(betas)), formula=y~x)
ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/04_jointEffectDist.png', width=11,height=8)

# mean effect per QTL
nlsfit= nlsfit(data.frame(r2MeanEffect$maf1012,r2MeanEffect$MeanAbsBeta), model=6)
a <- nlsfit$Parameters[row.names(nlsfit$Parameters) == 'coefficient a',]
b <- nlsfit$Parameters[row.names(nlsfit$Parameters) == 'coefficient b',]
ggplot(r2MeanEffect)+geom_point(alpha=.4, size=1, aes(x=maf1012,y=MeanAbsBeta,color=density))+#scale_color_viridis()+
    scale_x_continuous(name='peak marker minor allele frequency from 1012 yeast panel', breaks=seq(0,.5,.05)) +
    scale_y_continuous(name='average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1))
    #stat_function(fun=function(x) a*exp(b*x), colour = "red")+
    #stat_smooth(method='lm', aes(x=maf1012,y=MeanAbsBeta), formula=y~x)
ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/04_jointEffectDistQTLMean.png', width=11,height=8)


ggplot(r2MeanEffect, aes(x=(maf1012), y=MeanAbsBeta))+geom_point(alpha=.3, size=.75)+facet_wrap(~trait)+
    scale_x_continuous(name='peak marker minor allele frequency from 1012 yeast panel', breaks=seq(0,.5,.1)) +
    scale_y_continuous(name='average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1))+
    geom_smooth(method='lm', formula=y~x, size=.25)+theme_bw()
ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/04_perTraitQTLMean.png', width=11,height=11)



ggplot(r2MeanEffect, aes(x=(maf1012), y=MeanAbsBeta))+geom_point(alpha=.5, size=.75)+facet_wrap(~trait)+
    scale_x_log10(name='peak marker minor allele frequency from 1012 yeast panel', limits=c(.001,.75), breaks=c(0,0.001,0.01,0.25)) +
    scale_y_log10(name='average absolute effect, SD units', limits=c(.01,1))+geom_smooth(method='lm', formula=y~x, size=.25)+theme_bw()
ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/04_perTraitQTLMeanLOG10.png', width=11,height=11)




    #stat_function(fun=function(x) a*exp(b*x), colour = "red")+
    #stat_smooth(method='lm', aes(x=maf1012,y=MeanAbsBeta), formula=y~x)




ggplot(r2, aes(x=(maf1012), y=abs(betas)))+geom_point(alpha=.15, size=.75)+facet_wrap(~trait)+
    scale_x_continuous(name='minor allele frequency from 1012 yeast panel', limits=c(.001,.75)) +
    scale_y_continuous(name='absolute effect, SD units', limits=c(.0001,1))+geom_smooth(method='lm', formula=y~x)


df2=r2MeanEffect
mclustBIC(cbind(df2$MeanAbsBeta, df2$maf1012Fill))
test=mclustBIC(cbind(df2$MeanAbsBeta, df2$maf1012Fill))

#spikeAndSlab(df$MeanAbsBeta, r2MeanEffect$maf1012Fill)
ggplot(r2MeanEffect, aes(x=(maf1012), y=r2MeanEffect$MeanAbsBeta))+geom_point(alpha=.7, size=1)+facet_wrap(~trait)+
    scale_x_continuous(name='minor allele frequency from 1012 yeast panel', limits=c(.001,.5)) +
    scale_y_continuous(name='absolute effect, SD units', limits=c(.0001,1))+geom_smooth(method='lm', formula=y~x)






nlsfit= nlsfit(data.frame(r2$maf1012Fill,abs(r2$betas)), model=6)
a <- nlsfit$Parameters[row.names(nlsfit$Parameters) == 'coefficient a',]
b <- nlsfit$Parameters[row.names(nlsfit$Parameters) == 'coefficient b',]

ggplot(r2)+geom_point(alpha=.15, size=.75, aes(x=maf1012Fill,y=abs(betas),color=densityF))+scale_color_viridis()+
    scale_x_continuous(name='minor allele frequency from 1012 yeast panel') +
    scale_y_continuous(name='absolute effect, SD units')+
    stat_function(fun=function(x) a*exp(b*x), colour = "red")

ggplot(r2)+geom_point(alpha=.15, size=.75, aes(x=maf1012Fill,y=abs(betas),color=densityF))+scale_color_viridis()+
    scale_x_continuous(name='minor allele frequency from 1012 yeast panel') +
    scale_y_continuous(name='absolute effect, SD units')+
    stat_function(fun=function(x) a*exp(b*x), colour = "red")








ggplot(r2, aes(x=(maf1012), y=abs(betas), color=cross))+geom_point(alpha=.25, size=.75)+facet_wrap(~trait)
ggplot(r2, aes(x=(maf1012), y=abs(betas)))+geom_point(alpha=.15, size=.75)+facet_wrap(~trait)+
    scale_x_log10(name='minor allele frequency from 1012 yeast panel', limits=c(.001,.5)) +
    scale_y_log10(name='absolute effect, SD units', limits=c(.0001,1))+geom_smooth(method='lm', formula=y~x)

ggplot(r2, aes(x=(maf1012Fill), y=abs(betas)))+geom_point(alpha=.15, size=.75)+
    scale_x_log10(name='minor allele frequency from 1012 yeast panel', limits=c(.001,.5)) +
    scale_y_log10(name='absolute effect, SD units', limits=c(.0001,1))+geom_smooth(method='lm', formula=y~x)

library(MASS)
library(ggplot2)
library(viridis)
theme_set(theme_bw(base_size = 16))

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.


    stat_smooth(method='lm', formula=log(y)~x)


+facet_wrap(~trait)

#r=rbindlist(jointPeakEffects[[35]])

par(mfrow=c(3,1))
plot(jitter(test), abs(r2$betas), col='#00000022',cex=.75, xlab='# of crosses a variant is segregating in', ylab='sd units')
#r2$maf1012[is.na(r2$maf1012)]=0
plot(r2$maf1012, abs(r2$betas), col='#00000044',cex=.75, xlab='maf in 1012 panel', ylab='sd units')
plot(r2$maf1012, abs(r2$betas), ylim=c(0,0.4), col='#00000044',cex=.75, xlab='maf in 1012 panel', ylab='variance explained')


plot(r2$maf1012, abs(r2$vexp)




    cpeaks=cross.peaks[[cross.name]] #[grep('Mang', cross.peaks[[cross.name]]$gene),]
    cpeaks=separate(cpeaks, pmarker, c('chr', 'pos', 'ref', 'alt', 'index'), sep='_', convert=T, remove=F)
    cpeaks$marker=parents.list[[cross.name]]$marker.name[match(cpeaks$pmarker, parents.list[[cross.name]]$marker.name.n)]
    cpeaks.by.trait=split(cpeaks, cpeaks$trait)
    chromosomes=paste0('chr', as.roman(1:16))
    cvec=parents.list[[cross.name]]$chr
    gs.by.chr=list()
    vac=list()
    for(cc in chromosomes) {   
            gs.by.chr[[cc]]=g.s[,which(cvec %in% cc)]   
            vac[[cc]]=parents.list[[cross.name]][which(cvec %in% cc),]
    }
    
     for(tt in names(cpeaks.by.trait)) {
         print(tt)
         fmodel=cpeaks.by.trait[[tt]]
         bad.markers=duplicated(fmodel$marker)
         if(sum(bad.markers)>0) {
             fmodel=fmodel[-which(bad.markers),]
         }
         dff=data.frame(g.s[,fmodel$marker])
         qmodel=lm(mPhenos[,tt]~.-1, data=dff)
         yr=residuals(qmodel)
         aov.a = anova(qmodel)
         tssq  = sum(aov.a[,2])
         a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
         ps=drop1(qmodel, test='F')[-1,6]
         out=data.frame(
            trait=tt,
            fmodel[,c(3:7,13)],
            sig_covar=fmodel$marker,
            CI.l=parents.list[[cross.name]]$marker.name[match(fmodel$CI.l, parents.list[[cross.name]]$marker.name.n)],
            CI.r=parents.list[[cross.name]]$marker.name[match(fmodel$CI.r, parents.list[[cross.name]]$marker.name.n)],
            betas=as.vector(coef(qmodel)),
            vexp=a.effs,
            chr=fmodel$chr,
            p = ps,
            bootstrap_interval='', stringsAsFactors=F)


wilcox.test(sapply(rjq[names(which(sapply(rjq, function(x) sum(x$overlapQTL)>0)))], function(x) mean(abs(x$crossCount))),
            sapply(rjq[names(which(sapply(rjq, function(x) sum(x$overlapQTL)==0)))], function(x) mean(abs(x$crossCount))))
#wilcox.test(r2MeanEffect$crossCount, sapply(rjq[names(which(sapply(rjq, function(x) sum(x$overlapQTL)==0)))], function(x) mean(abs(x$maf1012))) ) 


#jpel2=list()
#for(i in 1:nrow(QTGsorted.resolved)) {
#    print(i)
#    tt=QTGsorted.resolved[i,]$trait
#    
#    cc=QTGsorted.resolved[i,]$cross
#    jpt=jointPeakEffects[[tt]][[cc]]     
#    chr2=tstrsplit(jpt$peaks ,'_', type.convert=T)[[1]]
#    ppos3=tstrsplit(jpt$peaks, '_', type.convert=T)[[2]]
#    jpti=makeGRangesFromDataFrame(data.frame(chr=chr2, start=ppos3, end=ppos3 ), ignore.strand=T)
#    cpti=QTGsorted.resolved[i,]$CI95
#    overlapQ=unique(queryHits(findOverlaps(jpti, cpti)))
#    if(length(overlapQ)==0) { 
#     cc=QTGsorted.resolved[i,]$cross.1
#     jpt=jointPeakEffects[[tt]][[cc]]     
#     chr2=tstrsplit(jpt$peaks ,'_', type.convert=T)[[1]]
#     ppos3=tstrsplit(jpt$peaks, '_', type.convert=T)[[2]]
#     jpti=makeGRangesFromDataFrame(data.frame(chr=chr2, start=ppos3, end=ppos3 ), ignore.strand=T)
#     cpti=QTGsorted.resolved[i,]$CI95.1
#     overlapQ=unique(queryHits(findOverlaps(jpti, cpti)))
#    } 
#    jpel2[[i]]=jpt[overlapQ,]
#
#
#}
#jpel2=rbindlist(jpel2)
#length(split(jpel2, paste(jpel2$trait, jpel2$peaks)))


