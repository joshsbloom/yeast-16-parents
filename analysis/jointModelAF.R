library(Hmisc)
library(tidyr)
library(viridis) 
library(easynls)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)
library(ggsci)

load('/data/rrv2/genotyping/RData/joint.cross.cnt.RData')
cross.count.lookup=stack(joint.cross.cnt)
#list by trait then cross
# use recoded genotypes such that direction of effect is - if reference allele increases, + if alternate allele increases

# removed from loop below
    # cross=cross.list[[cross.name]]
    #if(cross.name=='A') {       cross=subset(cross, ind=!grepl('A11', as.character(cross$pheno$id)))    }
    #snames = as.character(cross$pheno$id)
    ## use genotypes coded as reference vs non-reference 

    #g=pull.argmaxgeno(cross)
    # recode based on parental genotypes 
    #seg.pcode=recode.as.allele.number(t(g),parents.list[[cross.name]])
    #seg.recoded[[cross.name]]=seg.pcode
   
    # recode effects as BY =reference
    # are there fixed loci ?? (no)-------------------------------
    #g.af=apply(g,2,function(x) sum(x==1))
    #parents.list[[cross.name]]$fixed=(g.af==0 | g.af==nrow(g))
    #fixed.loci=which(parents.list[[cross.name]]$fixed)
    #if(length(fixed.loci)>0) {    g=g[,-fixed.loci] }
    #------------------------------------------------------------
    #g.r=g[,-which(duplicated(g, MARGIN=2))]
#    subPheno=lapply(NORMpheno, function(x) x[match(snames, names(x))])
#    mPheno  =sapply(subPheno, function(x) sapply(x, mean, na.rm=T))
#    mPheno=apply(mPheno,2, function(x) {x[is.na(x)]=mean(x, na.rm=T); return(x)})
#    mPhenos=scale(mPheno)


jointPeakEffects=list()
for(cross.name in crosses) {
    print(cross.name)
    g.s=scale(t(seg.recoded[[cross.name]]))
    #rename columns
    colnames(g.s)=parents.list[[cross.name]]$marker.name
    mPheno=pheno_extracted[[cross.name]]
    mPhenos  = scale(mPheno)

    for(tt in names(jPs)) {
         sig.joint.markers=parents.list[[cross.name]]$marker.name[ na.omit(match(jPs[[tt]]$fscan.markers[jPs[[tt]]$q<.05], parents.list[[cross.name]]$marker.name))]
         bad.markers=duplicated(sig.joint.markers)
         if(sum(bad.markers)>0) {
             sig.joint.markers=sig.joint.markers[-which(bad.markers)]
         }
         dff=data.frame(g.s[,sig.joint.markers])
         fC=findCorrelation(cor(dff), cutoff=.99)
         if(length(fC)>0) {
             sig.joint.markers=sig.joint.markers[-fC]
             dff=data.frame(g.s[,sig.joint.markers])
         }
         qmodel=lm(mPhenos[,tt]~.-1, data=dff)
         yr=residuals(qmodel)
         aov.a = anova(qmodel)
         tssq  = sum(aov.a[,2])
         a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
         ps=drop1(qmodel, test='F')[-1,6]
         jointPeakEffects[[tt]][[cross.name]]=data.frame(
            trait=tt,
            cross=cross.name,
            peaks=sig.joint.markers,
            betas=as.vector(coef(qmodel)),
            vexp=a.effs,
            p = ps,
            maf1012=iseq.freqs$maf1012[match(sig.joint.markers, iseq.freqs$marker)],
            # iseq.freqs[,8] is the allele frequency of the alternate allele 
            alt012=iseq.freqs[,8][match(sig.joint.markers, iseq.freqs$marker)],
            cross.cnt=cross.count.lookup[match(sig.joint.markers, rownames(cross.count.lookup)),'values'],
            stringsAsFactors=F)
          print(jointPeakEffects[[tt]][[cross.name]])

         }
}

#save(jointPeakEffects, file='/data/rrv2/genotyping/RData/jointPeakEffects_JS_variants.RData')
load('/data/rrv2/genotyping/RData/jointPeakEffects_JS_variants.RData')
r2=rbindlist(lapply(jointPeakEffects, rbindlist))
r2=r2[!(r2$trait %in% c("YPD;;2","YPD;;3")),]
r2$maf1012Fill=r2$maf1012
r2$crossCount=ifelse(r2$cross.cnt%%2==1, r2$cross.cnt+1, r2$cross.cnt)
r2$maf1012Fill[is.na(r2$maf1012Fill)]=c((r2$crossCount/2)/1012)[is.na(r2$maf1012Fill)]
#r2$density=get_density(r2$maf1012, abs(r2$betas))
#r2$densityF=get_density(r2$maf1012Fill, abs(r2$betas))
r2$absBeta=abs(r2$betas)
r2$ancestral=iseq.freqs$ancestral[match(r2$peaks, iseq.freqs$marker)]
r2$maf2=ifelse(r2$ancestral, r2$alt012, 1-r2$alt012)

# 0 = reference derived, alternate paradoxus (ancestral)
# 1 = reference paradoxus (ancestral), alternate cerevisiae


#r2MeanEffect=r2 %>% group_by_(.dots=c("trait","peaks")) %>% mutate(MeanAbsBeta=mean(absBeta)) %>% distinct(traits,peaks,maf1012,crossCount, MeanAbsBeta, ancestral, maf2)
r2MeanEffect=r2 %>% group_by_(.dots=c("trait","peaks")) %>% mutate(MeanAbsBeta=mean(betas)) %>% distinct(traits,peaks,maf1012,crossCount, MeanAbsBeta, ancestral, maf2)
r2VarEffect=r2 %>% group_by_(.dots=c("trait","peaks")) %>% mutate(VarAbsBeta=sd(betas)) %>% distinct(traits,peaks,maf1012,crossCount, VarAbsBeta, ancestral, maf2)

#r2MeanEffect$maf1012Fill=r2MeanEffect$maf1012
#r2MeanEffect$maf1012Fill[is.na(r2MeanEffect$maf1012Fill)]=c((r2MeanEffect$crossCount/2)/1012)[is.na(r2MeanEffect$maf1012Fill)]
r2MeanEffect$density=get_density(r2MeanEffect$maf1012, abs(r2MeanEffect$MeanAbsBeta))
#r2MeanEffect$densityF=get_density(r2MeanEffect$maf1012Fill, abs(r2MeanEffect$MeanAbsBeta))
#r2MeanEffect$ancestral=iseq.freqs$ancestral[match(r2MeanEffect$peaks, iseq.freqs$marker)]

r2MeanEffect$geneResolved=FALSE
r2MeanEffect$geneResolved.ORF=NA
r2MeanEffect$geneResolved.NAME=NA
r2MeanEffect$geneResolved.localFDR=NA

#lookup if variant is resolved in QTGsorted.resolved
#iterate through resolved peaks
load('/data/rrv2/genotyping/RData/FDR_QTGsorted.RData')
QTGsorted.resolved=QTGsorted[QTGsorted$localFDR<.2,]

for(i in 1:nrow(QTGsorted.resolved)){
    print(i)
    um=unique(c(names(QTGsorted.resolved[i,]$pCausal[[1]]), names(QTGsorted.resolved[i,]$pCausal.1[[1]])))
    ut=QTGsorted.resolved$trait[i]
    uorf=QTGsorted.resolved$ORF[i]
    uname=QTGsorted.resolved$NAME[i]
    ufdr=QTGsorted.resolved$localFDR[i]
    r2MeanEffect$geneResolved[which(r2MeanEffect==ut & r2MeanEffect$peaks %in% um)]=T
    r2MeanEffect$geneResolved.ORF[which(r2MeanEffect==ut & r2MeanEffect$peaks %in% um)]=uorf
    r2MeanEffect$geneResolved.NAME[which(r2MeanEffect==ut & r2MeanEffect$peaks %in% um)]=uname
    r2MeanEffect$geneResolved.localFDR[which(r2MeanEffect==ut & r2MeanEffect$peaks %in% um)]=ufdr
}

fisher.test(sapply(split(sign(r2MeanEffect$MeanAbsBeta), r2MeanEffect$maf2>.01), table))
#p-value = 3e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.198 1.664
#sample estimates: odds ratio  1.411 
#sr2=split(r2MeanEffect$MeanAbsBeta, r2MeanEffect$maf1012<.01)
#t.test(abs(sr2[[1]]),abs(sr2[[2]]))


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
uaf=ggplot(r2MeanEffect)+geom_point(alpha=.5, size=1, aes(x=maf2,y=MeanAbsBeta))+
     #scale_color_viridis()+
    scale_x_continuous(name='unfolded allele frequency spectrum (derived allele frequency)', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='average absolute effect, SD units', limits=c(-1,1), breaks=seq(-1,1,.1), expand=c(0,0))+theme_bw()
#ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/unfolded_spectrum.png', width=11,height=8)
ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/112918/unfolded_spectrum.png', width=7,height=5)

ggplot(r2MeanEffect)+geom_point(size=.5, aes(x=maf1012,y=abs(MeanAbsBeta),color=density))+
    scale_color_viridis(option = "inferno", direction=1,end=1)+
    scale_x_continuous(name='folded allele frequency spectrum (minor allele frequency)', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()+
    theme(panel.grid.major = element_line(colour = "#80808022"))+
    theme(panel.grid.minor = element_line(colour = "#80808022")) 

ggsave('/home/jbloom/Dropbox/Lab Meeting - Presentations/112918/folded_spectrum.png', width=7,height=5)
#rre=sapply(seq(0.1,.9,.1), function(x)
#      table(r2MeanEffect$maf1012>.01, abs(r2MeanEffect$MeanAbsBeta)>x))
#rownames(rre)=c('maf<.01_eff<T', 'maf<.01_eff>T', 'maf>.01_eff<T', 'maf<.01_eff>T')
#barplot(rre, beside=T, legend=rownames(rre), names.arg=seq(0.1,.9,.1), xlab='T')
#ggarrange(af,uaf,ncol=1, nrow=2, labels=c('a','b') )


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

