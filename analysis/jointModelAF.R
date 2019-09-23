library(Hmisc)
library(tidyr)
library(viridis) 
library(easynls)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)
library(ggsci)
library(caret)

# load extracted phenotypes (pheno_extracted)
load('/data/rrv2/genotyping/RData/extracted_average_phenotypes.RData')

# load recoded genotype data
# use recoded genotypes such that direction of effect is - if reference allele increases, + if alternate allele increases
load('/data/rrv2/genotyping/RData/FDR_seg.recoded.RData')

# output from genotyping/code/segregants_hmm.R
load('/data/rrv2/genotyping/RData/parents.list.RData')
parents.list=lapply(parents.list, function(x) {
                  z=x;
                  z$marker.name.n=paste0(z$marker.name, '_', seq(1:nrow(z)))
                  return(z) })

# load allele-frequency lookup table 
load('/data/rrv2/genotyping/RData/iseq.freqs.RData')


# additional helper functions
source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')
source('/data/rrv2/analysis/mapping_fx.R')

load('/data/rrv2/genotyping/RData/joint.cross.cnt.RData')
cross.count.lookup=stack(joint.cross.cnt)

# estimate joint QTL effects within each cross 
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

# Supplementary Table 3 ------------------------------------------------------------------------------------------
suptable3a=data.frame(r2MeanEffect)
suptable3b=data.frame(r2)
suptable3c=data.frame(cross.peaks.flat)
WriteXLS(c('suptable3a', 'suptable3b', 'suptable3c'), 
         "/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryTable3.xls", SheetNames=c('joint model (avg)',
                                                                                            'joint model (ind)',
                                                                                            'within-cross model' 
                                                                                            ))
#-----------------------------------------------------------------------------------------------------------------



############ some specific stats mentioned in the text -------------------------###################################

# number of QTL with effect less than 0.1 SD
sum(abs(r2MeanEffect$MeanAbsBeta)<.1)
# 2911 / 4552  = 64%

# number of QTL with lead variant common (MAF>1%)
sum(r2MeanEffect$maf1012>.01)
# 3568 / 4552  = 78%

# slice at relatively large effect (sd>.3)
fisher.test(t(table(abs(r2MeanEffect$MeanAbsBeta)>.3, r2MeanEffect$maf1012>.01)))
#          <.3   >.3
#  rare   839    145
#  common 3478    90
#
# odds ratio = 1/.1498 = 6.676


# median effect size
#median(abs(r2MeanEffect$MeanAbsBeta))
#0.08071

# slicing table of effect size and allele frequency -----------------------------------
# significant enrichment for variants with maf <1% having large effects
table(abs(r2MeanEffect$MeanAbsBeta)<.1,  r2MeanEffect$maf1012>.05)
fisher.test(table(abs(r2MeanEffect$MeanAbsBeta)<.1,  r2MeanEffect$maf1012>.01))

# signficant enrichment of derived allele freq lessthan 1% decreasing fitness 
fisher.test(sapply(split(sign(r2MeanEffect$MeanAbsBeta), r2MeanEffect$maf2>.01), table))

t(table(abs(r2MeanEffect$MeanAbsBeta)>.3, r2MeanEffect$maf1012>.01))
# right column is effect greater than 0.3
# bottom row is common variants
fisher.test(t(table(abs(r2MeanEffect$MeanAbsBeta)>.3, r2MeanEffect$maf1012>.01)))
#	Fisher's Exact Test for Count Data

#data:  
#p-value <2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.1126 0.1985
#sample estimates:
#odds ratio 
#    0.1498 
# enrichment is 1/.1498 = 6.676
# ----------------------------------------------------------------------------------------

# slicing ancestral allele analysis ------------------------------------------------------
ancestral.eff=(r2MeanEffect$MeanAbsBeta)[r2MeanEffect$maf2>.95]
derived.eff=(r2MeanEffect$MeanAbsBeta)[r2MeanEffect$maf2<.05]

# compare derived vs ancestral for effects greater than 0.1 SD
fisher.test(rbind(table(abs(derived.eff)>.1),
      table(abs(ancestral.eff)>.1)))
#data:  
#p-value = 9e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.4201 0.7590
#sample estimates:
#odds ratio 
#    0.5658 
# 1/0.5658 =  1.767

# compare sign of ancestral vs derived 
fisher.test(rbind(table(sign(derived.eff)),
      table(sign(ancestral.eff))))

#data:  rbind(table(sign(derived.eff)), table(sign(ancestral.eff)))
#p-value = 0.008
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.095 1.955
#sample estimates:
#odds ratio 
#     1.462 

#-----------------------------------------------------------------------------------------

#####################################################################################################################



# Figure 3 -------------------------------------------------------------------------------------------------------
# needs R 3.5 for ggplot stat_summary_bin breaks parameter to function properly
#load('E:/Dropbox/RR/PreviousVersions/testj.RData')

# from variance_components_by_AF.R
load('/data/rrv2/genotyping/RData/testj.RData')
# r2MeanEffect is defined above
#load('E:/Dropbox/RR/PreviousVersions/r2MeanEffect.RData')
library(ggpubr)
library(Hmisc)
library(ggplot2)

Fig3A=ggplot(testj, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    ylab(expression(Fraction~of~heritability~(h^2)))+
    scale_fill_manual(guide=F, values =c('lightgrey', 'lightblue'), labels=c( '> 0.01', '< 0.01')) +
    #guides(fill=guide_legend("minor allele frequency\n in 1,011 panel"))+ 
    scale_y_continuous(expand=c(0,0))+
    theme_bw()+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))+
    theme(axis.text.x=element_text(size=rel(1),color='#000000CC', angle=70,hjust=1),
          axis.text.y=element_text(size=rel(1),color='black'))
Fig3Blank=ggplot(r2MeanEffect)+geom_blank()+theme_classic()

Fig3B=ggplot(r2MeanEffect)+geom_point(alpha=.4, size=.75, aes(x=maf1012,y=abs(MeanAbsBeta)))+ #,color=density))+
    #scale_color_viridis(option = "inferno", direction=1,end=1)+
    scale_x_continuous(name='Minor allele frequency', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='Average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()+
    theme(panel.grid.major = element_line(colour = "#80808022"))+
    theme(panel.grid.minor = element_line(colour = "#80808022"))+
    theme(axis.text.x=element_text(size=rel(1),color='black'),
          axis.text.y=element_text(size=rel(1),color='black'))+
    stat_summary_bin(aes(x=maf1012,y=abs(MeanAbsBeta)), breaks=cut2(r2MeanEffect$maf1012,g=42, onlycuts=T),  col='red')

Fig3C=ggplot(r2MeanEffect)+geom_point(alpha=.4, size=.75, aes(x=maf2,y=MeanAbsBeta))+
     #scale_color_viridis()+
    scale_x_continuous(name='Derived allele frequency', breaks=seq(0,1,.1), expand=c(0.01,0)) +
    scale_y_continuous(name='Average effect, SD units', limits=c(-1,1), breaks=seq(-1,1,.1), expand=c(0,0))+theme_bw()
ggarrange(Fig3A,Fig3Blank,Fig3B, Fig3C, ncol=2, nrow=2, labels=c('A','','B', 'C'))
ggsave(file='/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/Figure3.pdf',width=11,height=11)
#------------------------------------------------------------------------------------------------------------------

#Supplementary Figure 5 ----------------------------
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
#------------------------------------------------------------------------------------


# Supplementary Figure 6  --------------------------------------------------------------------------------
sr2=r2MeanEffect
#sr2$ve=4*r2MeanEffect$maf1012*(1-r2MeanEffect$maf1012)*r2MeanEffect$MeanAbsBeta^2
#sr2$ve=abs(r2MeanEffect$MeanAbsBeta)
#sr2$ve=4*(r2MeanEffect$crossCount/32)*(1-(r2MeanEffect$crossCount/32))*r2MeanEffect$MeanAbsBeta^2
sr2$ve=2*(r2MeanEffect$crossCount/32)*(1-(r2MeanEffect$crossCount/32))*r2MeanEffect$MeanAbsBeta^2

sr2=split(sr2, sr2$trait)
names(sr2)=as.character(levels(WCV$Trait))
pdf('/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure6.pdf', width=15, height=20) 
par(mfrow=c(7,6))
for(i in 1:38) {
sr21=sr2[[i]]
sr21=sr21[order(sr21$maf1012, decreasing=F),]
sr21$cve=cumsum(sr21$ve)/sum(sr21$ve)
plot(sr21$maf1012, sr21$cve, main=names(sr2[i]), ylab='cummulative GVE', 
     xlab='MAF', xlim=c(0,.5), ylim=c(0,1), type='l', lwd=2)
#plot(density(log10(sr21$maf1012),sr21$ve, n=3))
abline(0,2, lty=2, col='grey')
#readline()
}
dev.off()
# --------------------------------------------------------------------------------------------------------



# compare QTL effect sizes for within-cross vs joint QTL mapping ----------------------------------
# specifically, find the QTL from the joint analysis that don't overlap the within-cross analysis 
# and investigate their effects 
# for each cross
jpel=list()
for(cc in names(cross.peaks)){
    # for each trait 
    print(cc)
     for(tt in unique(cross.peaks[[1]]$trait)[-c(37,38)]) {
         print(tt)
        cpt=cross.peaks[[cc]][cross.peaks[[cc]]$trait==tt,]
        jpt=jointPeakEffects[[tt]][[cc]]     
        if(nrow(cpt)==0) {
             jpt$overlapQTL=F
             #jpt$overlapQTL[overlapQ]=T
            jpel[[tt]][[cc]]=jpt
        
        } else {
        ppos1=tstrsplit(cpt$CI.l, '_', type.convert=T)[[2]]
        ppos2=tstrsplit(cpt$CI.r, '_', type.convert=T)[[2]]
        chr1=tstrsplit(cpt$CI.l, '_', type.convert=T)[[1]]
        chr2=tstrsplit(jpt$peaks ,'_', type.convert=T)[[1]]
        ppos3=tstrsplit(jpt$peaks, '_', type.convert=T)[[2]]



        cpti=makeGRangesFromDataFrame(data.frame(chr=chr1, start=ppos1, end=ppos2 ), ignore.strand=T)
        jpti=makeGRangesFromDataFrame(data.frame(chr=chr2, start=ppos3, end=ppos3 ), ignore.strand=T)
        overlapQ=unique(queryHits(findOverlaps(jpti, cpti)))
        jpt$overlapQTL=F
        jpt$overlapQTL[overlapQ]=T
        jpel[[tt]][[cc]]=jpt
        }
     }
}
rjpel=rbindlist(lapply(jpel, rbindlist))
rjpel=rjpel[!(rjpel$trait %in% c("YPD;;2","YPD;;3")),]


rjpel$crossCount=ifelse(rjpel$cross.cnt%%2==1, rjpel$cross.cnt+1, rjpel$cross.cnt)
#r2$density=get_density(r2$maf1012, abs(r2$betas))
#r2$densityF=get_density(r2$maf1012Fill, abs(r2$betas))
rjpel$absBeta=abs(rjpel$betas)
rjpel$ancestral=iseq.freqs$ancestral[match(rjpel$peaks, iseq.freqs$marker)]
rjpel$maf2=ifelse(rjpel$ancestral, rjpel$alt012, 1-rjpel$alt012)

# 0 = reference derived, alternate paradoxus (ancestral)
# 1 = reference paradoxus (ancestral), alternate cerevisiae

#r2MeanEffect=r2 %>% group_by_(.dots=c("trait","peaks")) %>% mutate(MeanAbsBeta=mean(absBeta)) %>% distinct(traits,peaks,maf1012,crossCount, MeanAbsBeta, ancestral, maf2)
rjMeanEffect=rjpel %>% group_by_(.dots=c("trait","peaks")) %>% mutate(MeanAbsBeta=mean(betas)) %>% distinct(trait,peaks,maf1012,crossCount, MeanAbsBeta, ancestral, maf2)
rjq=split(rjpel, paste(rjpel$trait,rjpel$peaks))
sapply(rjq, function(x) sum(x$overlapQTL)>0)

rjq[names(which(sapply(rjq, function(x) sum(x$overlapQTL)==0)))]
wilcox.test(sapply(rjq[names(which(sapply(rjq, function(x) sum(x$overlapQTL)>0)))], function(x) mean(abs(x$betas))),
            sapply(rjq[names(which(sapply(rjq, function(x) sum(x$overlapQTL)==0)))], function(x) mean(abs(x$betas))))
#W = 1e+06, p-value = 9e-05
#median=.08259 vs median 0.0711 
#------------------------------------------------------------------------------------------------------

