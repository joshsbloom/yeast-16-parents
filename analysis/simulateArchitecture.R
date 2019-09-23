library(caret)
library(regress)
library(qtl)
library(rrBLUP)
library(Matrix)
library(VariantAnnotation)
library(vcfR)
library(qtl)
library(qtl2)
library(ASMap)
library(data.table)
library(S4Vectors)
library(dplyr)
library(reshape2)
library(tidyr)
library(Rfast)
library(doMC)
library(Hmisc)
library(sisus)
library(qvalue)
library(ggplot2)
library(ggpubr)

source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')
source('/data/rrv2/analysis/mapping_fx.R')

#pheno_extracted
load('/data/rrv2/genotyping/RData/extracted_average_phenotypes.RData')

# output from genotyping/code/segregants_hmm.R
load('/data/rrv2/genotyping/RData/parents.list.RData')
parents.list=lapply(parents.list, function(x) {
                  z=x;
                  z$marker.name.n=paste0(z$marker.name, '_', seq(1:nrow(z)))
                  return(z) })

# allele frequency lookup from 1,011 yeast panel
load('/data/rrv2/genotyping/RData/iseq.freqs.RData')

# output from genotyping/code/segregants_hmm.R
load('/data/rrv2/genotyping/RData/cross.list.RData')

# push effects onto seg.recoded (coded as reference or not)-----------------------
set.seed(1000)
ifs.set=replicate(20, {
    ifs=iseq.freqs[sort(sample(1:nrow(iseq.freqs), 200)),]
    lt=ifs$maf1012<.01
    gt=ifs$maf1012>=.01
    elt=rnorm(sum(lt), 0, .3) #.3
    egt=rnorm(sum(gt), 0, .05)
    hist(abs(c(elt,egt)))
    ifs$eff=0
    ifs$eff[lt]=elt
    ifs$eff[gt]=egt
    return(ifs)
   }, simplify=F)
#plot(ifs$maf1012, abs(ifs$eff))
#save(ifs.set, file='/data/rrv2/genotyping/RData/ifs.set.RData')
#save(ifs.set, file='/data/rrv2/genotyping/RData/ifs.set.RData2')
#save(ifs.set, file='/data/rrv2/genotyping/RData/ifs.set.RData3')

h2=.6
sim.phenos=list()
for(cross.name in crosses) {
    ys=list()
    for(s in 1:20) {
        mm= match(ifs.set[[s]]$marker, parents.list[[cross.name]]$marker.name)
        if2=ifs.set[[s]][!is.na(mm),]
        if2$marker.name.n=parents.list[[cross.name]]$marker.name.n[mm[!is.na(mm)]]
        yg=sqrt(h2)*as.vector(if2$eff %*% seg.recoded[[cross.name]][if2$marker.name.n,])
        y=yg+rnorm(length(yg),mean=0,sd=sqrt((1-h2)/h2*var(yg)))
        #plot( (coef(lm(y~t(seg.recoded[[cross.name]][if2$marker.name.n,]))))[-1], if2$eff)
        #abline(0,1)
        ys[[s]]=y
     }
    ys=do.call('cbind', ys)
    colnames(ys)=as.roman(1:20)
    rownames(ys)=rownames(pheno_extracted[[cross.name]])
    sim.phenos[[cross.name]]=ys

}
#save(sim.phenos, file='/data/rrv2/genotyping/RData/sim.phenos.RData')
#save(sim.phenos, file='/data/rrv2/genotyping/RData/sim.phenos.RData2')
#save(sim.phenos, file='/data/rrv2/genotyping/RData/sim.phenos.RData3')

#----------------------------------------------------------------------------
cross.peaks.sim=list()
pheno.resids.sim=list()
for(cross.name in crosses) {
    print(cross.name)
    cross=cross.list[[cross.name]]
    # subset the BYxRM to ~10 96 well plates worth of segregants so statistical power is similar to other crosses 
    if(cross.name=='A') {       cross=subset(cross, ind=!grepl('A11', as.character(cross$pheno$id)))    }
    snames = as.character(cross$pheno$id)
    # get imputed genotypes from R/QTL
    g=pull.argmaxgeno(cross)
 
    # for joint analysis, recode alleles as reference or not (for within cross linkage mapping coding is one parent or the other)
    # recode based on parental genotypes 
    #seg.pcode=recode.as.allele.number(t(g),parents.list[[cross.name]])
    #seg.recoded[[cross.name]]=seg.pcode
   
    # are there fixed loci ?? (no)------------------------------
    g.af=apply(g,2,function(x) sum(x==1))
    parents.list[[cross.name]]$fixed=(g.af==0 | g.af==nrow(g))
    fixed.loci=which(parents.list[[cross.name]]$fixed)
    if(length(fixed.loci)>0) {    g=g[,-fixed.loci] }
    #------------------------------------------------------------
    
    # reduce the set of markers by pruning markers in complete LD
    g.r=g[,-which(duplicated(g, MARGIN=2))]
    # scale genotypes 
    g.s=scale(g.r)
    gall.s=scale(g)
    mPhenos=scale(sim.phenos[[cross.name]])

    cQTL=apply(mPhenos, 2, function(x) {set.seed(100); return(doTraitFDR(x, g.s, gall.s, FDR_thresh=.05, nperm=50, doLODdrop=F)) })
    ppc=rbindlist(cQTL, idcol='trait') #(do.call('rbind', cQTL)) #lapply(cQTL, function(x) do.call('rbind', x))))
    ppc$pmarker=ppc$fscan.marker
    cross.peaks.sim[[cross.name]]=ppc
    
     # get residuals per trait for model with fixed effects of significant qtl on other chromosomes and random effect on other chromosomes 
    chromResids=calcChromosomeResiduals(ppc[ppc$q<.05,], sim.phenos[[cross.name]], gall.s, g.s)
    pheno.resids.sim[[cross.name]]=chromResids
}
#save(cross.peaks.sim, file='/data/rrv2/genotyping/RData/cross.peaks.sim.RData')
#save(pheno.resids.sim, file='/data/rrv2/genotyping/RData/pheno.resids.sim.RData')

#save(cross.peaks.sim, file='/data/rrv2/genotyping/RData/cross.peaks.sim.RData2')
#save(pheno.resids.sim, file='/data/rrv2/genotyping/RData/pheno.resids.sim.RData2')

save(cross.peaks.sim, file='/data/rrv2/genotyping/RData/cross.peaks.sim.RData3')
save(pheno.resids.sim, file='/data/rrv2/genotyping/RData/pheno.resids.sim.RData3')

jointPeaksSim=mapJointQTLsJS_variants(n.perm=1e2, FDR_thresh=.05, parents.list, pheno.resids.sim, seg.recoded, iseq.freqs, filterJS=T)
#save(jointPeaksSim, file='/data/rrv2/genotyping/RData/jointPeaksSim.RData')
#save(jointPeaksSim, file='/data/rrv2/genotyping/RData/jointPeaksSim.RData2')
save(jointPeaksSim, file='/data/rrv2/genotyping/RData/jointPeaksSim.RData3')

jPsS=rbindlist(jointPeaksSim, idcol='chromosome')
jPsS=split(jPsS, jPsS$trait)

load('/data/rrv2/genotyping/RData/joint.cross.cnt.RData')
cross.count.lookup=stack(joint.cross.cnt)

jointPeakEffectsSim=list()
for(cross.name in crosses) {
    print(cross.name)
    g.s=scale(t(seg.recoded[[cross.name]]))
    #rename columns
    colnames(g.s)=parents.list[[cross.name]]$marker.name
    mPheno=sim.phenos[[cross.name]]
    mPhenos  = scale(mPheno)

    for(tt in names(jPsS)) {
         sig.joint.markers=parents.list[[cross.name]]$marker.name[na.omit(match(jPsS[[tt]]$fscan.markers[jPsS[[tt]]$q<.05], parents.list[[cross.name]]$marker.name))]
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
         jointPeakEffectsSim[[tt]][[cross.name]]=data.frame(
            trait=tt,
            cross=cross.name,
            peaks=sig.joint.markers,
            betas=as.vector(coef(qmodel)),
            vexp=a.effs,
            p = ps,
            maf1012=iseq.freqs$maf1012[match(sig.joint.markers, iseq.freqs$marker)],
            # iseq.freqs[,8] is the allele frequency of the alternate allele 
            #alt012=iseq.freqs[,8][match(sig.joint.markers, iseq.freqs$marker)],
            cross.cnt=cross.count.lookup[match(sig.joint.markers, rownames(cross.count.lookup)),'values'],
            stringsAsFactors=F)
          print(jointPeakEffectsSim[[tt]][[cross.name]])

         }
}
r2Sim=rbindlist(lapply(jointPeakEffectsSim, rbindlist))
#save(r2Sim, file ='/data/rrv2/genotyping/RData/r2Sim.RData')
#save(r2Sim, file ='/data/rrv2/genotyping/RData/r2Sim.RData2')
#save(r2Sim, file ='/data/rrv2/genotyping/RData/r2Sim.RData3')


#ifs.set

pdf(file='/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/SupplementaryFigure4.pdf', width=8, height=8)
par(mfrow=c(2,2), yaxs='i')
load('/data/rrv2/genotyping/RData/ifs.set.RData3')
load('/data/rrv2/genotyping/RData/r2Sim.RData3')

load('/data/rrv2/genotyping/RData/ifs.set.RData2')
load('/data/rrv2/genotyping/RData/r2Sim.RData2')


r2MeanEffectSim=r2Sim %>% group_by_(.dots=c("trait","peaks")) %>% mutate(MeanAbsBeta=mean(betas)) %>% distinct(trait,peaks,maf1012, MeanAbsBeta)
ifsdf=rbindlist(ifs.set)
#png(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/2019/031919/simulations_rr_mixture2.png', width=1024, height=512)

plot(ifsdf$maf1012, abs(ifsdf$eff)*.6,
     main='Simulated architecture 2', 
     ylim=c(0,1),
     ylab='Average absolute effect, SD units', 
     xlab='Minor allele frequency', col='#00000033')
#abline(v=0.025, col='red')

plot(r2MeanEffectSim$maf1012, abs(r2MeanEffectSim$MeanAbsBeta), 
     main='Detected architecture for simulation 2', 
     ylab='Average absolute effect, SD units',
     xlab='Minor allele frequency', ylim=c(0,1), col='#00000033')

#abline(v=0.025, col='red')
dev.off()

#ifs.set
#load('/data/rrv2/genotyping/RData/ifs.set.RData2')
#r2Sim
#load('/data/rrv2/genotyping/RData/r2Sim.RData2')

