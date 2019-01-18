#NORMpheno
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

load('/data/rr/Phenotyping/NORMpheno.RData')
load('/data/rrv2/genotyping/RData/cross.list.RData')
load('/data/rrv2/genotyping/RData/parents.list.RData')
source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')
source('/data/rrv2/analysis/mapping_fx.R')

#cross.peaks=list()
#pheno.resids=list()
#seg.recoded=list()
#cross.validated.r2=list()
#h2.vc.model=list()
#pheno.extracted=list()
#VC_A.list=list()

for(cross.name in crosses) {
    print(cross.name)
    cross=cross.list[[cross.name]]
    if(cross.name=='A') {       cross=subset(cross, ind=!grepl('A11', as.character(cross$pheno$id)))    }
    snames = as.character(cross$pheno$id)
    g=pull.argmaxgeno(cross)
 
    # recode based on parental genotypes 
    #seg.pcode=recode.as.allele.number(t(g),parents.list[[cross.name]])
    #seg.recoded[[cross.name]]=seg.pcode
   
    # are there fixed loci ?? (no)-------------------------------
    g.af=apply(g,2,function(x) sum(x==1))
    parents.list[[cross.name]]$fixed=(g.af==0 | g.af==nrow(g))
    fixed.loci=which(parents.list[[cross.name]]$fixed)
    if(length(fixed.loci)>0) {    g=g[,-fixed.loci] }
    #------------------------------------------------------------
    
    # reduce the set of markers by pruning markers in complete LD
    g.r=g[,-which(duplicated(g, MARGIN=2))]
    # scale 
    g.s=scale(g.r)
    gall.s=scale(g)

    # extract phenotype for genotyped strains 
    subPheno=lapply(NORMpheno, function(x) x[match(snames, names(x))])
    # 4NQO phenotyping failed, remove it
    mPheno  =sapply(subPheno[-1], function(x) sapply(x, mean, na.rm=T))

    mPheno=apply(mPheno,2, function(x) {x[is.na(x)]=mean(x, na.rm=T); return(x)})
    pheno_extracted[[cross.name]]=mPheno
    
    mPhenos=scale(mPheno)
    
    #calculate whole genome additive only model (fast)
    #VC_A.list[[cross.name]]=do_VC_additive_only_average(mPheno, g.s)

    # calculate whole genome additive only model with standard errors (slower)
    aVC=doAdditiveVC(mPhenos, g.s)
    h2.vc.model[[cross.name]]=aVC

    # cross-validated R^2 for additive QTL model
    cvr2=doTraitCV(mPhenos, g.s, gall.s)
    cross.validated.r2[[cross.name]]=cvr2
    
    # QTL mapping, within-cross only
    cQTL=apply(mPhenos, 2, function(x) {set.seed(100); return(doTraitFDR(x, g.s, gall.s, FDR_thresh=.2, nperm=1e4)) })
    ppc=rbindlist(cQTL, idcol='trait') #(do.call('rbind', cQTL)) #lapply(cQTL, function(x) do.call('rbind', x))))
    cross.peaks[[cross.name]]=ppc
    
     # get residuals per trait for model with fixed effects of significant qtl on other chromosomes and random effect on other chromosomes 
    chromResids=calcChromosomeResiduals(ppc[ppc$q<.05,], mPheno, gall.s, g.s)
    pheno.resids[[cross.name]]=chromResids
}
#save(pheno_extracted, file='/data/rrv2/genotyping/RData/extracted_average_phenotypes.RData')
#save(h2.vc.model, file= '/data/rrv2/genotyping/RData/h2_VC_model.RData')
#save(cross.validated.r2, file = '/data/rrv2/genotyping/RData/FDR_cross_validatedR2.RData')
#save(cross.peaks, file= '/data/rrv2/genotyping/RData/FDR_cross.peaks.RData')
#save(pheno.resids, file = '/data/rrv2/genotyping/RData/FDR_pheno.resids.RData')
#save(seg.recoded, file = '/data/rrv2/genotyping/RData/FDR_seg.recoded.RData')

load('/data/rrv2/genotyping/RData/h2_VC_model.RData')
load('/data/rrv2/genotyping/RData/FDR_cross_validatedR2.RData')
load('/data/rrv2/genotyping/RData/FDR_cross.peaks.RData')
load('/data/rrv2/genotyping/RData/FDR_pheno.resids.RData')
load('/data/rrv2/genotyping/RData/FDR_seg.recoded.RData')
load('/data/rrv2/genotyping/RData/extracted_average_phenotypes.RData')

# QTL per cross
# remove the results from the two additional YPD replicates 
sum(sapply(cross.peaks, function(x) nrow(x[x$q<.05 & !(x$trait %in% c("YPD;;2","YPD;;3")) ,])))
sapply(cross.peaks, function(x) (y=x[x$q<.05 & !(x$trait %in% c("YPD;;2","YPD;;3")) ,]))
qtl.per.cross.per.trait=unlist(sapply(cross.peaks, function(x) {
                                 y=x[x$q<.05 & !(x$trait %in% c("YPD;;2","YPD;;3")) ,];
                                 sapply(split(y, y$trait), nrow);
                                 }))

qcvr=(sapply(cross.validated.r2, function(x) rowMeans(x, na.rm=T)))[-c(37,38),]
h2s=(sapply(h2.vc.model, function(x) x[,1]))[-c(37,38),]
test=qvcr/h2s
test2=gather(data.frame(test), cross)
test2=data.frame(condition=rep(rownames(test), 16), test2)
#median(qcvr/h2s, na.rm=T)
#[1] 0.6834
# variance explained BY vs RM
plot(h2s, qcvr, xlim=c(0,1), ylim=c(0,1))
abline(0,1)
points(h2s[,2], qcvr[,2], col='red')
#hist(qcvr/h2s, breaks=100)


parents.list=lapply(parents.list, function(x) {
                  z=x;
                  z$marker.name.n=paste0(z$marker.name, '_', seq(1:nrow(z)))
                  return(z) })
#save(jointPeaks5, file='/data/rrv2/genotyping/RData/jointPeaks5.RData')


# JS panel allele frequency data
js.rr.overlap.allele.frequencies='/data/rrv2/1002genomes/isec_ouput/out.frq.mod'
sacCer3_CBS432_alignment.variants='/data/rrv2/spar_alignment/filt.snps'
sacCer3_CBS432_alignment.coords='/data/rrv2/spar_alignment/out.mcoords'
iseq.freqs=buildJS_variants_annotation_table(js.rr.overlap.allele.frequencies,sacCer3_CBS432_alignment.variants,sacCer3_CBS432_alignment.coords)
#save(iseq.freqs, file='/data/rrv2/genotyping/RData/iseq.freqs.RData')
load('/data/rrv2/genotyping/RData/iseq.freqs.RData')

# do multi-cross analysis using all called variants 
#jointPeaks5=mapJointQTLs(n.perm=1000, FDR.thresh=.05, parents.list, pheno.resids, seg.recoded) 
#load('/data/rrv2/genotyping/RData/jointPeaks5.RData')
#jP=rbindlist(jointPeaks5, idcol='chromosome')
#jPs=split(jP, jP$trait)

# do multi-cross analysis using JS markers only
#jointPeaksJS=mapJointQTLsJS_variants(n.perm=1e3, FDR_thresh=.05, parents.list, pheno.resids, seg.recoded, iseq.freqs, filterJS=T)
#save(jointPeaksJS, file='/data/rrv2/genotyping/RData/jointPeaksJS.RData')
load('/data/rrv2/genotyping/RData/jointPeaksJS.RData')
jP=rbindlist(jointPeaksJS, idcol='chromosome')
jPs=split(jP, jP$trait)
#sum(sapply(jPs[-c(37,38)], nrow))

# further functionalize this code for publication
#training.sets
set.seed(20)
strain.names=sapply(pheno_extracted, function(x) rownames(x))
strain.names.cv=lapply(strain.names, function(x) {
       cvg=cut(sample(1:length(x)),10)
       levels(cvg)=as.roman(1:10)
       cvg=as.character(cvg)
       return(cvg) })

library(openblasctl)
openblas_set_num_threads(48)
save.prefix='/data/rrv2/genotyping/RData/jointJSCV_'
#save.prefix='/data/home/jbloom/misc/jointJSCV_'
setme=as.character(as.roman(1:10))
for(setm in setme) {
    print(setm)
    prs=mapply(function(x,y){
                  lapply(x,function(z) z[which(!(y%in%setm)),-c(37,38)])
                  },x=pheno.resids,y=strain.names.cv, SIMPLIFY=F)
    srs=mapply(function(x,y){
                  x[,which(!(y%in%setm))]
                  },x=seg.recoded,y=strain.names.cv, SIMPLIFY=F)
    print(sapply(srs, ncol))
    temp=mapJointQTLsJS_variants(n.perm=3e2, FDR_thresh=.05, parents.list, prs, srs, iseq.freqs, filterJS=T)
    mjq=rbindlist(temp,idcol='chromosome')
    saveRDS(mjq, file=paste0(save.prefix,setm,'.RDS'))
    # get coefficients from training data, predict effects on data left out
}
cvjPs=list()
for(setm in setme) {
    print(setm)
    mjq=readRDS(paste0(save.prefix, setm, '.RDS'))
    cvjPs[[setm]]=split(mjq, mjq$trait)
}

# for each cross-validation set
jointPeaksCV_R2=list()
for(setm in setme) {
    print(setm)
    tcross=matrix(NA, length(names(cvjPs[[setm]])), length(crosses))
    rownames(tcross)=names(cvjPs[[setm]])
    colnames(tcross)=crosses
    for(cross.name in crosses) {
         print(cross.name)
         g.s=scale(t(seg.recoded[[cross.name]]))
         #rename columns
         colnames(g.s)=parents.list[[cross.name]]$marker.name
         mPheno=pheno_extracted[[cross.name]]
         mPhenos  = scale(mPheno)

          for(tt in names(cvjPs[[setm]])) {
             sig.joint.markers=parents.list[[cross.name]]$marker.name[ na.omit(match(cvjPs[[setm]][[tt]]$fscan.markers[cvjPs[[setm]][[tt]]$q<.05], parents.list[[cross.name]]$marker.name))]
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
             X2=dff
             yr=mPhenos[,tt]
             yr.train=yr
             # double check this
             yr.train[strain.names.cv[[cross.name]]==setm]=NA
             yr.test=yr
             yr.test[strain.names.cv[[cross.name]]!=setm]=NA    
             
             if(ncol(X2)==0) { next;}
             # fit model on training data
             fitme=(lm(yr.train~.-1,X2))
             if(is.null(dim(X2))){
                predicted=data.matrix(X2[strain.names.cv[[cross.name]]==setm])*coef(fitme)
             }else {
                 # estimate variance explained  
                predicted=data.matrix(X2[strain.names.cv[[cross.name]]==setm,])%*%coef(fitme)
            }
            tcross[tt,cross.name]=cor(yr.test[strain.names.cv[[cross.name]]==setm], predicted)^2
        }
    }
    jointPeaksCV_R2[[setm]]=tcross
}
#save(jointPeaksCV_R2, file='/data/rrv2/genotyping/RData/jointJSCVR2.RData')
aa=do.call('abind', c(jointPeaksCV_R2, along=3))
aam=apply(aa,c(1,2), mean)

cF=as.factor(rep(names(pheno_extracted), sapply(pheno_extracted, nrow)))
pe=do.call('rbind', (pheno_extracted))
cross.vars=sapply(split(data.frame(pe[,-c(37,38)]), cF), function(x) apply(x,2,var))
cross.cnts.table=table(cF)
aam=aam[,colnames(cross.vars)]
median(aam/h2s[,colnames(aam)])

plot(as.vector(h2s[,colnames(aam)]), as.vector(aam), ylim=c(0,1), xlim=c(0,1),xlab='additive heritability (whole genome)', ylab='CV QTL model variance explained')
abline(0,1)

wsum=rowSums(cross.vars*(as.vector(cross.cnts.table)))
wvar=rowSums(aam*cross.vars*(as.vector(cross.cnts.table)))/rowSums(cross.vars*(as.vector(cross.cnts.table)))
hvar=rowSums(h2s[,colnames(aam)]*cross.vars*(as.vector(cross.cnts.table)))/rowSums(cross.vars*(as.vector(cross.cnts.table)))


aweighted=rep(0,ncol(aam))
for(i in 1:nrow(aam)) {
    for(j in 1:ncol(aam)){
      aweighted[i]= 
      anorm3[i]=sum(sw[[i]]$additive*cross.cnts.table[sw[[i]]$cross]*cross.vars[i,sw[[i]]$cross])/sum(cross.cnts.table[sw[[i]]$cross]*cross.vars[i,sw[[i]]$cross])
    
}


# run code for variance component analysis












