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

# output from phenotyping/code/process_images.R 
load('/data/rr/Phenotyping/NORMpheno.RData')

# output from genotyping/code/segregants_hmm.R
load('/data/rrv2/genotyping/RData/cross.list.RData')

# output from genotyping/code/segregants_hmm.R
load('/data/rrv2/genotyping/RData/parents.list.RData')

# additional helper functions
source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')
source('/data/rrv2/analysis/mapping_fx.R')

# For Figure 1 see -------------------------------------
#
# BuildTree.R 
#
#--------------------------------------------------------

## initialize lists

# cross.peaks=list()
# pheno.resids=list()
# seg.recoded=list()
# cross.validated.r2=list()
# h2.vc.model=list()
# pheno.extracted=list()
# VC_A.list=list()

### --- within-cross QTL mapping analysis   --------------------------------------------------------------------------------------
# for each cross
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

    # extract phenotype for genotyped strains 
    subPheno=lapply(NORMpheno, function(x) x[match(snames, names(x))])
    # 4NQO phenotyping failed for technical reasons, remove it
    # calculate the average phenotype value for each segregant
    mPheno  =sapply(subPheno[-1], function(x) sapply(x, mean, na.rm=T))

    # there are typically only a handful of missing values, if any. If there are missing values replace with mean
    # 95% of traits and crosses have 0 missing the data, the rest typically have only one or two missing data points,
    # four specific trait cross combinations have more than 10 missing phenotypes
    mPheno=apply(mPheno,2, function(x) {x[is.na(x)]=mean(x, na.rm=T); return(x)})
    pheno_extracted[[cross.name]]=mPheno
    
    # scale phenotypes
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
# remove the results from the two additional YPD replicates (very similar to rep 1) 
sum(sapply(cross.peaks, function(x) nrow(x[x$q<.05 & !(x$trait %in% c("YPD;;2","YPD;;3")) ,])))
sapply(cross.peaks, function(x) (y=x[x$q<.05 & !(x$trait %in% c("YPD;;2","YPD;;3")) ,]))
qtl.per.cross.per.trait=unlist(sapply(cross.peaks, function(x) {
                                 y=x[x$q<.05 & !(x$trait %in% c("YPD;;2","YPD;;3")) ,];
                                 sapply(split(y, y$trait), nrow);
                                 }))
median(qtl.per.cross.per.trait[grepl('^A', names(qtl.per.cross.per.trait))])
median(qtl.per.cross.per.trait[!grepl('^A', names(qtl.per.cross.per.trait))])
### --------------------------------------- end within-cross QTL mapping analysis ------------------------------------------------



### ---Joint QTL mapping analysis -------------------------------------------------------------------------------------------------
#
#pre-processing
parents.list=lapply(parents.list, function(x) {
                  z=x;
                  z$marker.name.n=paste0(z$marker.name, '_', seq(1:nrow(z)))
                  return(z) })

# load JS 1,011 isolates allele frequency data into a structure
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
#save(jointPeaks5, file='/data/rrv2/genotyping/RData/jointPeaks5.RData')

# do multi-cross analysis using 1,011 panel variant data only
jointPeaksJS=mapJointQTLsJS_variants(n.perm=1e3, FDR_thresh=.05, parents.list, pheno.resids, seg.recoded, iseq.freqs, filterJS=T)
#save(jointPeaksJS, file='/data/rrv2/genotyping/RData/jointPeaksJS.RData')

### do multi-cross analysis with cross-validation to estimate qtl r^2 -------------------------------------------------------------
#
# jointModelCrossValidation.R
#
#----------------------------------------------------------------------------------------------------------------------------------
#
### ----end Joint QTL mapping analysis --------------------------------------------------------------------------------------------


### do within-cross variance component analysis -----------------------------------------------------------------------------------
#
# see variance_components_within_cross.R
#
#----------------------------------------------------------------------------------------------------------------------------------


### do joint variance component analysis (across the whole panel split by allele frequencies in 1,011 isolate collection ----------
#
# see variance_components_by_AF.R
#
#----------------------------------------------------------------------------------------------------------------------------------


### do causal gene identification and GO analysis ---------------------------------------------------------------------------------
# 
# see QTL_causality.R
# 
#----------------------------------------------------------------------------------------------------------------------------------


### do simulation analysis to investigate how mapping procedure captures allele-frequency to effect size coupling -----------------
#
# see simulateArchitecture.R
#
#----------------------------------------------------------------------------------------------------------------------------------


### investigate relationship between allele-frequency and effect size from joint mapping procedure  -------------------------------
#
# see jointModelAF.R
#
#----------------------------------------------------------------------------------------------------------------------------------


### Supplementary Figure 2 plots (QTL architecture with arrow plots)  -------------------------------------------------------------
#
# see additional_within_cross_models_and_plots.R
#
#----------------------------------------------------------------------------------------------------------------------------------
