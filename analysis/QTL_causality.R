library(S4Vectors)
library(Rfast)
load('/data/rr/Phenotyping/NORMpheno.RData')
load('/data/rrv2/genotyping/RData/cross.list.RData')
load('/data/rrv2/genotyping/RData/parents.list.RData')
source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')
# cross peaks 
load( '/data/rrv2/genotyping/RData/FDR_cross.peaks.RData')

library(data.table)
#source('/data/rrv2/genotyping/code/mapping_fx.R')
library(tidyr)

# add marker.name.n column
parents.list=lapply(parents.list, function(x) {
                  z=x;
                  z$marker.name.n=paste0(z$marker.name, '_', seq(1:nrow(z)))
                  return(z) })

causalSets.by.cross2=list()
nsim=500
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
    #g.r=g[,-which(duplicated(g, MARGIN=2))]
    g.s=scale(g)
    #rename columns
    colnames(g.s)=parents.list[[cross.name]]$marker.name
    subPheno=lapply(NORMpheno, function(x) x[match(snames, names(x))])
    mPheno  =sapply(subPheno, function(x) sapply(x, mean, na.rm=T))
    mPheno=apply(mPheno,2, function(x) {x[is.na(x)]=mean(x, na.rm=T); return(x)})
    mPhenos=scale(mPheno)

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
    causalSets=list()
    #tt=names(cpeaks.by.trait)[21]
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
        #relocalized_peak=cpeaks.by.trait[[tt]]$marker,
            #nbeta=as.vector(coef(qmodel)),

            out=DataFrame(out)
            pCausal=list()
            for(peak in 1:nrow(out)){pCausal[[peak]]=NA }
            out$pCausal=pCausal
            # only bother with large effect peaks  
            for(peak in 1:nrow(out)) {
                print(peak)
                if(out[peak,]$vexp<0.02) { next;}
                chr=out[peak,]$chr #snMMC3[peak,]$chrom
                vacc=vac[[chr]]
                gsc=gs.by.chr[[chr]] #(g.by.chr[[chr]][fmodel$segs.to.keep,])
                peak.marker=out[peak,]$sig_covar
                peak.index=match(peak.marker, colnames(gsc))
                peak.index.genome=match(peak.marker, colnames(g.s))

                # relocalize peak (regress out effects of other peaks and find position of max stat)
                #nr=as.vector(residuals(lm(mPhenos[,tt]~g.s[,cpeaks.by.trait[[tt]]$marker[-peak]]-1)))
                #new peak position
                #out[peak,]$relocalized_peak=names(which.max(fasterLOD(length(nr),scale(nr), gsc)[1,]))
                #peak.marker=out[peak,]$relocalized_peak
                #peak.index=match(peak.marker, colnames(gsc))
                #peak.index.genome=match(peak.marker, colnames(g.s))
                #if(peak.index.genome %in% cpeaks.by.trait[[tt]]$marker[-peak]  ) {
                #    out[peak,]$relocalized_peak=out[peak,]$sig_covar
                #    out[peak,]$nbeta=out[peak,]$beta
                #} else {
                #    out[peak,]$nbeta=coef(lm(mPhenos[,tt]~g.s[,c(peak.marker, cpeaks.by.trait[[tt]]$marker[-peak])]-1))[1]
                # }
                pL=findInterval(vacc$pos[peak.index]-25000, vacc$pos)
                if(pL==0) { pL=1}
                pR=findInterval(vacc$pos[peak.index]+25000, vacc$pos)
                new.coords=pL:pR
                peak.index=match(peak.index, new.coords)
                gsc=gsc[,new.coords]

                #gscm=apply(gsc, 2, function(m) out[peak,]$nbeta *m)
                gscm=out[peak,]$betas*gsc
                ysimt=replicate(nsim, sample(yr), simplify=F)
                ysimt=lapply(ysimt, function(x) x+gscm)
                stm=do.call('rbind',ysimt)
                #stm=apply(gsc, 2, function(m) {
                #        replicate(nsim, { sample(yr)+out[peak,]$nbeta*m })
                #})

                # n individuals * nsim   X nrow(gsc)
                stmt=matrix(as.vector(stm),nrow(gsc))
                r2=(crossprod(scale(stmt), gsc)/(nrow(gsc)-1))^2
                #mL=apply(r2,1,which.max)
                mL=rowMaxs(r2, value=F)
                # rows are simulations, #columns are markers
                # data is the index of the marker that explains the max variance
                maxLODs=matrix(mL,nsim,ncol(gsc))

                #pAA prob that A is the peak given A is causal
                #pAB prob that A is the peak given B is causal
                pAN=(apply(maxLODs,2, function(x) sum(x==peak.index)/nsim))
                pCausal=pAN/sum(pAN)
                names(pCausal)=colnames(gsc)
                pCausal=pCausal[pCausal>0]
                
                #pCausal1=pCausal
                #pCausal2=pCausal
                #p2s=pCausal2[names(pCausal2) %in% names(pCausal1)]
                #p1s=pCausal1[names(pCausal1) %in% names(pCausal2)]
                #(p1s*p2s)/sum(p1s*p2s)
                
                #pAA=sum(maxLODs[,peak.index]==peak.index)/nsim
                #pAB=sum(apply(maxLODs[,-peak.index],2, function(x) sum(x==peak.index)/nsim))
                #pCausal=pAA/(pAA+pAB)

                bootInt=quantile(maxLODs[, peak.index],c(.025,.05,.95,.975))
                bootIntMarkers=colnames(gsc)[round(bootInt)] 
                out[peak,]$bootstrap_interval=paste(bootIntMarkers, collapse=' ')
                out[peak,]$pCausal[[1]]=pCausal
                #plot(as.numeric(sapply(strsplit(names(pCausal), '_'), function(x) x[2])), pCausal, main=paste(as.character(out[peak,'sig_covar'])))
                print(bootInt)
             }
            causalSets[[tt]]=out
        }
causalSets.by.cross2[[cross.name]]=causalSets
}


#save(causalSets.by.cross2, file = '/data/rrv2/genotyping/RData/causalSets2.RData')

#Start here
#load('/data/rrv2/genotyping/RData/causalSets2.RData')

#save(causalSets.by.cross, file = '/data/rrv2/genotyping/RData/causalSets.RData')
#load('/data/rrv2/genotyping/RData/causalSets.RData')



#pCausal1=causalSets.by.cross[['375']][[21]][9,]$pCausal
#pCausal2=causalSets.by.cross[['A']][[21]][7,]$pCausal

#p2s=pCausal2[names(pCausal2) %in% names(pCausal1)]
#p1s=pCausal1[names(pCausal1) %in% names(pCausal2)]
#plot( (p1s*p2s)  / sum((p1s*p2s)*((1-p1s)*p2s)*(p1s*(1-p2s))*((1-p1s)*(1-p2s)))  )
#(p1s*p2s)/sum(p1s*p2s) 
crosses.to.parents=list(
     '375'=c("M22", "BYa"),
     'A'  =c("BYa", "RMx"),
     '376'=c("RMx", "YPS163a"),
     'B'  =c("YPS163a", "YJM145x"),
     '377'=c("YJM145x", "CLIB413a"),
     '393'=c("CLIB413a", "YJM978x"),
     '381'=c("YJM978x", "YJM454a"),
    '3008'=c("YJM454a", "YPS1009x"),
    '2999'=c("YPS1009x", "I14a"),
    '3000'=c("I14a", "Y10x"),
    '3001'=c("Y10x", "PW5a"),
    '3049'=c("PW5a", "273614xa"),
    '3003'=c("273614xa", "YJM981x"),
    '3004'=c("YJM981x", "CBS2888a"),
    '3043'=c("CBS2888a", "CLIB219x"),
    '3028'=c("CLIB219x", "M22")
    )
parents=unique(unlist(crosses.to.parents))
parents.to.crosses=sapply(parents, function(x) sapply(crosses.to.parents, function(y) x%in%y))
ptc=apply(parents.to.crosses, 2, function(x) rownames(parents.to.crosses)[x] )

library(VariantAnnotation)
library(GenomicFeatures)
library(ensemblVEP)
library('org.Sc.sgd.db')

txdb=makeTxDbFromGFF(file=paste0('/data/CRISPR_variant_engineering/rr_variant_oligos/',
                                 'reference/saccharomyces_cerevisiae.gff'))
tscl=as.list(txdb)
tscm=merge(tscl$transcripts, tscl$genes, by='tx_id')
library(org.Sc.sgd.db)
xx <- as.list(org.Sc.sgdGENENAME) #org.Sc.sgdCOMMON2ORF)
flatlist=sapply(xx, function(x) x)
flatlist[is.na(flatlist)]=''
gene.GR=GRanges(seqnames=tscm$tx_chrom, ranges=IRanges(start=tscm$tx_start, end=tscm$tx_end), 
             strand=tscm$tx_strand, ORF=tscm$gene_id)
gene.GR$NAME=as.character(flatlist[gene.GR$ORF])

library(BSgenome.Scerevisiae.UCSC.sacCer3)
# - strand
#promoters(gene.GR[61], upstream=1000, downstream=width(gene.GR[61])+250)
# + strand 
#promoters(gene.GR[1000], upstream=1000, downstream=width(gene.GR[1000])+250)

# more principled but creates overlaps between genes and non-genes 
# call gene as 1kb upstream + 250 bp downstream
#ge.list=list()
#for(g in 1:length(gene.GR)) {
#    print(g) 
#    ge.list[[g]]=promoters(gene.GR[g], upstream=1000, downstream=width(gene.GR[g])+250)
#}
#genes.extended=do.call(getMethod(c, "GenomicRanges"), ge.list) # IRanges::unlist(ge.list)

gsplit=split(gene.GR, as.vector(seqnames(gene.GR)))
gsplit=lapply(gsplit, function(x) {
    x=x[order(start(x)),]
    e=end(x)
    s=start(x)
    esd=((start(x[-1])-end(x))/2)
    esd[esd<0]=0
    e=e+(floor(esd))-1
    s[-1]=s[-1]-ceiling(esd[-length(esd)])
    start(x)=s
    end(x)=e
    return(x)
})
geneExtend.GR2=stack(GRangesList(gsplit))


causalSets.by.cross=causalSets.by.cross2

# flatten, add columns for confidence intervals left and right  from LOD score 
# replace with confidence intervals for left and right given ppcs 
# add flag to note if PPCs were calculated (vexp>0.02)
#y=lapply(crossOut, function(x) do.call('rbind', x))
QTL.table=list()

for(cross.name in crosses) {
    print(cross.name)
    QTL.table[[cross.name]]=do.call('rbind', causalSets.by.cross[[cross.name]])
    QTL.table[[cross.name]]$cross=cross.name
}
QTL.table=do.call('rbind', QTL.table)
CI=makeGRangesFromDataFrame(data.frame(start=tstrsplit(QTL.table$CI.l, '_', type.convert=T)[[2]],
           end=tstrsplit(QTL.table$CI.r, '_', type.convert=T)[[2]],
           chr=QTL.table$chr), ignore.strand=T)
QTL.table$large_effect=QTL.table$vexp>.02
#sum(QTL.table$large_effect)
QTL.table$CI95=CI
t95=t(sapply(QTL.table[QTL.table$large_effect,]$pCausal, function(x) {
       pos=tstrsplit(names(x),'_',type.convert=T)[[2]]
       dorder=order(x, decreasing=T)
       pdorder=pos[dorder]
       pCum=cumsum(x[dorder]) 
       if(pCum[1]>.95) {return(c(pdorder[1], pdorder[1])) 
       }else {
            return(range(pdorder[1:max(which(pCum<.95))]))
       } }))
t95df=data.frame(chr=QTL.table[QTL.table$large_effect,]$chr, t95)
names(t95df)[2:3]=c('start','end')
QTL.table[QTL.table$large_effect,]$CI95=makeGRangesFromDataFrame(t95df, ignore.strand=T)

# now iterate through QTL confidence intervals per trait
QTs=split(QTL.table, QTL.table$trait)
QTG.out=list()
for(tt in names(QTs)) {
    #tt=names(QTs)
    print(tt)
    QTsT=QTs[[tt]]
    print(tt)
    Ql=QTsT[QTsT$large_effect,]
    Ql.big=Ql[order(Ql$vexp, decreasing=T),]
    by.qtl=list()
    for(qoi in 1:length(Ql.big) ) { 
        ppos=tstrsplit(names(Ql.big[qoi,]$pCausal[[1]]), '_', type.convert=T)[[2]]
        gdf=makeGRangesFromDataFrame(data.frame(chr=Ql.big[qoi,]$chr, start=ppos, end=ppos ), ignore.strand=T)
        gdf$pCausal=as.vector(Ql.big[qoi,]$pCausal[[1]])

        overlaps = findOverlaps(gdf, geneExtend.GR2)
        signal = gdf$pCausal[queryHits(overlaps)]
        averagedSignal = aggregate(signal, list(subjectHits(overlaps)), sum)
        ndf=(data.frame(geneExtend.GR2[averagedSignal[,1],], pCausalSum=averagedSignal[,2]))
        by.gene=cbind(do.call('rbind', replicate(nrow(ndf),Ql.big[qoi,])),ndf)
        by.qtl[[as.character(qoi)]]=by.gene
        print(by.gene)
    }
    QTG.out[[tt]]=by.qtl
}

QTG.table=lapply(QTG.out,function(x) do.call('rbind', x))
QTG.table=do.call('rbind', QTG.table)



# shared parent analysis
QTGR.out=list()
for(tt in names(QTs)) {
    print(tt)
    #tt=names(QTs)
    QTsT=QTs[[tt]]
    print(tt)
    Ql=QTsT[QTsT$large_effect,]
    Ql.big=Ql[order(Ql$chr, Ql$vexp),]
    Ql.big$matched=F
    #Ql.big=Ql[order(Ql$vexp, decreasing=T),]
    osets=sapply(Ql.big$pCausal, names)
    oo=list()
    for(i in 1:length(osets)) {
        oo[[i]]=sapply(osets, function(x) osets[[i]] %in% x)
    }
    oQ=findOverlaps(Ql.big$CI95, Ql.big$CI95)
    sOq=split(subjectHits(oQ), queryHits(oQ))
    #Ql.big[sOq[[2]],][1,]$pCausal
    #Ql.big[sOq[[2]],][2,]$pCausal
    overlapping.QTL=as.vector(which(sapply(sOq, length)>1))
    Ql.big=Ql.big[overlapping.QTL,]
    osets=sapply(Ql.big$pCausal, names)
    oo=list()
    for(i in 1:length(osets)) {
        oo[[i]]=sapply(osets, function(x) osets[[i]] %in% x)
    }
    oQ=findOverlaps(Ql.big$CI95, Ql.big$CI95)
    sOq=split(subjectHits(oQ), queryHits(oQ))
    for(i in 1:length(sOq)){ sOq[[i]]=sOq[[i]][!(sOq[[i]] %in% i)] }

    by.qtl=list()
    for(qoi in 1:length(sOq) ) { 
        Q1=Ql.big[qoi,]
        if(Q1$matched) { next; } 
        Qset2=Ql.big[sOq[[qoi]],]
        Qset2=Qset2[!Qset2$matched,]
        #Qset2=Qset2[order(Qset2$vexp, decreasing=T),]
        if(nrow(Qset2)==0 ) { next; }
        qc1=Q1$cross
        qc2=Qset2$cross
        for(q2 in 1:length(qc2)) {
           mm= colSums(parents.to.crosses[c(qc1, qc2[q2]),])
           if(sum(mm==2)>0) {
                shared.parent=names(mm)[mm==2]
                Q2=Qset2[q2,]
                signme1=sign(Q1$betas)
                direction1=ifelse(min(grep(shared.parent, names(parents.list[[qc1]])))==7, signme1, signme1*-1)
                signme2=sign(Q2$betas)
                direction2=ifelse(min(grep(shared.parent, names(parents.list[[Q2$cross]])))==7, signme2, signme2*-1)
                if(direction1==-1 & direction2==-1) { effect='increase' } 
                else if (direction1==1 & direction2==1) {effect='decrease' }
                else {effect='discrepant'; next;}
           
                pCausal1=Q1$pCausal[[1]] 
                pCausal2=Q2$pCausal[[1]]
                p2s=pCausal2[names(pCausal2) %in% names(pCausal1)]
                p1s=pCausal1[names(pCausal1) %in% names(pCausal2)]
                p1s=p1s[unique(names(p1s))]
                p2s=p2s[unique(names(p2s))]
                if(length(p1s)==0 | length(p2s)==0) {next;}
                pjoint=(p1s*p2s)/sum(p1s*p2s)
           
                ppos=tstrsplit(names(pjoint), '_', type.convert=T)[[2]]
                gdf=makeGRangesFromDataFrame(data.frame(chr=Q1$chr, start=ppos, end=ppos ), ignore.strand=T)
                gdf$pCausal=as.vector(pjoint)

                overlaps = findOverlaps(gdf, geneExtend.GR2)
                signal = gdf$pCausal[queryHits(overlaps)]
                averagedSignal = aggregate(signal, list(subjectHits(overlaps)), sum)
                ndf=(data.frame(geneExtend.GR2[averagedSignal[,1],], pCausalSum=averagedSignal[,2]))
                Q1$maxPPC=max(Q1$pCausal[[1]])
                Q1$whichmaxPPC=names(which.max(Q1$pCausal[[1]]))
                Q2$maxPPC=max(Q2$pCausal[[1]])
                Q2$whichmaxPPC=names(which.max(Q2$pCausal[[1]]))

                by.gene=cbind(do.call('rbind', replicate(nrow(ndf),cbind(Q1,Q2))),ndf)
                by.gene$jointmaxPPC=max(pjoint)
                by.gene$whichjointmaxPPC=names(which.max(pjoint))
                by.gene$shared.parent=shared.parent
                by.gene$effect=effect
                
                by.qtl[[as.character(qoi)]]=by.gene
                Ql.big[c(qoi, sOq[[qoi]][q2]),]$matched=T
                print(by.gene)
                #readline()
           }

        }
    }
    QTGR.out[[tt]]=by.qtl
}   
QTGJ.table=lapply(QTGR.out,function(x) do.call('rbind', x))
QTGJ.table=do.call('rbind', QTGJ.table)
QTGsorted=QTGJ.table[order(QTGJ.table$pCausalSum, decreasing=T),]
QTGsorted=QTGsorted[!(QTGsorted$trait %in% c("YPD;;2","YPD;;3")),]
# Remove hits from the two extra ypd experiments
QTGsorted$localFDR=(cumsum(1-QTGsorted$pCausalSum))/(1:nrow(QTGsorted))
QTGsorted$localFDR[QTGsorted$localFDR<0]=0
#save(QTGsorted, file = '/data/rrv2/genotyping/RData/FDR_QTGsorted.RData')

QTGsorted.resolved=QTGsorted[QTGsorted$localFDR<.2,]
Qwrite=data.frame(QTGsorted.resolved[,-grep('pCausal\\.|pCausal$', colnames(QTGsorted.resolved))])
Qwrite=Qwrite[Qwrite$localFDR<.25,]
library(WriteXLS)
WriteXLS(Qwrite, '/data/rrv2/genotyping/RData/QTGresolved.xls')

# load 
library(gdata)
prev.mapped=read.xls('/data/rrv2/genotyping/RData/NIHMS544073-supplement-01.xls', pattern='Table S1')
pm=unique(as.character(prev.mapped[-1,1]))[-c(96,97)][-16]
#R> unique(QTGsorted.resolved$NAME)[which(unique(QTGsorted.resolved$NAME) %in% pm)]
# [1] "PCA1"  "RPI1"  "CYS4"  "HO"    "PHO84" "PDR5"  "GAL3"  "CAT5"  "IRA2" 
#[10] "IRA1"  "MKT1"  "END3"  "FLO11" "SAL1"  "CYR1"  "TAO3"  "RGA1"  "HAP1" 
#18 
# SWH1 (Wang et al.)
# TOR1 and WHI2 (Treusch et al.)
# ENA1 (Steinmetz)
# ENA5 (?)
# KRE33 (desai)
# PMR1 (us)
# MAl11 Lit


library(topGO)
myGene2GO.full.table=read.delim('/data/eQTL/reference/go_annotations_sgd.txt',  sep='\t', header=F, stringsAsFactors=F)
# data frame of gene then GO
myGene2GO=cbind( sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,5])
myGene2GO=na.omit(myGene2GO)
SYS2ORF=unique(cbind(sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,3]))
SYS2ORF.key=list()
SYS2ORF.key[SYS2ORF[,1]]=SYS2ORF[,2]
gene2GOList = lapply(unique(myGene2GO[,1]), function(x){myGene2GO[myGene2GO[,1] == x, 2]})
names(gene2GOList) = unique(myGene2GO[,1])


plotGOToTree <- function(GOdat, GOres, sigThres = 0.0005){
     # only plot if there are any significant GO terms (SEE ABOVE for 
    #"significance"; I am somewhat lenient here):
    #     # we need these extra lines because very small p-values are 
    #reported as a text string "< X", rather than a numeric
     toTest <- as.numeric(GenTable(GOdat, pVal = GOres)[1,6])
     if(is.na(toTest)){toTest <- 0.000000000000000000000000000001}
     if (toTest < sigThres){
         showSigOfNodes(GOdat, score(GOres), firstSigNodes = 
        length(score(GOres)[score(GOres) < sigThres]), useInfo = "all")
     }else{
         plot(1,1)
         legend("topright", legend="no significant GO categories", 
        cex=0.4, box.lty=0)
     }
}


#testGenes= factor(0+(colnames(t.tpm.matrix) %in% geneList)) # colnames(t.tpm.matrix)[(which(vcA[,1]>.7))]) #totest )
#names(testGenes)=colnames(t.tpm.matrix)
#testGeneSet=names(testGenes[testGenes==1])

allG=geneExtend.GR2$ORF[1:6663]
sigG=unique(QTGsorted[QTGsorted$localFDR<.2,]$ORF)

testGenes=factor(0+allG %in% sigG)
names(testGenes)=allG

GO.results=list()
for( thisOntology in c("BP", "MF", "CC") ) {
        GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=3)
        GOresult = runTest(GOData, algorithm="classic", statistic="fisher")
        pdf(file=paste0('/data/rrv2/Figures_and_Tables/GO/', thisOntology, '.pdf'), width=25, height=25)
        plotGOToTree(GOData, GOresult, sigThres = 0.00005)
        dev.off()
        gt=GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult)))
        
        genesinterms=genesInTerm(GOData, gt[,1])
        genes.enriched.list=lapply(genesinterms, function(x) x[x%in%names(testGenes[testGenes==1])])
        genes.enriched.list.simple=lapply(genes.enriched.list, function(x) as.character(SYS2ORF.key[x]))
        gt$result1=as.numeric(gt$result1)
        gt$fdr=p.adjust(gt$result1,method='fdr')

        gt$Genes=as.vector(sapply(genes.enriched.list.simple, paste, collapse=','))
        gt$GenesSystematic=   as.vector(sapply(genes.enriched.list, paste, collapse=','))
        GO.results[[thisOntology]]=gt
        write.table(gt, file=paste0('/data/rrv2/Figures_and_Tables/GO/', '_', thisOntology, '.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
    }

#GO.results[[1]][p.adjust(GO.results[[1]]$result1,method='fdr')<.05,]
#GO.results[[2]][p.adjust(GO.results[[2]]$result1,method='fdr')<.05,]
#GO.results[[3]][p.adjust(GO.results[[3]]$result1,method='fdr')<.05,]

# cAMP and glucose sensing
# plasma membrane
# PCA1 allelic heterogeneity












test.marker='chrXIV:467219_A/G' #chrI:203922_A/G'
test.marker=gsub(':|/', '_', test.marker)
# match against one of the two crosses

QTGsortedYPD=QTGsorted[QTGsorted$trait.1=='YPD;;1',]
otest=as.logical(sapply(QTGsortedYPD$pCausal, function(x) sum(grepl(test.marker, names(x)))))
QTGsortedYPD[otest,][1,]


table(QTGsorted[QTGsorted$localFDR<.05,]$shared.parent,QTGsorted[QTGsorted$localFDR<.05,]$effect)

Qwrite=data.frame(QTGsorted[,-grep('pCausal\\.|pCausal$', colnames(QTGsorted))])
Qwrite=Qwrite[Qwrite$localFDR<.25,]
library(WriteXLS)
WriteXLS(Qwrite, '~/Desktop/QTGs_030818.xls')

tgtable=t(table(QTGsorted[QTGsorted$localFDR<.25,]$shared.parent,QTGsorted[QTGsorted$localFDR<.25,]$effect))
tgtable=tgtable[,colnames(parents.to.crosses)]

barplot(tgtable, 
        beside=2, las=2, col=c('black', 'grey'), main='large effect QTL direction by shared parent')
legend("topright", legend=c("decreasing", "increasing"), fill=c("black", "grey"))



# visualizing PPCs       
trait='Cobalt_Chloride;2mM;2'
cross.name='375'
cross.name='A'
p1='375'
p2='A'

    cross=cross.list[[cross.name]]
    if(cross.name=='A') {       cross=subset(cross, ind=!grepl('A11', as.character(cross$pheno$id)))    }
    snames = as.character(cross$pheno$id)
    g=pull.argmaxgeno(cross)
 
    # recode based on parental genotypes 
    seg.pcode=recode.as.allele.number(t(g),parents.list[[cross.name]])
    #seg.recoded[[cross.name]]=seg.pcode
   
    # are there fixed loci ?? (no)-------------------------------
    g.af=apply(g,2,function(x) sum(x==1))
    parents.list[[cross.name]]$fixed=(g.af==0 | g.af==nrow(g))
    fixed.loci=which(parents.list[[cross.name]]$fixed)
    if(length(fixed.loci)>0) {    g=g[,-fixed.loci] }
    #------------------------------------------------------------
    g.r=g[,-which(duplicated(g, MARGIN=2))]
    g.s=scale(g.r)

    subPheno=lapply(NORMpheno, function(x) x[match(snames, names(x))])
    mPheno  =sapply(subPheno, function(x) sapply(x, mean, na.rm=T))

    mPheno=apply(mPheno,2, function(x) {x[is.na(x)]=mean(x, na.rm=T); return(x)})
    mPhenos=scale(mPheno)

    lods_B=fasterLOD(nrow(mPhenos), mPhenos, g.s)
    m_B=(tstrsplit(colnames(lods_B), '_', type.convert=T))

    lods_377=fasterLOD(nrow(mPhenos), mPhenos, g.s)
    m_377=(tstrsplit(colnames(lods_377), '_', type.convert=T))

#B=YPSxYJM
#377 = YJMxCLIB
QGset=QTGsorted[QTGsorted$shared.parent=='BYa',]

qtgi=34
qgic=QGset[qtgi,]$chr
qgit=QGset[qtgi,]$trait

m_B_c2=m_B[[1]]==qgic
m_377_c2=m_377[[1]]==qgic
tn=match(qgit, rownames(lods_B))

xxlim=c(QGset[qtgi,]$pos-50000, QGset[qtgi,]$pos+50000)

par(mfrow=c(2,1),mar = c(5,5,2,5))
plot(m_B[[2]][m_B_c2],lods_B[tn, m_B_c2], ylab=paste(p1,'LOD'), xlab='chr pos', xlim=xxlim,
     main=paste(qgit, qgic,QGset[qtgi,]$shared.parent,QGset[qtgi,]$efffect, QGset[qtgi,]$NAME), 
     sub=round(QGset[qtgi,]$pCausalSum,2))
abline(v=tstrsplit(QGset[qtgi,]$CI.l, '_', type.convert=T)[[2]])
abline(v=tstrsplit(QGset[qtgi,]$CI.r, '_', type.convert=T)[[2]])

segments(start(gene.GR)[as.vector(seqnames(gene.GR)==qgic)],
         0, end(gene.GR)[as.vector(seqnames(gene.GR)==qgic)],0)
par(new = T)
pB=QGset[qtgi,]$pCausal[[1]]
plot(as.numeric(tstrsplit(names(pB), '_',type.convert=T)[[2]]),
     pB, col='red', axes=F, ylab=NA, xlab=NA, xlim=xxlim)
axis(side = 4)
mtext(side = 4, line = 3, 'posterior probability of causality')

plot(m_377[[2]][m_377_c2],lods_377[tn, m_377_c2], ylab=paste(p2, 'LOD'), xlab='chr pos', xlim=xxlim)
abline(v=tstrsplit(QGset[qtgi,]$CI.l.1, '_', type.convert=T)[[2]])
abline(v=tstrsplit(QGset[qtgi,]$CI.r.1, '_', type.convert=T)[[2]])

segments(start(gene.GR)[as.vector(seqnames(gene.GR)==qgic)],
         0, end(gene.GR)[as.vector(seqnames(gene.GR)==qgic)],0)
p377=QGset[qtgi,]$pCausal.1[[1]]
par(new = T)
plot(as.numeric(tstrsplit(names(p377), '_',type.convert=T)[[2]]),
     p377, col='red', axes=F, ylab=NA, xlab=NA, xlim=xxlim)
axis(side = 4)
mtext(side = 4, line = 3, 'posterior probability of causality')



p377=QTGsorted[1,]$pCausal.1[[1]]
as.numeric(tstrsplit(names(p377), '_',type.convert=T)[[2]]),p377, col='red')






#for(cross.name in crosses) {
#    print(cross.name)
#    parent.vcf.file=   paste0('/data/rrv2/genotyping/parents/cross_vcf/', cross.name,'_w_svar.vcf.recode.vcf.gz')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#    
#    library(VIGoR)
#chrom='chrVII'
#plist.chr=sapply(parents.list, function(x) x[x$chr==chrom,])
#chr.markers=sapply(parents.list, function(x){ return(x$marker.name.n[x$chr==chrom]) } )
#seg.chr.phenos=lapply(pheno.resids, function(x) { x[[chrom]] })
#seg.chr.n=lapply(seg.chr.phenos, nrow)
#seg.chr.phenos.scaled=lapply(seg.chr.phenos,scale)
#
#seg.chr.genos=mapply(function(cmm,seg.r) { t(seg.r[cmm,])} , cmm=chr.markers, seg.r=seg.recoded)
#seg.chr.genos.scaled=lapply(seg.chr.genos, scale)
#cmm=sapply(parents.list, function(x){ return(x$marker.name[x$chr==chrom]) } )
##rename seg.chr.genos.scaled
#seg.chr.genos.scaled=mapply(function(cmm,seg.r) { colnames(seg.r)=cmm; return(seg.r);}, cmm=cmm, seg.r=seg.chr.genos.scaled)
#
#mpmarkers=melt(cmm)
#mpmarkers[,1]=as.character(mpmarkers[,1])
#names(mpmarkers)=c('marker', 'cross')
#
#mpmarkers=separate(mpmarkers, marker, c('chr', 'pos', 'ref', 'alt'), sep='_', convert=T, remove=F)
#mpmarkers=mpmarkers[order(mpmarkers$pos),]
#
#ccnt=split(mpmarkers$cross, mpmarkers$marker)
#crosses.per.marker=sapply(ccnt, function(x) unique(x))
#cross.cnt.per.marker=sapply(crosses.per.marker, length)
##seg.rare=names(cross.cnt.per.marker[cross.cnt.per.marker<3])
#
#mm=mpmarkers[match(names(crosses.per.marker), mpmarkers$marker),-6]
#mm=mm[order(mm$pos),]
#
## mapping cross marker index to joint marker table  
#mupos=lapply(plist.chr, function(x) match(x$marker.name, mm$marker))
## mapping cross marker to joint marker table name 
#mucpos = mapply(function(m, p) { p$marker.name[ match(mm$marker[m], p$marker.name)] }, m=mupos, p=plist.chr )
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#cc='A'
#tt=
#cpeaks=cross.peaks[[1]][grep('Mang', cross.peaks[[1]]$gene),]
#cpeaks=separate(cpeaks, pmarker, c('chr', 'pos', 'ref', 'alt', 'index'), sep='_', convert=T, remove=F)


