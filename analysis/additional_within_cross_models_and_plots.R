# plot composite mapping per cross 
peaksModel=list()
interactionPeaks=list()
marginalR=list()
for(cross.name in crosses) {
    print(cross.name)
    cross=cross.list[[cross.name]]
    if(cross.name=='A') {       cross=subset(cross, ind=!grepl('A11', as.character(cross$pheno$id)))    }
    snames = as.character(cross$pheno$id)
    subPheno=lapply(NORMpheno, function(x) x[match(snames, names(x))])
    mPheno  =sapply(subPheno, function(x) sapply(x, mean, na.rm=T))
    mPheno=apply(mPheno,2, function(x) {x[is.na(x)]=mean(x, na.rm=T); return(x)})
    g=pull.argmaxgeno(cross)
    # are there fixed loci ?? (no)-------------------------------
    #g.af=apply(g,2,function(x) sum(x==1))
    #parents.list[[cross.name]]$fixed=(g.af==0 | g.af==nrow(g))
    #fixed.loci=which(parents.list[[cross.name]]$fixed)
    #if(length(fixed.loci)>0) {    g=g[,-fixed.loci] }
    #------------------------------------------------------------
    #g.r=g[,-which(duplicated(g, MARGIN=2))]
    g.s=scale(g)
    #A=tcrossprod(g.s)/(ncol(g.s))
    mPheno=scale(mPheno)
    marginalR[[cross.name]]=(crossprod(mPheno,g.s)/(nrow(mPheno)-1))
    #.1581 = var exp ~.05
    #.2236 = var exp ~.16
    cps=cross.peaks[[cross.name]]
    cps=cps[cps$q<.05,]
    # remove 4NQO, YPD;;2 and YPD;;3 
    for(pheno in names(subPheno)[-c(1,38,39)]) {
            print(pheno)
            cpQTL=cps[cps$trait==pheno,]
            if(length(cpQTL$pmarker)!=0) {
                apeaks = unique(match(cpQTL$fscan.markers, colnames(g.s)))
                X=data.frame(g.s[,apeaks])
                fitme=lm(mPheno[,pheno]~.-1, data=X)
                aov.a = anova(fitme)
                tssq  = sum(aov.a[,2])
                a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
                coeffs=coefficients(fitme)
                cpQTL$var.exp=a.effs
                cpQTL$lm.coeff=as.vector(coeffs)
                cpQTL$chr=sapply(strsplit(cpQTL$pmarker, '_'), function(x) x[1])
                cpQTL$pos=as.numeric(sapply(strsplit(cpQTL$pmarker, '_'), function(x) x[2]))
                cpQTL$cross=cross.name
                names(cpQTL)[1]='trait'
                print(cpQTL)
                peaksModel[[cross.name]][[pheno]]=cpQTL
             }
             
            if(length(cpQTL$pmarker)>1) {
                qtl.combs=combn(apeaks,2)
                null=lm(mPheno[,pheno]~g.s[,apeaks]-1)
                int.coef1=rep(NA, ncol(qtl.combs))
                int.coef2=rep(NA, ncol(qtl.combs))
                int.coef=rep(NA, ncol(qtl.combs))
                int.pvalue=rep(NA, ncol(qtl.combs))
                
                 for(ist in 1:ncol(qtl.combs)){
                        full=lm(mPheno[,pheno]~g.s[,apeaks]+g.s[,qtl.combs[1,ist]]*g.s[,qtl.combs[2,ist]]-1)
                        int.pvalue[ist]=anova(null,full)$'Pr(>F)'[2]
                        int.coef1[ist]=coef(full)[paste0("g.s[, apeaks]",colnames(g.s)[qtl.combs[1,ist]])]
                        int.coef2[ist]=coef(full)[paste0("g.s[, apeaks]",colnames(g.s)[qtl.combs[2,ist]])]
                        int.coef[ist]=coef(full)[length(coef(full))] #anova(null,full)$'Pr(>F)'[2]
                    }
                tqc=t(qtl.combs)
                dfi=data.frame(m1=colnames(g.s)[tqc[,1]], m2=colnames(g.s)[tqc[,2]], int.coef1, int.coef2, int.coef, int.pvalue, stringsAsFactors=F)
                dfi$cross=cross.name
                dfi$chr1=sapply(strsplit(dfi$m1, '_'), function(x) x[1])
                dfi$chr2=sapply(strsplit(dfi$m2, '_'), function(x) x[1])
                dfi$pos1=as.numeric(sapply(strsplit(dfi$m1, '_'), function(x) x[2]))
                dfi$pos2=as.numeric(sapply(strsplit(dfi$m2, '_'), function(x) x[2]))
                dfi$trait=pheno
                interactionPeaks[[cross.name]][[pheno]]=dfi
                #interactions_per_trait[[pheno]]=dfi
            }
    }
}
#save(marginalR,file='/data/rrv2/genotyping/RData/FDR_marignalR.RData')
#save(peaksModel,file='/data/rrv2/genotyping/RData/FDR_wcPeaksModel.RData')
#save(interactionPeaks, file='/data/rrv2/genotyping/RData/FDR_wcInteractionPeaksModel.RData')
load('/data/rrv2/genotyping/RData/FDR_marignalR.RData')
load('/data/rrv2/genotyping/RData/FDR_wcPeaksModel.RData')
load('/data/rrv2/genotyping/RData/FDR_wcInteractionPeaksModel.RData')

cross.peaks.flat=do.call('rbind', lapply(peaksModel, function(y) { do.call('rbind', y)} )) #cross.peaks)
cross.peaks.flat$gcoord=gcoord.key[cross.peaks.flat$chr]+cross.peaks.flat$pos
interactionPeaks.flat=do.call('rbind', lapply(interactionPeaks, function(y){ do.call('rbind', y)}))
qs.int=qvalue(interactionPeaks.flat$int.pvalue, fdr.level=.1)
interactionPeaks.flat$significant=qs.int$qvalues<.1
intP=interactionPeaks.flat[interactionPeaks.flat$significant,]
intP$gcoord1=gcoord.key[intP$chr1]+intP$pos1
intP$gcoord2=gcoord.key[intP$chr2]+intP$pos2
intP=na.omit(intP)

ssi=split(intP, paste(intP$trait, intP$cross) )
hist(c(sapply(ssi, nrow), rep(0, (38*16)-length(ssi))))

glength=sum(unlist(chr.lengths))

#load('/data/rrv2/genotyping/RData/jointPeaks5.RData')
#jP=rbindlist(jointPeaks5, idcol='chromosome')
#jPs=split(jP, jP$trait)
jointPeaksFlat=rbindlist(jointPeaksJS, idcol='chromosome')
#data.frame(do.call('rbind', jointPeaks5), stringsAsFactors=F)
names(jointPeaksFlat)[1]='chr'
#jointPeaksFlat$chr=sapply(strsplit(jointPeaksFlat$marker, '_'), function(x) x[1])
jointPeaksFlat$pos=as.numeric(sapply(strsplit(jointPeaksFlat$fscan.markers, '_'), function(x) x[2]))
jointPeaksFlat$gpos=gcoord.key[jointPeaksFlat$chr ]+jointPeaksFlat$pos



utraits.orig=unique(cross.peaks.flat$trait)
utraits=utraits.orig
utraits[34]='YNB_ph8'
utraits[36]='YPD_15C'
utraits[33]='YNB_ph3'
utraits[10]='EtOH_Glu'
utraits[37]='YPD_37C'
utraits=gsub(';.*','', utraits)
utraits=gsub('_', ' ', utraits)

pdf(file=paste0('/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure2.pdf'), width=11, height=8)

for(piter in 1:length(utraits))  {
    png(file=paste0('/home/jbloom/Dropbox/RR/Figures and Tables/other formats/SuplementaryFigure2_', piter, '.png'), width=1100, height=800)

    pheno.iter=utraits.orig[piter]
#pdf(file=paste0('/home/jbloom/Dropbox/RR/Figures and Tables/', filename.clean(pheno.iter), '_joint.pdf'), width=1024, height=800)

pcnt=0
op <- par(mfrow = c(16,1),
          oma = c(5,8,5,.5) + 0.01,
          mar = c(0,4,0,0) + 0.01,
          xaxs='i',
          yaxs='i'
          )

joint.peaks.toplot=jointPeaksFlat[jointPeaksFlat$trait==pheno.iter,]
joint.peaks.toplot=joint.peaks.toplot[!duplicated(joint.peaks.toplot$gpos),]

parent.vec=c('M22', 'BY', 'RM', 'YPS163', 'YJM145', 'CLIB413', 'YJM978', 'YJM454',
             'YPS1009', 'I14', 'Y10', 'PW5', '273614N', 'YJM981', 'CBS2888', 'CLIB219')
glength=1.2e7
for(cross.iter in 1:length(crosses)){
    cross.name=crosses[cross.iter]
    jptlookup=joint.peaks.toplot[joint.peaks.toplot$fscan.markers %in% parents.list[[cross.name]]$marker.name,]
    cross.sub.p=cross.peaks.flat[cross.peaks.flat$trait==pheno.iter & cross.peaks.flat$cross==crosses[cross.iter],]
    cross.sub.pi=intP[intP$trait==pheno.iter & intP$cross==crosses[cross.iter],]
    mpM=marginalR[[crosses[cross.iter]]]
    mpM.marker=tstrsplit(colnames(mpM), '_', type.convert=T)
    mpM.gcoord=gcoord.key[mpM.marker[[1]]]+mpM.marker[[2]]
    if(nrow(cross.sub.p)>0) {
        #plot(0,0, type='n', xlim=c(0, glength), yaxt='n', ylab='', xaxt='n', yaxt='n', ylim=c(-1,1),cex.lab=1.5)
        plot(0,0, type='n', xlim=c(0, glength), yaxt='n', ylab='', xaxt='n', yaxt='n', ylim=c(0,1),cex.lab=1.5)
        abline(h=0)
        abline(v=jptlookup$gpos, col='lightgreen')
        axis(2,at=1,labels=parent.vec[cross.iter], cex.axis=1.5, las=2) 
        pcnt=pcnt+nrow(cross.sub.p)
        signme=sign(cross.sub.p$lm.coeff)
        if(min(grep(crosses.to.parents[[cross.iter]][1], names(parents.list[[cross.iter]])))==7) {
            signme=signme
            signR=-1*mpM[pheno.iter,]
        } else{
         signme=-1*signme
         signR=mpM[pheno.iter,]
        }
        # flip sign ... (if negative then point to strain that increases growth)
        signme=-1*signme
        splus=signme==1
        sminus=signme==-1
        #points(mpM.gcoord,  signR)
        
        cross.sub.p$lm.ceiling=rep(.1, nrow(cross.sub.p))
        cross.sub.p$lm.ceiling[cross.sub.p$var.exp>0]=.2
        cross.sub.p$lm.ceiling[cross.sub.p$var.exp>.04]=.5
        cross.sub.p$lm.ceiling[cross.sub.p$var.exp>.08]=.75
        cross.sub.p$lm.ceiling[cross.sub.p$var.exp>.25]=1
        if(sum(splus)>0){
            arrows(cross.sub.p$gcoord[splus],0, cross.sub.p$gcoord[splus], cross.sub.p$lm.ceiling[splus], code=2, length=.12, lwd=4, 
                   col=ifelse(cross.sub.p$lm.ceiling[splus]>.2, 'black', 'grey'))
        }
        if(sum(sminus)>0){
            arrows(cross.sub.p$gcoord[sminus],cross.sub.p$lm.ceiling[sminus], cross.sub.p$gcoord[sminus],0 , code=2, length=.12, lwd=4,
                   col=ifelse(cross.sub.p$lm.ceiling[sminus]>.2, 'black', 'grey')
                   )
        }
        abline(v=gcoord.key, lty=2, col='lightblue')
    } else {
             plot(0,0, type='n', xlim=c(0, max(glength)), ylim=c(0,1),  xaxt='n' ) #ylab=crosses[cross.iter] ,
             #abline(h=0, lty=3, col='grey')
             abline(v=cumsum(genome.chr.lengths), lty=2, col='lightblue')

        }
   # if(nrow(cross.sub.pi)>0) {
   #         peak.number=c(seq_along(cross.sub.pi[,1]), c(seq_along(cross.sub.pi[,2])))
   #         #peak.chr=c(cross.sub.pi$chr1, cross.sub.pi$chr2)
   #         #peak.pos=as.numeric(sapply(strsplit(sapply(strsplit(c(cross.sub.pi.sig[,1], cross.sub.pi.sig[,2]), ':'), function (x) x[2]), '_'), function(x)x[1]))
   #         peak.gpos=c(cross.sub.pi$gcoord1, cross.sub.pi$gcoord2)
   #         text(peak.gpos, (peak.number/max(peak.number))*.9, '*', col='red', cex=4)
   #       }

    if(cross.iter==16){ axis(1, at=gcoord.key, labels=names(gcoord.key), cex.axis=1.5)}
}
title(xlab='genomic position', ylab='', outer=TRUE, cex.lab=2, 
      main=paste(utraits[piter], '    ',  pcnt, 'total QTL       |  ', 
                 length(joint.peaks.toplot$gpos), 'joint QTL'
                 ))
dev.off()

}
   # dev.off()


