load('/data/rrv2/genotyping/RData/FDR_seg.recoded.RData')
load('/data/rrv2/genotyping/RData/extracted_average_phenotypes.RData')
load('/data/rrv2/genotyping/RData/jointPeakEffects_JS_variants.RData')
load('/data/rrv2/genotyping/RData/iseq.freqs.RData')
load('/data/rrv2/genotyping/RData/parents.list.RData')
parents.list=lapply(parents.list, function(x) {
                  z=x;
                  z$marker.name.n=paste0(z$marker.name, '_', seq(1:nrow(z)))
                  return(z) })


source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')

source('/data/rrv2/analysis/mapping_fx.R')

library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=48)

jPs=jointPeakEffects
# Interactions 
#jointPeakEffectsInt=list()
# might gain additional 
jointInteractionPeaks=list()
for(cross.name in crosses) {
    print(cross.name)
    g.s=t(seg.recoded[[cross.name]]) #scale(t(seg.recoded[[cross.name]]))
    #rename columns
    colnames(g.s)=parents.list[[cross.name]]$marker.name
    mPheno=pheno_extracted[[cross.name]]
    mPhenos  = scale(mPheno)
    
    jointInteractionPeaks[[cross.name]]=foreach(tt = names(jPs) ) %dopar% {
        print(tt)
        dfi=NULL
#        if(length(jointPeakEffects[[tt]][[cross.name]]$peaks)>1) {
#                apeaks=jointPeakEffects[[tt]][[cross.name]]$peaks
#                qtl.combs=combn(apeaks,2)
#                #null=lm(mPhenos[,tt]~g.s[,apeaks]-1)
#                int.coef1=rep(NA, ncol(qtl.combs))
#                int.coef2=rep(NA, ncol(qtl.combs))
#                int.coef=rep(NA, ncol(qtl.combs))
#                int.pvalue=rep(NA, ncol(qtl.combs))
#                
#                 for(ist in 1:ncol(qtl.combs)){
#                        #print(ist)
#                        full=lm(mPhenos[,tt]~g.s[,apeaks]+g.s[,qtl.combs[1,ist]]*g.s[,qtl.combs[2,ist]])
#                        int.pvalue[ist]=drop1(full, 'g.s[, qtl.combs[1, ist]]:g.s[, qtl.combs[2, ist]]', test='Chisq')[[5]][2]
#                        #anova(null,full)$'Pr(>F)'[2]
#                        coefs=coef(full)
#                        int.coef1[ist]=coefs[paste0("g.s[, apeaks]",qtl.combs[1,ist])]
#                        int.coef2[ist]=coefs[paste0("g.s[, apeaks]",qtl.combs[2,ist])]
#                        int.coef[ist]=coefs[length(coefs)] #anova(null,full)$'Pr(>F)'[2]
#                    }
#                tqc=t(qtl.combs)
#                dfi=data.frame(m1=tqc[,1], m2=tqc[,2], int.coef1, int.coef2, int.coef, int.pvalue, stringsAsFactors=F)
#                dfi$cross=cross.name
#                dfi$trait=tt
#                dfi$chr1=sapply(strsplit(dfi$m1, '_'), function(x) x[1])
#                dfi$pos1=as.numeric(sapply(strsplit(dfi$m1, '_'), function(x) x[2]))
#                dfi$pos1_maf1012=iseq.freqs$maf1012[match(dfi$m1, iseq.freqs$marker)]
#                dfi$pos1_alt012=iseq.freqs[,8][match(dfi$m1, iseq.freqs$marker)]
#                dfi$pos1_ancestral=iseq.freqs$ancestral[match(dfi$m1, iseq.freqs$marker)]
#                dfi$pos1_maf2=ifelse(dfi$pos1_ancestral, dfi$pos1_alt012, 1-dfi$pos1_alt012)
#
#                dfi$chr2=sapply(strsplit(dfi$m2, '_'), function(x) x[1])
#                dfi$pos2=as.numeric(sapply(strsplit(dfi$m2, '_'), function(x) x[2]))
#                dfi$pos2_maf1012=iseq.freqs$maf1012[match(dfi$m2, iseq.freqs$marker)]
#                dfi$pos2_alt012=iseq.freqs[,8][match(dfi$m2, iseq.freqs$marker)]
#                dfi$pos2_ancestral=iseq.freqs$ancestral[match(dfi$m2, iseq.freqs$marker)]
#                dfi$pos2_maf2=ifelse(dfi$pos2_ancestral, dfi$pos2_alt012, 1-dfi$pos2_alt012)
#                #jointInteractionPeaks[[cross.name]][[tt]]=dfi
#                #interactions_per_trait[[pheno]]=dfi
#        }
        dfii=NULL
        if(length(jointPeakEffects[[tt]][[cross.name]]$peaks)>2) {
                apeaks=jointPeakEffects[[tt]][[cross.name]]$peaks
                qtl.combs=combn(apeaks,3)
                tqc=t(qtl.combs)
                
                dfii=data.frame(m1=tqc[,1], m2=tqc[,2], m3=tqc[,3],stringsAsFactors=F)
                dfii$cross=cross.name
                dfii$trait=tt
                dfii$chr1=sapply(strsplit(dfii$m1, '_'), function(x) x[1])
                dfii$pos1=as.numeric(sapply(strsplit(dfii$m1, '_'), function(x) x[2]))
                dfii$chr2=sapply(strsplit(dfii$m2, '_'), function(x) x[1])
                dfii$pos2=as.numeric(sapply(strsplit(dfii$m2, '_'), function(x) x[2]))
                dfii$chr3=sapply(strsplit(dfii$m3, '_'), function(x) x[1])
                dfii$pos3=as.numeric(sapply(strsplit(dfii$m3, '_'), function(x) x[2]))
                
                qtl.combs=qtl.combs[,which(dfii$chr1!=dfii$chr2  &  dfii$chr1!=dfii$chr3)]
                dfii=dfii[which(dfii$chr1!=dfii$chr2  &  dfii$chr1!=dfii$chr3),]

                int.pvalue=rep(NA, ncol(qtl.combs))
                int.coef=rep(NA, ncol(qtl.combs))
                int.vexp=rep(NA, ncol(qtl.combs))
                null=lm(mPhenos[,tt]~g.s[,apeaks])
                for(ist in 1:ncol(qtl.combs)){
                    if(ist==ncol(qtl.combs)) {print(tt)}
                    full=lm(mPhenos[,tt]~g.s[,apeaks]+g.s[,qtl.combs[1,ist]]*g.s[,qtl.combs[2,ist]]*g.s[,qtl.combs[3,ist]])
                    #mm=(model.matrix(mPhenos[,tt]~g.s[,apeaks]+g.s[,qtl.combs[1,ist]]*g.s[,qtl.combs[2,ist]]*g.s[,qtl.combs[3,ist]]))
                    #full=lm.fit(mm, mPhenos[,tt])
                    int.vexp[ist]=anova(null,full)[2,4]/(nrow(mPhenos)-1)
                    coefs=coef(full)
                    int.pvalue[ist]=drop1(full, 'g.s[, qtl.combs[1, ist]]:g.s[, qtl.combs[2, ist]]:g.s[, qtl.combs[3, ist]]', test='Chisq')[[5]][2]
                    int.coef[ist]=coefs[length(coefs)] #anova(null,full)$'Pr(>F)'[2]
                }
                dfii$int.coef=int.coef
                dfii$int.pvalue=int.pvalue
                dfii$int.vexp=int.vexp
        }
        return(list(twoLocus=dfi, threeLocus=dfii))
    }
    names(jointInteractionPeaks[[cross.name]])=names(jPs)
}
#save(jointInteractionPeaks, file='/data/rrv2/genotyping/RData/jointInteractionPeaks.RData')
#save(jointInteractionPeaks, file='/data/rrv2/genotyping/RData/jointInteractionPeaksVE.RData')
load('/data/rrv2/genotyping/RData/jointInteractionPeaks.RData')
interactionPeaks.flat=rbindlist(lapply(jointInteractionPeaks, function(y){ rbindlist(lapply(y, function(z) z$twoLocus),idcol='trait') } ) , idcol='cross')
interactionPeaks.flat=interactionPeaks.flat[interactionPeaks.flat$trait!="YPD;;2" & interactionPeaks.flat$trait!="YPD;;3",]
plot(qvalue(interactionPeaks.flat$int.pvalue))

q1=(qvalue(interactionPeaks.flat$int.pvalue))
sum(q1$qvalue<.05, na.rm=T)

interactionPeaks.flat3D=rbindlist(lapply(jointInteractionPeaks, function(y){ rbindlist(lapply(y, function(z) z$threeLocus),idcol='trait') } ) , idcol='cross')
#whoops
interactionPeaks.flat3D=interactionPeaks.flat3D[interactionPeaks.flat3D$chr2!=interactionPeaks.flat3D$chr3,]

q3d=qvalue(interactionPeaks.flat3D$int.pvalue)
sum(q3d$qvalue<.05, na.rm=T)

plot(q3d)

ipSig=interactionPeaks.flat3D[which(q3d$qvalues<.05),]
ipSig[order(ipSig$int.pvalue, decreasing=F)[1:10],]
































#do.call('rbind', y$twoLocus)}))
#with(interactionPeaks.flat, {plot_ly(x=(pos1_maf1012), y=(pos2_maf1012), z=abs(int.coef), sizes=c(.5,1), type="scatter3d", mode="markers", color=abs(int.coef)) })
# p.coef=predict(lm(abs(int.coef)~pos1_maf1012*pos2_maf1012, data=interactionPeaks.flat))
#
#summary(lm(log10(abs(int.coef))~log10(pos1_maf1012)*log10(pos2_maf1012), data=interactionPeaks.flat))
#obs=summary(lm(log10(abs(int.coef))~log10(pos1_maf1012)*log10(pos2_maf1012), data=interactionPeaks.flat))$coefficients[4,3]
#eee=rep(NA, 1000)
#for(i in 1:1000) {
#    print(i)
#    eee[i]=summary(lm(log10(sample(abs(int.coef)))~log10(pos1_maf1012)*log10(pos2_maf1012), data=interactionPeaks.flat))$coefficients[4,3]
#
#}
#
# p.coef=predict(lm(log10(abs(int.coef))~log10(pos1_maf1012)*log10(pos2_maf1012), data=interactionPeaks.flat))
#
#with(interactionPeaks.flat, {plot_ly(x=(pos1_maf1012), y=(pos2_maf1012), z=10^p.coef, sizes=c(.5,1), type="scatter3d", mode="markers", color=abs(int.coef)) })
#
#
#
#qs.int=qvalue(interactionPeaks.flat$int.pvalue, fdr.level=.1)
#interactionPeaks.flat$significant=qs.int$qvalues<.1
#intP=interactionPeaks.flat[interactionPeaks.flat$significant,]
#intP$binned=cut(abs(intP$int.coef), breaks=c(0,0.1,0.2,0.3,0.4), include.lowest=F)
#intP=intP[!is.na(intP$int.coef),]
#intP$minAF=ifelse(intP$pos1_maf1012<intP$pos2_maf1012, intP$pos1_maf1012,intP$pos2_maf1012)
#intP$maxAF=ifelse(intP$pos1_maf1012>intP$pos2_maf1012, intP$pos1_maf1012,intP$pos2_maf1012)
#plot(intP$minAF, abs(intP$int.coef), pch=21, cex=.5)
#plot(intP$maxAF, abs(intP$int.coef), pch=21, cex=.5)
#
## sample from main effect QTL
#rint=cbind(sample(r2$maf1012), sample(r2$maf1012))
#r2s=split(r2, paste(r2$trait, r2$cross))
#intPs=split(intP, paste(intP$trait, intP$cross))
#
#expected=matrix(0,1000,3)
#for(i in 1:1000) {
#    print(i)
#    rsamp=list()
#    for(tt in names(intPs)) {
#        #print(tt)
#       rs=r2s[[tt]]
#       all.possible=combn(rs$maf1012, 2)
#       nint=nrow(intPs[[tt]])
#       rsamp[[tt]]=all.possible[,sample(1:ncol(all.possible), nint)]
#    }
#    rint=t(do.call('cbind', rsamp))
#    expected[i,]=c(sum(rint[,1]<.01 & rint[,2]<.01), sum(rint[,1]<.01 | rint[,2]<.01), sum(rint[,1]>.01 &  rint[,2]>.01))
#}
#
#E=expected/rowSums(expected)
#mE=melt(E)
#stripchart(mE$value*1236 ~ mE$Var2, vertical=T, method='jitter', group.names=c('both rare', 'one rare', 'both common' ), 
#           ylab='# of QTL-QTL interactions', pch=21, cex=.75 , col='#00000022')
#points(c(1,2,3), observed*1236, col='red',cex=3, pch=20)
#
#expected=expected/sum(expected)
#observed=c(sum(intP$pos1_maf1012<.01 & intP$pos2_maf1012<.01), 
#sum(intP$pos1_maf1012<.01 | intP$pos2_maf1012<.01 ),
#sum(intP$pos1_maf1012>.01 & intP$pos2_maf1012>.01) )
#observed=observed/sum(observed)
#
#rmin=apply(,1, min)
#fisher.test(rbind(table(rmin<.01), table(intP$minAF<.01)))
#ggplot(intP, aes(x=pos1_maf1012, y=pos2_maf1012))+geom_density_2d()+geom_point()
#x11()
#plot(intP$minAF, abs(intP$int.coef), pch=21, cex=.5)
#plot(intP$maxAF, abs(intP$int.coef), pch=21, cex=.5)
#
#plot(density(intP$minAF), col='red', lwd=2)
#points(density(rmin), color='black', type='l', lwd=2)
#legend('topright', c('expected lowest allele frequency', 'observed lowest allele frequency'), 
#       fill=c('black', 'red'))
#apply(cbind(sample(r2$maf1012), sample(r2$maf1012)),1, min)
#
#
#library(plot3D)
#scatter3D(intP$pos1_maf1012, intP$pos2_maf1012, abs(intP$int.coef))
#plot3d(intP$pos1_maf1012, intP$pos2_maf1012, abs(intP$int.coef), col=intP$binned)
#
#
#save(intP, file='~/Desktop/int.RData')
#
#,z=abs(int.coef), color=binned)
#
#
#stat_bin2d(bins=100)
#d+stat_density_2d(geom = "point", aes(size=stat(density)))
#
#
#, n = 20, contour = FALSE)
#
#    stat_density_2d(aes(fill = binned))
#
#    scale_color_continuous2(low="white",mid="blue",high="red", limits=c(0,.4))+
#    geom_point()
#+
#    geom_density_2d()
#
#
#
#geom_tile()+
#
#
#
#\stat_density_2d(aes(fill= ..density..),geom='tile')
#
#geom_contour()+facet_wrap(~trait)
#
#with(intP, {plot_ly(x=(pos1_maf1012), y=(pos2_maf1012), z=abs(int.coef), sizes=c(.5,1), type="scatter3d", mode="markers", color=abs(int.coef)) })


  # sig.joint.markers=parents.list[[cross.name]]$marker.name[ na.omit(match(jPs[[tt]]$fscan.markers[jPs[[tt]]$q<.05], parents.list[[cross.name]]$marker.name))]
        # bad.markers=duplicated(sig.joint.markers)
        # if(sum(bad.markers)>0) {
        #     sig.joint.markers=sig.joint.markers[-which(bad.markers)]
        # }
        # dff=data.frame(g.s[,sig.joint.markers])
        # fC=findCorrelation(cor(dff), cutoff=.99)
        # if(length(fC)>0) {
        #     sig.joint.markers=sig.joint.markers[-fC]
        #     dff=data.frame(g.s[,sig.joint.markers])
        # }
        # qmodel=lm(mPhenos[,tt]~.-1, data=dff)
        # yr=residuals(qmodel)
        # aov.a = anova(qmodel)
        # tssq  = sum(aov.a[,2])
        # a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
        # ps=drop1(qmodel, test='F')[-1,6]
        # jointPeakEffects[[tt]][[cross.name]]=data.frame(
        #    trait=tt,
        #    cross=cross.name,
        #    peaks=sig.joint.markers,
        #    betas=as.vector(coef(qmodel)),
        #    vexp=a.effs,
        #    p = ps,
        #    maf1012=iseq.freqs$maf1012[match(sig.joint.markers, iseq.freqs$marker)],
        #    # iseq.freqs[,8] is the allele frequency of the alternate allele 
        #    alt012=iseq.freqs[,8][match(sig.joint.markers, iseq.freqs$marker)],
        #    cross.cnt=cross.count.lookup[match(sig.joint.markers, rownames(cross.count.lookup)),'values'],
        #    stringsAsFactors=F)
        #   print(jointPeakEffects[[tt]][[cross.name]])


