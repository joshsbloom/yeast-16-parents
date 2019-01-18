# collapse phenotypes across the crosses
pe=do.call('rbind', (pheno_extracted))
# a factor object describing which segregant belongs to which cross
cF=as.factor(rep(names(pheno_extracted), sapply(pheno_extracted, nrow)))

# extract genotypes and process for joint variance component analysis 
joint.genotypes=list()
joint.genotypes.scaled=list()
joint.cross.cnt=list()
for(chrom in unique.chrs) {
    print(chrom)
    #chrom='chrI'
    plist.chr=sapply(parents.list, function(x) x[x$chr==chrom,])
    chr.markers=sapply(parents.list, function(x){ return(x$marker.name.n[x$chr==chrom]) } )
    seg.chr.phenos=lapply(pheno.resids, function(x) { x[[chrom]] })
    seg.chr.n=lapply(seg.chr.phenos, nrow)
    seg.chr.phenos.scaled=lapply(seg.chr.phenos, scale)

    cmm=sapply(parents.list, function(x){ return(x$marker.name[x$chr==chrom]) } )

    seg.chr.genos=mapply(function(cmm,seg.r) { t(seg.r[cmm,])}, cmm=chr.markers, seg.r=seg.recoded)
    seg.chr.genos=mapply(function(cmm,seg.r) { colnames(seg.r)=cmm; return(seg.r);}, cmm=cmm, seg.r=seg.chr.genos)

    seg.chr.genos.scaled=lapply(seg.chr.genos, scale)
    #rename seg.chr.genos.scaled

    mpmarkers=melt(cmm)
    mpmarkers[,1]=as.character(mpmarkers[,1])
    names(mpmarkers)=c('marker', 'cross')

    mpmarkers=separate(mpmarkers, marker, c('chr', 'pos', 'ref', 'alt'), sep='_', convert=T, remove=F)
    mpmarkers=mpmarkers[order(mpmarkers$pos),]

    ccnt=split(mpmarkers$cross, mpmarkers$marker)
    crosses.per.marker=sapply(ccnt, function(x) unique(x))
    cross.cnt.per.marker=sapply(crosses.per.marker, length)
    #seg.rare=names(cross.cnt.per.marker[cross.cnt.per.marker<3])

    mm=mpmarkers[match(names(crosses.per.marker), mpmarkers$marker),-6]
    mm=mm[order(mm$pos),]

    # mapping cross marker index to joint marker table  
    mupos=lapply(plist.chr, function(x) match(x$marker.name, mm$marker))
    # mapping cross marker to joint marker table name 
    mucpos = mapply(function(m, p) { p$marker.name[ match(mm$marker[m], p$marker.name)] }, m=mupos, p=plist.chr )

    big.mat=matrix(0, sum(sapply(seg.chr.phenos.scaled,nrow)), nrow(mm))
    rownames(big.mat)=as.vector(unlist((sapply(seg.chr.phenos.scaled,rownames))))
    colnames(big.mat)=mm$marker
    for(cc in names(seg.chr.genos.scaled)) {
        x=(seg.chr.genos[[cc]]*2)-1
        big.mat[rownames(x), colnames(x)]=x
    }
    bmr.scaled=apply(big.mat,2, function(x) {
                goodx=x!=0
                y=x[goodx]
                y=scale(y)
                x[goodx]=y
                return(x)
                })
    cross.cnt.per.marker=cross.cnt.per.marker[colnames(big.mat)]
    joint.genotypes[[chrom]]=big.mat
    joint.genotypes.scaled[[chrom]]=bmr.scaled
    joint.cross.cnt[[chrom]]=cross.cnt.per.marker
}
#save(joint.genotypes.scaled, file='/data/rrv2/genotyping/RData/joint.genotypes.scaled.RData')
#save(joint.cross.cnt, file='/data/rrv2/genotyping/RData/joint.cross.cnt.RData')
load('/data/rrv2/genotyping/RData/joint.genotypes.scaled.RData')
load('/data/rrv2/genotyping/RData/joint.cross.cnt.RData')

# whole genome additive covariance
#GG=lapply(joint.genotypes.scaled, tcrossprod)
#GG=Reduce('+', GG)/(sum(sapply(joint.cross.cnt, length))-1)
#GG=do.call(lapply(joint.genotypes.scaled, tcrossprod), sum)/(sum(sapply(joint.cross.cnt, length))-1)

#Construct private and non-private genetic covariance matrices based on data for 16 parents only 
rseg.cnt=lapply(joint.cross.cnt, function(x) which(x<3))
cseg.cnt=lapply(joint.cross.cnt, function(x) which(x>2))

# Gower's centered matrix 
#https://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-43-1#Sec2
# or equation 5 in this paper https://www.nature.com/articles/ng.548#methods

# two-component variance component analysis for variants partitioned by private or non-private from the 16 crosses--------------
    Grr=mapply(function(x,y) { tcrossprod(x[,y]) }, joint.genotypes.scaled, rseg.cnt, SIMPLIFY=FALSE)
    Gr=Grr[[1]]
    for(i in 2:length(Grr)){
        print(i)
        Gr=Gr+Grr[[i]]
    }
    Gr=Gr/(sum(sapply(rseg.cnt, length))-1)
    rm(Grr)
    Gr=Gr/mean(diag(Gr))

    Gcc=mapply(function(x,y) { tcrossprod(x[,y]) }, joint.genotypes.scaled, cseg.cnt, SIMPLIFY=FALSE)
    Gc=Gcc[[1]]
    for(i in 2:length(Gcc)){
        print(i)
        Gc=Gc+Gcc[[i]]
    }
    Gc=Gc/(sum(sapply(cseg.cnt, length))-1)
    rm(Gcc)
    Gc=Gc/mean(diag(Gc))

    #save(Gr,file='/data/rrv2/genotyping/RData/full_panel_private_additive_covariance_Gr.RData')
    #save(Gc,file='/data/rrv2/genotyping/RData/full_panel_non_private_additive_covariance_Gc.RData')

    joint.VCs=list()
    for(pheno in colnames(pe)) {
        joint.VCs[[pheno]]=extractVarCompResultsR(regress(pe[,pheno]~cF,~Gr+Gc, verbose=T))
    }
    #save(joint.VCs, file='/data/rrv2/genotyping/RData/joint.VCs.private_vs_nonprivate_16_parents_only.RData')
#---------------------------------------------------------------------------------------------------------------------------------
load('/data/rrv2/genotyping/RData/joint.VCs.private_vs_nonprivate_16_parents_only.RData')


# 7-component variance component analysis for variants partitioned by allele frequencies in Joseph's panel -------------------------
    load('/data/rrv2/genotyping/RData/iseq.freqs.RData')

    ibins=split(iseq.freqs, iseq.freqs$bins)
    RKHS.af=list()
    for(i in names(ibins)) {
        print(i)
        mcnt=sapply(joint.genotypes.scaled, function(x) sum(colnames(x) %in% ibins[[i]]$marker))
        print(sum(mcnt))
        jsub=lapply(joint.genotypes.scaled, function(x) tcrossprod(x[,colnames(x) %in% ibins[[i]]$marker]))
        jsub=Reduce('+',jsub)/(sum(mcnt)-1)
        jsub=jsub/mean(diag(jsub))
        RKHS.af[[i]]=jsub
    }
    # might be faster to recreate then load
    #save(RKHS.af,file='/data/rrv2/genotyping/RData/7_covariance_matrices_binned_by_JS_allele_frequency.RData')

    # some diagnostic visualizations
    jcc=data.frame(cnt=do.call('c', joint.cross.cnt))
    jcc$marker.name=do.call('c', sapply(joint.cross.cnt, names)
    imerge=merge(iseq.freqs, jcc, by.x='marker', by.y='marker.name')
    imerge$cross.count=2*round((imerge$cnt+.1)/2)

    # visualize allele frequency of variants in JS1012  in RR
    ggplot(imerge, aes(x=maf1012))+facet_wrap(~ cross.count)+geom_histogram(binwidth=.005)
    ggplot(imerge, aes(x=maf1012))+facet_wrap(~cross.count)+geom_histogram(binwidth=.005)+xlab('minor allele frequency of rr variants in 1012 yeast panel')+
        ggtitle('allele frequencies split by number of crosses a variant segregates in')+theme_bw()+ylim(0,5e3)
    #ggsave(file='/data/rrv2/genotyping/maf_by_number_of_crosses_variant_segregates_in.pdf')

    joint.VCs.af=list()
    for(pheno in colnames(pe)[21:41]) {
        print(pheno)
        joint.VCs.af[[pheno]]=extractVarCompResultsR(
                regress(pe[,pheno]~cF,~RKHS.af[[1]]+RKHS.af[[2]]+RKHS.af[[3]]+RKHS.af[[4]]+RKHS.af[[5]]+RKHS.af[[6]]+RKHS.af[[7]], 
                                           verbose=T, pos=rep(T,8),
                                           tol=.1, maxcyc=20) )
    }
    #save(joint.VCs.af, file='/data/rrv2/genotyping/RData/joint.VCs.7_JS_allele_frequency_bins.RData')
#------------------------------------------------------------------------------------------------------------------------------------------
load('/data/rrv2/genotyping/RData/joint.VCs.7_JS_allele_frequency_bins.RData')


# 2-component variance component analysis for variants partitioned by allele frequencies in Joseph's panel -------------------------
    ibins2=split(iseq.freqs, iseq.freqs$bins2)
    RKHS.af2=list()
    for(i in names(ibins2)) {
        print(i)
        mcnt=sapply(joint.genotypes.scaled, function(x) sum(colnames(x) %in% ibins2[[i]]$marker))
        print(sum(mcnt))
        jsub=lapply(joint.genotypes.scaled, function(x) tcrossprod(x[,colnames(x) %in% ibins2[[i]]$marker]))
        jsub=Reduce('+',jsub)/(sum(mcnt)-1)
        jsub=jsub/mean(diag(jsub))
        RKHS.af2[[i]]=jsub
    }
    names(RKHS.af2)=c('<.01', '>.01')
    
    joint.VCs.af2=list()
    for(pheno in colnames(pe)) {
        starthere=as.vector(c(sum(joint.VCs.af[[pheno]]$sigma[1:2]),
                    sum(joint.VCs.af[[pheno]]$sigma[3:7]), 
                        joint.VCs.af[[pheno]]$sigma[8]))
        if(starthere[1]<.01) { starthere[1]=starthere[2] }
        joint.VCs.af2[[pheno]]=extractVarCompResultsR(
                        regress(pe[,pheno]~cF,
                                ~RKHS.af2[[1]]+RKHS.af2[[2]], 
                                start=starthere,
                                pos=rep(T,3),
                                verbose=T))
    }
    #save(joint.VCs.af2, file='/data/rrv2/genotyping/RData/joint.VCs.2_JS_allele_frequency_bins.RData')

#---------------------------------------------------------------------------------------------------------------------------------



# 7-component variance component analysis for variants partitioned by allele frequencies in Joseph's panel, but only using private variants from RR panel
    rseg.markers=do.call('c', sapply(rseg.cnt, names))
    RKHS.af.rare=list()
    for(i in names(ibins)) {
        print(i)
        mcnt=sapply(joint.genotypes.scaled, function(x) sum(colnames(x) %in% ibins[[i]]$marker & colnames(x) %in% rseg.markers ))
        print(sum(mcnt))
        jsub=lapply(joint.genotypes.scaled, function(x) tcrossprod(x[,colnames(x) %in% ibins[[i]]$marker & colnames(x) %in% rseg.markers ]))
        jsub=Reduce('+',jsub)/(sum(mcnt)-1)
        jsub=jsub/mean(diag(jsub))
        RKHS.af.rare[[i]]=jsub
    }

    joint.VCs.af.private=list()
    for(pheno in colnames(pe)[38:41]) {
        print(pheno)
        #startvec=as.vector(joint.VCs.af[[pheno]]$sigma)
        #startvec=c(1+rep(sum(joint.VCs.af[[pheno]]$sigma[1:7])/7, 7), joint.VCs.af[[pheno]]$sigma[8])
        joint.VCs.af.private[[pheno]]=extractVarCompResultsR(
                                   regress(pe[,pheno]~cF,~RKHS.af.rare[[1]]+RKHS.af.rare[[2]]+RKHS.af.rare[[3]]+RKHS.af.rare[[4]]+RKHS.af.rare[[5]]+RKHS.af.rare[[6]]+RKHS.af.rare[[7]], 
                                           verbose=T, pos=rep(T,8),
                                           tol=.1, maxcyc=20) )
    }
    #save(joint.VCs.af.private, file='/data/rrv2/genotyping/RData/joint.VCs.7_JS_allele_frequency_bins_only_private_in_RR.RData')
load('/data/rrv2/genotyping/RData/joint.VCs.7_JS_allele_frequency_bins_only_private_in_RR.RData')
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

# 2-component variance component analysis for variants partitioned by allele frequencies in Joseph's panel, but only using private variants from RR panel
    rseg.markers=do.call('c', sapply(rseg.cnt, names))
    RKHS.af.rare2=list()
    for(i in names(ibins2)) {
        print(i)
        mcnt=sapply(joint.genotypes.scaled, function(x) sum(colnames(x) %in% ibins2[[i]]$marker & colnames(x) %in% rseg.markers ))
        print(sum(mcnt))
        jsub=lapply(joint.genotypes.scaled, function(x) tcrossprod(x[,colnames(x) %in% ibins2[[i]]$marker & colnames(x) %in% rseg.markers ]))
        jsub=Reduce('+',jsub)/(sum(mcnt)-1)
        jsub=jsub/mean(diag(jsub))
        RKHS.af.rare2[[i]]=jsub
    }
    names(RKHS.af.rare2)=c('<.01', '>.01')

    joint.VCs.af.private2=list()
    for(pheno in colnames(pe)) {
        print(pheno)
        #startvec=as.vector(joint.VCs.af[[pheno]]$sigma)
        #startvec=c(1+rep(sum(joint.VCs.af[[pheno]]$sigma[1:7])/7, 7), joint.VCs.af[[pheno]]$sigma[8])
        starthere=as.vector(joint.VCs.af2[[pheno]]$sigma[1:3])
        #if(starthere[1]<.01) { starthere[1]=starthere[2] }
        joint.VCs.af.private2[[pheno]]=extractVarCompResultsR(
                                   regress(pe[,pheno]~cF,~RKHS.af.rare2[[1]]+RKHS.af.rare2[[2]], 
                                           verbose=T, pos=rep(T,3), start=starthere)
                                            )
    }
    #save(joint.VCs.af.private2, file='/data/rrv2/genotyping/RData/joint.VCs.2_JS_allele_frequency_bins_only_private_in_RR.RData')
    #load('/data/rrv2/genotyping/RData/joint.VCs.2_JS_allele_frequency_bins_only_private_in_RR.RData')

# Visualize results from two-component variance component analysis for variants partitioned by private or non-private from the 16 crosses--------------

# private vs non-private within panel only     
JCV=extractVC_long_format(joint.VCs)
levels(JCV$Component)=c('E', 'Not-private', 'Private')

# split by JS panel allele frequencies (7-components)
JCV.af7=extractVC_long_format(joint.VCs.af)
levels(JCV.af7$Component)=c('E', rev(names(ibins))) # 'Not-private', 'Private')
#JCV.af$Trait=factor(JCV.af$Trait, levels=names(JC.af.sigma.norm[3,])[order(JC.af.sigma.norm[7,]+JC.af.sigma.norm[8,])])
JCV.af7=JCV.af7[!JCV.af7$Component=='E',]

# split by JS panel allele frequencies, private variants only (7-components)
JCV.af.private7=extractVC_long_format(joint.VCs.af.private)
levels(JCV.af.private7$Component)=c('E', rev(names(ibins))) # 'Not-private', 'Private')
#JCV.af.private$Trait=factor(JCV.af.private$Trait, levels=names(JC.af.private.sigma.norm[3,])[order(JC.af.sigma.norm[7,]+JC.af.sigma.norm[8,])])
JCV.af.private7=JCV.af.private7[!JCV.af.private7$Component=='E',]

# split by JS panel allele frequencies (2-components)
JCV.af2=extractVC_long_format(joint.VCs.af2)
levels(JCV.af2$Component)=c('E', 'Not-private', 'Private')
#levels(JCV.af2$Component)=c('E', rev(names(ibins2))) # 'Not-private', 'Private')
#JCV.af$Trait=factor(JCV.af$Trait, levels=names(JC.af.sigma.norm[3,])[order(JC.af.sigma.norm[7,]+JC.af.sigma.norm[8,])])
#JCV.af2=JCV.af2[!JCV.af2$Component=='E',]

JCV.af2.private=extractVC_long_format(joint.VCs.af.private2)
levels(JCV.af2.private$Component)=c('E', 'Not-private', 'Private')



# comparison of weighted sum of additive variance components and two-component joint model, private vs non private within cross:
cross.vars=sapply(split(data.frame(pe[,-c(37,38)]), cF), function(x) apply(x,2,var))
cross.cnts.table=table(cF)
Wc.additive.table=do.call('rbind', lapply(VC_A.list, function(x) x[-c(37,38),]))
Wc.additive.table$cross=rep(names(VC_A.list), sapply(VC_A.list, function(x) nrow(x[-c(37,38),])))
names(Wc.additive.table)[1]='additive'
sw=split(Wc.additive.table, Wc.additive.table$trait)
anorm3=rep(0,length(levels(JCV$Trait)))
for(i in 1:length(sw)) {
      anorm3[i]=sum(sw[[i]]$additive*cross.cnts.table[sw[[i]]$cross]*cross.vars[i,sw[[i]]$cross])/sum(cross.cnts.table[sw[[i]]$cross]*cross.vars[i,sw[[i]]$cross])
    
}
# evaluate whether joint estimate of heritability is inflated 
Joint.additive=sapply(split(JCV, JCV$Trait), function(x) sum(x$fraction_of_variance[2:3]))
Joint.additive.se=sapply(split(JCV, JCV$Trait), function(x) sqrt(sum(x$se^2)))

par(oma=c(2,2,2,2))
plot(anorm3, Joint.additive, xlim=c(0,1), ylim=c(0,1), ylab=bquote(paste(sigma[a]^2, ' from joint model')), 
    xlab=expression(Sigma(sigma[a[wc]]^2*sigma[p]^2*n)/Sigma(sigma[p]^2*n)), cex=1, pch=19)
    segments(anorm3, Joint.additive, anorm3, Joint.additive+Joint.additive.se,lty=1)
    segments(anorm3, Joint.additive, anorm3, Joint.additive-Joint.additive.se,lty=1)
    abline(0,1)
legend('bottomright', legend=paste(c('b=', 'x='), round(coef(lm(Joint.additive~anorm3)),3)), cex=2)

# compare estimate of variance coming from private variants only
# x-axis is from classifying variants as private or non-private from 16 panel only
# y-axis is from classifying variants from JS panel with MAF <.01 or not
plot(JCV$fraction_of_variance[JCV$Component=='Private'], JCV.af2$fraction_of_variance[JCV.af2$Component=='Private'],
        xlab='within-panel allele frequencies, private', 
        ylab='JS panel maf<.01', xlim=c(0,1), ylim=c(0,1), 
        main='variance due to rare variants')
abline(0,1)
# compare estimate of variance coming from private variants only
plot(JCV.af2.private$fraction_of_variance[JCV.af2.private=='Private'], JCV.af2$fraction_of_variance[JCV.af2$Component=='Private'],
        xlab='JS panel maf<.01 and only using private variants from within 16 panel', 
        ylab='JS panel maf<.01', xlim=c(0,1), ylim=c(0,1), 
        main='variance due to rare variants')
abline(0,1)

# compare ration of variance coming from private vs common variants only
plot(JCV.af2.private$fraction_of_variance[JCV.af2.private$Component=='Private']/JCV.af2.private$fraction_of_variance[JCV.af2.private$Component=='Not-private'], 
     JCV.af2$fraction_of_variance[JCV.af2$Component=='Private']/JCV.af2$fraction_of_variance[JCV.af2$Component=='Not-private'],
        xlab='JS panel maf<.01 and only using private variants from within 16 panel', 
        ylab='JS panel maf<.01', 
        main='ratio of private/non-private')
abline(0,1)

c1=ggplot(iseq.freqs, aes(x=maf1012,fill=maf1012<.01))+geom_histogram(binwidth=.0025)+
    scale_fill_manual(values =c('lightgrey', 'lightblue'))+
    xlab('minor allele frequency of segregating variants from 16 round-robin parents in 1011 yeast panel')+
    ylab('number of variants')+
    theme_bw()+guides(fill=FALSE)
ggsave(file='~/Dropbox/Lab Meeting - Presentations/112918/1011_af.png',width=7,height=5)

js.afs=fread('/data/rrv2/1002genomes/all.freq.frq')
js.afs$maf1012=ifelse(js.afs[[6]]<js.afs[[8]], js.afs[[6]], js.afs[[8]])
d1=ggplot(js.afs, aes(x=maf1012))+geom_histogram(binwidth=.001)+
           xlab('minor allele frequency of all biallelic variants in 1012 yeast panel')+ggtitle('allele frequencies of all bialleic variants from yeast panel')+
           theme_bw()
d2=ggplot(js.afs, aes(x=maf1012))+geom_histogram(binwidth=.001)+
           xlab('minor allele frequency of all biallelic variants in 1012 yeast panel')+ggtitle('Zoomed - allele frequencies of all bialleic variants from yeast panel')+
           theme_bw()+coord_cartesian(ylim=c(0, 1e4))

ggarrange(c1,d1, d2, ncol=1, nrow=3, labels=c('a','b', 'c') )

# how to reorder barplot 
#JCV$Trait=factor(JCV$Trait, levels=names(JC.sigma.norm[3,])[order(JC.sigma.norm[3,])])
# Visualize results from 2- component variance component analysis for variants partitioned by JS allele frequencies --------------

xx=(JCV.af2[JCV.af2$Component=='Private',])
JCV$Trait=factor(JCV$Trait, levels=as.character(xx$Trait[order(xx$fraction_of_variance)]))
JCV.af2$Trait=factor(JCV.af2$Trait, levels=as.character(xx$Trait[order(xx$fraction_of_variance)]))

aa=ggplot(JCV, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    scale_fill_manual(values =c('white', 'lightgrey', 'lightblue')) +theme_bw()+ylim(0,1)+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))+
    theme(axis.text.x=element_text(angle=70,hjust=1))
ggsave(file='~/Dropbox/Lab Meeting - Presentations/112918/vc_analysis_joint_within_panel.png',width=7,height=5)



bb=ggplot(JCV.af2, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    scale_fill_manual(values =c('white', 'lightgrey', 'lightblue')) +theme_bw()+ylim(0,1)+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))+
    theme(axis.text.x=element_text(angle=70,hjust=1))
ggsave(file='~/Dropbox/Lab Meeting - Presentations/112918/vc_analysis_joint.png',width=7,height=5)

b1=ggplot(JCV.af2.private, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    scale_fill_manual(values =c('white', 'lightgrey', 'lightblue')) +theme_bw()+ylim(0,1)+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))+
    theme(axis.text.x=element_text(angle=70,hjust=1))

ggarrange(aa,bb, ncol=1, nrow=2, heights=c(2,2),  labels=c('a','b'))

ggarrange(aa,bb,b1, ncol=1, nrow=3, heights=c(2,2,2),  labels=c('a','b','c'))

ggarrange(aa,bb,b1, ncol=1, nrow=3, heights=c(2,2,2),  labels=c('a','b','c'))
x=JCV.af2[JCV.af2$Component=='Private','fraction_of_variance']
y=JCV.af7[JCV.af7$Component=="[0.00000,0.00495)",'fraction_of_variance']+JCV.af7[JCV.af7$Component=="[0.00495,0.01088)",'fraction_of_variance']



# Visualize results from 7-component variance component analysis for variants partitioned by JS allele frequencies --------------
cc=ggplot(JCV.af7, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    theme_bw()+ylim(0,1)+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))+
    theme(axis.text.x=element_text(angle=70,hjust=1))+ scale_fill_manual(name = "allele frequency", values=tim.colors(7))+
    ggtitle('split by 1012 genomes allele frequencies')

# Visualize results from 7-component variance component analysis for variants partitioned by JS allele frequencies, private RR variants only --------------
dd=ggplot(JCV.af.private7, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    theme_bw()+ylim(0,1)+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))+
    theme(axis.text.x=element_text(angle=70,hjust=1))+ scale_fill_manual(name = "allele frequency", values=tim.colors(7))+
    ggtitle('split by 1012 genomes allele frequencies - only private variants')
ggarrange(cc,dd, ncol=1, nrow=2, heights=c(2,2),  labels=c('a','b'))


c1=ggplot(iseq.freqs, aes(x=maf1012,fill=iseq.freqs$bins))+geom_histogram(binwidth=.001)+scale_fill_manual(name='allele frequency', values=rev(tim.colors(7)))+
    xlab('minor allele frequency of rr16 variants in 1012 yeast panel')+ggtitle('allele frequencies of rr16 variants in yeast panel')+theme_bw()+guides(fill=FALSE)
js.afs=fread('/data/rrv2/1002genomes/all.freq.frq')
js.afs$maf1012=ifelse(js.afs[[6]]<js.afs[[8]], js.afs[[6]], js.afs[[8]])
d1=ggplot(js.afs, aes(x=maf1012))+geom_histogram(binwidth=.001)+
           xlab('minor allele frequency of all biallelic variants in 1012 yeast panel')+ggtitle('allele frequencies of all bialleic variants from yeast panel')+
           theme_bw()
d2=ggplot(js.afs, aes(x=maf1012))+geom_histogram(binwidth=.001)+
           xlab('minor allele frequency of all biallelic variants in 1012 yeast panel')+ggtitle('Zoomed - allele frequencies of all bialleic variants from yeast panel')+
           theme_bw()+coord_cartesian(ylim=c(0, 1e4))

ggarrange(c1,d1, d2, ncol=1, nrow=3, labels=c('a','b', 'c') )




