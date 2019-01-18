library(VariantAnnotation)
library(vcfR)
library(qtl)
library(qtl2)
library(ASMap)
library(S4Vectors)

unique.chrs=paste0('chr', as.roman(1:16))

chr.lengths=list(
chrI=    230218,
chrII=   813184,
chrIII=  316620,
chrIV=   1531933,
chrV=    576874,
chrVI=   270161,
chrVII=  1090940,
chrVIII= 562643,
chrIX=   439888,
chrX=    745751,
chrXI=   666816,
chrXII=  1078177,
chrXIII= 924431,
chrXIV=  784333,
chrXV=   1091291,
chrXVI=  948066)

gcoord.key=cumsum(c(0,chr.lengths))[-17]
names(gcoord.key)=names(chr.lengths)

# global variables
blocked.regions=list(
        chrI=rbind(c(1,23500),
                   c(219009,230218)),
        chrII=rbind(c(1,10000),
                    c(800006,813184)),
        chrIII=rbind(c(1,6000),
                     c(303503,316620)),
        chrIV=rbind(c(1,19000),
                    c(1522015,1531933)),
        chrV=rbind(c(1,10000),
                   c(19700,23500),
                   c(565005,576874)),
        chrVI=rbind(c(1,16800)),
        chrVII=rbind(c(1,25000),
                     c(1075993,1090940)),
        chrVIII=rbind(c(1,13002),
                      c(211006,217998),
                      c(525000,562543)),
        chrIX=rbind(c(1,28850),
                    c(425003,439888)),
        chrX=rbind(c(1,26000),
                   c(726009,745751)),
        chrXI=rbind(c(1,3399),
                    c(644362,666816)),
        chrXII=rbind(c(1,14000),
                     c(449999,490999),
                     c(1059002,1078177)),
        chrXIII=rbind(c(1,9000),
                      c(915002,924431)),
        chrXIV=rbind(c(1,18000),
                     c(772000,784333)),
        chrXV=rbind(c(1,33001),
                    c(1070002,1091291)),
        chrXVI=rbind(c(1,26400),
                     c(930004,948066)))

blocked=DataFrame(seqnames=rep(names(blocked.regions), sapply(blocked.regions, nrow)), do.call('rbind', blocked.regions))
names(blocked)[c(2,3)]=c('start', 'end')
blocked=makeGRangesFromDataFrame(blocked)


crosses=c('375', 'A', '376', 'B', '377', '393', '381', '3008', '2999', '3000' , '3001', '3049', '3003', '3004', '3043', '3028')
genome.chr.lengths=c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)	
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




recode.as.parental=function(seg.mat, parents) { 
    GT.cols=grep('GT', names(parents))
    z=rep(NA,nrow(seg.mat))
    names(z)=rownames(seg.mat)
    g.recode=apply(seg.mat, 2, function(x) {
        z[which(x==parents[,GT.cols[1]])]=1
        z[which(x==parents[,GT.cols[2]])]=2
        return(z)
    })
    return(g.recode)
}

# important, flip back (GT calls in vcf)
recode.as.allele.number=function(seg.mat, parents) {
    GT.cols=grep('GT', names(parents))
    z=rep(NA,nrow(seg.mat))
    names(z)=rownames(seg.mat)

    g.recode.back=apply(seg.mat, 2, function(x) {
        w1=which(x==1)
        w2=which(x==2)
        z[w1]=parents[,GT.cols[1]][w1]
        z[w2]=parents[,GT.cols[2]][w2]
       return(z)
    })
    return(g.recode.back)
}

#add in physical map as temporary substitute for genetic map 
fake_gmap = function(cross, delim=':|_', pfield=2, p.to.d=2200) {
    pmap=lapply(cross$geno, function(y) {
         z=as.numeric(sapply(strsplit(colnames(y$data), delim), function(x) x[pfield]))/p.to.d
         names(z)=colnames(y$data)
        #return(z)
        y$map=z  
        return(y)
    })
    cross$geno=pmap
    return(cross)
}

# function to impute marker genotypes to speed up genetic map estimation 
# modified from ASMap
quickEstJB=function (object, chr,  map.function = "kosambi",...)  {
    if (missing(chr)) { 
        chr <- names(nmar(object)) }
    imf <- switch(map.function, kosambi = imf.k, haldane = imf.h, morgan = imf.m, cf = imf.cf)
    # new code to build a temporary genetic map given physical map
    # depends on marker name being in the second field delimited as such
    nm <- nmar(object)
    for (i in chr) {
        temp <- subset(object, chr = i)
        nc   <- ncol(temp$geno[[i]]$data)
        tempa <- argmax.geno(temp, step = 0, map.function = map.function, ...)
        tempa$geno[[i]]$data <- tempa$geno[[i]]$argmax
        tempa$geno[[i]] <- tempa$geno[[i]][-3]
        esta <- est.rf(tempa)$rf
        era <- esta[cbind(2:nc, 1:(nc - 1))]
        if (class(object)[1] == "riself") 
            era <- (era/2)/(1 - era)
        object$geno[[i]]$map <- c(0, cumsum(imf(era)))
        names(object$geno[[i]]$map) <- dimnames(object$geno[[i]]$data)[[2]]
    }
    object
}

 
quickEstJB.mc=function (object, ep=.005, map.function = "kosambi",...)  {
    chrs <- names(nmar(object)) 
    imf <- switch(map.function, kosambi = imf.k, haldane = imf.h, morgan = imf.m, cf = imf.cf)
    # new code to build a temporary genetic map given physical map
    # depends on marker name being in the second field delimited as such
    nm <- nmar(object)
    mapout <- mclapply(chrs, function(i, ...){
        temp <- subset(object, chr = i)
        nc   <- ncol(temp$geno[[i]]$data)
        tempa <- argmax.geno(temp, step = 0, map.function = map.function, error.prob=ep)
        tempa$geno[[i]]$data <- tempa$geno[[i]]$argmax
        tempa$geno[[i]] <- tempa$geno[[i]][-3]
        esta <- est.rf(tempa)$rf
        era <- esta[cbind(2:nc, 1:(nc - 1))]
        if (class(object)[1] == "riself") { era <- (era/2)/(1 - era) }
        map=c(0, cumsum(imf(era)))
        names(map)=dimnames(object$geno[[i]]$data)[[2]]
        return(map)
        #object$geno[[i]]$map <- c(0, cumsum(imf(era)))
        #names(object$geno[[i]]$map) <- dimnames(object$geno[[i]]$data)[[2]]
        
    }, object, ep, mc.cores=16)
    names(mapout)=chrs
    class(mapout)='map'
    #return(mapout)
    return(replace.map(object, mapout))
}

# hack a r/qtl cross object without the annoying formatting and i/o 
# g.recode is markers by segregants
hackCrossObject=function(g.recode, parents, unique.chrs) {
    
    g.recode.df=data.frame(g.recode)
    rownames(g.recode.df)=paste0(parents$marker.name, '_', seq(1,nrow(parents)))

    g.split=lapply(split(g.recode.df, parents$chr), data.matrix)
    g.split=g.split[unique.chrs]

    g.split=lapply(g.split, function(x) list(data=t(x)) )
    g.split=lapply(g.split, function(x) { class(x)='A'; return(x); })
    pmat=data.frame(id=colnames(g.recode))

    jcross=list(geno=g.split, pheno=pmat)
    class(jcross)=c('riself', 'cross')
    attr(jcross, 'alleles')=c('A', 'B')

    return(jcross)

}

# build parent DataFrame object
buildParentsDF=function(parent.vcf.file, possible_het_ratio, parent.min.read.support, parent.gq.cutoff, blocked.regions) {

    pvcf=read.vcfR(parent.vcf.file) #, param, genome='sacCer3')
    p.gt=extract.gt(pvcf, element='GT', as.numeric=T)
    p.dp=extract.gt(pvcf, element='DP', as.numeric=T)
    p.gq=extract.gt(pvcf, element='GQ', as.numeric=T)
    p.ad=extract.gt(pvcf, element='AD')
    # reference counts
    p.ad1=masplit(p.ad, record=1, sort=0)
    # alt counts 
    p.ad2=masplit(p.ad, record=2, sort=0)

    maxcnt1=apply(cbind(p.ad1[,1], p.ad2[,1]), 1, max)
    mincnt1=apply(cbind(p.ad1[,1], p.ad2[,1]), 1, min)

    maxcnt2=apply(cbind(p.ad1[,2], p.ad2[,2]), 1, max)
    mincnt2=apply(cbind(p.ad1[,2], p.ad2[,2]), 1, min)

    # unreliable calls 
    p1_unreliable=(mincnt1/maxcnt1)>possible_het_ratio
    p2_unreliable=(mincnt2/maxcnt2)>possible_het_ratio

    parents=DataFrame(chr=getCHROM(pvcf), pos=getPOS(pvcf), qual=getQUAL(pvcf), id=getID(pvcf), ref=getREF(pvcf), alt=getALT(pvcf), GT=p.gt,DP=p.dp, GQ=p.gq,  AD1=p.ad1, AD2=p.ad2)
    #, P1_unreliable=p1_unreliable, P2_unreliable=p2_unreliable)
    is.structural.variant=grepl('simple|complex', rownames(parents))

    GT.cols=grep('GT', names(parents))
    DP.cols=grep('DP', names(parents))
    GQ.cols=grep('GQ', names(parents))

    bad.vars=is.na(parents[,GT.cols[1]]) | is.na(parents[,GT.cols[2]]) | (parents[,DP.cols[1]]<parent.min.read.support) |
          (parents[,DP.cols[2]]<parent.min.read.support)  | (parents[,GQ.cols[1]]<parent.gq.cutoff) | (parents[,GQ.cols[2]]<parent.gq.cutoff) |     p1_unreliable | p2_unreliable

    # add back the structural variants 
    bad.vars[is.structural.variant]=FALSE

    parents=parents[-which(bad.vars),]
    parents$marker.name=paste0(parents$chr, '_', parents$pos, '_', parents$ref, '_', parents$alt)
    # now reorder
    parents=do.call('rbind', split(parents, parents$chr)[names(blocked.regions)])
    return(parents)
}

# make GRanges to block out a window around structural variants
blockSVarRegions=function(parents, structural_variant_padding) {
    sv.annot=rownames(parents)[grep('simple|complex', parents$id)]
    sv.annot.mat=do.call('rbind', strsplit(sv.annot, ':|-|_'))
    sv.blocked.regions=makeGRangesFromDataFrame(DataFrame(seqnames=sv.annot.mat[,1], start=as.numeric(sv.annot.mat[,2]),end=as.numeric(sv.annot.mat[,3]) ))
    sv.blocked.regions=reduce(sv.blocked.regions)
    sv.blocked.regions=sv.blocked.regions+structual_variant_padding
    return(sv.blocked.regions)
}

#read segregant vcf and recode as parental allele codes
buildSegregantMatrix=function(segregant.vcf.file,segregant.quality.cutoff,segregant.gq.cutoff, blocked, sv.blocked.regions, parents, 
                              segregant.missing.marker.fraction,segregant.missing.by.marker.cutoff,segregant.afH, segregant.afL ) {
    # now segs may contain variants not present in parents
    svcf=read.vcfR(segregant.vcf.file) # , param, genome='sacCer3')
    segregants=DataFrame(chr=getCHROM(svcf), start=getPOS(svcf),end=getPOS(svcf)) #, qual=getQUAL(svcf), ref=getREF(svcf), alt=getALT(svcf))
    segregants=makeGRangesFromDataFrame(segregants)
    segregants$qual=getQUAL(svcf)
    segregants$ref =getREF(svcf)
    segregants$alt =getALT(svcf)

    segregants$marker.name=paste0( as.character(seqnames(segregants)), '_',start(segregants), '_', segregants$ref, '_', segregants$alt)
    s.gt=extract.gt(svcf, element='GT', as.numeric=T)
    s.gq=extract.gt(svcf, element='GQ', as.numeric=T)
    rownames(s.gt)=segregants$marker.name
    rownames(s.gq)=segregants$marker.name
        
    #apply filters
    s.gt[which(segregants$qual<segregant.quality.cutoff),]=NA
    s.gt[which(s.gq<segregant.gq.cutoff, arr.ind=T)]=NA
    s.gt[subjectHits(findOverlaps(blocked, segregants)),]=NA
    s.gt[subjectHits(findOverlaps(sv.blocked.regions, segregants)),]=NA

    # segregant marker name must exist in parent table
    s.gt=s.gt[which(segregants$marker.name %in% parents$marker.name),]

    # match all segregant markers to parent markers and fill in missing parent markers with NAs 
    s.gt.x=matrix(NA, nrow(parents), ncol(s.gt))
    rownames(s.gt.x)=parents$marker.name
    colnames(s.gt.x)=colnames(s.gt)
    s.gt.x[match(rownames(s.gt), rownames(s.gt.x)),]=s.gt
    # ------------------------------------------------------------------------------------------

    #p.ad=extract.gt(pvcf, element='AD', as.numeric=F)
    #p.ad1=masplit(p.ad, record=1)
    #p.ad2=masplit(p.ad, record=2)

    # now recode as parental alleles
    g.recode=recode.as.parental(s.gt.x, parents)
    #rle(sapply(strsplit(rownames(g.recode), '_'), function(x) x[1]))

    # stats per segregant
    #s.cnt1=apply(g.recode,2, function(x) sum(x==1, na.rm=T))
    #s.cnt2=apply(g.recode,2, function(x) sum(x==2,na.rm=T))
    s.cntNA=apply(g.recode,2, function(x) sum(is.na(x)))
    s.fracNA=s.cntNA/nrow(g.recode)
    # stats per marker
    m.cnt1=apply(g.recode,1, function(x) sum(x==1,na.rm=T))
    m.cnt2=apply(g.recode,1, function(x) sum(x==2,na.rm=T))
    m.cntNA=apply(g.recode,1, function(x) sum(is.na(x)))
    m.af = m.cnt1/(m.cnt1+m.cnt2)

    #remove segregants with extremely low coverage sequencing
    g.recode=g.recode[,-which(s.fracNA>segregant.missing.marker.fraction)]

    #NA out markers with too much 
    g.recode[which(m.cntNA>segregant.missing.by.marker.cutoff),]=NA
    g.recode[which(m.af>segregant.afH),]=NA
    g.recode[which(m.af<segregant.afL),]=NA
    return(g.recode)
}


filterPoorQualitySegs=function(jcross, error.prob, max_XO_after_preimp){
    # goal here is to remove segregants that are misbehaving
    # mimic ASMap quickEst -----------------------------------------------
    # slow, but way way faster than redoing the whole map reconstruction
    # jmap=quickEstJB(jcross, error.prob=.005)
    # jmap=quickEstJB.mc(jcross, ep=error.prob)

    #jmap=replace.map(jcross, jmap)
    jmap=argmax.geno(jcross, step=0, map.function='kosambi', error.prob=error.prob)
    jmap$geno = lapply(jmap$geno, function(x) { rownames(x$argmax)=rownames(x$data); return(x) } )

    # replace data with viterbi calls 
    jmap.am=jmap
    jmap.am$geno=lapply(jmap.am$geno, function(x) return(list(data=x$argmax, map=x$map)))
    jmap.am$geno=lapply(jmap.am$geno, function(x) { class(x)='A'; return(x); })


    #jmap.am.error.lod=calc.errorlod(jmap.am, error.prob=.005, map.function='kosambi')
    jmap.cXO=countXO(jmap.am, bychr=TRUE)
    median.cXO=apply(jmap.cXO,1,median)
    # set to #450 
    # remove individuals with too many crossovers
    #jmap=subset(jmap, ind=-which(jmap.cXO>max_XO_after_preimp))
    #jmap.am= subset(jmap.am, ind=-which(jmap.cXO>500))

    jcross=subset(jcross, ind=-which(median.cXO>max_XO_after_preimp))
    return(jcross)
}

hackScanOneObject=function(cross) {
    gtf=data.frame(chr=rep(names(cross$geno), sapply(cross$geno, function(x) ncol(x$data))), pos=unlist((pull.map(cross))))
    class(gtf)=c('scanone', 'data.frame')
    rownames(gtf)=gsub('^.*\\.', '', rownames(gtf))

    # add some allele frequency data
    g.xs=lapply(cross$geno, function(x) x$argmax)
    g.xs.af=lapply(g.xs, function(g.x) apply(g.x, 2, function(x) sum(x==1, na.rm=T)/(sum(!is.na(x))) ))
    g.xs.af.smooth=lapply(g.xs.af, function(g.x) predict(smooth.spline(g.x))$y)

    gtf$AF=do.call('c', g.xs.af)
    gtf$AF.smooth=do.call('c', g.xs.af.smooth)
    gtf$AF.resid=abs(gtf$AF-gtf$AF.smooth)

    return(gtf)
}
