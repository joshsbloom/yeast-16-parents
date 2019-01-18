source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')

library(foreach)
library(doMC)
registerDoMC(cores=70)


# filtering parameters -------------------------
parent.min.read.support=3
parent.gq.cutoff=98

segregant.quality.cutoff=100
segregant.gq.cutoff=30

segregant.missing.by.marker.cutoff=960*.9
segregant.missing.marker.fraction = .98
segregant.afH= .99
segregant.afL= .01
possible_het_ratio= 0.1
structual_variant_padding= 200
error.prob1=0.1
error.prob2=0.01 #075
max_XO_after_preimp=15
# -----------------------------------------------

cross.list=list()
parents.list=list()

for(cross.name in crosses) {

    print(cross.name)
    parent.vcf.file=   paste0('/data/rrv2/genotyping/parents/cross_vcf/', cross.name,'_w_svar.vcf.recode.vcf.gz')
    segregant.vcf.file=paste0('/data/rrv2/genotyping/segregants/vcfs/',cross.name, '.vcf.recode.vcf.gz')


    # build parent DataFrame variant stucture (structural variant names are in the id field)
    parents=buildParentsDF(parent.vcf.file, possible_het_ratio, parent.min.read.support, parent.gq.cutoff, blocked.regions)

    parents.list[[cross.name]]=parents

    #create an additional set of blocked variants based on called structural variants
    sv.blocked.regions=blockSVarRegions(parents, structural_variant_padding)

    g.recode=buildSegregantMatrix(segregant.vcf.file,segregant.quality.cutoff,segregant.gq.cutoff, blocked, sv.blocked.regions, parents, 
                                  segregant.missing.marker.fraction,segregant.missing.by.marker.cutoff,segregant.afH, segregant.afL )

    m.cnt1=apply(g.recode,1, function(x) sum(x==1,na.rm=T))
    m.cnt2=apply(g.recode,1, function(x) sum(x==2,na.rm=T))
    m.af = m.cnt1/(m.cnt1+m.cnt2)
    
    # careful, this will fail if first 96 from cross 3008 are filtered out in
    if(cross.name=='3008') { #3008_G1_96_vertical_flip
        snames = colnames(g.recode)
        pot.flip.ind=grep('G1', snames)
        snames[pot.flip.ind]= snames[rev(pot.flip.ind)]  
        colnames(g.recode)=snames
        g.recode=g.recode[,order(colnames(g.recode))]
    }
    if(cross.name=='B') {
        snames = colnames(g.recode)
        snames=gsub('^B_', 'B',snames)
        colnames(g.recode)=snames
    }

    #---------------------------------------------------------------------

    # build an r/qtl cross object
    cross=hackCrossObject(g.recode,parents,unique.chrs)

    # set physical map as genetic map to speed up map estimation and hmm for segregant genotype reconstruction
    cross = fake_gmap(cross)

    # remove segregants with too many crossovers
    cross =filterPoorQualitySegs(cross, error.prob1, max_XO_after_preimp)
    print(nind(cross))

    # re-calculate a map
    cross =quickEstJB.mc(cross, ep=error.prob2)

    #rerun imputation
    cross =  argmax.geno(cross, step=0, map.function='kosambi', error.prob=error.prob2)
    #cross =  calc.genoprob(cross, step=0, map.function='kosambi', error.prob=.0075)
    cross$geno = lapply(cross$geno, function(x) { rownames(x$argmax)=rownames(x$data); return(x) } )

    cross.list[[cross.name]]=cross
    # core analysis is done a this point


    # diagnostic visualizations 
    genotable = hackScanOneObject(cross)
    genotable$gpos=gcoord.key[parents$chr]+parents$pos
    
    pnames=names(parents)[grep('GT', names(parents))]
    pnames=gsub('GT.', '', pnames)

    # plot allele frequencies
    pdf(file=paste0('/data/rrv2/genotyping/plots/af_', cross.name, '.pdf'), width=11,height=7)
        plot(genotable$gpos, m.af, ylim=c(0,1), main=cross.name, xlab='genome pos', ylab=paste(pnames[1], '/', pnames[1], '+', pnames[2]) )
        points(genotable$gpos, genotable$AF, col='blue', ylim=c(0,1))
        abline(v=gcoord.key, lty=2)
    dev.off()
   
    # plot estimated genetic map
    pdf(file=paste0('/data/rrv2/genotyping/plots/genetic_map_', cross.name, '.pdf'), width=11,height=11)
        plotMap(cross, main=cross.name, sub=paste(pnames[1], 'X', pnames[2]) )
    dev.off()
    
}

setwd('/data/rrv2/genotyping/segregants/vcfs/')
for(cross.name in crosses) {
    print(cross.name)
    system(paste0('vcftools --depth  --gzvcf ', cross.name,'.vcf.recode.vcf.gz --out ', cross.name))
}
depth.list=list()
for(cross.name in crosses) {
    depth.list[[cross.name]]=read.delim(paste0(cross.name, '.idepth'), header=T, stringsAsFactors=F)

}
# tail is due to high coverage for BYxRM and YPSxYJM crosses
median(unlist(sapply(depth.list, function(x) x[,3]))
#2.272

       
  ## plot pairwise marker ld measurements
  #  rawmat = pull.geno(cross)
  #  rawmat=ifelse(rawmat==2, .5, -.5)

  #  segmat = pull.argmaxgeno(cross)
  #  segmat=ifelse(segmat==2, 1, -1)

  #  segcor=crossprod(scale(segmat))/(nrow(segmat)-1)
  #  png(file=paste0('/data/rrv2/genotyping/plots/marker_ld_', cross.name, '.png'), width=1920,height=1920)
  #   image.plot(segcor[seq(1,ncol(segmat),10), seq(1,ncol(segmat),10)]^2, main=paste0('LD_', cross.name), axes=F, sub=expression(r^2))
  #   axis(1, at=gcoord.key/c(max(gcoord.key)+chr.lengths[[16]]), labels =names(chr.lengths) )
  #   axis(2, at=gcoord.key/c(max(gcoord.key)+chr.lengths[[16]]), labels =names(chr.lengths) )
  #   abline(v= gcoord.key/c(max(gcoord.key)+chr.lengths[[16]]), col='white', lwd=1)
  #   abline(h= gcoord.key/c(max(gcoord.key)+chr.lengths[[16]]), col='white', lwd=1)
  #  dev.off()

 ## plot segregant genotype calls 
 #   dir.create(paste0('/data/rrv2/genotyping/plots/', cross.name, '/'))
 #   n=foreach(i=1:nrow(segmat)) %dopar% {
 #       png(file=paste0('/data/rrv2/genotyping/plots/', cross.name, '/', rownames(segmat)[i], '.png') , width=1920,height=600)
 #       plot(genotable$gpos,segmat[i,], col='#000000ff', type='h', main=rownames(segmat)[i], xlab='genome position', yaxt='n', xaxt='n', ylab='')
 #       axis(2, at=c(-1,1), labels=c(pnames[1], pnames[2]))
 #       axis(1, at=(gcoord.key+c(gcoord.key, max(gcoord.key)+chr.lengths[[16]])[-1])/2 , labels=names(gcoord.key))
 #       points(genotable$gpos,rawmat[i,], col='lightgrey', type='h')
 #       abline(v=gcoord.key, lty=1, col='blue', lwd=2)
 #       dev.off()
 #   }



#}

#save(parents.list, file='/data/rrv2/genotyping/RData/parents.list.RData')
#save(cross.list, file='/data/rrv2/genotyping/RData/cross.list.RData')






























# consistently uncertain markers
#jmap.cg=do.call('cbind' , lapply(jmap$geno, function(x) x$prob[,,1]))
#jmap.cg.b=ifelse(jmap.cg>.5, 1,2)
#jcga=apply(jmap.cg.b, 2, function(x) sum(x==1, na.rm=T)/(sum(!is.na(x))) )
# bad markers
#bms=apply(jmap.cg,2, function(x) sum(x>.01 & x<.99))

jmap.am=jmap
jmap.am$geno=lapply(jmap.am$geno, function(x) return(list(data=x$argmax, map=x$map)))
jmap.am$geno=lapply(jmap.am$geno, function(x) { class(x)='A'; return(x); })

#gtf=geno.table(jmap.am2)



gtf$AF=do.call('c', g.xs.af)
gtf$AF.smooth=do.call('c', g.xs.af.smooth)
gtf$AF.resid=abs(gtf$AF-gtf$AF.smooth)
gtf$rawAF=m.af
#gtf$cgAF=jcga
plot(m.af, ylim=c(0,1))
points(gtf$AF, ylim=c(0,1), col='blue')
#plot(jcga)






jraw=do.call('cbind', lapply(jmap$geno, function(x) x$data))
# extract imputed genotype calls 
jimp=do.call('cbind', lapply(jmap$geno, function(x) x$argmax))

gtm=data.frame(chr=rep(names(jmap.am$geno), sapply(jmap.am$geno, function(x) ncol(x$data))), pos=unlist((pull.map(jmap.am))))
class(gtm)=c('scanone', 'data.frame')

raw.geno=t(jraw)
raw.geno=raw.geno-1
imp.geno=t(jimp)-1
gtm.geno = data.frame(gtm, raw=raw.geno, imp=imp.geno)
class(gtm.geno)=c('scanone', 'data.frame')

i=229
par(mfrow=c(2,1))
plot(gtm.geno,lodcolumn=c(i), type='p',  
     bandcol='grey90', alternate.chrid=T, col=c('black'),
      main=colnames(raw.geno)[i]
     ) 
plot(gtm.geno,lodcolumn=c(i+ncol(raw.geno)), type='p', 
     bandcol='grey90', alternate.chrid=T, col=c('black') )













plot(gtf, lodcolumn=c(4,1))
#plot(gtf, lodcolumn=c(1))


bad.markers=rownames(gtf)[gtf$AF.resid>.375]
jmap3=jmap
jmap3$geno=lapply(jmap3$geno, function(x) {
        y=x$data
        bic=bad.markers %in% colnames(x$data)
        if(sum(bic)>0){
            bic.hits=bad.markers[bic]
            y[,bic.hits]=NA
        }
        x$data =y 
        return(x)
        }
)
jmap3=quickEstJB(jmap3, error.prob=.005)
jmap3=argmax.geno(jmap3, step=0, map.function='kosambi', error.prob=.005)
jmap3$geno = lapply(jmap3$geno, function(x) { rownames(x$argmax)=rownames(x$data); return(x) } )

g.xs=lapply(jmap3$geno, function(x) x$argmax)
g.xs.af=lapply(g.xs, function(g.x) apply(g.x, 2, function(x) sum(x==1, na.rm=T)/(sum(!is.na(x))) ))
g.xs.af.smooth=lapply(g.xs.af, function(g.x) predict(smooth.spline(g.x))$y)

gtf=data.frame(chr=rep(names(jmap3$geno), sapply(jmap3$geno, function(x) ncol(x$data))), pos=unlist((pull.map(jmap3))))
class(gtf)=c('scanone', 'data.frame')
rownames(gtf)=gsub('^.*\\.', '', rownames(gtf))
gtf$AF=do.call('c', g.xs.af)
jmap4=fake_gmap(jmap3,p.to.d=1000 )
gtf$pos=do.call('c', pull.map(jmap4))
plot(gtf, lodcolumn=1)

jmap.am3=jmap3
jmap.am3$geno=lapply(jmap.am3$geno, function(x) return(list(data=x$argmax, map=x$map)))
jmap.am3$geno=lapply(jmap.am3$geno, function(x) { class(x)='A'; return(x); })

system.time({
  jmap.final.map=est.map(jmap.am3, error.prob=.005, map.function='kosambi', n.cluster=16, tol=1e-4)
})





g.x=do.call('cbind',


            )
gtf$AF=apply(g.x,2, function(x) sum(x==1, na.rm=T)/(sum(!is.na(x))) )
gtfs=split(gtf, gtf$chr)
#



mcor=lapply(jmap.am2$geno, function(x) cor(x$data))




# remove segregants with too many apparent crossovers 
#jmap.lcXO=locateXO(jmap.am, 'chrII')

# hack a scanone object for plotting 
# extract raw genotype calls 
#, col=c('#00000011','#ff000011') )




jsub=subset(jmap.am, c('chrI', 'chrIV'))
jsub$geno=lapply(jsub$geno, function(x) { class(x)='A'; return(x); })

system.time({
ej=est.map(jsub, 'chrIV', error.prob=.005, map.function='kosambi', verbose=T) #,error.prob=0)
})

#jmap.probs=calc.genoprob(jmap, step=0, error.prob=.005, map.function='kosambi')

# relevant data strucu






















#jcross.s=subset(jcross, 'chrI')
#m1=est.map(jcross.s, tol=1e-4, error.prob=.005,verbose=T)
#jcross2=convert2cross2(jcross)

#--------------------------------------------------------

# downsample markers, estimate error rate
#   c6=subset(jcross, chr='chrVI')
#    mnames6=colnames(c6$geno[[1]]$data)
#   mnames.sub=mnames6[seq(1, length(mnames6),15)]
#   c6=pull.markers(c6, mnames.sub)
# for each map 
#   e6=est.map(c6, error.prob=.005, map.function='kosambi', verbose=T)


# use ASMap----------------------------------------------
# preimpute genotypes 

#jcross_preimpute=argmax.geno(jcross, step=0, map.function='kosambi', error.prob=.005)
#jcross_preimpute$geno=lapply(jcross_preimpute$geno, function(x) x$geno=x$argmax)

# timings with ASMap
#mstcross=mstmap.cross(jcross, anchor=T, bychr=T, miss.thresh=1, p.value=2, id='id', detectBadData=TRUE, trace=FALSE)
# mstmap changes chromosome order ... fix this 
#mstcross$geno=mstcross$geno[names(jcross$geno)]
#=========================================================
#g.x=do.call('cbind', lapply(jcross$geno, function(x) x$data))
#g.xaf=apply(g.x,2, function(x) sum(x==1, na.rm=T)/(sum(!is.na(x))) )




gcsvsr=data.frame(colnames(gdata), nchrom, Bc2)
names(gcsvsr)[1:2]=c('id', 'chrom')
genoCallsRQTL = paste('/data/rr/Segregants/' ,crosses[cc], '/' ,crosses[cc], '.rqtl_geno.csvsr', sep='')
write.table(gcsvsr, file=genoCallsRQTL, sep='\t' , row.names=FALSE, quote=FALSE)
readCross = function(g,p) {   return( 
       read.cross(format='csvsr', 
          genfile = g,
          phefile = p,
          na.strings='NA', genotypes=c('B','R'), alleles=c('B','R'), 
          sep='\t', estimate.map=FALSE, comment.char='')   ) 
        }
genoCallsRQTL = paste('/data/rr/Segregants/' ,crosses[cc], '/' ,crosses[cc], '.rqtl_geno.csvsr', sep='')
phenoRQTL = paste('/data/rr/Segregants/' ,crosses[cc], '/' ,crosses[cc], '.rqtl_pheno_norm.csvsr', sep='')
class(cross)[1]='riself'
if(crosses[[cc]]=='A'){
         filt.ind=c(screen_segs[[crosses[cc]]], rep(FALSE,96))
} else { filt.ind=screen_segs[[crosses[cc]]] }
cross=subset(cross,ind=filt.ind) 



# updated code to quickly reestimate genetic map
# needs some map as input





library(qtl)
library(qtl2)
library(ASMap)
load('/data/rr/Segregants/A/A.cross.RData')

# 1) create cross object
# 2) create a proxy for the genetic map using the physical map
# 3) impute missing genotypes given this map
# 4) create new genetic map
# 5) re-impute missing genotypes given this new genetic  map
#loads cross
#








c3=est.map(cross, error.prob=.005, verbose=T)
#c3=quickEstJB(cross, error.prob=.005, map.function='haldane')
c3=calc.genoprob(c2, step=0, error.prob=.005, map.function='kosambi')
c4=argmax.geno(c3, error.prob=.005, map.function='kosambi')



system.time({
c1=mstmap.cross(cross, anchor=T, bychr=T, miss.thresh=1, p.value=2, id='id', detectBadData=TRUE, trace=T, mvest.bc=T)
})
str(cross$geno[[3]])
str(c1$geno[[3]])

plot(c1$geno[[1]]$map, c2$geno[[1]]$map)


test2=quickEst(cc,chr=3, error.prob=.01)








test$geno[[1]]$map=seq(1,length(test$geno[[1]]$map))
test$geno[[3]]$map=seq(1,length(test$geno[[3]]$map))

test2=quickEst(test,chr=3, error.prob=.01)

test$geno[[1]]$map=NULL
test=quickEst(cross,chr=1)

test=pushCross(c1)

c3=quickEstJB(cross)
(c3$geno[[3]]$map)

temp=subset(cross, chr=3)
est=est.rf(temp)


c3=convert2cross2(cross)
str(c3)

pmap=lapply(c3$gmap, function(y) {
    z=as.numeric(sapply(strsplit(names(y), ':|_'), function(x) x[2]))/1e6
    names(z)=names(y)
    return(z)  
       })
c4=c3
c4$pmap=pmap
c5=c4
c5$gmap=c5$pmap


c6=est_map(c5,quiet=F, cores=6)      

c4$gmap=NULL
c5=calc_genoprob(c4)






























pall=geno(pvcf)
gall=geno(svcf)   #re


gc()

        g=rowRanges(svcf)
        ss.vcf=svcf[names(g) %in% names(rowRanges(cross.vcf.segregating))]

        gall=geno(ss.vcf)   #re
        g=rowRanges(ss.vcf) #re














library(ASMap)
library(vcfR)
rv=read.vcfR('/data/rrv2/genotyping/segregants/vcfs/375_chrI.vcf')
grv=extract.gt(rv)
grv=grv[,grep('375', colnames(grv))]

grv[grv=='2']=NA
grv[grv=='3']=NA
grv[grv=='4']=NA
grv[grv=='5']=NA
grv[grv=='6']=NA
grv[grv==0]='1'
grv[grv==1]='2'
grv[is.na(grv)]='U'
af=apply(grv, 1, function(x) sum(x=='A')/sum(!is.na(x)))
grv[af==0,]='U'
grvd=data.frame(grv, stringsAsFactors=F)
omap=mstmap.data.frame(grvd, p.value=2, detectBadData=T, return.imputed=T, miss.thresh=1.5, anchor=T)



grv=extract.gt(rv, as.numeric=T)
grv=grv[,grep('375', colnames(grv))]
grv[grv>1]=NA
grv=grv+1


om=list(geno=list(chrI=list(data=t(grv))), pheno=data.frame(Genotype=colnames(grv)))
class(om)=c('riself', 'cross')
omap=mstmap(om, chr='chrI', anchor=TRUE, detectBadData=T, return.imputed=T, miss.thresh=.99 ,p.value=2 )






, miss.thresh=100) # , chr=rep('1', nrow(grv)),
                       detectBadData =TRUE, return.imputed=TRUE, trace=TRUE, anchor=T)

                       ) # , miss.thresh=950)

