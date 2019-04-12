# Reconstruct neighbor-joining tree for Figure 1, also outputs crossing ring design for figure 1

#This file was produced by vcfisec.
#The command line was:	bcftools isec  -p isec_ouput /data/rrv2/genotyping/parents/parents.vcf.gz /data/rrv2/1002genomes/1011MatrixRename.vcf.gz

#Using the following file names:
#isec_ouput/0000.vcf	for records private to	/data/rrv2/genotyping/parents/parents.vcf.gz
#isec_ouput/0001.vcf	for records private to	/data/rrv2/1002genomes/1011MatrixRename.vcf.gz
#isec_ouput/0002.vcf	for records from /data/rrv2/genotyping/parents/parents.vcf.gz shared by both	/data/rrv2/genotyping/parents/parents.vcf.gz /data/rrv2/1002genomes/1011MatrixRename.vcf.gz
#isec_ouput/0003.vcf	for records from /data/rrv2/1002genomes/1011MatrixRename.vcf.gz shared by both	/data/rrv2/genotyping/parents/parents.vcf.gz /data/rrv2/1002genomes/1011MatrixRename.vcf.gz

##vcftools --gzvcf /data/rrv2/1002genomes/1011MatrixRename.vcf.gz --freq --out all.freq  --min-alleles 2 --max-alleles 2
##bcftools gtcheck -g 0002.vcf.gz 0003.vcf.gz > testQuery

library(vcfR)
parents=read.vcfR('0002.vcf.gz', verbose=T)
panel=read.vcfR('0003.vcf.gz', verbose=T)
par.gt=extract.gt(parents,as.numeric=T)
pan.gt=extract.gt(panel, as.numeric=T)
par.gt.s=par.gt[rownames(pan.gt),]

#rownames(pan.gt) %in% rownames(par.gt)
pgtc=cor(pan.gt, use='pairwise.complete.obs')
pgtcd=1-abs(pgtc^2)
findBest=cor(par.gt.s, pan.gt, use='pairwise.complete.obs')
best.match=colnames(findBest)[apply(findBest,1, which.max)]
colnames(par.gt.s)
#save(par.gt.s,file='/data/rrv2/1002genomes/pargts.RData')
#save(best.match,file='/data/rrv2/1002genomes/bestMatch.RData')

#library(fields)
library(ape)
library(igraph)
library(RColorBrewer)
load('/data/rrv2/1002genomes/pargts.RData')
load('/data/rrv2/1002genomes/bestMatch.RData')

#edge.col=rep('red', 81), tip.col=rep('red', 81))
source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')
# to load crosses.to.parents list
RR=crosses.to.parents
#annotate origin
locations=c('clinical', 'lab','wine', 'oak', 'clinical', 'fermenting rice','oak', 'clinical', 'oak','wine','fruit', 'palm wine', 'clinical', 'clinical', 'soil', 'wine')
strain.loc.color=data.frame(cbind(do.call('rbind', RR)[,1], locations), stringsAsFactors=F)
colorscheme=brewer.pal(length(unique(locations)), 'Set3')#tim.colors(length(unique(locations)))
strain.loc.color$color=colorscheme[match(as.numeric(as.factor(strain.loc.color[,2])), 1:length(colorscheme))]


dissim=read.delim('/data/rrv2/1002genomes/1011DistanceMatrixBasedOnSNPs.tab', header=T, stringsAsFactors=F)
rownames(dissim)=dissim$STD
dissim=dissim[,-1]
bd=bionj(data.matrix(dissim))
bd$edge.length[bd$edge.length>.5]=.5
# truncate the chinese strains for plotting
# relabel strains
matched.strains=match(best.match, bd$tip.label)
bd$tip.label[matched.strains]=paste0(colnames(par.gt.s), ' (', best.match, ')')
bd$tip.label.renamed[matched.strains]=colnames(par.gt.s)

tcol=strain.loc.color[match(bd$tip.label.renamed, strain.loc.color[,1]),3]
tcol[is.na(tcol)]='#00000000'

# Figure 1 --------------------------------------------------------------------------------------------------
# basically this, but slightly tweaked in inkscape
pdf(file='/home/jbloom/Dropbox/RR/Figures and Tables/Figure1.pdf', width=14, height=8)
par(mfrow=c(1,2), oma=c(1,1,1,1))
bdp=plot(bd, type='unrooted', 
     edge.lty=1,
     node.pos=1,
     show.tip.label=T,
     tip.color=tcol,
     no.margin=F,
     font=2,
     label.offset=.15,
     align.tip.label=T,
     lab4ut='axial'
)
hh=locator()
for(i in 1:3){
    segments(hh$x[i], hh$y[i], hh$x[i+1], hh$y[i+1])
}
stc=strain.loc.color[!duplicated(strain.loc.color[,2]),]
legend('topleft', legend=stc[,2], fill=stc$color)

g=graph.ring(16,circular=T)
plot(g, edge.label= names(RR),
      vertex.label=do.call('rbind', RR)[,1],
      vertex.label.dist=0,
      vertex.color=strain.loc.color$color, 
      vertex.label.color='black',
      edge.label.color= 'black') 
dev.off()
#------------------------------------------------------------------------------------------------------------------
 
