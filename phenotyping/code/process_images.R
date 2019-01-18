######## <PART 1>  ###################################################################################################################
# Image processing code V3 ... greatly simplified from V2 using new EBImage functions and kmeans for image segmentation
require(EBImage)
require(gdata)
library(locfit)
options(browser='chromium-browser')
library(foreach)
library(doMC)
library(Matrix)
library(plyr)
registerDoMC(cores=6)

# a couple of global variables -------------------------------------------
plate.types  <- list('96'=c(8,12), '384'=c(16,24), '1536'=c(32,48) )
dimensions = c(5184,3456)
corners = data.frame(X=c(525,4618,525,4624),
                     Y=c(388,400,3060,3064))
#-----------------------------------------------------------------------
fx.file1='/data/rrv2/phenotyping/code/rr_fx.R'
source(fx.file1)

# base paths should be aware of cross status 
crosses=c('375', 'A', '376', 'B', '377', '393', '381', '3008', '2999', '3000' , '3001', '3049', '3003', '3004', '3043', '3028')

base.paths=c('/data/rr/Phenotyping/RR_Round1/')
base.paths=c('/data/rr/Phenotyping/RR_Round2/')
base.paths=c('/data/rr/Phenotyping/RR_Round3/')

image.paths=paste(base.paths, 'Images/', sep='')
sapply(image.paths, dir.create)
layout.paths=paste(base.paths, 'keys/', sep='')
sapply(layout.paths, dir.create)
key.files = paste(layout.paths, 'Key.csv', sep='')
Results.files = paste(base.paths, 'Results.RData', sep='')
Results.processed.files=paste(base.paths, 'Results.processed.RData', sep='')

out.dirs=paste(base.paths, 'out/', sep='')
sapply(out.dirs, dir.create)


for(cc in 1:length(base.paths)){
    image.path = image.paths[cc]
    #/data/CRISPR_LOH/phenotyping/Images/'
    layout.path = layout.paths[cc]
    #/data/CRISPR_LOH/phenotyping/keys/'
    key.file =  key.files[cc]
    #/data/CRISPR_LOH/phenotyping/Key.xls'
    Results.file =Results.files[cc]
    #/data/CRISPR_LOH/phenotyping/Results.RData'
    outDir=out.dirs[cc]
    #/data/CRISPR_LOH/phenotyping/out/'
    #dir.create(outDir)

    key=read.delim(key.file, stringsAsFactors=F, header=T, sep='\t')
    # Run image processing
    Results=foreach(i=1:nrow(key)) %dopar% {  processImages(i, key, image.path, layout.path, outDir, plate.types[['384']], corners) }
    plate.names=apply(key, 1, paste, collapse='::' )
    names(Results)=plate.names
    #save(Results, file=Results.file)
 }

#writeTXTfiles(outDir, Results)
####### STOPPING POINT <PART 1>  ###########################################################################################################

#load('/data/rr/Phenotyping/RR_Round1/Results.RData')
#load('/data/rr/Phenotyping/RR_Round2/Results.RData')
#load('/data/rr/Phenotyping/RR_Round3/Results.RData')

library(jpeg)
library(raster)

bad.colonies=list()
for(i in 1:length(Results)) {
    #rj=readJPEG(paste0('/data/rr/Phenotyping/RR_Round1/out/', i, '.jpeg'), native=T)
    #rj=readJPEG(paste0('/data/rr/Phenotyping/RR_Round3/out/', i, '.jpeg'), native=T)
    #rj=readJPEG(paste0('/data/rr/Phenotyping/RR_Round1/out/', i, '.jpeg'), native=T)
    xy=Results[[8]][,c('m.cx', 'm.cy')]
    plot(1, type='n', xlim=c(1,5184/4), ylim=c(1,3456/4), main=paste(i,names(Results)[i] ))
    rasterImage(rj, 1,1,5184/4,3456/4, interpolate=FALSE)
    l=locator()
    if(!is.null(l)) {
    l.x=l$x*4
    l.y=(3456/4-l$y)*4

    bad.spots=sapply(1:length(l.x), function(ll) {
        xd=abs(xy[1]-l.x[ll])
        yd=abs(xy[2]-l.y[ll])
        which.min(sqrt(xd^2+yd^2))  })
    bad.spots=unique(bad.spots)
    bad.colonies[[ names(Results)[i]]]=bad.spots
    } else {
    bad.colonies[[ names(Results)[i]]]=NULL
    }
}

#save(bad.colonies, file='/data/rr/Phenotyping/RR_Round1/bad_colonies.RData')
#save(bad.colonies, file='/data/rr/Phenotyping/RR_Round2/bad_colonies.RData')
#save(bad.colonies, file='/data/rr/Phenotyping/RR_Round3/bad_colonies.RData')

# now read in and filter results based on manually removed colonies
# also concatenate Results lists
load('/data/rr/Phenotyping/RR_Round1/Results.RData')
# #Results
load('/data/rr/Phenotyping/RR_Round1/bad_colonies.RData')
# #bad.colonies
for( n in names(bad.colonies)) {
    bc=bad.colonies[[n]]
    Results[[n]][bc,3:ncol(Results[[n]])]=NA
}
#Bad plates
Results[[862]][,3:ncol(Results[[n]])]=NA
Results[[378]][,3:ncol(Results[[n]])]=NA

#These guys are obviously flipped 
#tmp1=Results[[325]]; ntmp1=names(Results)[325]
#tmp2=Results[[326]]; ntmp2=names(Results)[326]
#Results[[325]]=tmp2; names(Results)[325]=ntmp2
#Results[[326]]=tmp1; names(Results)[326]=ntmp1
#
#tmp1=Results[[633]]; ntmp1=names(Results)[633]
#tmp2=Results[[634]]; ntmp2=names(Results)[634]
#Results[[633]]=tmp2; names(Results)[633]=ntmp2
#Results[[634]]=tmp1; names(Results)[634]=ntmp1

Results[[325]][,3:ncol(Results[[n]])]=NA
Results[[326]][,3:ncol(Results[[n]])]=NA

Results[[633]][,3:ncol(Results[[n]])]=NA
Results[[634]][,3:ncol(Results[[n]])]=NA

Results1=Results
Key1=reconstructKey(Results1)
Key1$PermutationGroup=as.factor(1)
names(Results1)=apply(Key1, 1, paste, collapse='::')

load('/data/rr/Phenotyping/RR_Round2/Results.RData')
# #Results
load('/data/rr/Phenotyping/RR_Round2/bad_colonies.RData')
# #bad.colonies
for( n in names(bad.colonies)) {
    bc=bad.colonies[[n]]
    Results[[n]][bc,3:ncol(Results[[n]])]=NA
}
Results[[668]][,3:ncol(Results[[n]])]=NA

Results2=Results
Key2=reconstructKey(Results2)
Key2$PermutationGroup=as.factor(2)
names(Results2)=apply(Key2, 1, paste, collapse='::')

load('/data/rr/Phenotyping/RR_Round3/Results.RData')
# #Results
load('/data/rr/Phenotyping/RR_Round3/bad_colonies.RData')
# #bad.colonies
for( n in names(bad.colonies)) {
    bc=bad.colonies[[n]]
    Results[[n]][bc,3:ncol(Results[[n]])]=NA
}
Results3=Results
Key3=reconstructKey(Results3)
Key3$PermutationGroup=as.factor(3)
names(Results3)=apply(Key3, 1, paste, collapse='::')

#Add back batch information before concatenating 
Results=(c(Results1, Results2, Results3))
key=rbind(Key1[,1:6], Key2[,1:6])
key=rbind(key, Key3[,1:6])

SEG_Results.df= getResults.df(Results, key, c(1:length(Results)))
SEG_Results.df$index=1:nrow(SEG_Results.df)
more.bad.colonies=c(1187334,915925,800697,776356,1202684,1202320,370182,836210,943691,331553,458420,534324,2281,16825,846409,968997,1014931,206569,392425,37609,3781,3817,46057,460340,536244,840599,1015696,1180048,1201160,1024885, 1020661,1185020,1058296)
SEG_Results.df[more.bad.colonies,2:12]=NA
#save(SEG_Results.df, file='/data/rr/Phenotyping/SEG_Results_df_manual_filtered.RData')

# main 
load('/data/rr/Phenotyping/SEG_Results_df_manual_filtered.RData')

foi='s.radius.mean'
SEGpheno=split(SEG_Results.df[,foi], SEG_Results.df$condition)
SEGpheno.names=split(SEG_Results.df$strain, SEG_Results.df$condition)
SEGpheno.index=split(SEG_Results.df$index, SEG_Results.df$condition)
SEGpheno.batch=split(SEG_Results.df$batch, SEG_Results.df$condition)
SEGpheno.layout=split(SEG_Results.df$layout, as.character(SEG_Results.df$condition))

Spheno=mapply(function(x,n){split(x,n)}, SEGpheno, SEGpheno.names, SIMPLIFY=F)
Spheno2=mapply(function(x,n){split(x,n)}, SEGpheno.index, SEGpheno.names, SIMPLIFY=F)
Spheno3=mapply(function(x,n){split(x,n)}, SEGpheno.layout, SEGpheno.names, SIMPLIFY=F)
#Spheno2=mapply(function(x,n){split(x,n)}, SEGpheno.index, SEGpheno.names, SIMPLIFY=F)

xc=lapply(Spheno, function(x) do.call('rbind', x) )
#36 38 39
for(i in 1:41) {
   # png(file=paste0('/data/rr/Phenotyping/plots/filtered_radius/', i, '.png'), width=1080, height=1080)
    plot(xc[[i]][,1], xc[[i]][,2], xlab='radius rep1', ylab='radius rep2', main=names(Spheno)[i], sub=paste0('R=', round(cor(xc[[i]][,1], xc[[i]][,2], use='pairwise.complete.obs'),3)))
  #  xy=xy.coords(xc[[i]][,1], xc[[i]][,2])
  #  mid=identify(xy)
  #  print(Spheno2[[i]][mid])
  #  print(Spheno[[i]][mid])
    readline()
   # dev.off()
    #readline()
}

# for each batch ----------------------------------------------------------------------------------------------------
b=1
# for each layout
lay="2999_B.txt"
condition = "Galactose;;1"
control  = "YPD;;1"

sind.p= which( SEG_Results.df$condition==condition &
SEG_Results.df$batch==b &
SEG_Results.df$layout==lay)

sind.c=which( SEG_Results.df$condition==control&
SEG_Results.df$batch==b &
SEG_Results.df$layout==lay)
#average control plate

sind.c=which( SEG_Results.df$condition==control&
SEG_Results.df$batch==b &
SEG_Results.df$layout==lay)

cmean=sapply(split(SEG_Results.df[sind.c,]$s.radius.mean, SEG_Results.df[sind.c,]$strain), mean, na.rm=T)
pmean=sapply(split(SEG_Results.df[sind.p,]$s.radius.mean, SEG_Results.df[sind.p,]$strain), mean, na.rm=T)
plot(pmean,cmean)
#--------------------------------------------------------------------------------------------------------------------


#fix 3008

# note, with column order 180 rotation is reversing of index order
#plate.in='3008_G1'
#template96=matrix(1:96, 8,12)
#template96=matrix(paste(plate.in, sprintf("%02d", template96), sep='_'), 8, 12)
#rotate <- function(x) t(apply(x, 2, rev))
#m96.flip = rotate(rotate(template96))
#inorder=as.vector(template96)
#oflip=as.vector(m96.flip)
#vertical.flip.ind=match(inorder, oflip)


#mapping!!!!
library(Matrix)
library(rrBLUP)
#source('~/Dropbox/Public/CRISPR_LOH/calcMM.R')
#stats.galore=list()
#load(pheno.file)
#for A and B
#for( cc in c(2,4) ) {

#asn=names(Spheno[[1]])
crosses=c('375', 'A', '376', 'B', '377', '393', '381', '3008', '2999', '3000' , '3001', '3049', '3003', '3004', '3043', '3028')


# now append information about control plate growth
load('/data/rr/Phenotyping/working/SEG_Results_df_manual_filtered.RData')
cross.id=sapply(strsplit(as.character(SEG_Results.df$layout), '_'), function(x) x[1])
SEG_Results.df$cross.id=cross.id

scl=dlply(SEG_Results.df, c("cross.id", "condition", "layout"))
n.measures=sapply(scl, function(x) nrow(x))
tech.rep.plates=names(which(n.measures>384))
for(tnp in tech.rep.plates) {
   print(tnp)
   spn=split(scl[[tnp]], scl[[tnp]]$strain)
   one.rep=do.call('rbind', lapply(spn, function(x) x[1,]))
   avg.vals=do.call('rbind', lapply(spn, function(x) apply(x[,2:12], 2, mean, na.rm=T)))
   if(identical(rownames(avg.vals), one.rep$strain) ) {
        one.rep[,2:12]=avg.vals
   } else {print(paste(tnp, 'fatal' ))}
   scl[[tnp]]=one.rep
}
scl.groups=do.call('rbind', lapply(scl, function(x) { y=(c(x[1,c('cross.id', 'condition', 'batch', 'layout')])) 
    y$layout=as.character(y$layout)
    return(as.character(y))
}))
colnames(scl.groups)=c("cross.id", "condition", "batch", "layout")
rownames(scl.groups)=names(scl)

for(n in 1:nrow(scl.groups)) {
    print(n)
    scg=scl.groups[n,]
    if(!grepl('YNB', scg[2]) ) {
        print('YPD')
        control.plate= which(scl.groups[,1]==scg[1] & grepl('YPD;;', scl.groups[,2]) & scl.groups[,3]==scg[3] & scl.groups[,4]==scg[4])
        print(control.plate)
    } else  {
        print('YNB')
        control.plate= which(scl.groups[,1]==scg[1] & grepl('YNB;;', scl.groups[,2]) & scl.groups[,3]==scg[3] & scl.groups[,4]==scg[4])
        print(control.plate)
    }
    control.df=scl[[control.plate]][,2:12]
    colnames(control.df)=paste0('ctrl.', colnames(control.df))
    scl[[n]]=cbind(scl[[n]], control.df)    
    
    norm.pheno=matrix(NA, 384, 3)
    colnames(norm.pheno)=c('norm.s.area', 'norm.s.perimeter', 'norm.s.radius.mean')
    rownames(norm.pheno)=1:384

    # precompute corrections for plate effects 
    if(grepl('YNB;;', scg[2]) |grepl('YPD;;', scg[2]) | sum(is.na(scl[[n]]$s.area))==384 ) {
        norm.pheno[,1]=scl[[n]]$s.area
        norm.pheno[,2]=scl[[n]]$s.perimeter
        norm.pheno[,3]=scl[[n]]$s.radius.mean
    } else {
        r1=(residuals(lm(scl[[n]]$s.area~scl[[n]]$ctrl.s.area)))
        norm.pheno[match(names(r1), rownames(norm.pheno)), 1]=as.vector(r1)
        r2=(residuals(lm(scl[[n]]$s.perimeter~scl[[n]]$ctrl.s.perimeter)))
        norm.pheno[match(names(r2), rownames(norm.pheno)), 2]=as.vector(r2)
        r3=(residuals(lm(scl[[n]]$s.radius.mean~scl[[n]]$ctrl.s.radius.mean)))
        norm.pheno[match(names(r3), rownames(norm.pheno)), 3]=as.vector(r3)
   }
    rownames(norm.pheno)=NULL
    scl[[n]]=cbind(scl[[n]], norm.pheno)    
}
seg_phenos_norm = rbind.fill.matrix(scl)
#seg_phenos_norm=data.frame(seg_phenos_norm)
write.table(seg_phenos_norm, file='/data/rr/Phenotyping/seg_phenos_norm.txt', row.names=F, sep='\t', quote=F)

seg_phenos=read.delim('/data/rr/Phenotyping/seg_phenos_norm.txt', header=T, stringsAsFactors=F, sep='\t')


foi='s.radius.mean'
SEGpheno=split(as.numeric(seg_phenos[,foi]), seg_phenos[,'condition'])
SEGpheno.names=split(seg_phenos[,'strain'], seg_phenos[,'condition'])
RAWpheno=mapply(function(x,n){split(x,n)}, SEGpheno, SEGpheno.names, SIMPLIFY=F)

foi='norm.s.radius.mean'
SEGpheno=split(as.numeric(seg_phenos[,foi]), seg_phenos[,'condition'])
SEGpheno.names=split(seg_phenos[,'strain'], seg_phenos[,'condition'])
NORMpheno=mapply(function(x,n){split(x,n)}, SEGpheno, SEGpheno.names, SIMPLIFY=F)
save(NORMpheno, file='/data/rr/Phenotyping/NORMpheno.RData')

