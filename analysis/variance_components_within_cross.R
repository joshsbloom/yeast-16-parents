library(qtl)
library(regress)
library(ggplot2)
library(cowplot)
library(reshape)
library(data.table)
library(Matrix)

# within cross  variance components analysis ---------------------------------------------------------------------------


# accessory functions
source('/data/rrv2/analysis/mapping_fx.R')

# output from segregants_hmm.R
load('/data/rrv2/genotyping/RData/cross.list.RData')
# output from segregants_hmm.R
load('/data/rrv2/genotyping/RData/parents.list.RData')
# phenotypes
load('/data/rr/Phenotyping/NORMpheno.RData')
# QTL mapping results
load('/data/rrv2/genotyping/RData/FDR_cross.peaks.RData')

source('/data/rrv2/genotyping/code/segregants_hmm_fx.R')
# list to hold output
cross.VCs=list()
cross.VCs2=list()
for(cross.name in crosses[1:16]) {
    print(cross.name)
    cross=cross.list[[cross.name]]
    if(cross.name=='A') {       cross=subset(cross, ind=!grepl('A11', as.character(cross$pheno$id)))    }
    snames = as.character(cross$pheno$id)
    g=pull.argmaxgeno(cross)
 
    # are there fixed loci ?? (no)-------------------------------
    g.af=apply(g,2,function(x) sum(x==1))
    parents.list[[cross.name]]$fixed=(g.af==0 | g.af==nrow(g))
    fixed.loci=which(parents.list[[cross.name]]$fixed)
    if(length(fixed.loci)>0) {    g=g[,-fixed.loci] }
    #------------------------------------------------------------
    g.r=g[,-which(duplicated(g, MARGIN=2))]
    g.s=scale(g.r)
    
    A=tcrossprod(g.s)/(ncol(g.s))
    cps=cross.peaks[[cross.name]]
    subPheno=lapply(NORMpheno, function(x) x[match(snames, names(x))])
    for(pheno in names(subPheno)[-1]) {
            print(pheno)
            cpQTL=cps[cps$trait==pheno,]
            # only QTL with FDR<5%
            cpQTL=cpQTL[cpQTL$q<.05,]
            pproc=processPhenos_for_MM(NORMpheno[[pheno]], snames)
            unames=unique(names(pproc$y))
            # for maltose cross 3001 ??? what happened here?
            if(length(cpQTL$pmarker)==0) {
                cpQTL=data.frame(pmarker=sample(colnames(g.s),2))
            }
            Q=(tcrossprod(g.s[unames,cpQTL$pmarker])/length(cpQTL$pmarker))
            Z=as.matrix(pproc$Z)

            ZZt=Z%*%t(Z)
            ZQZt=Z%*%Q%*%t(Z)
            ZAZt=Z%*%A[unames,unames]%*%t(Z)
            ZQQZt=Z%*%(Q*Q)%*%t(Z)
            ZQAZt=Z%*%(A[unames,unames]*Q)%*%t(Z)
            ZAAZt=Z%*%(A[unames,unames]*A[unames,unames])%*%t(Z)
            # ZQAAZt=as.matrix(pproc$Z%*%(A*A*Q)%*%t(pproc$Z))
            # ZQQAZt=as.matrix(pproc$Z%*%(Q*Q*A)%*%t(pproc$Z))
            # ZQQQZt=as.matrix(pproc$Z%*%(Q*Q*Q)%*%t(pproc$Z))
            # ZAAAZt=as.matrix(pproc$Z%*%(A*A*A)%*%t(pproc$Z))
            tryCatch({
                cross.VCs[[cross.name]][[pheno]]=extractVarCompResultsR(regress(pproc$y~1,~ZZt+ZQZt+ZAZt+ZQQZt+ZQAZt+ZAAZt, verbose=T,  pos=rep(T, 12)))
                cross.VCs2[[cross.name]][[pheno]]=extractVarCompResultsR(regress(pproc$y~1,~ZZt+ZQZt+ZAZt+ZAAZt, verbose=T,  pos=rep(T, 5)))

            },error=function(e) {
                cross.VCs[[cross.name]][[pheno]]=extractVarCompResultsR(regress(scale(pproc$y)~1,~ZZt+ZQZt+ZAZt+ZQQZt+ZQAZt+ZAAZt, verbose=T))
                cross.VCs2[[cross.name]][[pheno]]=extractVarCompResultsR(regress(scale(pproc$y)~1,~ZZt+ZQZt+ZAZt+ZAAZt, verbose=T))
            })
   }
}

# currently finished B
#save(cross.VCs, file='/data/rrv2/genotyping/RData/withinCrossVCs_FDR.RData')
#save(cross.VCs2, file='/data/rrv2/genotyping/RData/withinCrossVCs2_FDR.RData')
load(  '/data/rrv2/genotyping/RData/withinCrossVCs_FDR.RData')

#save(cross.VCs, file='/data/rrv2/genotyping/RData/withinCrossVCs.RData')
#load(  '/data/rrv2/genotyping/RData/withinCrossVCs.RData')

# extract components and normalize
# remove the repeat data for YPD::2 and YPD;;3 (37, and 38)
WC.sigma=sapply(cross.VCs, function(x) sapply(x[-c(37,38)], function(y) y$sigma[rev(c(2:6,1,7))] ))
WC.sigma.se=sapply(cross.VCs, function(x) sapply(x[-c(37,38)], function(y) sqrt(diag(y$sigma.cov)[rev(c(2:6,1,7))])))
WC.sigma.nf=sapply(WC.sigma, colSums)
WC.sigma.norm=mapply(function(x,y) {
            z=t(t(x)/y)
            #z=z[rev(c(2:6,1,7)),]
            return(z)
            }, x=WC.sigma, y=WC.sigma.nf)


WC.sigma.norm.pos=lapply(WC.sigma.norm, function(x) {
                    apply(x, 2, function(x) rev(cumsum(rev(x))))
            })
WC.sigma.se.norm=mapply(function(x,y) {
            z=t(t(x)/y)
            #z=z[rev(c(2:6,1,7)),]
            return(z)
            }, x=WC.sigma.se, y=WC.sigma.nf)


# histogram of heritability explained by additive QTL
Q=unlist(sapply(WC.sigma.norm, function(x) x[7,]))
Q_A=unlist(sapply(WC.sigma.norm, function(x) x[7,]+x[6,]))
Q_AA=unlist(sapply(WC.sigma.norm, function(x) (x[3,]+x[4,]+x[5,])))

# additive explained
hist(unlist(sapply(WC.sigma.norm, function(x) x[7,]/(x[6,]+x[7,]))), breaks=100)
hist(unlist(sapply(WC.sigma.norm, function(x) x[5,]/(x[3,]+x[4,]+x[5,]))), breaks=100)

#library(gplots)
#plotCI(Q_A,Q, uiw=Q_A.se, err='x', ylim=c(0,1), xlim=c(0,1))
#plotCI(Q_A,Q, Q.se, err='y', add=T)
plot(Q_A,Q)
abline(0,1)
BYxRM.logical=grepl('^A', names(Q_A))
WCggplot=data.frame(trait=as.vector(do.call('c', sapply(WC.sigma.norm, function(x) names(x[1,])))),
               cross=rep(names(WC.sigma.norm), sapply(WC.sigma.norm, ncol)),
               additive=Q_A, non_additive=Q_AA, additive_plus_non_additive=Q_A+Q_AA)


# Figure 2-------------------------------------------------------------------------------------------------------------
h2sm=melt(h2s)
colnames(h2sm)=c('trait', 'cross', 'av')
qvcrm=melt(qcvr)
colnames(qvcrm)=c('trait', 'cross', 'pve')
h2sm=data.frame(h2sm, qvcrm[,3])
colnames(h2sm)[4]='pve'
#Additive variance (whole genome)
#Phenotypic variance explained by QTL
h2sm$pvej=melt(aam)[,3]
main_plot=ggplot(data=h2sm, aes(x=av, y=pvej))+geom_point(alpha=.5)+
    geom_point(data=h2sm[h2sm$cross=='A',], aes(x=av, y=pvej), size=1.5, colour='red')+
    scale_y_continuous(limits=c(0,1), expand=c(0.01,0), position='right')+
    scale_x_continuous(limits=c(0,1), expand=c(0.01,0))+geom_abline(intercept=0, slope=1)+
    xlab('Additive variance (whole genome)')+
    ylab('Phenotypic variance explained by QTL')+theme_classic()+
    theme(axis.text.x=element_text(size=rel(1.2),color='black'),
          axis.text.y=element_text(size=rel(1.2),color='black')
          )
inset_plot=ggplot(WCggplot, aes(x=non_additive/additive))+geom_histogram(binwidth=.01)+
    scale_x_continuous(limits=c(0,1), expand=c(0,0))+
    xlab(expression(frac('Pairwise interaction variance','Additive variance')))+
    ylab('Count')+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.text.x=element_text(size=rel(.8),color='black'),
          axis.text.y=element_text(size=rel(.8),color='black'),
          axis.title.x=element_text(size=rel(.8), color='black'),
          axis.title.y=element_text(size=rel(.8), color='black'),
          )+
    theme(panel.background = element_rect(fill = "transparent"))
 ggdraw() +
    draw_plot(main_plot)+
    draw_plot(inset_plot, x=.01, y=.52, width=.45, height=.45)
ggsave('/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/Figure2.pdf', width=8.5,height=8.5)
#---------------------------------------------------------------------------------------------------------------------------






############################################################################################################################3
# reviewer comment to fit single component epistasis model 
WC.sigma2=sapply(cross.VCs2, function(x) sapply(x[-c(37,38)], function(y) y$sigma[rev(c(2:4,1,5))] ))
WC.sigma2.se=sapply(cross.VCs2, function(x) sapply(x[-c(37,38)], function(y) sqrt(diag(y$sigma.cov)[rev(c(2:4,1,5))])))
WC.sigma2.nf=sapply(WC.sigma2, colSums)
WC.sigma2.norm=mapply(function(x,y) {
            z=t(t(x)/y)
            #z=z[rev(c(2:6,1,7)),]
            return(z)
            }, x=WC.sigma2, y=WC.sigma2.nf)


WC.sigma2.norm.pos=lapply(WC.sigma2.norm, function(x) {
                    apply(x, 2, function(x) rev(cumsum(rev(x))))
            })
WC.sigma2.se.norm=mapply(function(x,y) {
            z=t(t(x)/y)
            #z=z[rev(c(2:6,1,7)),]
            return(z)
            }, x=WC.sigma2.se, y=WC.sigma2.nf)


# histogram of heritability explained by additive QTL
Q2=unlist(sapply(WC.sigma2.norm, function(x) x[5,]))
Q_A2=unlist(sapply(WC.sigma2.norm, function(x) x[5,]+x[4,]))
Q_AA2=unlist(sapply(WC.sigma2.norm, function(x) (x[3,])))

WCggplot$tc=paste0(WCggplot$trait, WCggplot$cross, sep="_")
WCggplot2$tc=paste0(WCggplot2$trait, WCggplot2$cross, sep="_")

twomodels.comp=merge(WCggplot, WCggplot2, by='tc')
png(file='/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/one_vs_three_QQ.png', width=512, height=512)
plot(twomodels.comp$non_additive.x, twomodels.comp$non_additive.y, 
     xlim=c(0,.5), ylim=c(0,.5),
     xlab=" Q ○ Q  + Q ○ A  +  A ○ A", ylab="A ○ A"
)
abline(0,1)
dev.off()

WCggplot2=data.frame(trait=as.vector(do.call('c', sapply(WC.sigma2.norm, function(x) names(x[1,])))),
               cross=rep(names(WC.sigma2.norm), sapply(WC.sigma2.norm, ncol)),
               additive=Q_A2, non_additive=Q_AA2, additive_plus_non_additive=Q_A2+Q_AA2)


WCV2=melt(WC.sigma2.norm)
WCV2.se=melt(WC.sigma2.se.norm)
ypos2=melt(WC.sigma2.norm.pos)
names(WCV2)=c('Component', 'Trait', 'fraction_of_variance', 'Cross')
WCV2$se=WCV2.se[,3]
WCV2$ypos=ypos2[,3]
#WCV=WCV[WCV$Trait!='4NQO;0.15ug/mL;3',]
WCV2$Cross=factor(WCV2$Cross, levels=crosses)
#WCV$se[as.character(WCV$Component)=='ZQZt']=NA
WCV2$Component=plyr::revalue(WCV2$Component, c("In"="Residual", "ZZt"="Repeatibility", "ZAAZt"="A ○ A",  "ZAZt"="A", "ZQZt"="Q"))
WCV2=WCV2[WCV2$Trait!='YPD;;2',]
WCV2=WCV2[WCV2$Trait!='YPD;;3',]
WCV2$Trait=droplevels(WCV2$Trait)

levels(WCV2$Trait)[34]='YNB_ph8'
levels(WCV2$Trait)[36]='YPD_15C'
levels(WCV2$Trait)[33]='YNB_ph3'
levels(WCV2$Trait)[10]='EtOH_Glu'
levels(WCV2$Trait)[37]='YPD_37C'
levels(WCV2$Trait)=gsub(';.*','', levels(WCV2$Trait))
levels(WCV2$Trait)=gsub('_', ' ', levels(WCV2$Trait))

sup.table.2c=WCV2
sup.table.2c=sup.table.2c[,c(2,4,1,3,5)]
WriteXLS(c('sup.table.2c'), 
         "/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/SupplementaryTable2_2.xls", SheetNames=c('multiple_components, 2'))
######################################################################################################################################################





















#png(file='/data/rrv2/Figures_and_Tables/additive_vs_interaction_scatter.png', width=768, height=768)
#plot(Q_A+Q_AA, Q_A,
#     pch=20,
#     xlim=c(0,1),
#     ylim=c(0,1), 
#     xlab='Additive + Interaction variance (whole genome)',
#     ylab='Additive variance (whole genome)' )
#points(c(Q_A+Q_AA)[BYxRM.logical], Q_A[BYxRM.logical], col='red', pch=20)
#abline(0,1)
#par(fig=c(.075,.5,.5,1),new=T)
#hist(Q_AA/Q_A, breaks=100, xlim=c(0,1), main='', xlab='(Interaction variance)/(Additive variance)')
#dev.off()
##plot(Q_A, QAAA, pch=20, xlim=c(0,1), ylim=c(0,1))
##points(c(Q_A)[BYxRM.logical], QAAA[BYxRM.logical], col='red', pch=20)
##abline(0,1)
##par(fig=c(.075,.5,.5,1),new=T)
##hist(QAAA/(QAAA+Q_A), breaks=100)

# Supplementary Figure 1 ---------------------------------------------------------------------------------------------------------------------------------------
WCV=melt(WC.sigma.norm)
WCV.se=melt(WC.sigma.se.norm)
ypos=melt(WC.sigma.norm.pos)
names(WCV)=c('Component', 'Trait', 'fraction_of_variance', 'Cross')
WCV$se=WCV.se[,3]
WCV$ypos=ypos[,3]
#WCV=WCV[WCV$Trait!='4NQO;0.15ug/mL;3',]
WCV$Cross=factor(WCV$Cross, levels=crosses)
#WCV$se[as.character(WCV$Component)=='ZQZt']=NA
WCV$Component=plyr::revalue(WCV$Component, c("In"="Residual", "ZZt"="Repeatibility", "ZAAZt"="A ○ A", "ZQAZt"="Q ○ A", "ZQQZt"="Q ○ Q", "ZAZt"="A", "ZQZt"="Q"))
WCV=WCV[WCV$Trait!='YPD;;2',]
WCV=WCV[WCV$Trait!='YPD;;3',]
WCV$Trait=droplevels(WCV$Trait)

levels(WCV$Trait)[34]='YNB_ph8'
levels(WCV$Trait)[36]='YPD_15C'
levels(WCV$Trait)[33]='YNB_ph3'
levels(WCV$Trait)[10]='EtOH_Glu'
levels(WCV$Trait)[37]='YPD_37C'
levels(WCV$Trait)=gsub(';.*','', levels(WCV$Trait))
levels(WCV$Trait)=gsub('_', ' ', levels(WCV$Trait))

sup.table.2b=WCV

colorsR=rev(colorRampPalette(brewer.pal(7, "Paired"))(7))
colorsR[1]='white'
colorsR[2]='mistyrose2' #colorsR[3]
colorsR[3:5]=c('darkseagreen', 'lawngreen' ,'forestgreen')

label_cross=unlist(lapply(crosses.to.parents, function(x) paste(x, collapse= '   x   ')))
png(file='/data/rrv2/Figures_and_Tables/Supp_withinCrossVC.png', width=1920, height=1080)
#pdf(file='/home/jbloom/Desktop/111417/withinCrossVC.pdf', height=8.5*2, width=11*2)
ggplot(WCV, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    facet_wrap(~Cross, ncol=4, labeller=as_labeller(label_cross))+coord_cartesian(ylim=c(0,1))+ylab('fraction of phenotypic variance')+
    scale_fill_manual(values =colorsR) +theme_bw()+theme(axis.text.x=element_text(angle=70,hjust=1))+
    theme(legend.key = element_rect(colour = "black"))+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))
ggsave('/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure1.png', width=18, height=10)
#ggsave('/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure1.png', width=18, height=10)

dev.off()                   

WCV1=WCV[WCV$Cross %in% c("375",  "A",    "376",  "B"),]
WCV2=WCV[WCV$Cross %in% c( "377",  "393",  "381",  "3008"),]
WCV3=WCV[WCV$Cross %in% c("2999", "3000","3001", "3049"),]
WCV4=WCV[WCV$Cross %in% c("3003", "3004", "3043", "3028"),]
pdf(file='/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure1.pdf', width=8.5, height=5.5)
ggplot(WCV1, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    facet_wrap(~Cross, ncol=2, labeller=as_labeller(label_cross))+coord_cartesian(ylim=c(0,1))+ylab('Fraction of phenotypic variance')+
    scale_fill_manual(values =colorsR) +theme_bw()+theme(axis.text.x=element_text(angle=70,hjust=1))+
    theme(legend.key = element_rect(colour = "black"))+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))
ggsave(file='/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure1A.png', width=8.5, height=5.5)

ggplot(WCV2, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    facet_wrap(~Cross, ncol=2, labeller=as_labeller(label_cross))+coord_cartesian(ylim=c(0,1))+ylab('Fraction of phenotypic variance')+
    scale_fill_manual(values =colorsR) +theme_bw()+theme(axis.text.x=element_text(angle=70,hjust=1))+
    theme(legend.key = element_rect(colour = "black"))+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))
ggsave(file='/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure1B.png', width=8.5, height=5.5)

ggplot(WCV3, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    facet_wrap(~Cross, ncol=2, labeller=as_labeller(label_cross))+coord_cartesian(ylim=c(0,1))+ylab('Fraction of phenotypic variance')+
    scale_fill_manual(values =colorsR) +theme_bw()+theme(axis.text.x=element_text(angle=70,hjust=1))+
    theme(legend.key = element_rect(colour = "black"))+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))
ggsave(file='/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure1C.png', width=8.5, height=5.5)

ggplot(WCV4, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    facet_wrap(~Cross, ncol=2, labeller=as_labeller(label_cross))+coord_cartesian(ylim=c(0,1))+ylab('Fraction of phenotypic variance')+
    scale_fill_manual(values =colorsR) +theme_bw()+theme(axis.text.x=element_text(angle=70,hjust=1))+
    theme(legend.key = element_rect(colour = "black"))+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))
ggsave(file='/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure1D.png', width=8.5, height=5.5)

dev.off()
#levels(WCV$Cross)
# [1] "375"  "A"    "376"  "B"    "377"  "393"  "381"  "3008" "2999" "3000"
#[11] "3001" "3049" "3003" "3004" "3043" "3028"

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------



#cF=as.factor(rep(names(pheno_extracted), sapply(pheno_extracted, nrow)))
#pe=do.call('rbind', (pheno_extracted))
#cross.vars=sapply(split(data.frame(pe[,-c(37,38)]), cF), function(x) apply(x,2,var))
#cross.cnts.table=table(cF)
#aam=aam[,colnames(cross.vars)]
#median(aam/h2s[,colnames(aam)])

#plot(as.vector(h2s[,colnames(aam)]), as.vector(aam), ylim=c(0,1), xlim=c(0,1),xlab='additive heritability (whole genome)', ylab='CV QTL model variance explained')
#abline(0,1)

#wsum=rowSums(cross.vars*(as.vector(cross.cnts.table)))
#wvar=rowSums(aam*cross.vars*(as.vector(cross.cnts.table)))/rowSums(cross.vars*(as.vector(cross.cnts.table)))
#hvar=rowSums(h2s[,colnames(aam)]*cross.vars*(as.vector(cross.cnts.table)))/rowSums(cross.vars*(as.vector(cross.cnts.table)))

#aweighted=rep(0,ncol(aam))
#for(i in 1:nrow(aam)) {
#    for(j in 1:ncol(aam)){
#      aweighted[i]= 
#      anorm3[i]=sum(sw[[i]]$additive*cross.cnts.table[sw[[i]]$cross]*cross.vars[i,sw[[i]]$cross])/sum(cross.cnts.table[sw[[i]]$cross]*cross.vars[i,sw[[i]]$cross])
#    
#}
# evaluate results for within-cross QTL model total variance explained
# run code for variance component analysis




# Supplementary Table 2 ------------------------------------------------------------------------------------- 
#load('/data/rrv2/genotyping/RData/FDR_cross_validatedR2.RData')
qcvr=(sapply(cross.validated.r2, function(x) rowMeans(x, na.rm=T)))[-c(37,38),]
#load('/data/rrv2/genotyping/RData/h2_VC_model.RData')
h2s=(sapply(h2.vc.model, function(x) x[,1]))[-c(37,38),]
test=qcvr/h2s
test2=gather(data.frame(test), cross)
test2=data.frame(condition=rep(rownames(test), 16), test2)
#median(qcvr/h2s, na.rm=T)
#[1] 0.6834

# load cross validated joint r^2 results
load('/data/rrv2/genotyping/RData/jointJSCVR2.RData')
aa=do.call('abind', c(jointPeaksCV_R2, along=3))
aam=apply(aa,c(1,2), mean)

  
# for table 
sup.table.2a=rbindlist(lapply(h2.vc.model, function(x) {y=data.frame(x); y$trait=rownames(x); return(y)}), idcol='cross')
sup.table.2a=sup.table.2a[sup.table.2a$trait!='YPD;;2',]
sup.table.2a=sup.table.2a[sup.table.2a$trait!='YPD;;3',]
sup.table.2a=data.frame(sup.table.2a, within_cross_variance_explained_by_qtl=h2sm[,4],  joint_variance_explained_by_qtl=melt(aam)[,3])
names(sup.table.2a)[1:6]=c('cross', 'additive_heritability', 'residual_error', 'additive_heritability_standard_error', 'residual_error_standard_error', 'trait')
sup.table.2a=sup.table.2a[,c(6,1,7,8,2,4,3,5)]                           

sup.table.2b=sup.table.2b[,c(2,4,1,3,5)]
WriteXLS(c('sup.table.2a', 'sup.table.2b'), 
         "/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryTable2.xls", SheetNames=c('additive_only', 'multiple_components'))
#-------------------------------------------------------------------------------------------------------------------
