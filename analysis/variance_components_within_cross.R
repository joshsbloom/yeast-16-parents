# within cross analysis 
cross.VCs=list()
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
            },error=function(e) {
                cross.VCs[[cross.name]][[pheno]]=extractVarCompResultsR(regress(scale(pproc$y)~1,~ZZt+ZQZt+ZAZt+ZQQZt+ZQAZt+ZAAZt, verbose=T))
            })
   }
}
# currently finished B
#save(cross.VCs, file='/data/rrv2/genotyping/RData/withinCrossVCs_FDR.RData')
load(  '/data/rrv2/genotyping/RData/withinCrossVCs_FDR.RData')

#save(cross.VCs, file='/data/rrv2/genotyping/RData/withinCrossVCs.RData')
#load(  '/data/rrv2/genotyping/RData/withinCrossVCs.RData')

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

library(gplots)
#plotCI(Q_A,Q, uiw=Q_A.se, err='x', ylim=c(0,1), xlim=c(0,1))
#plotCI(Q_A,Q, Q.se, err='y', add=T)
plot(Q_A,Q)
abline(0,1)
BYxRM.logical=grepl('^A', names(Q_A))
WCggplot=data.frame(trait=as.vector(do.call('c', sapply(WC.sigma.norm, function(x) names(x[1,])))),
               cross=rep(names(WC.sigma.norm), sapply(WC.sigma.norm, ncol)),
               additive=Q_A, non_additive=Q_AA, additive_plus_non_additive=Q_A+Q_AA)
library(cowplot)
#library(ggextra)
main_plot=ggplot(WCggplot, aes(x=additive_plus_non_additive, y=additive,color=cross=='A',size=cross=='A'))+scale_size_discrete(guide = FALSE) +
    geom_point()+ geom_abline(intercept=0,slope=1)+labs(color="Cross is BYxRM:")
#ggMarginal(main_plot,type='histogram')
inset_plot=ggplot(WCggplot, aes(x=non_additive/additive))+geom_histogram()+xlim(c(0,1))
ggdraw() +
    draw_plot(main_plot)+
    draw_plot(inset_plot, x=.075, y=.6, width=.4, height=.4)

png(file='/data/rrv2/Figures_and_Tables/additive_vs_interaction_scatter.png', width=768, height=768)
plot(Q_A+Q_AA, Q_A,
     pch=20,
     xlim=c(0,1),
     ylim=c(0,1), 
     xlab='Additive + Interaction variance (whole genome)',
     ylab='Additive variance (whole genome)' )
points(c(Q_A+Q_AA)[BYxRM.logical], Q_A[BYxRM.logical], col='red', pch=20)
abline(0,1)
par(fig=c(.075,.5,.5,1),new=T)
hist(Q_AA/Q_A, breaks=100, xlim=c(0,1), main='', xlab='(Interaction variance)/(Additive variance)')
dev.off()
#plot(Q_A, QAAA, pch=20, xlim=c(0,1), ylim=c(0,1))
#points(c(Q_A)[BYxRM.logical], QAAA[BYxRM.logical], col='red', pch=20)
#abline(0,1)
#par(fig=c(.075,.5,.5,1),new=T)
#hist(QAAA/(QAAA+Q_A), breaks=100)





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

colorsR=rev(colorRampPalette(brewer.pal(7, "Paired"))(7))
colorsR[1]='white'
colorsR[2]='mistyrose2' #colorsR[3]
colorsR[3:5]=c('darkseagreen', 'lawngreen' ,'forestgreen')

label_cross=unlist(lapply(crosses.to.parents, function(x) paste(x, collapse= '   x   ')))
png(file='/data/rrv2/Figures_and_Tables/Supp_withinCrossVC.png', width=1920, height=1080)
#pdf(file='/home/jbloom/Desktop/111417/withinCrossVC.pdf', height=8.5*2, width=11*2)
ggplot(WCV, aes(x=Trait,y=fraction_of_variance, fill=Component))+geom_bar(stat="identity")+
    facet_wrap(~Cross, ncol=4, labeller=as_labeller(label_cross))+ylim(-.25,1.25)+
    scale_fill_manual(values =colorsR) +theme_bw()+theme(axis.text.x=element_text(angle=70,hjust=1))+
    theme(legend.key = element_rect(colour = "black"))+
    geom_errorbar(color='grey10', position=position_dodge(width=.5), aes(ymax=ypos+se, ymin=ypos-se, width=0))
dev.off()                   








