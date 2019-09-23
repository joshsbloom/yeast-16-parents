library(openblasctl)
openblas_set_num_threads(48)

load('/data/rrv2/genotyping/RData/jointPeaksJS.RData')
jP=rbindlist(jointPeaksJS, idcol='chromosome')
jPs=split(jP, jP$trait)
#sum(sapply(jPs[-c(37,38)], nrow))

# further functionalize this code for publication
#training.sets
set.seed(20)

# load pheno_extracted
load('/data/rrv2/genotyping/RData/extracted_average_phenotypes.RData')
strain.names=sapply(pheno_extracted, function(x) rownames(x))
strain.names.cv=lapply(strain.names, function(x) {
       cvg=cut(sample(1:length(x)),10)
       levels(cvg)=as.roman(1:10)
       cvg=as.character(cvg)
       return(cvg) })

# output from segregants_hmm.R
load('/data/rrv2/genotyping/RData/parents.list.RData')
# load iseq.freqs
load('/data/rrv2/genotyping/RData/iseq.freqs.RData')

save.prefix='/data/rrv2/genotyping/RData/jointJSCV_'
#save.prefix='/data/home/jbloom/misc/jointJSCV_'
setme=as.character(as.roman(1:10))
for(setm in setme) {
    print(setm)
    prs=mapply(function(x,y){
                  lapply(x,function(z) z[which(!(y%in%setm)),-c(37,38)])
                  },x=pheno.resids,y=strain.names.cv, SIMPLIFY=F)
    srs=mapply(function(x,y){
                  x[,which(!(y%in%setm))]
                  },x=seg.recoded,y=strain.names.cv, SIMPLIFY=F)
    print(sapply(srs, ncol))
    temp=mapJointQTLsJS_variants(n.perm=3e2, FDR_thresh=.05, parents.list, prs, srs, iseq.freqs, filterJS=T)
    mjq=rbindlist(temp,idcol='chromosome')
    saveRDS(mjq, file=paste0(save.prefix,setm,'.RDS'))
    # get coefficients from training data, predict effects on data left out
}
cvjPs=list()
for(setm in setme) {
    print(setm)
    mjq=readRDS(paste0(save.prefix, setm, '.RDS'))
    cvjPs[[setm]]=split(mjq, mjq$trait)
}

# for each cross-validation set
jointPeaksCV_R2=list()
jointPeakEffectsUnbiased=list()
for(setm in setme[10]) {
    print(setm)
    tcross=matrix(NA, length(names(cvjPs[[setm]])), length(crosses))
    rownames(tcross)=names(cvjPs[[setm]])
    colnames(tcross)=crosses
    for(cross.name in crosses) {
         print(cross.name)
         g.s=scale(t(seg.recoded[[cross.name]]))
         #rename columns
         colnames(g.s)=parents.list[[cross.name]]$marker.name
         mPheno=pheno_extracted[[cross.name]]
         mPhenos  = scale(mPheno)

          for(tt in names(cvjPs[[setm]])) {
             sig.joint.markers=parents.list[[cross.name]]$marker.name[ na.omit(match(cvjPs[[setm]][[tt]]$fscan.markers[cvjPs[[setm]][[tt]]$q<.05], parents.list[[cross.name]]$marker.name))]
             bad.markers=duplicated(sig.joint.markers)
             if(sum(bad.markers)>0) {
                 sig.joint.markers=sig.joint.markers[-which(bad.markers)]
             }
             dff=data.frame(g.s[,sig.joint.markers])
             fC=findCorrelation(cor(dff), cutoff=.99)
             if(length(fC)>0) {
                 sig.joint.markers=sig.joint.markers[-fC]
                 dff=data.frame(g.s[,sig.joint.markers])
             }
             X2=dff
             yr=mPhenos[,tt]
             yr.train=yr
             # double check this
             yr.train[strain.names.cv[[cross.name]]==setm]=NA
             yr.test=yr
             yr.test[strain.names.cv[[cross.name]]!=setm]=NA    
             
             if(ncol(X2)==0) { next;}
             # fit model on training data
             fitme=(lm(yr.train~.-1,X2))
             if(is.null(dim(X2))){
                predicted=data.matrix(X2[strain.names.cv[[cross.name]]==setm])*coef(fitme)
             }else {
                 # estimate variance explained  
                predicted=data.matrix(X2[strain.names.cv[[cross.name]]==setm,])%*%coef(fitme)
            }
            tcross[tt,cross.name]=cor(yr.test[strain.names.cv[[cross.name]]==setm], predicted)^2
            
            fitme=(lm(yr.test~.-1,X2))

            #output unbiased estimates of variant effects 
            if(setm == 'X') {
            jointPeakEffectsUnbiased[[tt]][[cross.name]]=data.frame(
                trait=tt,
                cross=cross.name,
                peaks=sig.joint.markers,
                betas=as.vector(coef(fitme)),
                #vexp=a.effs,
                #p = ps,
                maf1012=iseq.freqs$maf1012[match(sig.joint.markers, iseq.freqs$marker)],
                # iseq.freqs[,8] is the allele frequency of the alternate allele 
                alt012=iseq.freqs[,8][match(sig.joint.markers, iseq.freqs$marker)],
                cross.cnt=cross.count.lookup[match(sig.joint.markers, rownames(cross.count.lookup)),'values'],
                stringsAsFactors=F)
             print(jointPeakEffectsUnbiased[[tt]][[cross.name]])
            }
        }
    }
    jointPeaksCV_R2[[setm]]=tcross
}
#save(jointPeaksCV_R2, file='/data/rrv2/genotyping/RData/jointJSCVR2.RData')

# make plot of MAF in JS panel vs effect size within each cross -----------------------------------------------------------------------

library(data.table)
r3=rbindlist(lapply(jointPeakEffectsUnbiased, rbindlist))
r3=r3[!(r3$trait %in% c("YPD;;2","YPD;;3")),]
r3$maf1012Fill=r3$maf1012
r3$crossCount=ifelse(r3$cross.cnt%%2==1, r3$cross.cnt+1, r3$cross.cnt)
r3$maf1012Fill[is.na(r3$maf1012Fill)]=c((r3$crossCount/2)/1012)[is.na(r3$maf1012Fill)]
#r3$density=get_density(r3$maf1012, abs(r3$betas))
#r3$densityF=get_density(r3$maf1012Fill, abs(r3$betas))
r3$absBeta=abs(r3$betas)
r3$ancestral=iseq.freqs$ancestral[match(r3$peaks, iseq.freqs$marker)]
r3$maf2=ifelse(r3$ancestral, r3$alt012, 1-r3$alt012)

r3MeanEffect=r3 %>% group_by_(.dots=c("trait","peaks")) %>% 
    mutate(MeanAbsBeta=mean(betas)) %>% distinct(traits,peaks,maf1012,crossCount, MeanAbsBeta, ancestral, maf2)
#r3VarEffect=r3 %>% group_by_(.dots=c("trait","peaks")) %>% mutate(VarAbsBeta=sd(betas)) %>% distinct(traits,peaks,maf1012,crossCount, VarAbsBeta, ancestral, maf2)


# 0 = reference derived, alternate paradoxus (ancestral)
# 1 = reference paradoxus (ancestral), alternate cerevisiae
library(ggplot2)
sfigUnbiased1=ggplot(r3MeanEffect)+geom_point(alpha=.4, size=.75, aes(x=maf1012,y=abs(MeanAbsBeta)))+ #,color=density))+
    #scale_color_viridis(option = "inferno", direction=1,end=1)+
    scale_x_continuous(name='Minor allele frequency', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='Average absolute effect, SD units', limits=c(0,1), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()+
    theme(panel.grid.major = element_line(colour = "#80808022"))+
    theme(panel.grid.minor = element_line(colour = "#80808022"))+
    theme(axis.text.x=element_text(size=rel(1),color='black'),
          axis.text.y=element_text(size=rel(1),color='black'))+
    stat_summary_bin(aes(x=maf1012,y=abs(MeanAbsBeta)), breaks=cut2(r3MeanEffect$maf1012,g=42, onlycuts=T),  col='red')

sfigUnbiased2=ggplot(r3MeanEffect)+geom_point(alpha=.4, size=.75, aes(x=maf2,y=MeanAbsBeta))+
     #scale_color_viridis()+
    scale_x_continuous(name='Derived allele frequency', breaks=seq(0,1,.1), expand=c(0.01,0)) +
    scale_y_continuous(name='Average effect, SD units', limits=c(-1,1), breaks=seq(-1,1,.1), expand=c(0,0))+theme_bw()
ggarrange(sfigUnbiased1,sfigUnbiased2,ncol=2, nrow=1, labels=c('A','B'))
sfigUnbiased1
#ggsave(file='/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/ReviewerComment4_1.png',width=11,height=5.5)
ggsave(file='/home/jbloom/Dropbox/Manuscripts/RR/elife/Figure3-figure_supplement_4.pdf', width=5.5, height=5.5)
#----------------------------------------------------------------------------------------------------------------------------------------






# see jointModelAF.R 
sr2=r2MeanEffect
sr2$af=(r2MeanEffect$crossCount/32)
sr2$ve2=2*(r2MeanEffect$crossCount/32)*(1-(r2MeanEffect$crossCount/32))*r2MeanEffect$MeanAbsBeta^2
# variance explained within 
sr2$veobs=2*(r2MeanEffect$maf1012)*(1-(r2MeanEffect$maf1012))*r2MeanEffect$MeanAbsBeta^2
sr2$veexp=2*(r2MeanEffect$maf1012)*(1-(r2MeanEffect$maf1012))*median((r2MeanEffect$MeanAbsBeta^2)) #[r2MeanEffect$maf1012>.05])
sr2=sr2[order(sr2$maf1012),]
sr2$cveRR=cumsum(sr2$ve2)/sum(sr2$ve2)
sr2$cveobs=cumsum(sr2$veobs)/sum(sr2$veobs)

plot(sr2$maf1012, sr2$cveRR)
points(sr2$maf1012, sr2$cveobs, col='red')

vexplainedRR=ggplot(sr2)+geom_jitter(alpha=.4, size=.75, aes(x=af,y=ve2))+
     #scale_color_viridis()+
    scale_x_continuous(name='Minor allele frequency', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='Variance explained', limits=c(0,.1), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()
   # stat_summary_bin(aes(x=af,y=ve2), breaks=cut2(r3MeanEffect$af,g=12, onlycuts=T),  col='red')

vexplainedJS=ggplot(sr2)+geom_point(alpha=.4, size=.75, aes(x=maf1012,y=ve3))+
     #scale_color_viridis()+
    scale_x_continuous(name='Minor allele frequency', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='Variance explained', limits=c(0,.005), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()+
    stat_summary_bin(aes(x=maf1012,y=ve3), breaks=cut2(r3MeanEffect$maf1012,g=42, onlycuts=T),  col='red')+
    stat_summary_bin(aes(x=maf1012,y=ve4), breaks=cut2(r3MeanEffect$maf1012,g=42, onlycuts=T),  col='blue')+


ggarrange(vexplainedRR,vexplainedJS,ncol=2, nrow=1, labels=c('A','B'))
ggsave(file='/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/ReviewerComment4_2.png',width=11,height=5.5)

test=sr2
test=test[order(test$maf1012),]
plot(test$maf1012, cumsum(test$ve3)/sum(test$ve3), col='red', ylab='Cummulative variance explained', xlab='Minor Allele Frequency')
points(test$maf1012, cumsum(test$ve4)/sum(test$ve4), col='blue')


sr4=split(sr2, sr2$trait)
#names(sr2)=as.character(levels(WCV$Trait))
#pdf('/home/jbloom/Dropbox/RR/Figures and Tables/SupplementaryFigure6.pdf', width=15, height=20) 
par(mfrow=c(7,6))
for(i in 1:38) {
sr21=sr4[[i]]
sr21=sr21[order(sr21$maf1012, decreasing=F),]
sr21$cve=cumsum(sr21$ve3)/sum(sr21$ve3)
plot(sr21$maf1012, sr21$cve, main=names(sr4)[i], ylab='cummulative GVE', 
     xlab='MAF', xlim=c(0,.5), ylim=c(0,1), type='l', lwd=2)
#plot(density(log10(sr21$maf1012),sr21$ve, n=3))
abline(0,2, lty=2, col='grey')
#readline()
}


sr3=r3MeanEffect
# variance explained within panel
sr3$ve2=2*(r3MeanEffect$crossCount/32)*(1-(r3MeanEffect$crossCount/32))*r3MeanEffect$MeanAbsBeta^2
# variance explained within 
sr3$ve3=2*(r3MeanEffect$maf1012)*(1-(r3MeanEffect$maf1012))*r3MeanEffect$MeanAbsBeta^2


vexplainedJSUnbiased=ggplot(sr3)+geom_point(alpha=.4, size=.75, aes(x=maf1012,y=ve3))+
     #scale_color_viridis()+
    scale_x_continuous(name='Minor allele frequency', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='Variance explained', limits=c(0,.1), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()+
    stat_summary_bin(aes(x=maf1012,y=ve3), breaks=cut2(r3MeanEffect$maf1012,g=42, onlycuts=T),  col='red')

vexplainedRRUnbiased=ggplot(sr3)+geom_point(alpha=.4, size=.75, aes(x=maf1012,y=ve2))+
     #scale_color_viridis()+
    scale_x_continuous(name='Minor allele frequency', breaks=seq(0,1,.05), expand=c(0.01,0)) +
    scale_y_continuous(name='Variance explained', limits=c(0,.1), breaks=seq(0,1,.1), expand=c(0,0))+theme_bw()+
    stat_summary_bin(aes(x=maf1012,y=ve2), breaks=cut2(r3MeanEffect$maf1012,g=42, onlycuts=T),  col='red')

ggarrange(vexplainedRRUnbiased,vexplainedJSUnbiased,ncol=2, nrow=1, labels=c('A','B'))
ggsave(file='/home/jbloom/Dropbox/Manuscripts/RR/Figures and Tables/ReviewerComment4_3_unbiased.png',width=11,height=5.5)
    
