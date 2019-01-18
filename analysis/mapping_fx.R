fasterLOD=function(n.pheno, pheno.s,gdata.s, betas=FALSE, sdx=1, pheno=NULL){
   r=crossprod(pheno.s, gdata.s)/(n.pheno-1)
   LOD=(-n.pheno*log(1-r^2))/(2*log(10))
   if(betas==FALSE) {
       return(LOD)
   } else {
      # beta=r*apply(cbind(pheno),2, sd,na.rm=T)/sdx
       return(list(r=r, LOD=LOD))
   }
}

fasterR2=function(n.pheno, pheno.s,gdata.s, betas=FALSE, sdx=1, pheno=NULL){
   r=(crossprod(pheno.s, gdata.s)/(n.pheno-1))
   tt=(r/sqrt(1-r^2))*sqrt(n.pheno-2)
   return(2*pt(-abs(tt),df=n.pheno-2))
   #LOD=(-n.pheno*log(1-r^2))/(2*log(10))
   #if(betas==FALSE) {
   #    return(LOD)
   #} else {
   #   # beta=r*apply(cbind(pheno),2, sd,na.rm=T)/sdx
   #    return(list(r=r, LOD=LOD))
   #}
}

processPhenos_for_MM=function(s, all.strain.names){
    #s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    s=s[match(all.strain.names, names(s))]
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'
    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------
    y=ny[!is.na(ny)]
    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
   
    #for constructing Strain Variance component
    Strain      = diag(strain.cnt)
    Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 
    strains.with.phenos=match(unique.sn, all.strain.names)
    n.strains=length(strains.with.phenos)
    return(list(y=y, Z=Z, Strain=Strain, strains.with.phenos=strains.with.phenos, n.strains=n.strains))
}  

avg.nvec=function(y, un=NULL){
    if(is.null(un)) {un=names(y) }
    y.avg.l=(by(y, un, mean, na.rm=T))
    nvec=names(y.avg.l)
    y.avg.l=as.vector(y.avg.l)
    names(y.avg.l)=nvec
    return(y.avg.l)
}


m.S=function (y, K = NULL, bounds = c(1e-09, 1e+09), theta=NULL, Q=NULL, X=NULL ) 
{
    n <- length(y)
    y <- matrix(y, n, 1)
    if(is.null(X) ) {  p <- 1    } else { p = ncol(X) }
    Z <- diag(n)
    m <- ncol(Z)
       
    omega <- crossprod(Q, y)
    omega.sq <- omega^2
    
    f.REML <- function(lambda, n.p, theta, omega.sq) {
        n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
    }
    soln <- optimize(f.REML, interval = bounds, n - p, theta,  omega.sq)
    lambda.opt <- soln$minimum
    
    df <- n - p
    Vu.opt <- sum(omega.sq/(theta + lambda.opt))/df
    Ve.opt <- lambda.opt * Vu.opt
    VCs=c(Vu.opt, Ve.opt)
    return(VCs)
}
#---------------------------------------------------------------------------------------------------------------

doEigenA_forMM=function(pheno.scaled,A ,X=NULL ) {
        n=nrow(pheno.scaled)
        if(is.null(X) ) {  X = matrix(rep(1, n), n, 1); p=1 } else {p=ncol(X) }
        XtX = crossprod(X, X)
        XtXinv = solve(XtX)
        S = diag(n) - tcrossprod(X %*% XtXinv, X)
        SHbS = S %*% A %*% S
        SHbS.system = eigen(SHbS, symmetric = TRUE)
        theta = SHbS.system$values[1:(n - p)] 
        Q = SHbS.system$vectors[, 1:(n - p)]
        return(list(theta=theta, Q=Q))
}
calc.BLUPS= function(G,Z,Vinv,y,X,B ){    G%*%crossprod(Z,Vinv)%*%(y- X%*%B)    }


find.background.QTL = function (t.tpm.matrix, pheno.scaled, gdata, gdata.scaled ,bck.thresh=3.5,
                                chromosomes=NULL, covariates=NULL) {
    if(is.null(covariates)) { covariates=matrix(1, nrow(t.tpm.matrix), 1) }

    if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
    cvec=(do.call('rbind', strsplit(colnames(gdata), '_'))[,1])

    # Iteration 1  (find background QTL effects)
    LODS=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
    # lists to hold LODs, markers, and scaled markers by chromosome----------------
    LODS.by.chr=list()
    gdata.s.by.chr=list()
    gdata.by.chr=list()
    for(cc in chromosomes) {   
            LODS.by.chr[[cc]]=LODS[,which(cvec %in% cc)]   
            gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
            gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
    }
    #-----------------------------------------------------------------------------
     
    # max LOD per chromosome
    max.lod.chr=lapply(LODS.by.chr, function(x) { apply(x,1,max) })
    # index of max LOD per chromosome
    max.lod.chr.ind=lapply(LODS.by.chr, function(x) { apply(x,1,which.max) })

    # pulls out the marker index (by chromosome) of the peak marker 
    pmarker.chr.ind =  mapply(function(x,y,z) { 
        h=y[x>bck.thresh]
        i=colnames(z)[h]
        names(i)=names(h)
        return(i)
        },  max.lod.chr, max.lod.chr.ind,gdata.by.chr)

    pmarker.chr.ind.df=data.frame(chr=rep(names(pmarker.chr.ind), sapply(pmarker.chr.ind, length)),
             gene=as.vector(do.call('c', lapply(pmarker.chr.ind, function(x) names(x)))),
             marker=as.vector(do.call('c', lapply(pmarker.chr.ind, function(x) (x)))), stringsAsFactors=F)
    pbg=split(pmarker.chr.ind.df$marker, pmarker.chr.ind.df$gene)

    # Iteration 2--------------------------------------------------------------------------------------------
    plist=lapply(colnames(t.tpm.matrix), function(i) { 
         if(!(i %in% names(pbg)) ) {  return(pheno.scaled[,i]) } 
         else{ residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,pbg[[i]]]) )  }
    })
    p.resid1=do.call('cbind', plist)
    colnames(p.resid1)=colnames(pheno.scaled)
    pheno.scaled=scale(p.resid1)
    LODS1=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
    LODS1.by.chr=list()
    for(cc in chromosomes) { LODS1.by.chr[[cc]]=LODS1[,which(cvec %in% cc)]   }
    max.lod1.chr = lapply(LODS1.by.chr, function(x) { apply(x,1,max) })
    max.lod1.chr.ind = lapply(LODS1.by.chr, function(x) { apply(x,1,which.max) })
    pmarker1.chr.ind = mapply(function(x,y,z) { 
        h=y[x>bck.thresh]
        i=colnames(z)[h]
        names(i)=names(h)
        return(i)
        },  max.lod1.chr, max.lod1.chr.ind,gdata.by.chr)
    pmarker1.chr.ind.df=data.frame(chr=rep(names(pmarker1.chr.ind), sapply(pmarker1.chr.ind, length)),
             gene=as.vector(do.call('c', lapply(pmarker1.chr.ind, function(x) names(x)))),
             marker=as.vector(do.call('c', lapply(pmarker1.chr.ind, function(x) (x)))), stringsAsFactors=F)
    pbg1=split(pmarker1.chr.ind.df$marker, pmarker1.chr.ind.df$gene)
    pbg.merge=pbg
    for(i in names(pbg1)){    pbg.merge[[i]]=c(pbg.merge[[i]], pbg1[[i]]) }

    # Iteration 3--------------------------------------------------------------------------------------------
    plist=lapply(colnames(t.tpm.matrix), function(i) { 
         if(!(i %in% names(pbg.merge)) ) {  return(pheno.scaled[,i])} 
         else{ residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,pbg.merge[[i]]]) )    }
    })
    p.resid2=do.call('cbind', plist)
    colnames(p.resid2)=colnames(pheno.scaled)
    pheno.scaled=scale(p.resid2)
    LODS2=fasterLOD(nrow(pheno.scaled),pheno.scaled,gdata.scaled)
    LODS2.by.chr=list()
    for(cc in chromosomes) {    LODS2.by.chr[[cc]]=LODS2[,which(cvec %in% cc)]   }
    max.lod2.chr=lapply(LODS2.by.chr, function(x) { apply(x,1,max) })
    max.lod2.chr.ind=lapply(LODS2.by.chr, function(x) { apply(x,1,which.max) })
    pmarker2.chr.ind =  mapply(function(x,y,z) { 
        h=y[x>bck.thresh]
        i=colnames(z)[h]
        names(i)=names(h)
        return(i)
        },  max.lod2.chr, max.lod2.chr.ind,gdata.by.chr)
    pmarker2.chr.ind.df=data.frame(chr=rep(names(pmarker2.chr.ind), sapply(pmarker2.chr.ind, length)),
             gene=as.vector(do.call('c', lapply(pmarker2.chr.ind, function(x) names(x)))),
             marker=as.vector(do.call('c', lapply(pmarker2.chr.ind, function(x) (x)))), stringsAsFactors=F)
    pbg2=split(pmarker2.chr.ind.df$marker, pmarker2.chr.ind.df$gene)

    pbg.final=pbg.merge
    for(i in names(pbg2)){    pbg.final[[i]]=c(pbg.final[[i]], pbg2[[i]]) }

    pm.chr.ind.df=rbind(pmarker.chr.ind.df, pmarker1.chr.ind.df, pmarker2.chr.ind.df)
    background.QTL=list()
    for(cc in chromosomes) {
        background.QTL[[cc]]=split(pm.chr.ind.df$marker[pm.chr.ind.df$chr!=cc], pm.chr.ind.df$gene[pm.chr.ind.df$chr!=cc]) 
    }
    return(background.QTL)
    }





mapQTL=function( background.QTL, t.tpm.matrix, pheno.scaled, 
                gdata, gdata.scaled, 
                chromosomes=NULL,
                n.perm=1000, FDR.thresh=.05, Arandom=TRUE, covariates=NULL) {

    if(is.null(covariates)) { covariates=matrix(1, nrow(t.tpm.matrix), 1) }
    if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
    cvec=(do.call('rbind', strsplit(colnames(gdata), '_'))[,1])

    gdata.s.by.chr=list()
    gdata.by.chr=list()
    for(cc in chromosomes) {   
            gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec %in% cc)]   
            gdata.by.chr[[cc]]=gdata[,which(cvec %in% cc)]   
    }

    peakList=list()
    #LODmatrix=list()
    for(cc in chromosomes) {  
         print(cc) 
         plist=lapply(colnames(t.tpm.matrix), function(i) { 
              if(!(i %in% names(background.QTL[[cc]])) ) { return(  residuals(lm(t.tpm.matrix[,i]~covariates)) ) }
              else{  residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,background.QTL[[cc]][[i]]]) ) }
         })
         preal=do.call('cbind', plist)
         colnames(preal)=colnames(t.tpm.matrix)
         
         presid=preal

         # remove residual additive effect
         if(Arandom) {
            cloco=match(cc, names(gdata.by.chr))
            clocoM=do.call('cbind', gdata.s.by.chr[-cloco] )
            Aloco=tcrossprod(clocoM)/ncol(clocoM)
            eigA=doEigenA_forMM(presid,Aloco)
            svdAloco=svd(Aloco)
            pb=txtProgressBar(min=1, max=ncol(presid), style=3)
            for(tp in 1:ncol(presid)){
                 setTxtProgressBar(pb,tp)
                rr=m.S(presid[,tp], K=Aloco,  theta=eigA$theta, Q=eigA$Q)
                ###W=solve(rr[1]*Aloco+rr[2]*diag(nrow(presid)))
                W=svdAloco$u %*% tcrossprod(diag(1/((svdAloco$d*rr[1])+(rr[2]))), svdAloco$v)
                if(rr[1]>0) {
                      blups=calc.BLUPS(rr[1]*Aloco,diag(nrow(presid)),W,presid[,tp],matrix(1,nrow(presid),1),0 )[,1]
                      presid[,tp]=as.vector(presid[,tp] - blups)
                }
            }
            rm(W)
            close(pb)
         }
         
         preal=presid
         presid=scale(presid)
         
         #pscale=pheno.scaled
         obsLOD = fasterLOD(nrow(preal),presid,gdata.s.by.chr[[cc]],betas=TRUE,pheno=preal)
         findingPeaks=TRUE
         jump=1

        while(findingPeaks==TRUE) {
            max.obsLOD=apply(obsLOD$LOD,1,max)   
            max.ind.obsLOD=apply(obsLOD$LOD,1,which.max)   

            permLOD=matrix(0, length(max.obsLOD), n.perm)
           
           # print('doing permutations')
           # pb =txtProgressBar(min = 1, max =n.perm, style = 3)
            for(i in 1:n.perm) {
           #      setTxtProgressBar(pb, i)
                 permLOD[,i]=apply(fasterLOD(nrow(presid),presid[sample(nrow(presid)),], gdata.s.by.chr[[cc]]),1,max)
            }
           #close(pb)
           
           obsPcnt = sapply(seq(1.5, 9, .05), function(thresh) { sum(max.obsLOD>thresh) }   )
           names(obsPcnt) = seq(1.5, 9, .05)
           if(sum(obsPcnt)<3) {break}
           # expected number of QTL peaks with LOD greater than threshold
           expPcnt = sapply(seq(1.5, 9, .05),  
                             function(thresh) { 
                                    #print(thresh); 
                                    mean(apply(permLOD, 2, function(ll) {sum(ll>thresh) }) )
                                } )
          names(expPcnt) = seq(1.5, 9, .05)
          pFDR = expPcnt/obsPcnt
          pFDR[is.na(pFDR)]=0
          pFDR[!is.finite(pFDR)]=0
          #to make sure this is monotonic
          pFDR = rev(cummax(rev(pFDR)))
          fdrFX=approxfun(pFDR, seq(1.5,9,.05))
          thresh=fdrFX(FDR.thresh)
          print(paste('FDR Thresh 5%', thresh)) 
          peak.gene.ind=max.ind.obsLOD[which(max.obsLOD>thresh)]

          if(length(peak.gene.ind)>0) { 
              findingPeaks=TRUE 
              pldf=data.frame(gene=names(peak.gene.ind),
                              pmarker=colnames(gdata.by.chr[[cc]])[peak.gene.ind], 
                              pcind=as.vector(peak.gene.ind),
                              stringsAsFactors=F)
              pldf$r=obsLOD$r[cbind(pldf$gene, pldf$pmarker)]
              pldf$LOD=obsLOD$LOD[cbind(pldf$gene, pldf$pmarker)]
              lodInt=matrix(NA,nrow(pldf),2)
              #Lmat=obsLOD$LOD[pldf$gene,]
              for(peak in 1:nrow(pldf)){
                 int.range=range(which(obsLOD$LOD[pldf[peak,]$gene,]>pldf[peak,]$LOD-1.5)) 
                 lodInt[peak,]=colnames(gdata.by.chr[[cc]])[int.range]
                 #if(is.null(nrow(Lmat))) {Lmat[!(c(1:length(Lmat)) %in% seq(int.range[1], int.range[2]))]=0 }
                 #else {Lmat[peak, !(c(1:ncol(Lmat)) %in% seq(int.range[1], int.range[2]))]=0 }
              }
              pldf$CI.l=lodInt[,1]
              pldf$CI.r=lodInt[,2]
              peakList[[cc]][[jump]]=pldf
              
              #LODmatrix[[cc]][[jump]]=Lmat
              presid.tmp=matrix(0, nrow(t.tpm.matrix), nrow(pldf))
              rownames(presid.tmp)=rownames(t.tpm.matrix)
              colnames(presid.tmp)=pldf$gene
             
              pldfm = do.call('rbind', peakList[[cc]])
              found.q.covs=split(pldfm$pmarker, pldfm$gene)
              
              for(i in pldf$gene){
                  # modified 041017
                  #if(!(i %in% names(background.QTL[[cc]])) ) {
                  #    presid.tmp[,i]=scale(residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,found.q.covs[[i]]]))) }
                  #else{
                  #   presid.tmp[,i]=scale(residuals(lm(t.tpm.matrix[,i]~covariates+gdata[,background.QTL[[cc]][[i]]]+gdata[,found.q.covs[[i]]])) ) }
                  presid.tmp[,i]=scale(residuals(lm(presid[,i]~gdata[,found.q.covs[[i]]])))
              }
              obsLOD = fasterLOD(nrow(presid.tmp),presid.tmp,gdata.s.by.chr[[cc]],betas=TRUE, pheno=presid.tmp)
              presid=presid.tmp
              jump=jump+1 
          } else {findingPeaks=FALSE} 
    }
    }
    return(peakList)
}



# calculate chromosome residuals
calcChromosomeResiduals=function(per_cross_peaks, pmatrix, gdata, gdata.scaled, chromosomes=NULL, covariates=NULL) {

    #pmatrix=mPheno
    #per_cross_peaks=cross.peaks[['A']]
    if(is.null(covariates)) { covariates=matrix(1, nrow(pmatrix), 1) }
    if(is.null(chromosomes)) { chromosomes=paste0('chr', as.roman(1:16)) } 
    cvec1=(do.call('rbind', strsplit(colnames(gdata), '_'))[,1])
    cvec2=(do.call('rbind', strsplit(colnames(gdata.scaled), '_'))[,1])

    gdata.s.by.chr=list()
    gdata.by.chr=list()
    for(cc in chromosomes) {   
            gdata.s.by.chr[[cc]]=gdata.scaled[,which(cvec2 %in% cc)]   
            gdata.by.chr[[cc]]=gdata[,which(cvec1 %in% cc)]   
    }
   
     # might as well ignore segregants that never have a phenotype
     #gdata=gdata[rownames(gdata) %in% all.strain.names.S,]
     #all.strain.names.S=all.strain.names.S[all.strain.names.S %in% rownames(gdata)]
     #print(identical(rownames(gdata), snames))

     # build a data structure like background.QTL that has a a list of each chromosome and the detected QTL not found on that chromosome 
     detected.QTL=list()
        for(crc in chromosomes) { 
            per_cross_peaks$chr=sapply(strsplit(per_cross_peaks$pmarker, '_'), function(x) x[1]) 
            detected.QTL[[crc]]= split(per_cross_peaks$pmarker[per_cross_peaks$chr!=crc], per_cross_peaks$trait[per_cross_peaks$chr!=crc]) 
        }
        

      print('calculating strain BLUP residuals')
      # this could be parallelized 
      blup.resid=list()
      for(chc in chromosomes) {  
             print(chc) 
             plist=lapply(colnames(pmatrix), function(i) { 
                  if(!(i %in% names(detected.QTL[[chc]])) ) {
                    return( pmatrix[,i]) }   else{  residuals(lm(pmatrix[,i]~gdata[,detected.QTL[[chc]][[i]]]) ) }
             })
             preal=do.call('cbind', plist)
             colnames(preal)=colnames(pmatrix)
             
            presid=preal
            # avoid unecessary scaling
            # presid=scale(preal)
            cloco=match(chc, names(gdata.by.chr))
            clocoM=do.call('cbind', gdata.s.by.chr[-cloco] )
            Aloco=tcrossprod(clocoM)/ncol(clocoM)

            # Eskin trick to speed up variance component calculation
            eigA=doEigenA_forMM(presid,Aloco)
            svdAloco=svd(Aloco)

            # calculate mixed model, one term for additive variance  -------------------------------------------
            pb=txtProgressBar(min=1, max=ncol(pmatrix), style=3)
            for(tp in 1:ncol(pmatrix)){
                setTxtProgressBar(pb,tp)
                #print(tp)
                # this was pheno.scaled[,tp] ... that can't be right???
                rr=m.S(presid[,tp], K=Aloco,  theta=eigA$theta, Q=eigA$Q)
                W=svdAloco$u %*% tcrossprod(diag(1/((svdAloco$d*rr[1])+(rr[2]))), svdAloco$v)
                if(rr[1]>0) {
                    blups=calc.BLUPS(rr[1]*Aloco,diag(nrow(presid)),W,presid[,tp],matrix(1,nrow(presid),1),0 )[,1]
                    # presid[,tp]=as.vector(scale(presid[,tp]-blups)[,1])
                    presid[,tp]=as.vector(presid[,tp] - blups)

                }
            }
            rm(W)
            blup.resid[[chc]]=presid
            close(pb)
        }
    return(blup.resid)

    }

mapJointQTLs=function(n.perm=100, FDR.thresh=.05, parents.list, pheno.resids, seg.recoded) {
    jointPeaks=list()
    jointLODs=list()

    for( chrom in unique.chrs) {
        plist.chr=sapply(parents.list, function(x) x[x$chr==chrom,])
        chr.markers=sapply(parents.list, function(x){ return(x$marker.name.n[x$chr==chrom]) } )
        seg.chr.phenos=lapply(pheno.resids, function(x) { x[[chrom]] })
        seg.chr.n=lapply(seg.chr.phenos, nrow)
        seg.chr.phenos.scaled=lapply(seg.chr.phenos,scale)

        seg.chr.genos=mapply(function(cmm,seg.r) { t(seg.r[cmm,])} , cmm=chr.markers, seg.r=seg.recoded)
        seg.chr.genos.scaled=lapply(seg.chr.genos, scale)
        cmm=sapply(parents.list, function(x){ return(x$marker.name[x$chr==chrom]) } )
        #rename seg.chr.genos.scaled
        seg.chr.genos.scaled=mapply(function(cmm,seg.r) { colnames(seg.r)=cmm; return(seg.r);}, cmm=cmm, seg.r=seg.chr.genos.scaled)

        mpmarkers=melt(cmm)
        mpmarkers[,1]=as.character(mpmarkers[,1])
        names(mpmarkers)=c('marker', 'cross')

       mpmarkers=separate(mpmarkers, marker, c('chr', 'pos', 'ref', 'alt'), sep='_', convert=T, remove=F)
       mpmarkers=mpmarkers[order(mpmarkers$pos),]

       ccnt=split(mpmarkers$cross, mpmarkers$marker)
       crosses.per.marker=sapply(ccnt, function(x) unique(x))
       cross.cnt.per.marker=sapply(crosses.per.marker, length)
       seg.rare=names(cross.cnt.per.marker[cross.cnt.per.marker<3])

       mm=mpmarkers[match(names(crosses.per.marker), mpmarkers$marker),-6]
       mm=mm[order(mm$pos),]

       # mapping cross marker index to joint marker table  
       mupos=lapply(plist.chr, function(x) match(x$marker.name, mm$marker))
       # mapping cross marker to joint marker table name 
       mucpos = mapply(function(m, p) { p$marker.name[ match(mm$marker[m], p$marker.name)] }, m=mupos, p=plist.chr )
       peakList=list()
       presid=seg.chr.phenos.scaled

       jointLOD=mapply(function(n,p,g) {
                             fasterLOD(n,p,g)
                             }, n=seg.chr.n, p=presid, g=seg.chr.genos.scaled)
       jL=matrix(0,41, nrow(mm))
       rownames(jL)=colnames(seg.chr.phenos[[1]])
       colnames(jL)=mm$marker
       upos=colnames(jL)

       for(cc in crosses){  
          censor.me=jointLOD[[cc]]
          jL[,mupos[[cc]]]=jL[,mupos[[cc]]]+censor.me 
       }
    
       findingPeaks=TRUE
       jump=1
       while(findingPeaks==TRUE) {
            max.obsLOD=apply(jL, 1, max)
            max.obsLOD.indG=apply(jL, 1, which.max)
            max.obsLOD.uposG=upos[apply(jL, 1, which.max)]
            max.ind.obsLOD=sapply(mucpos, function(x) sapply(max.obsLOD.uposG, function(y) which(x==y)))
            rownames(max.ind.obsLOD)=rownames(jL)
            permLOD=matrix(0, length(max.obsLOD), n.perm)
            pb =txtProgressBar(min = 1, max =n.perm, style = 3)
            for(i in 1:n.perm) {
                setTxtProgressBar(pb, i)
                jointLODp=mapply(function(n,p,g) {
                              fasterLOD(n,p[sample(n),],g)
                              }, n=seg.chr.n, p=presid, g=seg.chr.genos.scaled)
                jLp=matrix(0,ncol(presid[[1]]), length(upos))
                rownames(jLp)=colnames(presid[[1]])
                for(cc in crosses){  jLp[,mupos[[cc]]]=jLp[,mupos[[cc]]]+jointLODp[[cc]] }
                permLOD[,i]=apply(jLp, 1, max)
            }
            close(pb)

            obsPcnt = sapply(seq(1.5, 9, .05), function(thresh) { sum(max.obsLOD>thresh) }   )
                       names(obsPcnt) = seq(1.5, 9, .05)
                       if(sum(obsPcnt)<3) {break}
                       # expected number of QTL peaks with LOD greater than threshold
            expPcnt = sapply(seq(1.5, 9, .05),  
                                         function(thresh) { 
                                                #print(thresh); 
                                                mean(apply(permLOD, 2, function(ll) {sum(ll>thresh) }) )
                                            } )
            names(expPcnt) = seq(1.5, 9, .05)
            #plot(seq(1.5, 9, .05), expPcnt/obsPcnt)
            pFDR = expPcnt/obsPcnt
            pFDR[is.na(pFDR)]=0
            pFDR[!is.finite(pFDR)]=0
                      #to make sure this is monotonic
            pFDR = rev(cummax(rev(pFDR)))
            fdrFX=approxfun(pFDR, seq(1.5,9,.05))
            thresh=fdrFX(FDR.thresh)
            print(paste('FDR Thresh 5%', thresh)) 

            peak.gene.ind=max.ind.obsLOD[which(max.obsLOD>thresh),]
            if(is.null(dim(peak.gene.ind)) ) {
                peak.gene.ind=matrix(peak.gene.ind,1)
                rownames(peak.gene.ind)=names(which(max.obsLOD>thresh))
                colnames(peak.gene.ind)=colnames(max.ind.obsLOD)
                keep.traits=names(which(max.obsLOD>thresh))
            } else {    keep.traits=names(apply(peak.gene.ind, 1, function(x) length(unlist(x)))>0) }
            print(length(keep.traits))
            if(nrow(peak.gene.ind)>0 & (length(keep.traits)>0)) {
                jointLODs[[chrom]][[jump]]=jointLOD
                #tr=1
                #plot(upos, jL[tr,], col='red', main=paste(rownames(jL)[tr], jump), ylab='LOD', sub=max.obsLOD.uposG[tr])
                #for(i in 1:16) {
                #    points(upos[mupos[[i]]], jointLOD[[i]][tr,] )#, ylim=c(0,210))               
                # 
                # }
                # abline(v=max.obsLOD.uposG[tr], col='red')
                
                  findingPeaks=TRUE 
                  peak.gene.markers=matrix(NA, nrow(peak.gene.ind), ncol(peak.gene.ind))
                  rownames(peak.gene.markers)=rownames(peak.gene.ind)
                  colnames(peak.gene.markers)=colnames(peak.gene.ind)

                  for(cc in crosses) {
                     if(nrow(peak.gene.ind)==1) {
                        peak.gene.markers[,cc]=mucpos[[cc]][unlist(peak.gene.ind[,cc])[1]]
                     } else {
                     peak.ind=na.omit(sapply(peak.gene.ind[,cc], function(x) x[1]))
                        if(length(peak.ind)>0) {
                             peak.gene.markers[names(peak.ind),cc]=mucpos[[cc]][peak.ind]
                        }
                     }
                  }
                
                  peak.gene.markers=melt(peak.gene.markers)
                  peak.gene.markers=data.frame(na.omit(apply(peak.gene.markers,2, as.character)), stringsAsFactors=F)
                  names(peak.gene.markers)=c('trait', 'cross', 'marker')
                  peak.gene.markers$jointPeak=paste0(chrom, '_', jump)
                  peakList[[jump]]=peak.gene.markers
                
                 pldfm = do.call('rbind', peakList)
                 found.t.covs=split(pldfm, pldfm$trait)

                 # at this point remove traits for which no linkages were found in the current step
                 if(length(keep.traits)==1){
                        presid.tmp=lapply(presid, function(x) { z=matrix(x[,keep.traits], 1);
                                          rownames(z)=keep.traits; 
                                          colnames(z)=rownames(x);
                                          return(t(z));})
                 }else {        
                     presid.tmp=lapply(presid, function(x) x[,keep.traits])
                 }
                  for(tt in keep.traits) {
                       found.q.covs = split(found.t.covs[[tt]], found.t.covs[[tt]]$cross)
                       for(cc in names(found.q.covs)) {
                            presid.tmp[[cc]][,tt]=scale(residuals(lm(presid.tmp[[cc]][,tt]~seg.chr.genos.scaled[[cc]][,found.q.covs[[cc]]$marker])))
                       }
                 } 
                jointLOD=mapply(function(n,p,g) {
                              fasterLOD(n,p,g)
                              }, n=seg.chr.n, p=presid.tmp, g=seg.chr.genos.scaled)
                jL=matrix(0,length(keep.traits), length(upos))
                rownames(jL)=colnames(presid.tmp[[1]])
                for(cc in crosses){ 
                    censor.me=jointLOD[[cc]]
                    #censor.me[censor.me<3]=0
                    jL[,mupos[[cc]]]=jL[,mupos[[cc]]]+censor.me} 
                if(nrow(jL)==1) {findingPeaks=FALSE}
                presid=presid.tmp
                jump=jump+1 
            } else {findingPeaks=FALSE } 
        }
    jointPeaks[[chrom]]=do.call('rbind', peakList)
    }
    return(list(jointPeaks=jointPeaks, jointLODs=jointLODs))
}



doLODdrop=function(trait, genos.full, f.found) {
    ys=trait
    gs=genos.full
    nsegs=length(ys)
    #print(f.found)
    registerDoMC(cores=length(f.found))
    located=c()
    if(length(f.found)>1){
        located=foreach(j=1:length(f.found), .combine='rbind') %dopar% { 
             # in 1:nrow(zf5)) { 
            nm=lm(ys~gs[,f.found[-j]]-1)
            nllik=logLik(nm)/(log(10))
            coi=strsplit(f.found[j], '_')[[1]][1]
            gcoi=gs[,grep(paste0(coi,'_'), colnames(gs))]
            mnames=colnames(gcoi)
            LOD=rep(0, ncol(gcoi))
            for(g in 1:ncol(gcoi)){
                #if(g%%100==0) {print(g)}
                LOD[g]=(logLik(lm(ys~gs[,f.found[-j]]+gcoi[,g]-1))/log(10))-nllik
            }
           return(data.frame(LOD=max(LOD), pmarker=mnames[which.max(LOD)],
                             CI.l=mnames[min(which(LOD>max(LOD)-1.5))],
                             CI.r=mnames[max(which(LOD>max(LOD)-1.5))], stringsAsFactors=F))
    } }

  if(length(f.found)==1){
            nm=lm(ys~1)
            nllik=logLik(nm)/(log(10))
            coi=strsplit(f.found, '_')[[1]][1]
            gcoi=gs[,grep(paste0(coi,'_'), colnames(gs))]
            mnames=colnames(gcoi)
            LOD=rep(0, ncol(gcoi))
            for(g in 1:ncol(gcoi)){
                #if(g%%100==0) {print(g)}
                LOD[g]=(logLik(lm(ys~gcoi[,g]-1))/log(10))-nllik
            }
            located=(data.frame(LOD=max(LOD), pmarker=mnames[which.max(LOD)],
                             CI.l=mnames[min(which(LOD>max(LOD)-1.5))],
                             CI.r=mnames[max(which(LOD>max(LOD)-1.5))], stringsAsFactors=F))
   } 
    return(located)
}


# forward stepwise procedure with FDR control
# G'sell 2013 procedure to detect QTL per trait
doTraitFDR=function(trait, genos, genos.full, FDR_thresh=.05, nperm=1e4, doLODdrop=T) {
    f.found=c()
    p.found=c()
    q.found=c()
    m.found=c()

    n=length(trait)
    L= (crossprod(trait,genos)/(n-1))^2 
    mLi=which.max(L)
    mL=max(L)
    
    yperm=replicate(nperm, sample(trait))
    nullD=(crossprod(yperm,genos)/(n-1))^2
    
    permMax=rowMaxs(nullD,value=T)
    pNull=1-ecdf(permMax)(mL)
    if(pNull==0) {pNull=1/nperm}
    
    step=1
    
    repeat{
       p.temp=c(p.found, pNull)
       q=-mean(log(1-p.temp))
       if(q>FDR_thresh) {break;}
       p.found=c(p.found, pNull)
       q.found=c(q.found, q)
       m.found=c(m.found, colnames(genos)[mLi])
       f.found=c(f.found, mLi)
       #print(paste('step=', step, 'max index=', colnames(genos)[mLi], 'max r^2=', mL, 'pnull=', pNull, 'fdr=', q))
       yr=scale(residuals(lm(trait~genos[,f.found]) ))
       L=(crossprod(yr,genos)/(n-1))^2 
       mLi=which.max(L)
       mL=max(L)
       yperm=replicate(nperm, sample(yr))
       nullD=(crossprod(yperm,genos)/(n-1))^2
       permMax=rowMaxs(nullD, value=T) 
       pNull=1-ecdf(permMax)(mL)
       if(pNull==0) {pNull=1/nperm}
       step=step+1
   }
   results=data.frame(fscan.markers=m.found, index=f.found, p=p.found, q=q.found, stringsAsFactors=F) 
   if(doLODdrop) {
       drops=doLODdrop(trait, genos.full, results$fscan.markers)
       results=cbind(results,drops)
   }
   return(results)
}

# crossValidation to estimate variance explained
# 10-fold cross validation
doTraitCV=function(mPhenos, g.s, gall.s) {
    cvg=cut(sample(1:nrow(mPhenos)),10)
    levels(cvg)=as.roman(1:10)
    cvVE=list()
    for(cv in levels(cvg)) {
         print(cv)
         cQTLcv=apply(mPhenos[cvg!=cv,], 2, function(x) {
                          set.seed(100); 
                          return(doTraitFDR(scale(x), scale(g.s[cvg!=cv,]), scale(gall.s[cvg!=cv,]), FDR_thresh=.05, nperm=5e2, doLODdrop=F)) 
                          })
         predicted.r2=rep(NA, ncol(mPhenos))
         names(predicted.r2)=colnames(mPhenos)
         for(j in 1:ncol(mPhenos)){
             yr=mPhenos[,j]
             yr.train=yr
             yr.train[cvg==cv]=NA
             yr.test=yr
             yr.test[cvg!=cv]=NA    
             X2=data.frame(g.s[,unique(cQTLcv[[j]]$index)])
             if(ncol(X2)==0) { next;}
             # fit model on training data
             fitme=(lm(yr.train~.-1,X2))
             if(is.null(dim(X2))){
                predicted=data.matrix(X2[cvg==cv])*coef(fitme)
             }else {
                 # estimate variance explained  
                predicted=data.matrix(X2[cvg==cv,])%*%coef(fitme)
            }
            predicted.r2[j]=cor(yr.test[cvg==cv], predicted)^2
         }
         cvVE[[cv]]=predicted.r2
    }
    return(do.call('cbind', cvVE))
}

#test=fs(g.s, mPhenos[,1], maxsteps=20, verbose=T, normalize=F, intercept=F)
#testF=fsInf(test, verbose=T)
#testP=forwardStop(testF, alpha=.05)

extractVarCompResultsR = function(r) {list(sigma=r$sigma, sigma.cov=r$sigma.cov, #W=r$W,
                                          Bhat=as.vector(r$Beta), llik=as.vector(r$llik)) }

#additive only variance component model using averaged phenotype per segregant
# includes standard errors
doAdditiveVC=function(mPhenos, g.s) {
    A=tcrossprod(g.s)/ncol(g.s)
    aVC=(apply(mPhenos, 2, function(y) {
              extractVarCompResultsR(regress(y~1,~A, verbose=T)) }))
    A.sigma=sapply(aVC, function(y) y$sigma)
    #cross.VCs, function(x) sapply(x, function(y) y$sigma[rev(c(2:6,1,7))]))
    A.sigma.se=sapply(aVC, function(y) sqrt(diag(y$sigma.cov)))
    nf=colSums(A.sigma)
    A.sigma=t(t(A.sigma)/nf)
    A.sigma.se=t(t(A.sigma.se)/nf)
    results=cbind(t(A.sigma), t(A.sigma.se))
    colnames(results)[3:4]=paste0( colnames(results)[3:4], '.SE')
    return(results)
}
#faster version but no standard errors
do_VC_additive_only_average=function(mPheno, g.s) {
    A=tcrossprod(g.s)/(ncol(g.s))
    eigA=doEigenA_forMM(mPheno,A)
    svdAloco=svd(A)
    VC_A=matrix(0,2, ncol(mPheno))
    colnames(VC_A)=colnames(mPheno)
    rownames(VC_A)=c('A', 'E')
    for(tp in 1:ncol(mPheno)){
        VC_A[,tp]=m.S(mPheno[,tp], K=A,  theta=eigA$theta, Q=eigA$Q)
    }
    VC_A=data.frame(t(VC_A)/colSums(VC_A))
    VC_A$trait=rownames(VC_A)
    return(VC_A)
}
# joint QTLs 5 (estimate joint, fit within, repeat)
# fixed effect per joint detected marker, FDR control for model selection
mapJointQTLs3=function(n.perm=1e4, FDR_thresh=.2, parents.list, pheno.resids, seg.recoded) {
    FDR_thresh=.05
    n.perm=1e3
    jointPeaks5=list()
    for( chrom in unique.chrs) {
        print(chrom)
        plist.chr=sapply(parents.list, function(x) x[x$chr==chrom,])
        chr.markers=sapply(parents.list, function(x){ return(x$marker.name.n[x$chr==chrom]) } )
        seg.chr.phenos=lapply(pheno.resids, function(x) { x[[chrom]] })
        seg.chr.n=lapply(seg.chr.phenos, nrow)
        seg.chr.phenos.scaled=lapply(seg.chr.phenos,scale)
       
        cmm=sapply(parents.list, function(x){ return(x$marker.name[x$chr==chrom]) } )

        seg.chr.genos=mapply(function(cmm,seg.r) { t(seg.r[cmm,])} , cmm=chr.markers, seg.r=seg.recoded)
        seg.chr.genos=mapply(function(cmm,seg.r) { colnames(seg.r)=cmm; return(seg.r);}, cmm=cmm, seg.r=seg.chr.genos)

        seg.chr.genos.scaled=lapply(seg.chr.genos, scale)
        #rename seg.chr.genos.scaled
        #seg.chr.genos.scaled=mapply(function(cmm,seg.r) { colnames(seg.r)=cmm; return(seg.r);}, cmm=cmm, seg.r=seg.chr.genos.scaled)

        mpmarkers=melt(cmm)
        mpmarkers[,1]=as.character(mpmarkers[,1])
        names(mpmarkers)=c('marker', 'cross')

       mpmarkers=separate(mpmarkers, marker, c('chr', 'pos', 'ref', 'alt'), sep='_', convert=T, remove=F)
       mpmarkers=mpmarkers[order(mpmarkers$pos),]

       ccnt=split(mpmarkers$cross, mpmarkers$marker)
       crosses.per.marker=sapply(ccnt, function(x) unique(x))
       cross.cnt.per.marker=sapply(crosses.per.marker, length)
       seg.rare=names(cross.cnt.per.marker[cross.cnt.per.marker<3])

       mm=mpmarkers[match(names(crosses.per.marker), mpmarkers$marker),-6]
       mm=mm[order(mm$pos),]

       # mapping cross marker index to joint marker table  
       mupos=lapply(plist.chr, function(x) match(x$marker.name, mm$marker))
       # mapping cross marker to joint marker table name 
       mucpos = mapply(function(m, p) { p$marker.name[ match(mm$marker[m], p$marker.name)] }, m=mupos, p=plist.chr )
       peakList=list()
       presid=seg.chr.phenos.scaled
       upos=mm$marker

      
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
        bm2=big.mat
        bm2[bm2==0]=NA

       traits=colnames(seg.chr.phenos.scaled[[1]])
       jL.list=list()

       nczero=apply(bmr.scaled,2, function(x) sum(x!=0))
       
       for(trait in 1:length(traits)) {
      
       print(traits[trait])
       f.found=c()
       p.found=c()
       q.found=c()
       m.found=c()
       p=do.call('c', lapply(presid, function(x) x[,trait]))
       crossF=as.factor(rep(names(presid), sapply(presid, nrow)))
       
       jL=-log10(fasterR2(nczero,p,bmr.scaled)[1,])
       #jL2=cor(p, bm2, use='pairwise.complete.obs')[1,]^2
       mL=max(jL) 
       mLi=which.max(jL)
       
       plot(jL, main=paste(chrom, traits[trait]))
       step=1
       abline(v=mLi, col=step)
       
       #structured permutations
       set.seed(5)
       yp=do.call('rbind', mapply(function(n,y) { replicate(n.perm, {y[sample(n),trait]}) }, n=seg.chr.n, y=presid))
       tN=t(matrix(rep(nczero, n.perm),length(nczero)))
       jLp=-log10(fasterR2(tN, yp,bmr.scaled))
       #jLp=cor(yp,big.mat, use='pairwise.complete.obs')^2
       jLpm=rowMaxs(jLp, value=T) 
       #jLp2=cor(yp[,1:10], bm2, use='pairwise.complete.obs')[1,]^2

       pNull=1-ecdf(jLpm)(mL[1])
       if(pNull==0) {pNull=1/n.perm}
        
        repeat{
           p.temp=c(p.found, pNull)
           q=-mean(log(1-p.temp))
           if(q>FDR_thresh) {break;}
           p.found=c(p.found, pNull)
           q.found=c(q.found, q)
           m.found=c(m.found, upos[mLi])
           f.found=c(f.found, mLi)
           print(paste('step=', step, 'max index=', upos[mLi],
                       'peakCrosses=', cross.cnt.per.marker[[names(mLi)]],
                        '1012 af=', iseq.freqs$maf1012[match(names(mLi), iseq.freqs$marker)],
                       'max r2=', mL, 'pnull=', pNull, 'fdr=', q))
           px=lapply(presid,function(x) x[,trait])
           presid2=mapply(function(p,g){
                kmarkers=upos[f.found] %in% colnames(g)
                if(sum(kmarkers>0)) {
                    return( as.vector(scale(residuals(lm(p~g[,upos[f.found][kmarkers]])))))
                 }      else {      return(scale(p))              }
            }, p=px, g=seg.chr.genos.scaled)
            pr2=do.call('c', presid2)

           jL=-log10(fasterR2(nczero,pr2,bmr.scaled)[1,])
            
            mL=max(jL) 
            mLi=which.max(jL)
            abline(v=mLi, col=step+1)
         
            yp=do.call('rbind', mapply(function(n,y) { replicate(n.perm, {y[sample(n)]}) }, n=seg.chr.n, y=presid2))
            jLp=-log10(fasterR2(tN,yp ,bmr.scaled))
            jLpm=rowMaxs(jLp, value=T) 
            pNull=1-ecdf(jLpm)(mL[1])
            if(pNull==0) {pNull=1/n.perm}
            step=step+1
        }
       results=data.frame(fscan.markers=m.found, index=f.found, p=p.found, q=q.found, stringsAsFactors=F) 
         jL.list[[traits[trait]]]=results                         
      }
       jointPeaks5[[chrom]]=rbindlist(jL.list, idcol='trait')
    }
 return(jointPeaks5)
}


# only JS strain variants 
# fixed effect per joint marker, FDR control for model selection
mapJointQTLsJS_variants=function(n.perm=1e3, FDR_thresh=.05, parents.list, pheno.resids, seg.recoded, iseq.freqs=NULL, filterJS=T) {
    jointPeaksJS=list()
    #jointLODs=list()
    # assume we're starting with residualized traits per chromosome 
    for( chrom in  names(pheno.resids[[1]])) {
        print(chrom)
        plist.chr=sapply(parents.list, function(x) x[x$chr==chrom,])
        chr.markers=sapply(parents.list, function(x){ return(x$marker.name.n[x$chr==chrom]) } )
        seg.chr.phenos=lapply(pheno.resids, function(x) { x[[chrom]] })
        seg.chr.n=lapply(seg.chr.phenos, nrow)
        seg.chr.phenos.scaled=lapply(seg.chr.phenos,scale)
       
        cmm=sapply(parents.list, function(x){ return(x$marker.name[x$chr==chrom]) } )

        seg.chr.genos=mapply(function(cmm,seg.r) { t(seg.r[cmm,])} , cmm=chr.markers, seg.r=seg.recoded)
        seg.chr.genos=mapply(function(cmm,seg.r) { colnames(seg.r)=cmm; return(seg.r);}, cmm=cmm, seg.r=seg.chr.genos)
        
        if(filterJS) {
        seg.chr.genos=lapply(seg.chr.genos, function(g) {
            return(g[,which(colnames(g) %in% iseq.freqs$marker)])
            })
        }
        seg.chr.genos.scaled=lapply(seg.chr.genos, scale)
        #rename seg.chr.genos.scaled
        #seg.chr.genos.scaled=mapply(function(cmm,seg.r) { colnames(seg.r)=cmm; return(seg.r);}, cmm=cmm, seg.r=seg.chr.genos.scaled)
        # now update cmm
        cmm=sapply(seg.chr.genos, function(x){ return(colnames(x)) }) 
        mpmarkers=melt(cmm)
        mpmarkers[,1]=as.character(mpmarkers[,1])
        names(mpmarkers)=c('marker', 'cross')

       mpmarkers=separate(mpmarkers, marker, c('chr', 'pos', 'ref', 'alt'), sep='_', convert=T, remove=F)
       mpmarkers=mpmarkers[order(mpmarkers$pos),]

       ccnt=split(mpmarkers$cross, mpmarkers$marker)
       crosses.per.marker=sapply(ccnt, function(x) unique(x))
       cross.cnt.per.marker=sapply(crosses.per.marker, length)
       seg.rare=names(cross.cnt.per.marker[cross.cnt.per.marker<3])

       mm=mpmarkers[match(names(crosses.per.marker), mpmarkers$marker),-6]
       mm=mm[order(mm$pos),]

       # mapping cross marker index to joint marker table  
       mupos=lapply(plist.chr, function(x) match(x$marker.name, mm$marker))
       # mapping cross marker to joint marker table name 
       mucpos = mapply(function(m, p) { p$marker.name[ match(mm$marker[m], p$marker.name)] }, m=mupos, p=plist.chr )
       peakList=list()
       presid=seg.chr.phenos.scaled
       upos=mm$marker

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
        bm2=big.mat
        bm2[bm2==0]=NA

       traits=colnames(seg.chr.phenos.scaled[[1]])
       jL.list=list()

       nczero=apply(bmr.scaled,2, function(x) sum(x!=0))
       
       for(trait in 1:length(traits)) {
      
       print(traits[trait])
       f.found=c()
       p.found=c()
       q.found=c()
       m.found=c()
       p=do.call('c', lapply(presid, function(x) x[,trait]))
       crossF=as.factor(rep(names(presid), sapply(presid, nrow)))
       
       jL=-log10(fasterR2(nczero,p,bmr.scaled)[1,])
       #jL2=cor(p, bm2, use='pairwise.complete.obs')[1,]^2
       mL=max(jL) 
       mLi=which.max(jL)
       
       #plot(jL, main=paste(chrom, traits[trait]))
       step=1
       #abline(v=mLi, col=step)
       
       #structured permutations
       set.seed(5)
       yp=do.call('rbind', mapply(function(n,y) { replicate(n.perm, {y[sample(n),trait]}) }, n=seg.chr.n, y=presid))
       tN=t(matrix(rep(nczero, n.perm),length(nczero)))
       jLp=-log10(fasterR2(tN, yp,bmr.scaled))
       #jLp=cor(yp,big.mat, use='pairwise.complete.obs')^2
       jLpm=rowMaxs(jLp, value=T) 
       #jLp2=cor(yp[,1:10], bm2, use='pairwise.complete.obs')[1,]^2

       pNull=1-ecdf(jLpm)(mL[1])
       if(pNull==0) {pNull=1/n.perm}
        
        repeat{
           p.temp=c(p.found, pNull)
           q=-mean(log(1-p.temp))
           if(q>FDR_thresh) {break;}
           p.found=c(p.found, pNull)
           q.found=c(q.found, q)
           m.found=c(m.found, upos[mLi])
           f.found=c(f.found, mLi)
           # live updates of progress
           print(paste('step=', step, 'max index=', upos[mLi],
                       'peakCrosses=', cross.cnt.per.marker[[names(mLi)]],
                        '1012 af=', iseq.freqs$maf1012[match(names(mLi), iseq.freqs$marker)],
                       'max r2=', mL, 'pnull=', pNull, 'fdr=', q))
           px=lapply(presid,function(x) x[,trait])
           presid2=mapply(function(p,g){
                kmarkers=upos[f.found] %in% colnames(g)
                if(sum(kmarkers>0)) {
                    return( as.vector(scale(residuals(lm(p~g[,upos[f.found][kmarkers]])))))
                 }      else {      return(scale(p))              }
            }, p=px, g=seg.chr.genos.scaled)
            pr2=do.call('c', presid2)

           jL=-log10(fasterR2(nczero,pr2,bmr.scaled)[1,])
            
            mL=max(jL) 
            mLi=which.max(jL)
            #abline(v=mLi, col=step+1)
         
            yp=do.call('rbind', mapply(function(n,y) { replicate(n.perm, {y[sample(n)]}) }, n=seg.chr.n, y=presid2))
            jLp=-log10(fasterR2(tN,yp ,bmr.scaled))
            jLpm=rowMaxs(jLp, value=T) 
            pNull=1-ecdf(jLpm)(mL[1])
            if(pNull==0) {pNull=1/n.perm}
            step=step+1
        }
       results=data.frame(fscan.markers=m.found, index=f.found, p=p.found, q=q.found, stringsAsFactors=F) 
         jL.list[[traits[trait]]]=results                         
      }
       jointPeaksJS[[chrom]]=rbindlist(jL.list, idcol='trait')
    }
    return(jointPeaksJS)
}


# for plotting allele frequency vs effect size 
get_density <- function(x, y, n = 100) {
    x1= x[which(!is.na(x) & !is.na(y))]
    y1= y[which(!is.na(x) & !is.na(y))]
  dens <- MASS::kde2d(x =x1 , y = y1, n = n)
  ix <- findInterval(x1, dens$x)
  iy <- findInterval(y1, dens$y)
  ii <- cbind(ix, iy)
  zo=dens$z[ii]
  z=rep(NA,length(x))
  z[which(!is.na(x) &!is.na(y))]=zo
  return(z)
}


# annotate JS variants as ancestral or not based on alignment of sacCer3 to paradoxus 
buildJS_variants_annotation_table=function(js.rr.overlap.allele.frequencies,sacCer3_CBS432_alignment.variants,sacCer3_CBS432_alignment.coords) {
    
    iseq.freqs=read.delim(js.rr.overlap.allele.frequencies, header=F, stringsAsFactors=F)
    iseq.freqs$maf1012=ifelse(iseq.freqs[,6]<iseq.freqs[,8], iseq.freqs[,6], iseq.freqs[,8])
    iseq.freqs$marker=paste(iseq.freqs[,1], iseq.freqs[,2], iseq.freqs[,5], iseq.freqs[,7], sep='_')
    ifm=sapply(strsplit(iseq.freqs$marker, '_'), function(x) paste0(x[1], '_' , x[2] , '_'))

    palign=read.delim(sacCer3_CBS432_alignment.variants,header=T, sep='\t', stringsAsFactors=F)
    palign$refC=sapply(strsplit(palign$refChr, ':'), function(x) x[2])
    rtext=gsub('chr', '', palign$refC)
    palign$chr=paste0('chr', as.roman(as.numeric(rtext)))
    palign$merged=paste0(palign$chr, '_', palign$P1, '_')

    pr=paste0(palign$chr, '_', palign$P1, '_', palign$REF)
    pa=paste0(palign$chr, '_', palign$P1, '_', palign$ALT)

    #reference_iseqfreqs
    ifmr=sapply(strsplit(iseq.freqs$marker, '_'), function(x) paste0(x[1], '_' , x[2] , '_', x[3]))
    #alt iseq_freqs
    ifma=sapply(strsplit(iseq.freqs$marker, '_'), function(x) paste0(x[1], '_' , x[2] , '_', x[4]))
    #sum(!is.na(match(ifmr, pr)))
    #sum(!is.na(match(ifma, pa)))
    #sum(!is.na(match(ifmr, pa)))
    #sum(!is.na(match(ifma, pr)))

    ifrmatches=cbind(match(ifmr, pr), match(ifma, pa),match(ifma, pr))

    palign.alignments=read.delim(sacCer3_CBS432_alignment.coords,header=F, sep='\t', stringsAsFactors=F)
    palign.alignments$refC=sapply(strsplit(palign.alignments[,12], ':'), function(x) x[2])
    rtext=gsub('chr', '', palign.alignments$refC)
    palign.alignments$chr=paste0('chr', as.roman(as.numeric(rtext)))
    palign.alignments=palign.alignments[palign.alignments$chr %in% unique.chrs,]
    mummer.alignned=GRanges(seqnames=palign.alignments$chr, ranges=IRanges(start=palign.alignments[,1], end=palign.alignments[,2]), strand='*')
    js.snps=GRanges(seqnames=iseq.freqs[,1], ranges=IRanges(start=iseq.freqs[,2], end=iseq.freqs[,2]), strand='*')
    js.snps.aligned=unique(queryHits(findOverlaps(js.snps, mummer.alignned)))
    ifrmatches=cbind(ifrmatches, rep(NA, nrow(ifrmatches)))
    ifrmatches[js.snps.aligned,4]=1
    
    rdap=which(!is.na(ifrmatches[,1]) & !is.na(ifrmatches[,2]))
    rpad=which((!(!is.na(ifrmatches[,1]) & !is.na(ifrmatches[,2])) | !is.na(ifrmatches[,3])) & ifrmatches[,4]==1  )
    # break into even sized (equal number of markers per bin) bins 
    iseq.freqs$bins=cut2(iseq.freqs$maf1012, g=7)
    iseq.freqs$bins2=iseq.freqs$maf1012>.01
    iseq.freqs$ancestral=NA
    iseq.freqs$ancestral[rdap]=0
    iseq.freqs$ancestral[rpad]=1
    return(iseq.freqs) 
}

# by default remove indices 37 and 38 corresponding to the two extra YPD replicate experiments
extractVC_long_format=function(joint.VCs, remove.ind=-c(37,38)) {
    JC.sigma=sapply(joint.VCs[-c(37,38)], function(x) x$sigma)
    JC.sigma.se=sapply(joint.VCs[-c(37,38)], function(x) sqrt(diag(x$sigma.cov)))
    JC.sigma.nf=colSums(JC.sigma)
    JC.sigma.norm=t(t(JC.sigma[rev(1:nrow(JC.sigma)),])/JC.sigma.nf)
    JC.sigma.se.norm=t(t(JC.sigma.se[rev(1:nrow(JC.sigma)),])/JC.sigma.nf)
    JC.pos =apply(JC.sigma.norm,2,function(x) rev(cumsum(rev(x))))
    JCV=melt( JC.sigma.norm) 
    colnames(JCV)=c('Component', 'Trait', 'fraction_of_variance')
    JCV.se=melt(JC.sigma.se.norm) 
    JCV$se=JCV.se[,3]
    JCV$ypos=melt(JC.pos)[,3]
    #levels(JCV$Component)=c('E', 'Not-private', 'Private')
    return(JCV)
}




