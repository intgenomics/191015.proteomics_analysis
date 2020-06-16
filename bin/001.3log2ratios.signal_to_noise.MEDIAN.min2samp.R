


####============================================================================================================
##  signed v unsigned GTEX7 beta10 ----------------------
#╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╔╦╦╗
options(stringsAsFactors=F);library(colorout);rm(list=ls());ls()#╚═╣║ ╚╣║¯\_(•_•)_/¯║╚╣╔╣╔╣║║║║╚╣
options(menu.graphics=F)  #═╦╩╬╝╚╬╬╚╗═╚╩╗═╩╠╬╝╔╚╗╬╚║╣══╣╦╬╬╗╠╔╗╔╣╣╝╝╣╠╔╠╚╔╔═╦╩╬╣╦╣╔╚╬╦╣╩╬╚╩╗╣╝╚╠╣
library('adds') #╔╦╣═╩╚╣║╔╔╣╦═║║╔╗║╔╚╔╣╩╚╚╦╣║╩╔╦║║ ╚╩╣╚╚╣║╣╚╩╔╦╩╚╦╚╩╣╬╝╚╗╔╝╬╚╝ ╔╣═╦╦╦╩╠╔╠╗╔╝╚═╗╩║
# devtools::install_github("ks4471/addR") #╗╣╠═╩╠╣╠╬═╚╬╗╩╩═║╚╝╝╣╠╗╗╠╔║╩╬╠╝╣╬╔╬╬╚╦╝╔╗╩╠╚╝╠═╝╝╦╔═╚╠╝╣║
#╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩╩═╝

pdat=read.delim('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/dtb/Identified_proteinGroups_for_Kirill.txt')
# zdat=read.delim('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/dtb/Identified_proteinGroups_for_Kirill.txt')
  str(pdat)
  Table(pdat$label=='')
# 1 FALSE  3314    0.986
# 2  TRUE    46    0.014

# colnames(pdat)=gsub('Lam.vs.DCA','DCA.vs.Lam',colnames(pdat))


##  check mapping - remove duplicated or empty gene names
plog=unique(pdat[,c('Protein.IDs','label',colnames(pdat)[grepl('Log2[.]',colnames(pdat))])])
##  unassigned and missing 'labels'
  get.duplicates(plog,'Protein.IDs')
  get.duplicates(plog,'label')
  str(plog[plog$label=='',])

# pdf('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/img/missing.log2ratios.pdf',height=12,width=12)
#   is.missing(plog,T)
# dev.off()
#   ------------------   of total samples 3360 ,  1086  ( 32 %) remaining if "complete.cases()"  ------------------
#                                                      count percent
# Log2.Ratio.H.L.normalized.Lam.vs.DCA_Rep1             1952   0.581
# Log2.Ratio.H.L.normalized.Central.vs.Peripheral_Rep1  1753   0.522
# Log2.Ratio.H.L.normalized.Lam.vs.DCA_Rep2             1057   0.315
# Log2.Ratio.H.L.normalized.SNA.vs.Sham_Rep2            1049   0.312
# Log2.Ratio.H.L.normalized.SNA.vs.Sham_Rep3            1017   0.303
# Log2.Ratio.H.L.normalized.Central.vs.Peripheral_Rep3  1006   0.299
# Log2.Ratio.H.L.normalized.Lam.vs.DCA_Rep3              985   0.293
# Log2.Ratio.H.L.normalized.Central.vs.Peripheral_Rep2   936   0.279
# Log2.Ratio.H.L.normalized.SNA.vs.Sham_Rep1             935   0.278
# label                                                    0   0.000



##  collate unique names for all detected genes (any experiment / sample)
ugen=unique(plog$label)
  str(ugen)
ugen=ugen[ugen!='']
  str(ugen)
bkg=ugen
  overlap(ugen,plog$label)






##  for qc purposes - collage all genes that have info in at least one replicate for a given experiment
condi=c(
  'SNA.vs.Sham'
  ,'DCA.vs.Lam'
  ,'Central.vs.Peripheral'
)
gnams=list()
for(icon in condi){
# is.missing(plog[apply(plog[,grepl(icon,colnames(plog))],1,function(x){sum(!is.na(x))==0}),grepl(icon,colnames(plog))])      ##  confirm all missing selected only..
# is.missing(plog[apply(plog[,grepl(icon,colnames(plog))],1,function(x){sum(!is.na(x))==0}),grepl(icon,colnames(plog))])
Table(apply(plog[,grepl(icon,colnames(plog))],1,function(x){sum(!is.na(x))>0}))
#
  gnams[[icon]]=unique(plog$label[apply(plog[,grepl(icon,colnames(plog))],1,function(x){sum(!is.na(x))>0})])
  gnams[[icon]]=gnams[[icon]][gnams[[icon]]!='']
}

  dummy=overlap(ugen,unlist(unique(gnams)))

  # is.missing(plog[plog$label%in%dummy$ina,],T)
  ##  confirm all missing except id and Label.





##  merge genes with multiple datapoints (ie 2 peptides mapped to same gene)
mranks=matrix('',nrow=length(ugen),ncol=ncol(plog)-2)
  colnames(mranks)=colnames(plog)[grepl('Log2[.]',colnames(plog))]
  rownames(mranks)=ugen

# igen='Ahnak2'
for(igen in ugen){
  mranks[igen,]=apply(plog[plog$label==igen,grepl('Log2[.]',colnames(plog))],2,median,na.rm=T)

}

  mranks=make.numeric(mranks)
  overlap(rownames(mranks),ugen)
  overlap(rownames(mranks),unlist(unique(gnams)))




##  check for proteins with no expression in any replicate / experiment
  # rownames(mranks)[(apply(mranks,1,function(x){sum(is.na(x))==length(x)}))]
# mranks[c('Nasp','Fbxo7','Dr1'),]
# is.missing(mranks[(apply(mranks,1,function(x){sum(is.na(x))==length(x)})),])


# plog[plog$label=='Sptbn5',]
# mranks['Sptbn5',]



# write.file(rownames(mranks[(apply(mranks,1,function(x){sum(is.na(x))==length(x)})),]),file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/txt/genes.all_missing.txt',row.names=F,col.names=F)



# pdf('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/img/missing.log2ratios.unique_genes.pdf',height=12,width=12)
#   is.missing(holder,T)
# dev.off()



prnk=list()
norn=list()
mrank=list()
igen='1700066M21Rik'
##  merge log2 ranks across 3 replicates into one
# icon=condi[2]
for(icon in condi){
   k=lprogr(icon,condi)

##  rm data where only one sample is available for the experiment
  holder=mranks[,grepl(icon,colnames(mranks))]
    # is.missing(holder,T)
  holder=holder[apply(holder,1,function(x){sum(!is.na(x))>1}),]
    # is.missing(holder,T)
##  median ranks across 3 replicates
  prnk[[icon]]=as.data.frame(apply(holder,1,median,na.rm=T))
    colnames(prnk[[icon]])=icon


  pdf(paste0('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/img/pretty/hist_norm_MEDIAN.log2ratios.',icon,'.pdf'),height=9,width=12)
      hist.norm(prnk[[icon]][!is.na(prnk[[icon]])],main=icon,breaks=300)
      abline(v=0,col='darkblue',lty=3,lwd=3)
    if(k==1){abline(v=prnk[[icon]]['Prkaa2',],col='darkred')}
      hist.dens(prnk[[icon]][!is.na(prnk[[icon]])],main=icon)
      abline(v=0,col='darkblue',lty=3,lwd=3)
    if(k==1){abline(v=prnk[[icon]]['Prkaa2',],col='cornflowerblue')}

##  deconvolution on median ranks - separate out the 3 distributions (manual cutoff currently)
##   the below code is located here for plotting inside .pdf
##   - alternatively - parameter scan to maximise the normality of the central distribution using shapiro.test() - may be sub-optimal for numerous possibilities, unless just scan across a set of pre-defined cutoffs based on manual observations..
  if(icon=='SNA.vs.Sham'){norn[[icon]]=deconv(make.numeric(as.data.frame(prnk[[icon]][!is.na(prnk[[icon]])])),T,mth=c(0.2,0.8),breaks=100,ylim_dens=c(0,2.1),ylim=c(0,2.1))}
  if(icon=='DCA.vs.Lam'){norn[[icon]]=deconv(make.numeric(as.data.frame(prnk[[icon]][!is.na(prnk[[icon]])])),T,mth=c(0.015,0.995),breaks=100,ylim_dens=c(0,2.7),ylim=c(0,2.7))}
  if(icon=='Central.vs.Peripheral'){norn[[icon]]=deconv(make.numeric(as.data.frame(prnk[[icon]][!is.na(prnk[[icon]])])),T,mth=c(0.03,0.94),breaks=100,ylim_dens=c(0,0.8),ylim=c(0,0.9))}

# deconv(make.numeric(as.data.frame(prnk[[icon]][!is.na(prnk[[icon]])])),T,mth=c(0.03,0.94),breaks=100,ylim_dens=c(0,0.8))
  
  cat('\t\t',sum(is.na(prnk[[icon]])),'\n')

  dev.off()

  if(k==1){   mrank=prnk[[icon]]}
  if(k>1){    mrank=rmerge(mrank,prnk[[icon]],all=T)}
  
}

##  trying to make histograms that overlap by frequency rather than density, neither work well
# exp_mat=make.numeric(as.data.frame(prnk[[icon]][!is.na(prnk[[icon]])]))
# mth=c(0.2,0.8)

# pdat_hist=hist(exp_mat,col=rgb(0,0,0,alpha=0.5),prob=F,breaks=100)
# pdat_hist=hist(upper,col=rgb(0,1,0,alpha=0.5),prob=F,breaks=50,add=T)
# pdat_hist=hist(lower,col=rgb(0,0,1,alpha=0.5),prob=F,breaks=100,add=T)
# lines(density(upper, na.rm=T), lty=1, lwd=2,col='darkred')
# lines(density(lower, na.rm=T), lty=1, lwd=2,col='dodgerblue')


  Head(mrank)
  overlap(rownames(mrank),ugen)
  overlap(rownames(mrank),unlist(unique(gnams)))
  # is.missing(mrank,T)

##  calculate pvalue statistics - how far each gene ratio is away from the null and 'noise-only' distribution


# igen=rownames(mrank)[apply(mrank,1,function(x){sum(is.na(x))})==3]#[7]
# mrank[igen,,drop=F]
# median(mranks[igen,,drop=F])

# mrank[igen,,drop=F]
# median(mranks[igen,,drop=F])


genp=list()
geng=list()
igen='Pex5l'
geno=list()

for(icon in names(prnk)){
   k=lprogr(icon,names(prnk))
  prnk[[icon]]=as.data.frame(prnk[[icon]])
  genp[[icon]]$mlog2r=prnk[[icon]][!is.na(prnk[[icon]]),,drop=F]
    colnames(genp[[icon]]$mlog2r)=paste0(icon,'.log2(ratio)')

  # genp[[icon]]$mlog2rs=scale(genp[[icon]]$mlog2r, center = TRUE, scale = TRUE)
  #     colnames(genp[[icon]]$mlog2rs)=paste0(icon,'.log2(ratio).normalised')

holder=list()
tolder=list()
    for(igen in rownames(genp[[icon]]$mlog2r)){

      if(genp[[icon]]$mlog2r[igen,]<0){
        holder[[igen]]=pnorm(genp[[icon]]$mlog2r[igen,],mean = mean(genp[[icon]]$mlog2r[,1]), sd = sd(genp[[icon]]$mlog2r[,1]),lower.tail=T)
        tolder[[igen]]=pnorm(genp[[icon]]$mlog2r[igen,],mean = mean(norn[[icon]]$middle[,1]), sd = sd(norn[[icon]]$middle[,1]),lower.tail=T)
      }

      if(genp[[icon]]$mlog2r[igen,]>0){
        holder[[igen]]=pnorm(genp[[icon]]$mlog2r[igen,],mean = mean(genp[[icon]]$mlog2r[,1]), sd = sd(genp[[icon]]$mlog2r[,1]),lower.tail=F)
        tolder[[igen]]=pnorm(genp[[icon]]$mlog2r[igen,],mean = mean(norn[[icon]]$middle[,1]), sd = sd(norn[[icon]]$middle[,1]),lower.tail=F)
      }
    }


  # genp[[icon]]$pvalo=t(as.data.frame(holder))
  genp[[icon]]$pvaln=t(as.data.frame(tolder))
    # colnames(genp[[icon]]$pvalo)=paste0(icon,'.pvalo')    ##  original  -- full distribution
    colnames(genp[[icon]]$pvaln)=paste0(icon,'.pvaln')    ##  new       -- relative to the 'mid' distribution

  # genp[[icon]]$fdro=as.data.frame(p.adjust(genp[[icon]]$pvalo,method='fdr'))
        # colnames(genp[[icon]]$fdro)=paste0(icon,'.fdro')
        # rownames(genp[[icon]]$fdro)=rownames(genp[[icon]]$pvalo)

  genp[[icon]]$fdrn=as.data.frame(p.adjust(genp[[icon]]$pvaln,method='fdr'))
        colnames(genp[[icon]]$fdrn)=paste0(icon,'.fdrn')
        rownames(genp[[icon]]$fdrn)=rownames(genp[[icon]]$pvaln)

  if(k==1){geno=as.data.frame(genp[[icon]])}
  if(k>1){geno=rmerge(geno,genp[[icon]],all=T)}
  

  geng[[icon]]$bkg=rownames(genp[[icon]]$fdrn)
  geng[[icon]]$sigen=rownames(genp[[icon]]$fdrn)[(genp[[icon]]$fdrn)<0.05]

    Table(genp[[icon]]$pvaln<0.05)
}




  str(genp)
  str(geng)
  str(geno)

  overlap(rownames(geno),ugen)
  overlap(rownames(mranks),ugen)
  str(union(rownames(genp$SNA.vs.Sham$mlog2r),union(rownames(genp$Lam.vs.DCA$mlog2r),rownames(genp$Central.vs.Peripheral$mlog2r))))

# plog[plog$label%in%c('Fap',  'Hddc2',  'Hk2',  'Itgb2',  'Scpep1'),]
# mranks[c('Fap',  'Hddc2',  'Hk2',  'Itgb2',  'Scpep1'),]
# mrank[c('Fap',  'Hddc2',  'Hk2',  'Itgb2',  'Scpep1'),]


#   Head(genp[[icon]]$pvalo)
#     Table(genp[[icon]]$pvalo<0.05)
#     Table(genp[[icon]]$fdro<0.05)
#   Head(genp[[icon]]$pvaln)
#     Table(genp[[icon]]$pvaln<0.05)
#     Table(genp[[icon]]$fdrn<0.05)
#     Table(genp[[icon]]$fdrn<0.05)
#   Head(genp[[icon]]$mlog2r)



  overlap(rownames(geno),gnams$sep)
  overlap(ugen,unlist(unique(gnams)))



##  check using log2 info - at least one data point present to validate the new ranks & pvals

str(geno)

for(icon in condi){
  lprogr(icon,condi,T)

  holder=geno[,grepl(icon,colnames(geno))]
  holder=holder[complete.cases(holder),]
  overlap(gnams[[icon]],rownames(holder))

}

plog[plog$label%in%c('Acp1',  'Anxa6',  'Crmp1',  'Hnrnpa3',  'Mocs2' ),(colnames(plog)=='label' | grepl(icon,colnames(plog)))]
geno[c('Acp1',  'Anxa6',  'Crmp1',  'Hnrnpa3',  'Mocs2' ),]



plog[plog$label=='Sptbn5',]
mranks['Sptbn5',]
geno['Sptbn5',]     ##  check - is excluded (single value)
geno['Sept3',]




 



##  numbers significantly detected 
  Table(geno$DCA.vs.Lam.fdrn<0.05)
  Table(geno$SNA.vs.Sham.fdrn<0.05)
  Table(geno$SNA.vs.Sham.fdrn<0.05 & geno$SNA.vs.Sham.log2.ratio. <0)
  Table(geno$Central.vs.Peripheral.fdrn<0.05)

# dummy=write.file(geno,file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/txt/deconvolution_tests_3_cond.MEDIAN.min2samp.txt')

# save(geno,condi,mranks,file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/dtb/deconvolution.results.MEDIAN.Rdata')


####===================================================================================================
##   log2 plots   -------------------------------------------


# icon=condi[3]
for(icon in condi){
   lprogr(icon,condi,T)
  holder=mranks[,grepl(icon,colnames(mranks))]
  dummy=geno[,grepl('fdr|log2[.]ratio',colnames(geno))&grepl(icon,colnames(geno)),drop=F]
  cutoff_log2r=min(abs(dummy[dummy[grepl('fdr',colnames(dummy))]<0.05,grepl('log2[.]ratio',colnames(dummy)) ]),na.rm=T)
  holder=rmerge(holder, dummy )

  holder=holder[apply(holder,1,function(x){sum(is.na(x))<3}),]
    str(holder)

  minmaxx=c(floor(min(holder[,!grepl('fdr|log2[.]ratio',colnames(holder))],na.rm=T)),ceiling(max(holder[,!grepl('fdr|log2[.]ratio',colnames(holder))],na.rm=T)))
  minmaxy=c(0,ceiling(-log10(min(
    holder[,grepl('fdr',colnames(holder))]
    ,na.rm=T))))

  min(abs(holder[,4]),na.rm=T)   ## median log2 ratio is 4th column "by construction" (sort of)

  pdf(paste0('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/img/log10_pval.',icon,'.pdf'),height=12,width=7)
    for(irep in 1:3){
      pdat=holder[,c(irep,ncol(holder))]
      pdat=pdat[complete.cases(pdat),]

      psig=pdat[pdat[,2]<=0.05 & abs(pdat[,1])>=cutoff_log2r,]
        summary(psig)
      potr=pdat[!rownames(pdat)%in%rownames(psig),]
        summary(potr)


      if(irep==1){
        plot(psig[,1], -log10(psig[,2]), pch = 16, col=rgb(0.4,0,0,alpha=0.8), frame.plot = F, xlim=minmaxx, ylim=minmaxy,xlab='log2 Ratio Heavy/Light',ylab='-log10 FDR',main=paste0(icon))
        points(potr[,1], -log10(potr[,2]), pch = 16, col=rgb(0,0,0,alpha=0.5))

      }
      if(irep!=1){
        points(psig[,1], -log10(psig[,2]), pch = 16, col=rgb(0.4,0,0,alpha=0.8))
        points(potr[,1], -log10(potr[,2]), pch = 16, col=rgb(0,0,0,alpha=0.5))
      }
    }
  dev.off()

}




































































































