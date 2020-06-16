

####============================================================================================================
##   ----------------------
#╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╔╦╦╗
options(stringsAsFactors=F);library(colorout);rm(list=ls());ls()#╚═╣║ ╚╣║¯\_(•_•)_/¯║╚╣╔╣╔╣║║║║╚╣
options(menu.graphics=F)  #═╦╩╬╝╚╬╬╚╗═╚╩╗═╩╠╬╝╔╚╗╬╚║╣══╣╦╬╬╗╠╔╗╔╣╣╝╝╣╠╔╠╚╔╔═╦╩╬╣╦╣╔╚╬╦╣╩╬╚╩╗╣╝╚╠╣
library('adds') #╔╦╣═╩╚╣║╔╔╣╦═║║╔╗║╔╚╔╣╩╚╚╦╣║╩╔╦║║ ╚╩╣╚╚╣║╣╚╩╔╦╩╚╦╚╩╣╬╝╚╗╔╝╬╚╝ ╔╣═╦╦╦╩╠╔╠╗╔╝╚═╗╩║
# devtools::install_github("ks4471/addR") #╗╣╠═╩╠╣╠╬═╚╬╗╩╩═║╚╝╝╣╠╗╗╠╔║╩╬╠╝╣╬╔╬╬╚╦╝╔╗╩╠╚╝╠═╝╝╦╔═╚╠╝╣║
#╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩╩═╝


Load('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/dtb/proteo.targetd_pulldown.sna.sham_OFF.Rdata')

proteo=holder
holder=lm.mat(proteo)

pdf('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/img/plot.pulldown_sna_sham.rep1_vs_rep2_OFF.pdf')
  Plot(proteo$sna.log2Ratio.fwd,proteo$sna.log2Ratio.rev,line45deg=T)
  Plot(proteo$sham.log2Ratio.fwd,proteo$sham.log2Ratio.rev,line45deg=T)
dev.off()


	Head(proteo)
proteo=proteo[apply(proteo,1,function(x){sum(is.na(x))<4}),]
	Head(proteo)

proteo$sna.mean=apply(proteo[,1:2],1,mean,na.rm=T)
proteo$sham.mean=apply(proteo[,3:4],1,mean,na.rm=T)
	# is.missing(holder,T)
##  keep 1 only or if expressed in forward and reverse ?


pdat=make.numeric(proteo)
pdat[is.na(pdat)]=0
# pdf('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/img/heat.mean_replicates.log2_ratios.pdf',height=45,width=9)
#   Heat(pdat,values='cor')
# dev.off()

pdf('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/img/heat.mean_only.log2_ratios.pdf',height=53,width=9)
  Heat(pdat[,grepl('mean',colnames(pdat))],values='cor')
dev.off()





##  sig detected deconvolution approach
mranks=proteo[,grepl('mean',colnames(proteo))]
	str(mranks)

##Rpl10 - ribosomal complex proteins - not expected to change after treatment - both in the "UP" tail
contr=rownames(proteo)[rownames(proteo)%in%c('Cpsf3','Rpl10','Rpl36a')]
contr=unique(c(contr,rownames(proteo)[grepl('Rpl',rownames(proteo))]))
  mranks[contr,]


##  heatmap of log2 ratio to illustrate 'differential' levels
##  venn diagram with genes uniquely
##    -  significanly detected in one experiment and not another
##    -  levels of differential activity - one but not the other - high expression && changing between two states ..
##    -  differences in log2 ratios... ?!  <<<<<**************


condat=mranks[,1,drop=F]
condat=make.numeric(as.data.frame(condat[!is.na(condat),1,drop=F]))
	Head(condat)
norn=deconv(condat,do_plots=F,mth=c(0.3,0.7),breaks=50)
  # dev.off()

# holder=list()
tolder=list()
    for(igen in rownames(condat)){

      if(condat[igen,]<0){
        # holder[[igen]]=pnorm(genp[[icon]]$mlog2r[igen,],mean = mean(genp[[icon]]$mlog2r[,1]), sd = sd(genp[[icon]]$mlog2r[,1]),lower.tail=T)
        tolder[[igen]]=pnorm(condat[igen,],mean = mean(norn$middle[,1]), sd = sd(norn$middle[,1]),lower.tail=T)
      }

      if(condat[igen,]>0){
        # holder[[igen]]=pnorm(genp[[icon]]$mlog2r[igen,],mean = mean(genp[[icon]]$mlog2r[,1]), sd = sd(genp[[icon]]$mlog2r[,1]),lower.tail=F)
        tolder[[igen]]=pnorm(condat[igen,],mean = mean(norn$middle[,1]), sd = sd(norn$middle[,1]),lower.tail=F)
      }
}
tolder=as.data.frame(unlist(tolder))
	colnames(tolder)='sna.pval'
tolder$sna.fdr=p.adjust(tolder$sna.pval,method='fdr')
proteo=rmerge(tolder,proteo)




rm(horn,tolder)
condat=mranks[,2,drop=F]
condat=make.numeric(as.data.frame(condat[!is.na(condat),1,drop=F]))
	Head(condat)
norn=deconv(condat,T,mth=c(0.3,0.7),breaks=50,ylim=c(0,0.4))
  dev.off()


##  detection above noise pvalues
# holder=list()
tolder=list()
    for(igen in rownames(condat)){

      if(condat[igen,]<0){
        # holder[[igen]]=pnorm(genp[[icon]]$mlog2r[igen,],mean = mean(genp[[icon]]$mlog2r[,1]), sd = sd(genp[[icon]]$mlog2r[,1]),lower.tail=T)
        tolder[[igen]]=pnorm(condat[igen,],mean = mean(norn$middle[,1]), sd = sd(norn$middle[,1]),lower.tail=T)
      }

      if(condat[igen,]>0){
        # holder[[igen]]=pnorm(genp[[icon]]$mlog2r[igen,],mean = mean(genp[[icon]]$mlog2r[,1]), sd = sd(genp[[icon]]$mlog2r[,1]),lower.tail=F)
        tolder[[igen]]=pnorm(condat[igen,],mean = mean(norn$middle[,1]), sd = sd(norn$middle[,1]),lower.tail=F)
      }
}
tolder=as.data.frame(unlist(tolder))
	colnames(tolder)='sham.pval'
tolder$sham.fdr=p.adjust(tolder$sham.pval,method='fdr')
proteo=rmerge(tolder,proteo)

# write.file(proteo,file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/txt/sna.sham.pulldown.csv',sep=',')

##  due to the nature of the puldown vs own IgG
##   + only significantly "UP-regulated"
##   + care only about proteins that are UP in one but not in other
##  significanly detected in one but not the other condition



glis1=proteo[proteo$sna.fdr<0.05 & proteo$sna.mean>0,c('sna.fdr','sna.mean')]
glis1=glis1[complete.cases(glis1),]
  # overlap(rownames(glis1),rownames(humpty)[(humpty[,'sna.fdr']<0.05 & humpty[,'sna.mean']>0)])

glis2=proteo[proteo$sham.fdr<0.05 & proteo$sham.mean>0,c('sham.fdr','sham.mean')]
glis2=glis2[complete.cases(glis2),]




geno=overlap(rownames(glis1),rownames(glis2))
# Venn(geno)

write.file(proteo[geno$ina,],file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/txt/sna_only_sig_up.pulldown.csv',sep=',')
write.file(proteo[geno$inb,],file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/txt/sham_only_sig_up.pulldown.csv',sep=',')
write.file(proteo[c(geno$ina,geno$inb),],file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/txt/sham_sna.sig_one_only_heatm.pulldown.csv',sep=',')
write.file(proteo[c(geno$inter),],file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/txt/sham_sna.sig_both_heatm.pulldown.csv',sep=',')



save(proteo,file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/dtb/proteo.targetd_pulldown.sna_sham.proteo.Rdata')
save(proteo,file='~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/dtb/proteo.targetd_pulldown.sna_sham.proteo_OFF.Rdata')






####============================================================================================================
##   ----------------------
options(stringsAsFactors=F);library(colorout);rm(list=ls());ls()
options(menu.graphics=F)
library('adds')
# devtools::install_github("ks4471/addR")

Load('~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/dtb/proteo.targetd_pulldown.sna_sham.proteo_OFF.Rdata')
str(proteo)

Library('metap')
# allmetap      Carry out all or some of the methods
# beckerp       Example data
# invchisq      Combine p values using inverse chi squared method
# invt          Combine p values using inverse t method
# logitp        Combine p values using logit method
# meanp         Combine p values by the mean p method
# meanz         Combine p values using mean z method
# metap         package Meta-Analysis of Significance Values
# schweder      Schweder and Spjotvoll plot
# sumlog        Combine p-values by the sum of logs (Fisher's) method
# sump          Combine p-values using the sum of p (Edgington's) method
# sumz          Combine p-values using the sum of z (Stouffer's) method
# two2one       Convert two-sided p-values to one-sided
# votep         Combine p-values by the vote counting method
# wilkinsonp    Combine p-values using Wilkinson's method


proteo$FisherP=NA
for(igen in rownames(proteo)){
  pvals=proteo[igen,grepl('pval',colnames(proteo))]
  if(sum(is.na(pvals))==0){
    proteo[igen,'FisherP']=sumlog(pvals)$p
  }
}

##  a few p-values are missing - likely cos only 1 observation ?
fishp=proteo[!is.na(proteo$FisherP),'FisherP',drop=F]
fishp$FisherFDR=p.adjust(fishp$FisherP,method='fdr')

proteo=proteo[,colnames(proteo)!='FisherP']
proteo=rmerge(proteo,fishp)

proteo$sna.log2Ratio.rev=proteo$sna.log2Ratio.rev * -1
proteo$sham.log2Ratio.rev=proteo$sham.log2Ratio.rev * -1

str(proteo)


Table(proteo$sham.fdr<0.05)
Table(proteo$sna.fdr<0.05)


# write.file(proteo,'~/Dropbox/PROJ/icl/giovanni.llaria/191015.proteomics_analysis/out/txt/sna.sham.pulldown.FisherP_OFF.txt')















