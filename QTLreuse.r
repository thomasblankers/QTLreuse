



### konpar ####

# konpar has been coded the wrong way, L. paranigra alleles are B, L. kona alles are A

konpar<-read.cross(format="csv", file="konpar_Comprehensive_CrossFile.csv", estimate.map=FALSE)
konpar_jitter<-jittermap(konpar)
konpar_jitter<-sim.geno(konpar_jitter, step=1, n.draws=1000, error.prob=0.001)
konpar_jitter<-calc.genoprob(konpar_jitter, step=1, error.prob=0.001)

## scanone

qtl_konpar_imp<-scanone(konpar_jitter,method="imp", pheno.col="pr")
qtl_konpar_hk<-scanone(konpar_jitter,method="hk", pheno.col="pr")
qtl_konpar_ehk<-scanone(konpar_jitter,method="ehk", pheno.col="pr")
qtl_konpar_em<-scanone(konpar_jitter,method="em", pheno.col="pr")

qtl_konpar_imp_1000perm<-scanone(konpar_jitter,method="imp", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE) # this takes about 2 minutes
qtl_konpar_hk_1000perm<-scanone(konpar_jitter,method="hk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
qtl_konpar_em_1000perm<-scanone(konpar_jitter,method="em", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
#qtl_konpar_ehk_1000perm<-scanone(konpar_jitter,method="ehk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)


## cim

#cim_konpar_imp<-cim(konpar_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=10)
cim_konpar_hk<-cim(konpar_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20)
#cim_konpar_ehk<-cim(konpar_jitter, pheno.col="pr", method="ehk", n.marcovar=3, window=10)
#cim_konpar_em<-cim(konpar_jitter, pheno.col="pr", method="em", n.marcovar=3, window=10)

cim_konpar_hk_1000perm<-cim(konpar_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20, n.perm=1000)

cairo_pdf("cim_konpar_hk_5_20.pdf", width=8, height=3.5)
plot(cim_konpar_hk, ylim=c(0,90));abline(h=summary(cim_konpar_hk_1000perm)[1,1])
dev.off()

qtl_konpar<-c("1","2","3","4","5","6","X")

## scantwo
twoqtl_konpar_imp<-scantwo(konpar_jitter, pheno.col="pr", method="imp", chr=qtl_konpar, clean.output=TRUE)
twoqtl_konpar_hk<-scantwo(konpar_jitter, pheno.col="pr", method="hk", chr=qtl_konpar, clean.output=TRUE)

plot(twoqtl_konpar_imp, lower="cond-int", upper="cond-add", c("3","5"))
plot(twoqtl_konpar_hk, lower="cond-int")

plot(twoqtl_konpar_imp, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_konpar_imp, chr=2, lower="cond-int", upper="cond-add")
# potential second qtl on chr 2
plot(twoqtl_konpar_imp, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_konpar_imp, chr=4, lower="cond-int", upper="cond-add")
plot(twoqtl_konpar_imp, chr=5, lower="cond-int", upper="cond-add")
# potential second qtl on chr 5
plot(twoqtl_konpar_imp, chr=6, lower="cond-int", upper="cond-add")
plot(twoqtl_konpar_imp, chr="X", lower="cond-int", upper="cond-add")

# permutate the scantwo function for LOD  penalties
twoqtl_konpar_hk_1000perm<-scantwo(konpar_jitter, pheno.col="pr", method="hk", chr=qtl_konpar, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 1,5 hours
summary(twoqtl_konpar_hk, threshold=as.matrix(as.data.frame(summary(twoqtl_konpar_hk_10000perm, 0.05)[1:5])[1,]))
penalties=calc.penalties(twoqtl_konpar_hk_1000perm, alpha=0.05)

    main    heavy    light 
3.464517 5.361266 3.197938 

max(qtl_konpar_hk)
               chr  pos  lod
S004383_170239   5 32.5 40.3

konpar_init_qtl<-makeqtl(konpar_jitter,chr=5, pos=32.5)
summary(fitqtl(konpar_jitter, qtl=konpar_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1))
konpar_init_qtl_add <- addqtl(konpar_jitter, qtl=konpar_init_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1)

konpar_expand_qtl<-addtoqtl(konpar_jitter,konpar_init_qtl, 4, 48)
summary(fitqtl(konpar_jitter, qtl=konpar_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
konpar_expand_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand_qtl)
summary(fitqtl(konpar_jitter, qtl=konpar_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
konpar_expand_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand_qtl, method="imp", formula=  y ~ Q1 + Q2)
konpar_expand_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2)


konpar_expand2_qtl<-addtoqtl(konpar_jitter,konpar_expand_qtl, 2, 57.6)
summary(fitqtl(konpar_jitter, qtl=konpar_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
konpar_expand2_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand2_qtl)
summary(fitqtl(konpar_jitter, qtl=konpar_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
konpar_expand2_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand2_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3) # Q1:Q3 and Q2:Q3 are both potentially significant, forward and backward selection shows that adding both violates light penalty, Q1:Q3 increases the LOD more so keep that one
konpar_expand2_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand2_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q1:Q3)

konpar_expand3_qtl<-addtoqtl(konpar_jitter,konpar_expand2_qtl, 1, 79.1)
summary(fitqtl(konpar_jitter, qtl=konpar_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q1:Q3))
konpar_expand3_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand3_qtl)
summary(fitqtl(konpar_jitter, qtl=konpar_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 +  Q1:Q3))
konpar_expand3_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand3_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 +  Q1:Q3) # now Q2:Q3 also is significant
konpar_expand3_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand3_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 +  Q1:Q3 + Q2:Q3)


konpar_expand4_qtl<-addtoqtl(konpar_jitter,konpar_expand3_qtl, 3, 60.7)
summary(fitqtl(konpar_jitter, qtl=konpar_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3))
konpar_expand4_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand4_qtl)
summary(fitqtl(konpar_jitter, qtl=konpar_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3))
konpar_expand4_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand4_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3)
konpar_expand4_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand4_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3 + Q1:Q5)
summary(fitqtl(konpar_jitter, qtl=konpar_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q1:Q3 + Q2:Q3 + Q1:Q5))

konpar_expand5_qtl<-addtoqtl(konpar_jitter,konpar_expand4_qtl, "X", 60)
summary(fitqtl(konpar_jitter, qtl=konpar_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q3 + Q2:Q3 + Q1:Q5))
konpar_expand5_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand5_qtl)
summary(fitqtl(konpar_jitter, qtl=konpar_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q3 + Q2:Q3 + Q1:Q5))
konpar_expand5_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand5_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q3 + Q2:Q3 + Q1:Q5)
konpar_expand5_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand5_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q1:Q3 + Q2:Q3 + Q1:Q5)


konpar_expand5b_qtl<-addtoqtl(konpar_jitter,konpar_expand5_qtl, 4, 50)
summary(fitqtl(konpar_jitter, qtl=konpar_expand5b_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5))
konpar_expand5b_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand5b_qtl)
summary(fitqtl(konpar_jitter, qtl=konpar_expand5b_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5))
konpar_expand5b_qtl_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand5b_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5)
konpar_expand5b_qtl_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand5b_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5)


konpar_expand6_qtl<-addtoqtl(konpar_jitter,konpar_expand5_qtl, 5, 57.6)
summary(fitqtl(konpar_jitter, qtl=konpar_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5))
konpar_expand6_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand6_qtl)
summary(fitqtl(konpar_jitter, qtl=konpar_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5))
konpar_expand6_qtl_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand6_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5)
konpar_expand6_qtl_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand6_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q1:Q3 + Q2:Q3 + Q1:Q5)

konpar_expand7_qtl<-addtoqtl(konpar_jitter,konpar_expand6_qtl, 2, 12.7)
summary(fitqtl(konpar_jitter, qtl=konpar_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5))
konpar_expand7_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand7_qtl)
summary(fitqtl(konpar_jitter, qtl=konpar_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q3 + Q2:Q3 + Q1:Q5)) # only Q2:Q3 is significant
summary(fitqtl(konpar_jitter, qtl=konpar_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q2 + Q2:Q3)) # Q1:Q2 is also significant
konpar_expand7_qtl_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand7_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q2:Q3) 
konpar_expand7_qtl_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand7_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q8 + Q2:Q3 + Q1:Q2)

cairo_pdf("konpar_fitqtl.pdf", width=8, height=3.5)
plotLodProfile(konpar_expand7_qtl,showallchr=TRUE, qtl.labels=FALSE); abline(h=3.46)
dev.off()

konpar_expand8_qtl<-addtoqtl(konpar_jitter,konpar_expand7_qtl, 6, 58.0)
summary(fitqtl(konpar_jitter, qtl=konpar_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3))
# strangely enough, refine here messes things up, lower LOD score, effect size, etc... So only refine QTL 1, which changes nothing (but is necessary for downstream stuff)
konpar_expand8_qtl<-refineqtl(konpar_jitter, pheno.col="pr", method="imp", qtl = konpar_expand8_qtl, formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3)
summary(fitqtl(konpar_jitter, qtl=konpar_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3)) 
konpar_expand8_qtl_qtl_addint<-addint(konpar_jitter, pheno.col="pr",  qtl=konpar_expand8_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3)  # Q3:Q7 is also significant
#konpar_expand8_qtl_qtl_add <- addqtl(konpar_jitter, qtl=konpar_expand8_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 	+ Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3)
summary(fitqtl(konpar_jitter, qtl=konpar_expand8_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q2 + Q2:Q3 + Q3:Q7))

cairo_pdf("konpar_fitqtl_new.pdf", width=8, height=3.5)
plotLodProfile(konpar_expand8_qtl,showallchr=TRUE, qtl.labels=FALSE); abline(h=3.46)
dev.off()

konpar_QTLbayesintervals<-data.frame(LG=character(), locus=character(), position=numeric())
CIsize_konpar<-c()
for(i in 1:9) {
	temp_df<-data.frame(LG= bayesint(konpar_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1], locus=names(pull.map(konpar_jitter)[[as.character(bayesint(konpar_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]),position=matrix(pull.map(konpar_jitter)[[as.character(bayesint(konpar_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]))

	from=min(bayesint(konpar_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	to=max(bayesint(konpar_expand8_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	CIsize_konpar<-c(CIsize_konpar,to-from)

	konpar_QTLbayesintervals<-rbind(konpar_QTLbayesintervals,temp_df[max(which(temp_df$position<=from)):min(which(temp_df$position>=to)),])
	}
	
konpar_QTLbayesintervals$scaffold<-gsub("_.*","",konpar_QTLbayesintervals$locus)
CIsize_konpar<-CIsize_konpar[c(2,8,5,4,3,1,7,9,6)]
write.table(konpar_QTLbayesintervals,"konpar_QTLbayesintervals.txt", quote=FALSE, row.names=FALSE, sep="\t")



fitdistr(konpar_QTLeffects$est, "exponential") # fit exponential distribution to data
cairo_pdf("konpar_effectsize_distr.pdf", width=3, height=3)
plot(density(rexp(10000, rate=fitdistr(konpar_QTLeffects$est, "exponential")$estimate)),xlim=c(0,1))
hist(konpar_QTLeffects$est, prob=TRUE, add=TRUE)
dev.off()

true_loci_prukoh<-otto_jones_ci(D=0.73,M=0.071595,nd=9,alpha=0.05,amin=0.018894,res=4,res2=2,max.loci=100) # for otto_jones_ci() script see https://github.com/thomasblankers/inter-island-mating-behavior-genetics-Laupala/blob/master/QTL_power_detect.r

### prukoh ###
# prukoh has been coded the right way, L. kohalensis alleles are A, L. pruna alleles are B

prukoh_newLG<-read.cross(format="csv", file="PruKoh_Comprehensive_CrossFile_newLG1.csv", estimate.map=FALSE)
prukoh_newLG_jitter<-jittermap(prukoh_newLG)
prukoh_newLG_jitter<-sim.geno(prukoh_newLG_jitter, step=1, n.draws=1000, error.prob=0.001)
prukoh_newLG_jitter<-calc.genoprob(prukoh_newLG_jitter, step=1, error.prob=0.001)

## scanone

qtl_prukoh_newLG_imp<-scanone(prukoh_newLG_jitter,method="imp", pheno.col="pr")
qtl_prukoh_newLG_hk<-scanone(prukoh_newLG_jitter,method="hk", pheno.col="pr")
#qtl_prukoh_newLG_ehk<-scanone(prukoh_newLG_jitter,method="ehk", pheno.col="pr")
#qtl_prukoh_newLG_em<-scanone(prukoh_newLG_jitter,method="em", pheno.col="pr")

qtl_prukoh_newLG_imp_1000perm<-scanone(prukoh_newLG_jitter,method="imp", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE) # this takes about 2 minutes
qtl_prukoh_newLG_hk_1000perm<-scanone(prukoh_newLG_jitter,method="hk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
#qtl_prukoh_newLG_em_1000perm<-scanone(prukoh_newLG_jitter,method="em", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
#qtl_prukoh_newLG_ehk_1000perm<-scanone(prukoh_newLG_jitter,method="ehk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)


## cim 

#cim_prukoh_newLG_imp<-cim(prukoh_newLG_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=10)
cim_prukoh_newLG_hk<-cim(prukoh_newLG_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20)
#cim_prukoh_newLG_ehk<-cim(prukoh_newLG_jitter, pheno.col="pr", method="ehk", n.marcovar=3, window=10)
#cim_prukoh_newLG_em<-cim(prukoh_newLG_jitter, pheno.col="pr", method="em", n.marcovar=3, window=10)

cim_prukoh_newLG_hk_1000perm<-cim(prukoh_newLG_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20, n.perm=1000)

cairo_pdf("cim_prukoh_newLG_hk_5_20.pdf", width=8, height=3.5)
plot(cim_prukoh_newLG_hk, ylim=c(0,35));abline(h=summary(cim_prukoh_newLG_hk_1000perm)[1,1])
dev.off()

qtl_prukoh_newLG<-c("1","2","3","4","5","6","X")

## scantwo
twoqtl_prukoh_newLG_imp<-scantwo(prukoh_newLG_jitter, pheno.col="pr", method="imp", chr=qtl_prukoh_newLG, clean.output=TRUE)
twoqtl_prukoh_newLG_hk<-scantwo(prukoh_newLG_jitter, pheno.col="pr", method="hk", chr=qtl_prukoh_newLG, clean.output=TRUE)

plot(twoqtl_prukoh_newLG_imp, lower="cond-int") # strong evidence for interaction chr1:chr2
plot(twoqtl_prukoh_newLG_hk, lower="cond-int")

plot(twoqtl_prukoh_newLG_hk, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_newLG_hk, chr=2, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_newLG_hk, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_newLG_hk, chr=4, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_newLG_hk, chr=5, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_newLG_hk, chr=6, lower="cond-int", upper="cond-add")
plot(twoqtl_prukoh_newLG_hk, chr="X", lower="cond-int", upper="cond-add")

# permutate the scantwo function for LOD  penalties
twoqtl_prukoh_newLG_hk_1000perm<-scantwo(prukoh_newLG_jitter, pheno.col="pr", method="hk", chr=qtl_prukoh_newLG, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes
summary(twoqtl_prukoh_newLG_hk, threshold=as.matrix(as.data.frame(summary(twoqtl_prukoh_newLG_hk_10000perm, 0.05)[1:5])[1,]))
penalties=calc.penalties(twoqtl_prukoh_newLG_hk_1000perm, alpha=0.05)


    main    heavy    light 
3.389887 5.461719 3.136901

max(qtl_prukoh_newLG_hk)

prukoh_newLG_init_qtl<-makeqtl(prukoh_newLG_jitter,chr=2, pos=71)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1))
prukoh_newLG_init_qtl_add <- addqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_init_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1)

prukoh_newLG_expand_qtl<-addtoqtl(prukoh_newLG_jitter,prukoh_newLG_init_qtl, 1, 46)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
prukoh_newLG_expand_qtl<-refineqtl(prukoh_newLG_jitter, pheno.col="pr", method="imp", qtl = prukoh_newLG_expand_qtl)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2))
prukoh_newLG_expand_qtl_addint<-addint(prukoh_newLG_jitter, pheno.col="pr",  qtl=prukoh_newLG_expand_qtl, method="imp", formula=  y ~ Q1 + Q2)
prukoh_newLG_expand_qtl_add <- addqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2)


prukoh_newLG_expand2_qtl<-addtoqtl(prukoh_newLG_jitter,prukoh_newLG_expand_qtl, 6, 123)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
prukoh_newLG_expand2_qtl<-refineqtl(prukoh_newLG_jitter, pheno.col="pr", method="imp", qtl = prukoh_newLG_expand2_qtl)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3))
prukoh_newLG_expand2_qtl_addint<-addint(prukoh_newLG_jitter, pheno.col="pr",  qtl=prukoh_newLG_expand2_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3)
prukoh_newLG_expand2_qtl_qtl_add <- addqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand2_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3)

prukoh_newLG_expand3_qtl<-addtoqtl(prukoh_newLG_jitter,prukoh_newLG_expand2_qtl, 5, 87)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4))
prukoh_newLG_expand3_qtl<-refineqtl(prukoh_newLG_jitter, pheno.col="pr", method="imp", qtl = prukoh_newLG_expand3_qtl)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4))
prukoh_newLG_expand3_qtl_addint<-addint(prukoh_newLG_jitter, pheno.col="pr",  qtl=prukoh_newLG_expand3_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4)
prukoh_newLG_expand3_qtl_add <- addqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand3_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4)


prukoh_newLG_expand4_qtl<-addtoqtl(prukoh_newLG_jitter,prukoh_newLG_expand3_qtl, "X", 12)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))
prukoh_newLG_expand4_qtl<-refineqtl(prukoh_newLG_jitter, pheno.col="pr", method="imp", qtl = prukoh_newLG_expand4_qtl)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5))
prukoh_newLG_expand4_qtl_addint<-addint(prukoh_newLG_jitter, pheno.col="pr",  qtl=prukoh_newLG_expand4_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5)
prukoh_newLG_expand4_qtl_add <- addqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand4_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5)

prukoh_newLG_expand5_qtl<-addtoqtl(prukoh_newLG_jitter,prukoh_newLG_expand4_qtl, 4, 85)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
prukoh_newLG_expand5_qtl<-refineqtl(prukoh_newLG_jitter, pheno.col="pr", method="imp", qtl = prukoh_newLG_expand5_qtl)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
prukoh_newLG_expand5_qtl_addint<-addint(prukoh_newLG_jitter, pheno.col="pr",  qtl=prukoh_newLG_expand5_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)
prukoh_newLG_expand5_qtl_add <- addqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand5_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)

prukoh_newLG_expand6_qtl<-addtoqtl(prukoh_newLG_jitter,prukoh_newLG_expand5_qtl, 3,79)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7))
prukoh_newLG_expand6_qtl<-refineqtl(prukoh_newLG_jitter, pheno.col="pr", method="imp", qtl = prukoh_newLG_expand6_qtl)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7))
prukoh_newLG_expand6_qtl_addint<-addint(prukoh_newLG_jitter, pheno.col="pr",  qtl=prukoh_newLG_expand6_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7)
prukoh_newLG_expand6_qtl_add <- addqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand6_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7)

prukoh_newLG_expand7_qtl<-addtoqtl(prukoh_newLG_jitter,prukoh_newLG_expand6_qtl, 2, 140)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))
prukoh_newLG_expand7_qtl<-refineqtl(prukoh_newLG_jitter, pheno.col="pr", method="imp", qtl = prukoh_newLG_expand7_qtl)
summary(fitqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand7_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8))
prukoh_newLG_expand7_qtl_addint<-addint(prukoh_newLG_jitter, pheno.col="pr",  qtl=prukoh_newLG_expand7_qtl, method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8)
prukoh_newLG_expand7_qtl_add <- addqtl(prukoh_newLG_jitter, qtl=prukoh_newLG_expand7_qtl, pheno.col="pr", method="imp", formula=  y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8)


prukoh_QTLbayesintervals<-data.frame(LG=character(), locus=character(), position=numeric())
CIsize_prukoh<-c()
for(i in 1:8) {
	temp_df<-data.frame(LG= bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1],locus=names(pull.map(prukoh_newLG_jitter)[[as.character(bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]),position=matrix(pull.map(prukoh_newLG_jitter)[[as.character(bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$chr[1])]]))

	from=min(bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	to=max(bayesint(prukoh_newLG_expand7_qtl,qtl.index=i,expandtomarkers=TRUE)$pos)
	
	CIsize_prukoh=c(CIsize_prukoh,to-from)

	prukoh_QTLbayesintervals<-rbind(prukoh_QTLbayesintervals,temp_df[max(which(temp_df$position<=from)):min(which(temp_df$position>=to)),])
	}
CIsize_prukoh<-CIsize_prukoh[c(1,7,8,2,6,4,3,5)]
prukoh_QTLbayesintervals$scaffold<-gsub("_.*","",prukoh_QTLbayesintervals$locus)

write.table(prukoh_QTLbayesintervals,"prukoh_QTLbayesintervals.txt", quote=FALSE, row.names=FALSE, sep="\t")



cairo_pdf("PruKoh_fitqtl.pdf", width=8, height=3.5)
plotLodProfile(prukoh_newLG_expand7_qtl,showallchr=TRUE, qtl.labels=FALSE); abline(h=3.38)
dev.off()

fitdistr(prukoh_QTLeffects$est, "exponential") # fit exponential distribution to data
cairo_pdf("prukoh_effectsize_distr.pdf", width=3, height=3)
plot(density(rexp(10000, rate=fitdistr(prukoh_QTLeffects$est, "exponential")$estimate)),xlim=c(0,1))
hist(prukoh_QTLeffects$est, prob=FALSE, add=TRUE, breaks=4)
dev.off()


fitdistr(cereuk_QTLeffects$est, "gamma") # fit gamma distribution to data
cairo_pdf("cereuk_effectsize_distr.pdf", width=3, height=3)
plot(density(rgamma(10000,shape=3.12, rate=37.13)),xlim=c(0,1))
hist(cereuk_QTLeffects$est, prob=TRUE, breaks=4, add=TRUE)
dev.off()

true_loci_prukoh<-otto_jones_ci(D=0.82,M=0.091081925,nd=8,alpha=0.05,amin=0.0479672,res=4,res2=2,max.loci=100)

######### probabilities of observed QTL overlap

qtl_overlap<-function(n_qtl=15,n_1=8,n_2=9,n_shared=6,res=1e5) {
	overlapping<-c()
	for(i in 1:res) {
		qtl1<-sample(c(1:n_qtl),n_1)
		qtl2<-sample(c(1:n_qtl),n_2)
		overlapping<-c(overlapping,length(which(qtl1 %in% qtl2)))
				}
	probability=(which(quantile(overlapping, seq(0,1,1/res))>=n_shared)[1]-1)/res
	return(1-probability)		
	}

qtl_overlap3<-function(n_qtl=12,n_1=8,n_2=9, n_3=7, n_shared=5, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) {
	overlapping<-c()
	unique_1<-c()
	unique_2<-c()
	unique_3<-c()
	for(i in 1:res) {
		qtl1<-sample(c(1:n_qtl),n_1)
		qtl2<-sample(c(1:n_qtl),n_2)
		qtl3<-sample(c(1:n_qtl),n_3)
		share12<-qtl1[which(qtl1 %in% qtl2)]
		overlapping<-c(overlapping,length(which(share12 %in% qtl3)))
		unique_1<-c(unique_1,length(qtl1)-length(which(qtl1 %in% qtl2 | qtl1 %in% qtl3)))
		unique_2<-c(unique_2,length(qtl2)-length(which(qtl2 %in% qtl1 | qtl2 %in% qtl3)))
		unique_3<-c(unique_3,length(qtl3)-length(which(qtl3 %in% qtl1 | qtl3 %in% qtl2)))
		}
	probability_sharing=1-(which(quantile(overlapping, seq(0,1,1/res))>=n_shared)[1]-1)/res
	probability_unique1=(max(which(quantile(unique_1, seq(0,1,1/res))<=n_unique1))-1)/res
	probability_unique2=(max(which(quantile(unique_2, seq(0,1,1/res))<=n_unique2))-1)/res
	probability_unique3=(max(which(quantile(unique_3, seq(0,1,1/res))<=n_unique3))-1)/res
	out=list(probability_sharing,probability_unique1,probability_unique2,probability_unique3)
	return(out)		
	}

# average konpar n_qtl estimate, 12 QTL:
qtl_overlap(n_qtl=12,n_1=9,n_2=8,n_shared=6,res=1e5)
#P = 0.74
qtl_overlap(n_qtl=12,n_1=9,n_2=7,n_shared=6,res=1e5)
#P = 0.36
qtl_overlap(n_qtl=12,n_1=7,n_2=8,n_shared=6,res=1e5)
#P = 0.15

# average konpar n_qtl estimate, 16 QTL:
qtl_overlap(n_qtl=16,n_1=9,n_2=8,n_shared=6,res=1e5)
#P = 0.
qtl_overlap(n_qtl=16,n_1=9,n_2=7,n_shared=6,res=1e5)
#P = 0.
qtl_overlap(n_qtl=16,n_1=7,n_2=8,n_shared=6,res=1e5)
#P = 0.

# average PruKoh n_qtl estimate, 21 QTL:
qtl_overlap(n_qtl=21,n_1=9,n_2=8,n_shared=6,res=1e5)
#P = 0.028
qtl_overlap(n_qtl=21,n_1=9,n_2=7,n_shared=6,res=1e5)
#P = 0.009
qtl_overlap(n_qtl=21,n_1=7,n_2=8,n_shared=6,res=1e5)
#P = 0.004


# minimum triple and unique:
qtl_overlap3(n_qtl=12,n_1=8,n_2=9, n_3=7, n_shared=5, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 


qtl_overlap3(n_qtl=12,n_1=8,n_2=9, n_3=7, n_shared=4, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 

# mean triple and unique:
qtl_overlap3(n_qtl=16,n_1=8,n_2=9, n_3=7, n_shared=5, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 


qtl_overlap3(n_qtl=16,n_1=8,n_2=9, n_3=7, n_shared=4, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 

# maximum triple and unique:
qtl_overlap3(n_qtl=21,n_1=8,n_2=9, n_3=7, n_shared=5, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 


qtl_overlap3(n_qtl=21,n_1=8,n_2=9, n_3=7, n_shared=4, n_unique1=1, n_unique2=2, n_unique3=0, res=1e5) 


# average CerEuk estimate, 15
#qtl_overlap(n_qtl=15,n_1=9,n_2=7,n_shared=6,res=1e5)
#P = 0.084	

#qtl_overlap(n_qtl=15,n_1=7,n_2=8,n_shared=6,res=1e5)
#P = 0.032


## overlap in the sense of randomly distributed windows:

average_windowsize=(11+29+41+19+40+26+70+8+14+8+46+14+1+15+2+4)/16 # based on konpar and PruKoh QTL, 95 Bayesian CI

LG_length_konpar<-c();for(i in 1:8) { LG_length_konpar<-c(LG_length_konpar,max(pull.map(konpar)[[i]]))}; sum(LG_length_konpar) # total map length konpar
LG_length_prukoh<-c();for(i in 1:8) { LG_length_prukoh<-c(LG_length_prukoh,max(pull.map(prukoh)[[i]]))}; sum(LG_length_prukoh) # total map length PruKoh

average_LG_length=((sum(LG_length_konpar)/8)+(sum(LG_length_prukoh)/8))/2

qtl_overlap_randomwindows<-function(chr_size=90,window_size=22,qtl_number=9,shared_qtl=6,nboot=100000) {
	overlapping<-c()
	for(n in 1:nboot) {
		sample_1_start<-sample(1:chr_size,1)
		sample_2_start<-sample(1:chr_size,1)
		sample1<-seq(1:chr_size)[sample_1_start:(sample_1_start+window_size)]
		sample2<-seq(1:chr_size)[sample_2_start:(sample_2_start+window_size)]
		if(length(which(sample1 %in% sample2)) > 1) { overlapping<-c(overlapping,n) }
		}
	single_prob=length(overlapping)/nboot
	iterative_prob=(single_prob^shared_qtl)*((1-single_prob)^(qtl_number-shared_qtl)) * ncol(combn(qtl_number, shared_qtl))
	output<-list()
	output$iterative_prob<-iterative_prob
	output$overlap_prob<-single_prob
	return(output)
	}

qtl_overlap_randomwindows(chr_size=average_LG_length,window_size=average_windowsize,qtl_number=9,shared_qtl=6,nboot=100000)
0.02463777
qtl_overlap_randomwindows(chr_size=average_LG_length,window_size=average_windowsize,qtl_number=8,shared_qtl=6,nboot=100000)
0.01198767
qtl_overlap_randomwindows(chr_size=average_LG_length,window_size=average_windowsize,qtl_number=7,shared_qtl=6,nboot=100000)
0.004445523


qtl_overlap_randomwindows(chr_size=150,window_size=10,qtl_number=8,shared_qtl=6,nboot=100)


# simulate generally representive linkage map
L <- c(150, 150, 150, 100, 100, 75, 75, 150)
simmap <- sim.map(L, L/2+1, eq.spacing=FALSE, include.x=TRUE) # randomly distributed markers on average every 2 cM

countOverlappingWindows_chain<-c()
nsims=1000
library(dplyr)
pb <- progress_estimated(nsims) # set progress bar

while(length(countOverlappingWindows_chain) < nsims) {
	
	for( k in 1:2) {
		# we want to simulate 16 QTL of exponentially distributed effect size
		# will take the average effect size across all detected and undetected QTL = 0.047 and convert that to 16 QTL of pure additive effects
		effect_sizes<-rexp(16,(1/0.047))
		percentage_explained<-effect_sizes/1.6 # this converts the average effect sizes to the fraction of the average phenotypic distance in the crosses explained
		model_effects <- 2 * sqrt(percentage_explained / (1 - 2 * percentage_explained)) # rQTL models phenotypes as normally distributed with variance = 1. This converts the fraction explained to that phenotypic scale (see rQTL manual for more detail)

		# randomly draw markers that are associated with the QTL, QTL cannot be within 30 cM

		simmap_df=data.frame(unlist(simmap)) # convert linkage map to data frame
		simmap_df$marker=rownames(simmap_df)
		simmap_df$LG=substr(simmap_df$marker,start=1,stop=1)
		colnames(simmap_df)[1]="position"
		
		# randomly place QTLs
		QTLs<-c()
		drawmap<-simmap_df
		drawmap$marker=rownames(drawmap)
		drawmap[,1]=drawmap[,1]+as.numeric(gsub("X","8",substr(rownames(drawmap),start=1, stop=1)))*500 # raise every position by a factor based on the LG (where X is 8) so that marker positions are continuous non-overlapping
		for(i in 1:16) {
			QTLs<-c(QTLs,sample(drawmap[,"marker"],1))
			exclude_positions<-c(drawmap[which(drawmap[,"marker"] == QTLs[i]),1]-30,drawmap[which(drawmap[,"marker"] == QTLs[i]),1]+30)
			drawmap<-drawmap[c(which(drawmap[,1] < exclude_positions[1]),which(drawmap[,1] > exclude_positions[2])),]
			}

		# create QTL model from effect sizes and randomly drawn QTL locations
		simmap_df$LG<-gsub("X",8,simmap_df$LG)
		mymodel=cbind(as.numeric(simmap_df[QTLs,"LG"]),as.numeric(simmap_df[QTLs,"position"]),model_effects,rep(0,16))

		# create QTL object with 200 F2s
		simcross <- sim.cross(simmap, type="f2", n.ind=200, model=mymodel)
		
		# perform QTL scane
		simcross<-calc.genoprob(simcross)
		#simqtl_hk<-scanone(simcross, method="hk") # single interval mapping
		#simqtl_cim<-cim(simcross, method="hk", window=20,n.marcovar=3) # composite interval mapping
		#simqtl_mqm<-mqmscan(simcross,model="additive", cofactors=mqmautocofactors(simcross)) # multiple QTL mapping
		
		simqtl_stepwise<-stepwiseqtl(simcross,method="hk",model="normal",additive.only=TRUE,penalties=c(3,5,3), verbose=FALSE) # fit a stepwise MQM
		# calculate bayesian confidence intervals around QTL
		assign(paste0("bayesian_intervals_",k),bayesint(simqtl_stepwise,qtl.index=1))
		for( i in 2:simqtl_stepwise$n.qtl) {
			assign(paste0("bayesian_intervals_",k),rbind(get(paste0("bayesian_intervals_",k)),bayesint(simqtl_stepwise,qtl.index=i)))
			}
		}
	# check for overlap between confidence intervals	
	library(IRanges)
	shared_LGs_1_2<-bayesian_intervals_1[which(bayesian_intervals_1$chr %in% bayesian_intervals_2$chr),] # retain only shared LGs
	shared_LGs_2_1<-bayesian_intervals_2[which(bayesian_intervals_2$chr %in% bayesian_intervals_1$chr),]
	shared_LGs_1_2$chr<-factor(shared_LGs_1_2$chr) # reset factor levels
	shared_LGs_2_1$chr<-factor(shared_LGs_2_1$chr)

	countOverlappingWindows=0
	for( chr in levels(shared_LGs_1_2$chr)) {
		temp_1_2<-shared_LGs_1_2[which(shared_LGs_1_2$chr == chr),] # only bayesian confidence interval(s) on the n-th linkage group
		temp_2_1<-shared_LGs_2_1[which(shared_LGs_2_1$chr == chr),]
		for( i in 1:(nrow(temp_1_2)/3)) { # there may be more than one QTL per linkage group, so check groups of three (begin, peak, end of interval) positions
			assign(paste0("range_1_2_",i),IRanges(temp_1_2[(i*3)-2,"pos"],temp_1_2[i*3,"pos"]))
			for( j in 1:(nrow(temp_2_1)/3)) {
				assign(paste0("range_2_1_",j),IRanges(temp_2_1[(j*3)-2,"pos"],temp_2_1[j*3,"pos"]))
				if(get(paste0("range_1_2_",i)) %over% get(paste0("range_2_1_",j)) == TRUE) {countOverlappingWindows=countOverlappingWindows+1} # count if there is overlap between the windows in the first and second QTL scan
				}
			}
		}
		
	countOverlappingWindows_chain<-c(countOverlappingWindows_chain,countOverlappingWindows)
	pb$pause(0.01)$tick()$print() # update progress bar
	}



############################# match qtl effect sizes
library(ggplot2)

QTLeffects_sharing<-read.delim("6sppQTL_effectsizes_sharing.txt")

#Fig 3A:
# cross factor order: cereuk konpar parkoh prukoh
pps<-c(1.66,1.46,3.01,1.64)
detectedQTLnumber<-c(7,9,8,8)
totalQTLnumber<-c(17,13,14,20)
avg_effect<-aggregate(na.omit(QTLeffects_sharing)[,"effect"], by=list(na.omit(QTLeffects_sharing)[,"cross"]), FUN=mean)$x
lm_phenodist_effect<-lm(avg_effect~pps)

Residuals:
        1         2         3         4 
-0.001028 -0.005206 -0.000669  0.006903 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)  
(Intercept)  0.01696    0.01014   1.673   0.2363  
pps          0.04099    0.00497   8.247   0.0144 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.006175 on 2 degrees of freedom
Multiple R-squared:  0.9714,    Adjusted R-squared:  0.9572 
F-statistic: 68.01 on 1 and 2 DF,  p-value: 0.01439

Response: avg_effect
          Df     Sum Sq    Mean Sq F value  Pr(>F)  
pps        1 0.00259336 0.00259336  68.014 0.01439 *
Residuals  2 0.00007626 0.00003813 

phenodist_effect_plot<-ggplot(data.frame(avg_effect=avg_effect,pps=pps),aes(x=avg_effect,y=pps)) +
geom_point() +
geom_smooth(method="lm",fullrange = TRUE) +
scale_x_continuous(limits=c(0,0.15)) + 
scale_y_continuous(limits=c(0,3.5)) + 
theme(legend.position="none")

lm_phenodist_detectedQTLnumber<-lm(detectedQTLnumber~pps)

Residuals:
       1        2        3        4 
-1.03660  0.93749  0.13831 -0.03919 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)  
(Intercept)   8.2517     1.6308   5.060   0.0369 *
pps          -0.1296     0.7996  -0.162   0.8862  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.9935 on 2 degrees of freedom
Multiple R-squared:  0.01296,   Adjusted R-squared:  -0.4806 
F-statistic: 0.02625 on 1 and 2 DF,  p-value: 0.8862

lm_phenodist_totalQTLnumber<-lm(totalQTLnumber~pps)

Residuals:
      1       2       3       4 
 0.6010 -3.6814 -0.4925  3.5728 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)  
(Intercept)   18.743      6.022   3.112   0.0896 .
pps           -1.412      2.953  -0.478   0.6797  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 3.669 on 2 degrees of freedom
Multiple R-squared:  0.1026,    Adjusted R-squared:  -0.3461 
F-statistic: 0.2287 on 1 and 2 DF,  p-value: 0.6797

# Fig 3B:

QTLeffects_cereuk_prukoh<-na.omit(merge(QTLeffects_sharing[which(QTLeffects_sharing$cross=="cereuk"),],QTLeffects_sharing[which(QTLeffects_sharing$cross=="prukoh"),],"QTL",sort=FALSE))
colnames(QTLeffects_cereuk_prukoh)[c(4,7)]<-c("effect_cereuk","effect_prukoh")
QTLeffects_cereuk_konpar<-na.omit(merge(QTLeffects_sharing[which(QTLeffects_sharing$cross=="cereuk"),],QTLeffects_sharing[which(QTLeffects_sharing$cross=="konpar"),],"QTL",sort=FALSE))
colnames(QTLeffects_cereuk_konpar)[c(4,7)]<-c("effect_cereuk","effect_konpar")
QTLeffects_konpar_prukoh<-na.omit(merge(QTLeffects_sharing[which(QTLeffects_sharing$cross=="konpar"),],QTLeffects_sharing[which(QTLeffects_sharing$cross=="prukoh"),],"QTL",sort=FALSE))
colnames(QTLeffects_konpar_prukoh)[c(4,7)]<-c("effect_konpar","effect_prukoh")

c(0,0.225)
QTL_effects_plot<-ggplot() +
 geom_point(data=QTLeffects_cereuk_konpar,aes(x=effect_cereuk,y=effect_konpar), color="blue") +
 geom_smooth(data=QTLeffects_cereuk_konpar,aes(x=effect_cereuk,y=effect_konpar), method="lm",color="blue", fullrange=TRUE) +
 geom_point(data=QTLeffects_konpar_prukoh,aes(x=effect_konpar,y=effect_prukoh), color="red") +
 geom_smooth(data=QTLeffects_konpar_prukoh,aes(x=effect_konpar,y=effect_prukoh), method="lm",color="red", fullrange=TRUE) + 
 geom_point(data=QTLeffects_cereuk_prukoh,aes(x=effect_cereuk,y=effect_prukoh), color="black") +
 geom_smooth(data=QTLeffects_cereuk_prukoh,aes(x=effect_cereuk,y=effect_prukoh), method="lm",color="black", fullrange=TRUE) +
 xlab("effect size (pps) cereuk (black/blue) / konpar (red)") +
 ylab("effect size (pps) prukoh (black/red) / konpar (blue)") +
 theme(legend.position="none")

# Fig 3 C:
QTLeffects_sharing_plot<-ggplot(QTLeffects_sharing) + 
geom_point(aes(x=factor(sharing,levels=c("unique","double","triple")), y=effect, color=as.factor(cross))) +
#geom_smooth(data=data.frame(na.omit(QTLeffects_sharing)),aes(x=sharing, y=effect), method="lm") +
theme(legend.position="none")

abline(lm(effect~sharing, data=QTLeffects_sharing[which(QTLeffects_sharing$cross != "parkoh"),])
lm_size_sharing<-lm(effect~cross+sharing, data=QTLeffects_sharing[which(QTLeffects_sharing$cross != "parkoh"),])


Residuals:
      Min        1Q    Median        3Q       Max 
-0.068642 -0.034331 -0.009349  0.015487  0.149176 

Coefficients:
               Estimate Std. Error t value Pr(>|t|)  
(Intercept)    0.074914   0.028328   2.645    0.016 *
crosskonpar   -0.006880   0.029357  -0.234    0.817  
crossprukoh    0.008244   0.030042   0.274    0.787  
sharingtriple  0.023686   0.026510   0.893    0.383  
sharingunique -0.031348   0.035623  -0.880    0.390  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05788 on 19 degrees of freedom
  (12 observations deleted due to missingness)
Multiple R-squared:  0.1488,    Adjusted R-squared:  -0.03036 
F-statistic: 0.8306 on 4 and 19 DF,  p-value: 0.5222

Response: effect
          Df   Sum Sq   Mean Sq F value Pr(>F)
cross      2 0.001659 0.0008296  0.2476 0.7832
sharing    2 0.009472 0.0047362  1.4135 0.2677
Residuals 19 0.063662 0.0033506 

library(gridExtra)
pdf("6sppQTL_Fig3_QTLeffects_sharing.pdf", width=4, height=9.5)
grid.arrange(phenodist_effect_plot,QTL_effects_plot,QTLeffects_sharing_plot, ncol=1)
dev.off()

library(lme4)
lme_size_sharing<-lmer(effect~sharing + (1|cross), data=QTLeffects_sharing[which(QTLeffects_sharing$cross != "parkoh"),])
lme_null<-lmer(effect~ (1|cross), data=QTLeffects_sharing[which(QTLeffects_sharing$cross != "parkoh"),])
anova(lme_size_sharing,lme_null)

Models:
lme_null: effect ~ (1 | cross)
lme_size_sharing: effect ~ sharing + (1 | cross)
                 Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
lme_null          3 -64.397 -60.863 35.198  -70.397                         
lme_size_sharing  5 -63.906 -58.016 36.953  -73.906 3.5092      2      0.173



