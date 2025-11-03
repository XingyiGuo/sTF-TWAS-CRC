##################
# conda activate tftwas_coloc
# 80G
##################

library(data.table)
library(lme4)
library(dplyr)
library(renv)
library(parallel)
renv::restore() 
args <- commandArgs(trailingOnly=TRUE)

set.seed(123)
KB<-100000
glm_res<-matrix(integer(0), nrow=900, ncol=6) %>% as.data.frame()
names(glm_res) <- c("Beta1", "CistromeID", "Slope", "SD", "Tvalue", "Pvalue")
wk="../../data"
snps=as.character(args[1])

###load sim_ys
y_sim=data.frame(fread(paste0(wk, "/simulation_values_res/sim_y_values/sim_y_null.psam"), header = TRUE))
tmp=colnames(y_sim)[9:length(colnames(y_sim))][1:1000]
sim_num = sub("^X", "", tmp)

###load annots
annots=data.frame()
for(chr in 1:22){
	annot_per_chr=fread(paste0(wk,"/84Tracks/annots_sim",snps,"/",chr,".annot.txt"))
	annots=rbind(annots, annot_per_chr)
}
annots=as.data.frame(annots)


count=1
beta1="NULL"
for(y_track in sim_num){
	print(y_track)
	gwas_file=paste0(wk,"/simulation_values_res/sim_y_values/sim_y_null_glm/null.",y_track,".glm.logistic.hybrid")
	###load GWAS SS
	if(file.exists(gwas_file)){
		gwas_ss=fread(gwas_file)
		names(gwas_ss)[names(gwas_ss)=="ID"]<-"SNP"
		gwas_ss <- gwas_ss[!is.na(Z_STAT) & `#CHROM` >= 1 & `#CHROM` <= 22]  #Z_STAT  P
		# gwas_ss$A2 = ifelse(gwas_ss$A1==gwas_ss$ALT, gwas_ss$REF, gwas_ss$ALT)
		# gwas_ss$N=489
		# ldsc_gwas_ss=gwas_ss[,c("SNP", "N", "A1", "A2", "Z_STAT", "P")]
		# colnames(ldsc_gwas_ss) <- c("SNP", "N", "A1", "A2", "Z", "P")
		# fwrite(ldsc_gwas_ss, paste0(wk,"/simulation_values_res/sim_y_values/sim_y_null_glm/null.",y_track,".glm.logistic.hybrid.ldsc.txt"), sep="\t")
		
		###merge GWAS SS and Annots
		
		# As the y is not based on genotype, we can use any TF genotype. Here we use "SIM1"
		gwas_ss_annots=left_join(gwas_ss, annots[, c("SNP", "SIM1")], by = "SNP")
		gwas_ss_annots$loci <- paste0(gwas_ss_annots[["#CHROM"]],'_',floor(gwas_ss_annots[["POS"]]/KB)) # CHR and POS are from annots files

		###Run glm and output
		out <- summary(lmer(I(gwas_ss_annots$Z_STAT^2) ~ gwas_ss_annots[["SIM1"]]+(1|gwas_ss_annots$loci),control = lmerControl(calc.derivs = FALSE)))
		p_ <- 2 * (1 - pnorm(abs(out$coef[2,1])))  		
		print(c(beta1, y_track, out$coef[2,1], out$coef[2,2], out$coef[2,3], p_))
		glm_res[count,]=c(beta1, y_track, out$coef[2,1], out$coef[2,2], out$coef[2,3], p_)
		count=count+1
	}
}
fwrite(as.data.frame(glm_res), paste0(wk, "/simulation_values_res/lmer_res/TFs_",snps,"_SNPs_y_null.lmer.csv"))
