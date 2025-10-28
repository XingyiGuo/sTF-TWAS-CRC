#########################################################
####  Computer system requirements                    ###
####  R package:  data.table, lme4                    ###
####  Unix operating system environment               ###
####  20 GB memory for per_geno, 120G for whole geno  ### 
#########################################################


### Senario1： Based on Genotype
# conda activate tftwas
library(data.table) # read in large data matrix fast
library(lme4)       # for mixed model
library(dplyr)
library(renv)
renv::restore()
args = commandArgs(trailingOnly = TRUE)

### Load TF annotation matrix
TFs_annots<-as.data.frame(fread(paste0("../../data/84Tracks/annots/1_22.annot.txt")))
# remove duplicate variants 
TFs_annots <- TFs_annots[!duplicated(TFs_annots$SNP),]
TFs_annots<-TFs_annots[TFs_annots$CHR>=1 & TFs_annots$CHR<=22,]

### Select causal SNPs and make causal SNP annotation files
set.seed(20251022) 
n_causal=as.numeric(args[1])
n_round=50
beta_1=1
h2=as.numeric(args[2]) # 0.01, 0.05, 0.1, 0.25, 0.5
SNP_causal_pool = TFs_annots[TFs_annots$ANNO==1,]$SNP
sim_y_matrix = matrix(integer(0), nrow = n_round , ncol= 489) %>% as.data.frame()

### Generate annot files (TF occupied SNPs as 1 and others as 0)
annot_dir <- paste0("./84Tracks/annots_sim", n_causal, "/")
annot_list <- vector("list", 22)
for (chr in 1:22) {
  file_path <- paste0(annot_dir, chr, ".annot.txt")
  if (file.exists(file_path)) {
    message("Reading existing file: ", file_path)
    annot_list[[chr]] <- fread(file_path)
  } else {
    message("Generating new annotation for chr ", chr)
    DT <- as.data.table(TFs_annots[, c("CHR", "BP", "SNP", "CM")])
    simNames <- paste0("SIM", seq_len(n_round))
    causal_mat <- matrix(0L, nrow = nrow(DT), ncol = n_round)
    causal_list <- replicate(n_round, sample(SNP_causal_pool, n_causal), simplify = FALSE)
    snp_index <- match(unlist(causal_list), DT$SNP)
    round_index <- rep(seq_along(causal_list), each = n_causal)
    causal_mat[cbind(snp_index, round_index)] <- 1L
    causal_dt <- as.data.table(causal_mat)
    setnames(causal_dt, paste0("SIM", seq_len(n_round)))
    SNP_causal_annots <- cbind(DT, causal_dt)
    SNP_causal_annots_per_chr <- SNP_causal_annots[SNP_causal_annots$CHR == as.character(chr), ]
    fwrite(SNP_causal_annots_per_chr, file_path, sep = "\t")
    annot_list[[chr]] <- SNP_causal_annots_per_chr
  }
}
SNP_causal_annots <- rbindlist(annot_list)
simNames <- colnames(SNP_causal_annots)[5:length(colnames(SNP_causal_annots))]

for(sim_i in 1:n_round){
	simName=paste0("SIM", sim_i)
	print(simName)
	# >> make ldscores absed on the simulated annotations
	### Simulate Y based on causal SNPs
	SNP_causal_geno=data.frame()
	for(chr in 1:22){
		SNP_causal_chr <- SNP_causal_annots[get(simName) == 1 & CHR == as.character(chr),SNP]
		print(paste0("Num of annotated SNPs is for ", chr, " is ",length(SNP_causal_chr)))
		chr_geno=data.frame(fread(paste0("./ref/Alkesgroup/1000G_EUR_Phase3_plink_hg38/processed/",chr,".nodup.sorted.1000G_hg38.traw")))
		SNP_causal_chr_geno = chr_geno[chr_geno$SNP %in% SNP_causal_chr, ]
		rownames(SNP_causal_chr_geno) <- SNP_causal_chr_geno$SNP
		SNP_causal_chr_geno_clean <- SNP_causal_chr_geno[, 7:dim(SNP_causal_chr_geno)[2]] #select only genotype
		SNP_causal_geno=rbind(SNP_causal_geno, SNP_causal_chr_geno_clean)
		rm(chr_geno)
	}
	print(paste0("Total num of annotated SNPs is ",dim(SNP_causal_geno)[1]))
	SNP_causal_geno_r0a1 <- 2 - SNP_causal_geno  ##Chek the coding for ref and alt from *.traw and *.tped file
	tmp = sub("_(.*)", "", colnames(SNP_causal_geno_r0a1))
	colnames(SNP_causal_geno_r0a1) <- tmp
	sample_size= dim(SNP_causal_geno_r0a1)[2]
	y_geno <-colSums(beta_1 * SNP_causal_geno_r0a1)
	y_err <- rnorm(sample_size, mean=0, sd=sqrt((1-h2)/h2*var(y_geno)))
	y_sim = y_geno + y_err
	#y_sim=y_geno
	y_sim_cc <- ifelse(y_sim >= median(y_sim), 1, 2) #1 for case, 2 for control
	sim_y_matrix[sim_i, ]= y_sim_cc
}
sim_y_matrix_t = transpose(sim_y_matrix)
sim_y_df_t <- as.data.frame(sim_y_matrix_t)  #sample N x simulate times t
sim_y_df_t["IID"]=names(y_sim_cc)  # Add IID column
sim_y_df_t_IID = setcolorder(sim_y_df_t, c("IID", setdiff(names(sim_y_df_t), "IID"))) # reorder col
colnames(sim_y_df_t_IID) <- c("IID", simNames)
fwrite(sim_y_df_t_IID,paste0("../../data/simulation_values_res/sim_y_values/beta1_1_h2_",h2,"_SNPs_",n_causal,".sim_y.txt"), sep="\t", row.names = FALSE, col.names = TRUE)


### Senario2： Not based on Genotype
### Generate a 100 x 489 matrix of 0/1 sampled from Binomial(1, 0.5)
sim_y_null_matrix <- matrix(
  rbinom(489 * 100, size = 1, prob = 0.5), 
  nrow = 489, 
  ncol = 100
) %>% as.data.frame()
sim_y_null_matrix=sim_y_null_matrix+1 #1 for case, 2 for control defined by PLINK
colnames(sim_y_null_matrix) <- paste0("SIM", 1:100)
prefix=fread("./y_prefix.txt")
sim_y_null_matrix_psam = cbind(prefix, sim_y_null_matrix)
fwrite(sim_y_null_matrix_psam,paste0("../../data/simulation_values_res/sim_y_values/sim_y_null.psam"), sep="\t", row.names = FALSE, col.names = TRUE)
