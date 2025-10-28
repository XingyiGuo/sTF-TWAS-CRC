library(data.table)
library(dplyr)
library(renv)
renv::restore() 
args = commandArgs(trailingOnly = TRUE)

snps=as.character(args[1])
wk="../../data/simulation_values_res/ldsc_res"
#load ldsc res
format_ldsc <- function(snps) {
  beta1="NULL"
  ldsc_fd=paste0(wk,"/y_null_",snps)
  files=dir(ldsc_fd)
  ldsc_res=matrix(integer(0), nrow=length(files), ncol=6) %>% as.data.frame()
  names(ldsc_res) <- c("Beta1", "CistromeID", "Coeff", "CoeffSE", "Zvalue", "Pvalue")
  count=1
  for(file in files){
    sim_y_track= sub("\\..*$", "", file) 
    lines<-readLines(paste0(wk,"/y_null_",snps,"/", file))
    Cates<- Coefs <-Coefs_se <- NULL
    for(line in lines){
		if (startsWith(line, "Total Observed scale h2:")) {
		  Coefs <- as.numeric(sub(".*h2:\\s*([-0-9.eE]+)\\s*\\(.*", "\\1", line))
		  Coefs_se <- as.numeric(sub(".*\\(([-0-9.eE]+)\\).*", "\\1", line))
		  print(line)
		  print(Coefs)
		  print(Coefs_se)
		  z=Coefs/Coefs_se
		  p_=0
		  if(z<0){ p_ <- pnorm(z)}else{p_ <- 1 - pnorm(z)}
		  ldsc_res[count,]=c(beta1, as.character(sim_y_track), Coefs, Coefs_se, z, p_)
		}
    }
  count=count+1
  }
  return(data.frame(ldsc_res))
}
ldsc_res = format_ldsc(snps)
fwrite(as.data.frame(ldsc_res), paste0(wk, "/y_null_",snps,".ldsc.csv"))