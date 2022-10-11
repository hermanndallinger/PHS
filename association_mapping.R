#------------------------------------------------------------------------------#
# title: "QTL mapping and GWAS in PHS dataset, using sommer package"
# performed on unimputed marker matrix
# author: "hgd"
#------------------------------------------------------------------------------#

library(reshape2)
library(plyr)
library(dplyr)
library(data.table)

# devtools::install_version("sommer", "3.9.3", force = FALSE)
library(sommer)

# -------------------------------------------------------------------------------------------------
options(max.print=100)
rm(list=ls()); gc()

#------------------------------------------------------------------------------#
setwd(paste(dirname(rstudioapi::getSourceEditorContext()$path), "/", sep=""))


## read data and prepare dataframes --------------------------------------------
blups <- read.csv( file = "PHS-trial_means.csv")[,-1]
results <- read.csv( file = "PHS-trial_means_stats.csv")[,-1]

# minor allele frequency selection --------------------------------------------
af <- function(x) sum(x+1, na.rm= TRUE)/(2*length(na.omit(x)))


# read clean marker files ------------------------------------------------------

gmark <- fread(file = "Markers_imp.csv", data.table = FALSE)
rownames(gmark) <- gmark[,1]
gmark <- gmark[,-1]

gmap <- read.csv(file = "Map.csv", dec = ",")

# prefilter markers ------------------------------------------------------------

mafpre <- 0.02 # ~20 homoz. individuals in the whole population
#mafpre <- 0.01 # ~10 homoz. individuals in the whole population
gmap <- gmap[gmap$maf > mafpre,]
gmark <- gmark[,colnames(gmark) %in% unique(gmap$Loc)]

mafsel <- 0.01
gmark <- gmark[rownames(gmark) %in% blups$GEN,]
# ~ 5 homoz. individuals in the whole (1000 sub) population 
# ~ 2 homoz. individuals in a 500er subpopulation

# prepare resutls dataframes ---------------------------------------------------
GWres  <- data.frame()

results <- rbind(cbind(results, "gwmod" = "MLMsom"))
results$yt <- paste(results$year,results$trait, sep = "-")
results <- merge(results, expand.grid("yt" = results$yt, "PCs" = 0:10), by = "yt")
results$gwcomb <- paste(results$gwmod, "-PCs:",results$PCs, sep = "")

results <- data.frame(results[order(results$year, results$trait, 
																		as.numeric(results$PCs)),])

rownames(results) <- 1:nrow(results)

i <-39

# GWAS loop --------------------------------------------------------------------
for(i in 1:nrow(results)){
  tr <- results[i,]    
  
  cat("\n", paste(c(i, nrow(results),
                    tr[c("trait","year","gwcomb")],
                    sep = " - ")))
  
  pheno <- droplevels(blups[blups$trait %in% tr$trait & 
                              blups$year %in% tr$year &
                              blups$GEN %in% rownames(gmark) & # genotyped
                              !is.na(blups$value),]) # NA values
  
  pheno <- pheno[order(pheno$GEN),]
  pheno$GEN <- as.factor(pheno$GEN)
  pheno$value <- as.numeric(pheno$value)

  
  if (nrow(pheno) > 0) {
    
    # MLM -----------------------------------
    
    gm <- as.matrix(gmark[rownames(gmark) %in% as.character(pheno$GEN),])
    
    Am <- sommer::A.mat(gm, min.MAF=mafsel, n.core = no_cores)
    
    # fit the model using sommer
    sink("/dev/null") 
    mod <- NA
    
    try(expr =
          mod <- sommer::GWAS(data=pheno, fixed = value ~ 1, rcov=~units,
                              n.core =4,#   parallel::detectCores()-1,
                              verbose = TRUE, 
                              random= ~ vs(GEN, Gu = Am), M = gm,
                              gTerm = "u:GEN", min.MAF = mafsel, 
                              n.PC = tr$PCs, P3D =TRUE)
        , silent = TRUE)
    sink()

    if(length(mod)>1){
      mod1 <- merge(x = gmap[,1:3], y = t(mod$scores),
                    by.x = "Loc", by.y = 0, all.x = FALSE, all.y = TRUE)
      colnames(mod1) <- c("Loc","chr","pos","beta","logpval","Fstat","R2","R2s")
      mod1$af <- apply(gm[,mod1$Loc],2,af)
      mod1$maf <- ifelse(mod1$af < .5, mod1$af, 1-mod1$af)
      
      # Bonferroni corrected significant values
      logpcrit <- signif(-log10(0.05/nrow(mod1)),digits = 3)
      mod1$sigBon <- mod1$logpval>logpcrit
      
      # FDR corrected p values
      mod1$pfdr <- p.adjust(p = 10^-mod1$logpval, method = "fdr")
      mod1$sigFDR <- mod1$pfdr<0.05

      mod1[,c("logpval","Fstat","R2","R2s","beta","af","maf","pfdr")] <- 
        apply(mod1[,c("logpval","Fstat","R2","R2s","beta","af","maf","pfdr")], 2,
              function(x) signif(as.numeric(x),4))
      
      res <- data.frame(tr[c("trait","year",
                             "gwcomb","gwmod","PCs")],
                        "lines" = length(unique(pheno$GEN)), # genotypes with records
                        "markers" = ncol(mod$scores),
                        "logpcrit" = logpcrit,
                        "nQTL" = sum(mod1$sigBon, na.rm=TRUE),
                        "nFDR" = sum(mod1$sigFDR, na.rm=TRUE),
                        "lambda" =  median(qchisq(1-10^-mod1$logpval,1))/qchisq(0.5,1), #genomic inflation
                        mod1, row.names = NULL)
      GWres <- rbind.fill(GWres,res)
    }
  }
}

sink() 

remove(i, mod, pheno, gm, Am, tr, mafsel, mafpre)

GWp <- unique(GWres[,c("trait","year","gwcomb",
                       "Loc","beta","logpval","R2","R2s","maf",
                       "sigBon","pfdr","sigFDR")])

GWstat <- unique(GWres[,c("trait","year","gwcomb",
                          "gwmod","PCs","nQTL","nFDR","lambda",
                          "lines","markers","logpcrit")])

write.csv(x = GWp, file = "PHS-trial-GWAS_p.csv")
write.csv(x = GWstat, file = "PHS-trial-GWAS_stats.csv")

# END ########################################################################
