#### Script for performing Diff. Methylation Analysis using bayesian binomial mixed model
# implemented in the package brms

### Libraries ####
library(Rcpp)
library(rstan)
rstan_options(auto_write = TRUE) #avoid recompiling programs in rstan
library(brms)

### USER VALUES ###
#setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/")
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
#if (length(args)==0) {
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#}

# Input should be the path to the methylation data
#INPUT_M <- "/home/downeyam/Github/2017OAExp_Oysters/input_files/DNAm"
#INPUT_META <- "/home/downeyam/Github/2017OAExp_Oysters/input_files/metadata"
#OUTPUT <- "/home/downeyam/Github/2017OAExp_Oysters/results/DNAM_diffMethylation"
MINp <-  args[1]
#MINp <- 1
MAXp <-  args[2]
#MAXp <- 2
#S_INPUT <-  args[3]
#S_OUTPUT <-  args[4]
S_INPUT <- "processed_samples/05_countSummary"
S_OUTPUT <- "processed_samples/07_brmsSummary"

#sprintf(MINp)
#sprintf(MAXp)

# Hyphothesis testing parameters
alph = 0.0001 # Significance threshold for hypothesis testing
BC_COR <- 3 # Shrinks alph for basic comparisons based on number of comparisons
PC_COR <- 3 # Shrinks alph for planned comparisons based on number of comparisons

### Data ####
mC <- readRDS(paste0(S_INPUT,"/CG_unstranded_mC_geneOnly_5.RData"))
mC <- mC[,colnames(mC) != "17099"]
uC <- readRDS(paste0(S_INPUT,"/CG_unstranded_umC_geneOnly_5.RData"))
uC <- uC[,colnames(uC) != "17099"]
meta_locus <- readRDS(paste0(S_INPUT,"/CG_unstranded_summaryTable_geneOnly_5.RData"))
meta_samp <- readRDS("metadata/metadata_20190811.RData")
meta_samp$tank <- as.factor(meta_samp$tank)
meta_samp$shelf <- as.factor(meta_samp$shelf)
levels(meta_samp$SFV) <- c("1C","2C","1E","2E")
meta_samp <- meta_samp[meta_samp$ID != "17099",]

# unit test - these should both be true
identical(colnames(mC), colnames(uC))
identical(colnames(uC), meta_samp$ID)

#### Performing logistic regression (brms with binomial family) ####
## Three main explanatory variables  time and treatment + interaction
## plus tank as random effect (1|tank:shelf)

## Testing single loci
#loci_range <- 38337 # Good example of very diff. methylation for trt at day 09
#loci_range <- 1:10 # Test for starting at the first locus
#target <- 1

loci_range<-c(MINp:MAXp)

sprintf("Starting model...")
for(target in loci_range){
  #tryCatch({
    sprintf(paste0("Loci ",target," of ",nrow(mC)))
    Sys.sleep(0.01)
    temp <- data.frame(mC=as.numeric(unlist(mC[target,])),
                       size=as.numeric(c(unlist(mC[target,])+unlist(uC[target,]))),
                       levels=meta_samp$SFV,
                       shelf=meta_samp$shelf,
                       tank=meta_samp$tank,
                       time=meta_samp$Time,
                       trt=meta_samp$Treatment,
                       tankID=meta_samp$tankID)
    # ### Model ###
    m.out2 <- brm(mC | trials(size)~ time*trt + (1|tankID), data=temp,
                  family = binomial(), iter=10000, control = list(max_treedepth = 15, adapt_delta=0.99))
    #Model summary
    m.sum2 <- summary(m.out2)
    me <- marginal_effects(m.out2,robust=TRUE)
    me2 <- marginal_effects(m.out2,effects = "trt:time",plot=FALSE,ask=FALSE)

    #Save diagnostic plots
    png(paste0(S_OUTPUT,"/img/CpG_",meta_locus$ID[target],"_brmsEffectsPlot.png"))
    plot(m.out2)
    dev.off()
    # Save marginal effects
    png(paste0(S_OUTPUT,"/img/CpG_",meta_locus$ID[target],"_brmsMarginalEffectsPlot.png"))
    plot(me2)
    Sys.sleep(0.05)
    dev.off()

    ### Basic and planned contrasts ###
    ## Script for planned comparisons when using the proper time*treatment factor formula in glmer

    # Not simplified (with Intercept)
    # Basic Comparisons
    tp2 <- c("(Intercept + Intercept + time80)/2 > (Intercept + trt2800 + Intercept + time80 + trt2800 + time80:trt2800)/2", # Treatment C > E
             "(Intercept + Intercept + time80)/2 < (Intercept + trt2800 + Intercept + time80 + trt2800 + time80:trt2800)/2", # Treatment C < E
             "(Intercept + Intercept + trt2800)/2 > (Intercept + time80 + Intercept + time80 + trt2800 + time80:trt2800)/2", # Time 9 > 80
             "(Intercept + Intercept + trt2800)/2 < (Intercept + time80 + Intercept + time80 + trt2800 + time80:trt2800)/2", # Time 9 < 80
             "(time80:trt2800)=0") # Interaction
    h.out2 <-hypothesis(m.out2,tp2, class="b", alpha=alph/BC_COR)
    # Planned comparisons
    tp3 <- c(
      # All comparisons with Trt C Day 9
      "(Intercept) > (Intercept + trt2800)", # Trt C Day 09 > Trt E day 9
      "(Intercept) < (Intercept + trt2800)", # Trt C Day 09< Trt E day 9
      "(Intercept) > (Intercept + time80)", # Trt C Day 09 > Trt C day 80
      "(Intercept) < (Intercept + time80)", # Trt C Day 09 < Trt C day 80
      "(Intercept) > (Intercept + time80 + trt2800 + time80:trt2800)", # Trt C day 9 > Trt E day 80
      "(Intercept) < (Intercept + time80 + trt2800 + time80:trt2800)", # Trt C day 9 < Trt E day 80
      # All remaining comparisons with Trt E Day 09
      "(Intercept + trt2800) > (Intercept + time80)", # Trt E day 09 > Trt C day 80
      "(Intercept + trt2800) < (Intercept + time80)", # Trt E day 09 < Trt C day day 80
      "(Intercept + trt2800) > (Intercept + time80 + trt2800 + time80:trt2800)", # Trt E day 9 > Trt E day 80
      "(Intercept + trt2800) < (Intercept + time80 + trt2800 + time80:trt2800)", # Trt E day 9 < Trt E day 80
      # The remaining Day 80 comparison
      "(Intercept + time80) > (Intercept + time80 + trt2800 + time80:trt2800)", # Trt C Day 80 > E day 80
      "(Intercept + time80) < (Intercept + time80 + trt2800 + time80:trt2800)") # Trt C Day 80 < E day 80
    
    h.out3 <-hypothesis(m.out2,tp3, class="b", alpha=alph/PC_COR)
    # 
    # # Creating list for storage on first loop
    if(target == MINp){
      #Creates this once, provides some summary info about the model we use for all CpGs
      model_param <- list(formula=m.out2$formula,
                          family=m.out2$family,
                          model=m.out2$model,
                          prior=m.out2$prior,
                          random_effect=m.out2$ranef,
                          version=m.out2$version)
      #Creates matrices in list for the model stats
      summ_lab <- c("Intercept","time80","trt2800","time80:trt2800")
      model_summary_stats <- list(Estimate <- matrix(0,ncol=length(summ_lab),nrow=nrow(mC)),
      Est.Error = matrix(0,ncol=length(summ_lab),nrow=nrow(mC)),
      CI.Lower  = matrix(0,ncol=length(summ_lab),nrow=nrow(mC)),
      CI.Upper  = matrix(0,ncol=length(summ_lab),nrow=nrow(mC)),
      Rhat      = matrix(0,ncol=length(summ_lab),nrow=nrow(mC)),
      Bulk_ESS  = matrix(0,ncol=length(summ_lab),nrow=nrow(mC)),
      Tail_ESS  = matrix(0,ncol=length(summ_lab),nrow=nrow(mC)))
      
      # Add column names
      for(j in 1:length(model_summary_stats)){
        colnames(model_summary_stats[[j]]) <- summ_lab
      }
      
      # Create list of matrices for marginal effects info
      ME_lab <-  c("Time_09","Time_80","Trt_400","Trt_2800",
                   "Time:Trt_09_400","Time:Trt_09_2800",
                   "Time:Trt_80_400","Time:Trt_80_2800")
      m_effects <- list(
        marg_estimate=matrix(0,ncol=length(ME_lab),nrow=nrow(mC)),
        se_estimate=matrix(0,ncol=length(ME_lab),nrow=nrow(mC)),
        l_estimate=matrix(0,ncol=length(ME_lab),nrow=nrow(mC)),
        u_estimate=matrix(0,ncol=length(ME_lab),nrow=nrow(mC)))
      
      for(j in 1:length(m_effects)){
        colnames(m_effects[[j]]) <- ME_lab
      }
      # Create matrices for planned comparisons
      # Given the time to run, I decided to make this a single comparisons table, so that we could easily
      # access both the initial basic comparisons and planned comparisons in a single row.
      # This means we actually perform the hypothesis test and store all planned comparisons, but we can
      # easily subset those using the $significant for 'interaction' column.
      Hyp_lab <-c(
        "Trt_400>2800","Trt_400<2800",
        "Time_09>80","Time_09<80",
        "Trt:Time_Interaction",
        "C.09>E.09","C.09<E.09",
        "C.09>C.80","C.09<C.80",
        "C.09>E.80","C.09<E.80",
        "E.09>C.80","E.09<C.80",
        "E.09>E.80","E.09<E.80",
        "C.80>E.80","C.09<E.80")
      basic_comparisons <- list(Estimate = matrix(0,ncol=length(Hyp_lab),nrow=nrow(mC)),
                                Est.Error = matrix(0,ncol=length(Hyp_lab),nrow=nrow(mC)),
                                CI.Lower = matrix(0,ncol=length(Hyp_lab),nrow=nrow(mC)),
                                CI.Upper = matrix(0,ncol=length(Hyp_lab),nrow=nrow(mC)),
                                Evid.Ratio = matrix(0,ncol=length(Hyp_lab),nrow=nrow(mC)),
                                Post.Prob = matrix(0,ncol=length(Hyp_lab),nrow=nrow(mC)),
                                Significant = matrix("",ncol=length(Hyp_lab),nrow=nrow(mC)))
      # Add column names
      for(j in 1:length(basic_comparisons)){
        colnames(basic_comparisons[[j]]) <- Hyp_lab
      }
    }
    
    # Saving basic model stats
    for(j in 1:length(model_summary_stats)){
      model_summary_stats[[j]][target,] <- c(m.sum2$fixed[,(j)])
    }
    # Saving marginal effects
    m_effects$marg_estimate[target,] <-  c(me$time$estimate__,
                                           me$trt$estimate__,
                                           me$`time:trt`$estimate__)
    m_effects$se_estimate[target,] <-  c(me$time$se__,
                                         me$trt$se__,
                                         me$`time:trt`$se__)
    m_effects$l_estimate[target,] <-  c(me$time$lower__,
                                        me$trt$lower__,
                                        me$`time:trt`$lower__)
    m_effects$u_estimate[target,] <-  c(me$time$upper__,
                                        me$trt$upper__)
    # Saving planned comparisons
    #Changes Infinite values to 10000 (number of posterior samples)
    h.out2$hypothesis$Evid.Ratio[is.infinite(h.out2$hypothesis$Evid.Ratio)] <- 10000
    h.out3$hypothesis$Evid.Ratio[is.infinite(h.out3$hypothesis$Evid.Ratio)] <- 10000
    #Adding values to planned comparison list
    for(j in 1:length(basic_comparisons)){
      basic_comparisons[[j]][target,] <- c(h.out2$hypothesis[,(j+1)],h.out3$hypothesis[,j+1])
    }
    
    saveRDS(model_param,paste0(S_OUTPUT,"/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelParam.RData"))
    saveRDS(model_summary_stats,paste0(S_OUTPUT,"/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelSummary.RData"))
    saveRDS(m_effects,paste0(S_OUTPUT,"/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelMarginalEffects.RData"))
    saveRDS(basic_comparisons,paste0(S_OUTPUT,"/DNAm_gene_BRMS_",MINp,"_",MAXp,"_plannedComparisons.RData"))
#  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

sprintf("Saving outputs...")
saveRDS(model_param,paste0(S_OUTPUT,"/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelParam.RData"))
saveRDS(model_summary_stats,paste0(S_OUTPUT,"/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelSummary.RData"))
saveRDS(m_effects,paste0(S_OUTPUT,"/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelMarginalEffects.RData"))
saveRDS(basic_comparisons,paste0(S_OUTPUT,"/DNAm_gene_BRMS_",MINp,"_",MAXp,"_plannedComparisons.RData"))

# Alt hypotheses
# m.out1 <- brm(mC | trials(size)~ levels + (1|tankID), data=temp,
#             family = binomial(),iter=5000, control = list(max_treedepth = 15, adapt_delta=0.99))
# (m.sum1 <- summary(m.out1))
# tp1 <- c("(levels2C)/2 > (levels1E + levels2E)/2", # Treatment C > E
#         "(levels2C)/2 < (levels1E + levels2E)/2", # Treatment C < E
#         "(levels1E)/2 > (levels2C + levels2E)/2", # Time 9 > 80
#         "(levels1E)/2 < (levels2C + levels2E)/2", # Time 9 < 80
#         "(-levels1E) = (levels2C - levels2E)") # Interaction (this still doesn't seem correct)
# h.out1 <-hypothesis(m.out1,tp1, class="b", alpha=alph/5)
#
# tp4 <- c(
#   # All comparisons with Trt C Day 9
#   "(Intercept) > (Intercept + levels1E)", # Trt C Day 09 > Trt E day 9
#   "(Intercept) < (Intercept + levels1E)", # Trt C Day 09< Trt E day 9
#   "(Intercept) > (Intercept + levels2C)", # Trt C Day 09 > Trt C day 80
#   "(Intercept) < (Intercept + levels2C)", # Trt C Day 09 < Trt C day 80
#   "(Intercept) > (Intercept + levels2E)", # Trt C day 9 > Trt E day 80
#   "(Intercept) < (Intercept + levels2E)", # Trt C day 9 < Trt E day 80
#   # All remaining comparisons with Trt E Day 09
#   "(Intercept + levels1E) > (Intercept + levels2C)", # Trt E day 09 > Trt C day 80
#   "(Intercept + levels1E) < (Intercept + levels2C)", # Trt E day 09 < Trt C day day 80
#   "(Intercept + levels1E) > (Intercept + levels2E)", # Trt E day 9 > Trt E day 80
#   "(Intercept + levels1E) < (Intercept + levels2E)", # Trt E day 9 < Trt E day 80
#   # The remaining Day 80 comparison
#   "(Intercept + levels2C) > (Intercept + levels2E)", # Trt C Day 80 > E day 80
#   "(Intercept + levels2C) < (Intercept + levels2E)") # Trt C Day 80 < E day 80
# h.out4 <-hypothesis(m.out1,tp4, class="b", alpha=alph/PC_COR)

