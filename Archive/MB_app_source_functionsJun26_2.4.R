print(getwd())

#######################Packages############ 

# _____        _____ _  __          _____ ______  _____ 
#|  __ \ /\   / ____| |/ /    /\   / ____|  ____|/ ____|
#| |__) /  \ | |    | ' /    /  \ | |  __| |__  | (___  
#|  ___/ /\ \| |    |  <    / /\ \| | |_ |  __|  \___ \ 
#| |  / ____ \ |____| . \  / ____ \ |__| | |____ ____) |
#|_| /_/    \_\_____|_|\_\/_/    \_\_____|______|_____/ 
                                                        
                                                        

library(BiocManager)
options(repos = BiocManager::repositories())
# if(!require("devtools"))
#   install.packages("devtools")
# library(rsconnect)
# if(!require(devtools)) {
#   BiocManager::install("devtools")
#   library(devtools)
# }
# 
# if(!require(rsconnect)) {
#   BiocManager::install("rsconnect")
#   library(rsconnect)
# }


if(!require(minfiData)) {
  BiocManager::install("minfiData")
  library(minfiData)
}
message("minfidata done")
if(!require(sva)){
  BiocManager::install("sva")
  library(sva)
}
message("sva done")
if(!require(bumphunter)){
  BiocManager::install("bumphunter")
  library(bumphunter)
}
message("bumphunter done")
# if(!require(lumi)){
#   BiocManager::install("lumi")
#   library(lumi)
# }
# message("lumi done")
# 
# if(!require(minfi)){
#   BiocManager::install("minfi")
#   library(minfi)
# }
# 
# message("minfi done")
 if(!require(IlluminaHumanMethylationEPICmanifest)){
   BiocManager::install("IlluminaHumanMethylationEPICmanifest")
   library(IlluminaHumanMethylationEPICmanifest)
 }
# message("illumina1 done")
 if(!require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)){
   BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
   library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
 }
#message("manifest loaded")
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
message("shiny done")

if(!require(DT)){
  install.packages("DT")
  library(DT)
}

message("DT done")

if(!require(shinyWidgets)){
  install.packages("shinyWidgets")
  library(shinyWidgets)
}
message("shinyWidgets done")


if(!require(shinydashboard)){
  install.packages("shinydashboard")
  library(shinydashboard)
}
message("shinydashboard done")


if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
message("ggplot2 done")


if(!require(ggpubr)){
  install.packages("ggpubr")
  library(ggpubr)
}
message("ggpubr done")

if(!require(ggnewscale)){
  install.packages("ggnewscale")
  library(ggnewscale)
}
message("ggnewscale done")

if(!require(foreach)){
  install.packages("foreach")
  library(foreach)
}
message("foreach done")

if(!require(ggrepel)){
  install.packages("ggrepel")
  library(ggrepel)
}
message("ggrepel done")

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}
message("pheatmap done")

if(!require(shinycssloaders)){
  install.packages("shinycssloaders")
  library(shinycssloaders)
}
message("shinycssloaders done")

if(!require(shinybusy)){
  install.packages("shinybusy")
  library(shinybusy)
}
message("shinybusy done")

if(!require(waiter)){
  install.packages("waiter")
  library(waiter)
}
message("waiter done")

if(!require(tinytex)){
  install.packages('tinytex')
  library(tinytex)
}

if (!require(survival)) {
  install.packages('survival')
  library(survival)
}

# In your packages file
if (!require(rhandsontable)) {
  install.packages("rhandsontable")
  library(rhandsontable)
}


library(gridExtra)
library(grid)

####################################

message("packages loaded")
 readRDS(file = "./mb_app/G3_G4_subgroup.model.rds") -> G3_G4_sub
readRDS(file = "./mb_app/G3_G4_no_subgroup.model.rds") -> G3_G4_no_sub
 readRDS(file = "./mb_app/SHHn49.model.rds") -> SHH

 # ── Load pre-computed survival objects ───────────────────────────────────────
 
 # SHH
 .shh_df  <- readRDS("~/MB_Risk_App/mb_app/precomp_shh_df.rds")
 .shh_fit <- readRDS("~/MB_Risk_App/mb_app/precomp_shh_fit.rds")
 
 .shh_env <- new.env(parent = baseenv())
 .shh_env$surv.object.train.shh <- Surv(.shh_df$OS_Time, .shh_df$OS_Status)
 .shh_env$fixedMG               <- .shh_df$fixedMG
 .shh_env$ConsensusMYCN         <- .shh_df$ConsensusMYCN
 environment(.shh_fit$terms)    <- .shh_env
 
 # G3G4 late
 .g34late_df  <- readRDS("~/MB_Risk_App/mb_app/precomp_g34late_df.rds")
 .g34late_fit <- readRDS("~/MB_Risk_App/mb_app/precomp_g34late_fit.rds")
 
 .g34late_env <- new.env(parent = baseenv())
 .g34late_env$surv.object.train.late <- Surv(.g34late_df$OS_Time, .g34late_df$OS_Status)
 .g34late_env$fixedMG                <- .g34late_df$fixedMG
 environment(.g34late_fit$terms)     <- .g34late_env
 
 # G3G4 early
 .g34early_df  <- readRDS("~/MB_Risk_App/mb_app/precomp_g34early_df.rds")
 .g34early_fit <- readRDS("~/MB_Risk_App/mb_app/precomp_g34early_fit.rds")
 
 .g34early_df$ConsensusMYC <- factor(
   .g34early_df$ConsensusMYC,
   levels = .g34early_fit$xlevels[["strata(ConsensusMYC)"]]
 )
 
 .g34early_env <- new.env(parent = baseenv())
 .g34early_env$surv.object.train <- Surv(.g34early_df$OS_Time, .g34early_df$OS_Status)
 .g34early_env$fixedMG           <- .g34early_df$fixedMG
 .g34early_env$M._versus_M.      <- .g34early_df$M._versus_M.
 .g34early_env$ConsensusMYC      <- .g34early_df$ConsensusMYC
 .g34early_env$strata            <- survival::strata
 environment(.g34early_fit$terms) <- .g34early_env
 
 message("Survival objects loaded.")


 
 surv.object.train               <- Surv(.g34early_df$OS_Time, .g34early_df$OS_Status)
 df.pheno.mb.combined.grp3.grp4 <- .g34early_df
 
 env <- environment(.g34early_fit$terms)
 env$surv.object.train <- surv.object.train
 env$fixedMG           <- .g34early_df$fixedMG
 env$M._versus_M.      <- .g34early_df$M._versus_M.
 env$ConsensusMYC      <- .g34early_df$ConsensusMYC
 
 message("Survival objects loaded.")
 
 message("Survival objects loaded.")

beta2m <- function (beta) {
  m <- log2(beta/(1 - beta))
  return(m)
}

process_idats <- function(basenames){
  

  # Split into mixed or single test e.g. 450k & Epic or one or the other
  file.sizes <- file.info(paste0(basenames, "_Red.idat"))$size
  
  length(which(file.sizes> 1e7))>0 & length(which(file.sizes< 1e7))>0 -> mixed.test
  length(which(file.sizes> 1e7))>0 | length(which(file.sizes< 1e7))>0 -> single.test
  
  if(mixed.test) {
    #  as before where you process them seperately then combine
    idx.epic <-
      file.info(paste0(basenames, "_Red.idat"))$size > 1e7
    m450k <- read.metharray(basenames[!idx.epic],
                            verbose = TRUE,
                            force = TRUE) 
    mEpic <- read.metharray(basenames[idx.epic],
                            verbose = TRUE,
                            force = TRUE) 
    comb.rgSet <- combineArrays(m450k, 
                                mEpic,
                                outType = "IlluminaHumanMethylation450k")
    comb.mSet <- preprocessNoob(comb.rgSet, 
                                dyeMethod = "single")
    comb.detP <- detectionP(comb.rgSet)
    comb.gmrSet <- mapToGenome(ratioConvert(comb.mSet))
    comb.beta <- getBeta(comb.gmrSet)
    comb.MRT.pd <- pData(comb.gmrSet)
  }
  
  if((!mixed.test)&single.test){
    # only need the one
    single <- read.metharray(basenames,
                             verbose = TRUE,
                             force = TRUE)
    comb.mSet <- preprocessNoob(single, dyeMethod = "single")
    comb.detP <- detectionP(single)
    comb.gmrSet <- mapToGenome(ratioConvert(comb.mSet))
    comb.beta <- getBeta(comb.gmrSet)
    comb.MRT.pd <- pData(comb.gmrSet)
  }
  
  return(list(betas = comb.beta, pvals = comb.detP))
}

get_basenames <- function(dir){
  list.files(path = dir, pattern = "_Red.idat", full.names = T) -> temp.files
  return(gsub("_Red.idat","", temp.files))
}

extract.metagene <- function(index, weights, exp.matrix, scaling) {
  exp.matrix[index, , drop = FALSE] -> temp.exp
  apply(temp.exp, 2, function(x) {
    mean(x * weights)
  }) -> raw.metagene
  round(raw.metagene, digits = 3)
  as.numeric(scale(raw.metagene, center = scaling[1], scale = scaling[2])) -> Risk_Value
  round(Risk_Value, digits = 3)
  names(Risk_Value) <- colnames(exp.matrix)
  df <- as.data.frame(Risk_Value)
  return(df)
}

digits = 0:9
createRandString<- function() {
  v = c(sample(LETTERS, 5, replace = TRUE),
        sample(digits, 4, replace = TRUE),
        sample(LETTERS, 1, replace = TRUE))
  return(paste0(v,collapse = ""))
}

add_suffix <- function(number) {
  if (number %% 10 == 1 & number %% 100 != 11) {
    suffix <- "st"
  } else if (number %% 10 == 2 & number %% 100 != 12) {
    suffix <- "nd"
  } else if (number %% 10 == 3 & number %% 100 != 13) {
    suffix <- "rd"
  } else {
    suffix <- "th"
  }
  return(paste0(number, suffix))
}

render_survival_plot <- function(metagene, score, covs, indexRow) {
  if (metagene == "SHH") {
    generate_survival_figure_shh(score, covs$mycn, indexRow)
  } else if (metagene == "Group3/4 (Early)") {
    generate_survival_figure_G3_G4_sub(score, covs$myc, covs$mets, indexRow)
  } else {
    generate_survival_figure_G3_G4_no_sub(score, indexRow)
  }
}




make_covariate_df <- function(sample_names, metagene) {
  n <- length(sample_names)
  if (metagene == "SHH") {
    data.frame(
      sample = sample_names,
      mycn   = rep(0L, n),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else if (metagene == "Group3/4 (Early)") {
    data.frame(
      sample = sample_names,
      myc    = rep(0L, n),
      mets   = rep(0L, n),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else {
    # G3_G4_no_sub has no covariates
    data.frame(
      sample = sample_names,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
}



##############################################
 #   _____ _____ _    _  
 #  / ____/ ____| |  | | 
 # | (___| (___ | |__| | 
 #  \___ \\___ \|  __  | 
 #  ____) |___) | |  | | 
 # |_____/_____/|_|  |_| 
 #  
  
  

generate_figure_highlight_SHH <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}
   temp.df <- readRDS(file = "./mb_app/SHHn49.dist.rds")
   

  ggplot(aes(x = 1:nrow(temp.df), y = SHH), data = temp.df) +
    geom_line() +
    scale_shape_manual(values = c(1, 4,  3)) +
    scale_color_manual(values = c('#E69F00', '#999999', "white")) +
    ylab("SHH") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> b
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  temp.df$SHH < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  temp.df$SHH < new.sample.meta.score[i]
                )),
                xend = max(which(
                  temp.df$SHH < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(temp.df$SHH)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(temp.df$SHH) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor) 
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}




generate_survival_figure_shh <- function(new.sample.meta.score, mycn, indexRow) {
  message("generate_survival_figure_shh called - new version")
  if (is.null(indexRow)) indexRow <- 1
  
  df2 <- data.frame(surv = .shh_df$surv5)
  
  fit <- survfit(
    .shh_fit,
    newdata = data.frame(
      fixedMG       = as.numeric(new.sample.meta.score),
      ConsensusMYCN = as.numeric(mycn)
    )
  )
  new_surv <- as.numeric(summary(fit, time = 5)$surv)
  
  if (length(new_surv) == 0 || all(is.na(new_surv))) {
    message("Warning: no survival estimate at time=5 for SHH")
    new_surv <- rep(NA_real_, length(new.sample.meta.score))
  }
  
  df.lines.ver <- data.frame(
    x      = new_surv,
    colour = factor(
      ifelse(seq_along(new_surv) == indexRow, "highlight", "no.highlight"),
      levels = c("highlight", "no.highlight")
    )
  )
  
  ggplot(df2, aes(x = surv)) +
    geom_density(fill = "steelblue", color = "black", alpha = 0.5) +
    geom_vline(aes(xintercept = x, colour = colour),
               linetype = "dashed", linewidth = 0.5,
               data = df.lines.ver) +
    scale_colour_manual(values = c(highlight = "orange", no.highlight = "grey60")) +
    theme_classic() +
    xlab("Estimated 5-year Survival") +
    xlim(0, 1) +
    theme(legend.position = "none", text = element_text(size = 15))
}



###################################################
#   _______  _______          _______  _   ___          __    _  _______         _______  __   __  _______ 
#   |       ||       |        |       || | |   |        |  |  | ||       |       |       ||  | |  ||  _    |
#   |    ___||___    |        |    ___|| |_|   |        |   |_| ||   _   | ____  |  _____||  | |  || |_|   |
#   |   | __  ___|   |        |   | __ |       |        |       ||  | |  ||____| | |_____ |  |_|  ||       |
#   |   ||  ||___    |        |   ||  ||___    |        |  _    ||  |_|  |       |_____  ||       ||  _   | 
#   |   |_| | ___|   | _____  |   |_| |    |   | _____  | | |   ||       |        _____| ||       || |_|   |
#   |_______||_______||_____| |_______|    |___||_____| |_|  |__||_______|       |_______||_______||_______|
#   

  
generate_figure_highlight_G3_G4_no_sub <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}

  temp.df <- readRDS(file = "./mb_app/G3_G4_no_subgroup.dist.rds")

  ggplot(aes(x = 1:nrow(temp.df), y = G3_G4_no_sub), data = temp.df) +
    geom_line() +
    scale_shape_manual(values = c(1, 4,  3)) +
    scale_color_manual(values = c('#E69F00', '#999999', "white")) +
    ylab("Group3/4 (Late)") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> b
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  temp.df$G3_G4_no_sub < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  temp.df$G3_G4_no_sub < new.sample.meta.score[i]
                )),
                xend = max(which(
                  temp.df$G3_G4_no_sub < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(temp.df$G3_G4_no_sub)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(temp.df$G3_G4_no_sub) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor) 
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}




generate_survival_figure_G3_G4_no_sub <- function(new.sample.meta.score, indexRow) {
  if (is.null(indexRow)) indexRow <- 1
  
  df2 <- data.frame(surv = .g34late_df$surv10)
  
  # Extract surv explicitly as numeric vector
  fit <- survfit(.g34late_fit,
                 newdata = data.frame(
                   fixedMG = as.numeric(new.sample.meta.score)
                 ))
  
  surv_summary <- summary(fit, time = 10)
  new_surv     <- as.numeric(surv_summary$surv)
  
  # Guard against NA or failed extraction
  if (length(new_surv) == 0 || all(is.na(new_surv))) {
    message("Warning: no survival estimate at time=10, check training data range")
    new_surv <- rep(NA_real_, length(new.sample.meta.score))
  }
  
  df.lines.ver <- data.frame(
    x      = new_surv,
    colour = factor(
      ifelse(seq_along(new_surv) == indexRow, "highlight", "no.highlight"),
      levels = c("highlight", "no.highlight")
    )
  )
  
  ggplot(df2, aes(x = surv)) +
    geom_density(fill = "steelblue", color = "black", alpha = 0.5) +
    geom_vline(aes(xintercept = x, colour = colour),
               linetype = "dashed", linewidth = 0.5,
               data = df.lines.ver) +
    scale_colour_manual(values = c(highlight = "orange", no.highlight = "grey60")) +
    theme_classic() +
    xlab("Estimated 10-year Survival") +
    xlim(0, 1) +
    theme(legend.position = "none", text = element_text(size = 15))
}


###################################################
#    _______  _______          _______  _   ___         _______  __   __  _______ 
#   |       ||       |        |       || | |   |       |       ||  | |  ||  _    |
#   |    ___||___    |        |    ___|| |_|   | ____  |  _____||  | |  || |_|   |
#   |   | __  ___|   |        |   | __ |       ||____| | |_____ |  |_|  ||       |
#   |   ||  ||___    |        |   ||  ||___    |       |_____  ||       ||  _   | 
#   |   |_| | ___|   | _____  |   |_| |    |   |        _____| ||       || |_|   |
#   |_______||_______||_____| |_______|    |___|       |_______||_______||_______|

generate_figure_highlight_G3_G4_sub <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}

  temp.df <- readRDS(file = "./mb_app/G3_G4_subgroup.dist.rds")

  ggplot(aes(x = 1:nrow(temp.df), y = G3_G4_sub), data = temp.df) +
    geom_line() +
    scale_shape_manual(values = c(1, 4,  3)) +
    scale_color_manual(values = c('#E69F00', '#999999', "white")) +
    ylab("Group3/4 (Early)") +
    theme_minimal() +
    theme(
      #axis.title.x = element_blank(),
      #axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> b
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  temp.df$G3_G4_sub < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  temp.df$G3_G4_sub < new.sample.meta.score[i]
                )),
                xend = max(which(
                  temp.df$G3_G4_sub < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(temp.df$G3_G4_sub)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(temp.df$G3_G4_sub) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor) 
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}

generate_survival_figure_G3_G4_sub <- function(new.sample.meta.score, myc, mets_val, indexRow) {
  message("=== G3_G4_sub plot called ===")
  message("myc values: ", paste(myc, collapse=", "))
  message("mets values: ", paste(mets_val, collapse=", "))
  message("indexRow: ", indexRow)
  if (is.null(indexRow)) indexRow <- 1
  
  df2 <- data.frame(surv = .g34early_df$surv5)
  
  # Define new_rows here before using it
  new_rows <- data.frame(
    fixedMG      = as.numeric(new.sample.meta.score),
    M._versus_M. = as.numeric(mets_val),
    ConsensusMYC = factor(
      paste0("ConsensusMYC=", as.numeric(myc)),
      levels = .g34early_fit$xlevels[["strata(ConsensusMYC)"]]
    )
  )
  
  fit      <- survfit(.g34early_fit, newdata = new_rows)
  new_surv <- as.numeric(summary(fit, time = 5)$surv)
  
  if (length(new_surv) == 0 || all(is.na(new_surv))) {
    message("Warning: no survival estimate at time=5 for Group3/4 (Early)")
    new_surv <- rep(NA_real_, length(new.sample.meta.score))
  }
  
  df.lines.ver <- data.frame(
    x      = new_surv,
    colour = factor(
      ifelse(seq_along(new_surv) == indexRow, "highlight", "no.highlight"),
      levels = c("highlight", "no.highlight")
    )
  )
  
  ggplot(df2, aes(x = surv)) +
    geom_density(fill = "steelblue", color = "black", alpha = 0.5) +
    geom_vline(aes(xintercept = x, colour = colour),
               linetype = "dashed", linewidth = 0.5,
               data = df.lines.ver) +
    scale_colour_manual(values = c(highlight = "orange", no.highlight = "grey60")) +
    theme_classic() +
    xlab("Estimated 5-year Survival") +
    xlim(0, 1) +
    theme(legend.position = "none", text = element_text(size = 15))
}
