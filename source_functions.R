print(getwd())

#######################Packages############ 

# _____        _____ _  __          _____ ______  _____ 
#|  __ \ /\   / ____| |/ /    /\   / ____|  ____|/ ____|
#| |__) /  \ | |    | ' /    /  \ | |  __| |__  | (___  
#|  ___/ /\ \| |    |  <    / /\ \| | |_ |  __|  \___ \ 
#| |  / ____ \ |____| . \  / ____ \ |__| | |____ ____) |
#|_| /_/    \_\_____|_|\_\/_/    \_\_____|______|_____/ 
                                                        
                                                        
# --- Setup ---------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
options(repos = BiocManager::repositories())

# --- Helper functions ---------------------------------------------------

install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  message(pkg, " loaded")
}

bioc_install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  message(pkg, " loaded")
}

# --- Bioconductor packages ----------------------------------------------

bioc_pkgs <- c(
  "minfiData",
  "sva",
  "bumphunter",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)


# --- CRAN packages ------------------------------------------------------

cran_pkgs <- c(
  "shiny",
  "DT",
  "shinyWidgets",
  "shinydashboard",
  "ggplot2",
  "ggpubr",
  "foreach",
  "ggrepel",
  "waiter",
  "survival",
  "rhandsontable",
  "gridExtra",
  "grid"
)


invisible(lapply(bioc_pkgs, bioc_install_and_load))
invisible(lapply(cran_pkgs, install_and_load))


####################################

#source("./slim_mb.R")   # predict_survival / percentile_of (prediction-only payloads)

message("packages loaded")

# ── Classifiers (clean: CpG weights + scaling) ───────────────────────────────
G3_G4_sub    <- readRDS("./mb_app/slim/G3_G4_subgroup.model.rds")
G3_G4_no_sub <- readRDS("./mb_app/slim/G3_G4_no_subgroup.model.rds")
SHH          <- readRDS("./mb_app/slim/SHHn49.model.rds")

# ── Prediction-only survival payloads (no coxph fits, no per-patient frames) ──
.shh_slim      <- readRDS("./mb_app/slim/shh_slim.rds")
.g34early_slim <- readRDS("./mb_app/slim/g34early_slim.rds")
.g34late_slim  <- readRDS("./mb_app/slim/g34late_slim.rds")

message("Slim survival objects loaded.")

# ── One survival entry point for all three models ────────────────────────────
# Returns survival probability at the model's landmark (5 yr early/SHH, 10 yr late)
# for one or many samples. Mirrors exactly the newdata the original app built.
mb_surv_prob <- function(metagene, score, mycn = 0, myc = 0, mets = 0) {
  score <- as.numeric(score)
  if (metagene == "SHH") {
    nd <- data.frame(fixedMG = score, ConsensusMYCN = as.numeric(mycn))
    predict_survival(.shh_slim, nd, times = .shh_slim$landmark)$surv
  } else if (metagene == "Group3/4 (Early)") {
    nd <- data.frame(
      fixedMG      = score,
      M._versus_M. = as.numeric(mets),
      ConsensusMYC = factor(paste0("ConsensusMYC=", as.numeric(myc)),
                            levels = .g34early_slim$strata_levels))
    predict_survival(.g34early_slim, nd, times = .g34early_slim$landmark)$surv
  } else {  # Group3/4 (Late)
    nd <- data.frame(fixedMG = score)
    predict_survival(.g34late_slim, nd, times = .g34late_slim$landmark)$surv
  }
}

# predict_survival takes one row, so map over samples for vectors.
.mb_surv_vec <- function(metagene, score, mycn = NULL, myc = NULL, mets = NULL) {
  vapply(seq_along(score), function(i) mb_surv_prob(
    metagene, score[i],
    mycn = if (is.null(mycn)) 0 else mycn[i],
    myc  = if (is.null(myc))  0 else myc[i],
    mets = if (is.null(mets)) 0 else mets[i]), numeric(1))
}

# Survival-figure background from the stored density (no per-patient surv5/surv10).
.mb_surv_density_df <- function(slim) data.frame(x = slim$surv_density$x,
                                                 y = slim$surv_density$y)

.mb_density_plot <- function(bg, new_surv, indexRow, landmark) {
  if (length(new_surv) == 0 || all(is.na(new_surv)))
    new_surv <- rep(NA_real_, max(1, length(new_surv)))
  df.lines.ver <- data.frame(
    x = new_surv,
    colour = factor(ifelse(seq_along(new_surv) == indexRow, "highlight", "no.highlight"),
                    levels = c("highlight", "no.highlight")))
  ggplot(bg, aes(x = x, y = y)) +
    geom_area(fill = "steelblue", colour = "black", alpha = 0.5) +
    geom_vline(aes(xintercept = x, colour = colour), linetype = "dashed",
               linewidth = 0.5, data = df.lines.ver) +
    scale_colour_manual(values = c(highlight = "orange", no.highlight = "grey60")) +
    theme_classic() +
    xlab(paste0("Estimated ", landmark, "-year Survival")) +
    ylab("Density") +
    xlim(0, 1) +
    theme(legend.position = "none", text = element_text(size = 15))
}

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
#    _____ _    _ _    _ 
#   / ____| |  | | |  | |
#  | (___ | |__| | |__| |
#   \___ \|  __  |  __  |
#   ____) | |  | | |  | |
#  |_____/|_|  |_|_|  |_|
  
  
  
  

generate_figure_highlight_SHH <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}
   temp.df <- readRDS(file = "./mb_app/slim/SHHn49.dist.rds")
   

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
  if (is.null(indexRow)) indexRow <- 1
  bg       <- .mb_surv_density_df(.shh_slim)
  new_surv <- .mb_surv_vec("SHH", new.sample.meta.score, mycn = as.numeric(mycn))
  .mb_density_plot(bg, new_surv, indexRow, .shh_slim$landmark)
}



###################################################
#    _____ ____     ___  _     _            _______ ______ 
#   / ____|___ \   / / || |   | |        /\|__   __|  ____|
#  | |  __  __) | / /| || |_  | |       /  \  | |  | |__   
#  | | |_ ||__ < / / |__   _| | |      / /\ \ | |  |  __|  
#  | |__| |___) / /     | |   | |____ / ____ \| |  | |____ 
#   \_____|____/_/      |_|   |______/_/    \_\_|  |______|
  
  

  
generate_figure_highlight_G3_G4_no_sub <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}

  temp.df <- readRDS(file = "./mb_app/slim/G3_G4_no_subgroup.dist.rds")

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
  bg       <- .mb_surv_density_df(.g34late_slim)
  new_surv <- .mb_surv_vec("Group3/4 (Late)", new.sample.meta.score)
  .mb_density_plot(bg, new_surv, indexRow, .g34late_slim$landmark)
}


###################################################
#    _____ ____     ___  _     ______          _____  _  __     __
#   / ____|___ \   / / || |   |  ____|   /\   |  __ \| | \ \   / /
#  | |  __  __) | / /| || |_  | |__     /  \  | |__) | |  \ \_/ / 
#  | | |_ ||__ < / / |__   _| |  __|   / /\ \ |  _  /| |   \   /  
#  | |__| |___) / /     | |   | |____ / ____ \| | \ \| |____| |   
#   \_____|____/_/      |_|   |______/_/    \_\_|  \_\______|_|   
  
  

generate_figure_highlight_G3_G4_sub <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}

  temp.df <- readRDS(file = "./mb_app/slim/G3_G4_subgroup.dist.rds")

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
  if (is.null(indexRow)) indexRow <- 1
  bg       <- .mb_surv_density_df(.g34early_slim)
  new_surv <- .mb_surv_vec("Group3/4 (Early)", new.sample.meta.score,
                           myc = as.numeric(myc), mets = as.numeric(mets_val))
  .mb_density_plot(bg, new_surv, indexRow, .g34early_slim$landmark)
}


#for running without shiny

run_MB_risk_calculator <- function(
    idat_dir,
    metagene = c("SHH", "Group3/4 (Early)", "Group3/4 (Late)"),
    MYCN_Amplified = 0,
    MYC_Amplified = 0,
    Metastatic = 0
) {
  
  metagene <- match.arg(metagene)
  
  # Process IDATs
  temp.base <- get_basenames(idat_dir)
  temp.processed <- process_idats(temp.base)
  
  # Select model
  meta <- switch(
    metagene,
    "SHH" = SHH,
    "Group3/4 (Early)" = G3_G4_sub,
    "Group3/4 (Late)" = G3_G4_no_sub
  )
  
  # --------------------------------------------------
  # Validate covariates
  # --------------------------------------------------
  
  if (metagene == "SHH") {
    
    if (MYC_Amplified != 0 || Metastatic != 0) {
      stop(
        paste(
          "Group3/4 covariates were supplied for the SHH model.",
          "Use only MYCN_Amplified with metagene = 'SHH'."
        )
      )
    }
    
  }
  
  if (metagene == "Group3/4 (Early)") {
    
    if (MYCN_Amplified != 0) {
      stop(
        paste(
          "MYCN_Amplified is only valid for the SHH model.",
          "Remove MYCN_Amplified or set it to 0."
        )
      )
    }
    
  }
  
  if (metagene == "Group3/4 (Late)") {
    
    if (MYCN_Amplified != 0 ||
        MYC_Amplified != 0 ||
        Metastatic != 0) {
      
      stop(
        paste(
          "Group3/4 (Late) does not use any covariates.",
          "Set MYCN_Amplified = 0,",
          "MYC_Amplified = 0 and",
          "Metastatic = 0."
        )
      )
    }
    
  }
  
  # Calculate risk score
  test.res <- extract.metagene(
    as.character(meta[[1]]$genes),
    as.numeric(meta[[1]]$weights),
    beta2m(temp.processed$betas),
    as.numeric(meta[[2]])
  )
  
  test.res <- round(test.res, 3)
  
  
  cat("\n========================\n")
  cat("Model:", metagene, "\n")
  cat("========================\n\n")
  
  print(test.res)
  
  
  figure.input <- test.res$Risk_Value
  names(figure.input) <- rownames(test.res)
  
  # Risk plot
  risk_plot <- switch(
    metagene,
    "SHH" =
      generate_figure_highlight_SHH(figure.input, 1),
    
    "Group3/4 (Early)" =
      generate_figure_highlight_G3_G4_sub(figure.input, 1),
    
    "Group3/4 (Late)" =
      generate_figure_highlight_G3_G4_no_sub(figure.input, 1)
  )
  
  print(risk_plot)
  
  # Survival calculations
  if (metagene == "SHH") {
    surv_value <- mb_surv_prob("SHH", figure.input, mycn = MYCN_Amplified)
    surv_plot  <- generate_survival_figure_shh(figure.input, MYCN_Amplified, 1)
    timepoint  <- "5-year"
  } else if (metagene == "Group3/4 (Early)") {
    surv_value <- mb_surv_prob("Group3/4 (Early)", figure.input,
                               myc = MYC_Amplified, mets = Metastatic)
    surv_plot  <- generate_survival_figure_G3_G4_sub(figure.input, MYC_Amplified, Metastatic, 1)
    timepoint  <- "5-year"
  } else {
    surv_value <- mb_surv_prob("Group3/4 (Late)", figure.input)
    surv_plot  <- generate_survival_figure_G3_G4_no_sub(figure.input, 1)
    timepoint  <- "10-year"
  }
  
  print(surv_plot)
  
  cat(
    paste0(
      "\nEstimated ", timepoint,
      " survival (%)\n"
    )
  )
  
  print(round(100 * surv_value, 1))
  
  return(
    list(
      risk_scores = test.res,
      survival_probability = surv_value,
      risk_plot = risk_plot,
      survival_plot = surv_plot
    )
  )
}

##slim_mb addition
# slim_mb.R ─ prediction-only payloads for the MB_Risk_App survival models.
# No $y, $linear.predictors, $residuals, no captured environment, no per-patient
# rows. Baseline survival is evaluated on a fixed clinical grid so exact training
# death times are never shipped.

suppressMessages(library(survival))

# ── Step-function evaluation of a survival curve on an arbitrary grid ──────────
# Survival curves are right-continuous step functions: S(t) = surv at the largest
# event time <= t, and S(t)=1 before the first event. findInterval gives exactly
# that. This is what severs the link to native death times.
.step_surv <- function(time, surv, grid) {
  ord  <- order(time)
  time <- time[ord]; surv <- surv[ord]
  idx  <- findInterval(grid, time)            # 0 before first event time
  out  <- ifelse(idx == 0, 1, surv[pmax(idx, 1)])
  out
}

# ── Convert one coxph to a slim payload ───────────────────────────────────────
slim_cox <- function(fit, data = NULL, grid = NULL, step = 0.25, landmarks = NULL) {
  stopifnot(inherits(fit, "coxph"))
  
  svars <- .strata_vars(fit)            # 0 or 1 strata variables for these models
  if (length(svars) > 1)
    stop("multi-variable strata not handled; extend .strata_level mapping")
  
  # If a training df is supplied, rebuild a CLEAN terms env from it so survfit()
  # can reconstruct the baseline regardless of how the original env was captured.
  if (!is.null(data)) fit <- .rebuild_fit_env(fit, data)
  
  # Design-matrix terms: rebuild from the non-strata term labels only, so
  # model.matrix never references strata(); env reset so nothing rides along.
  full_labs <- attr(delete.response(terms(fit)), "term.labels")
  keep_labs <- full_labs[!grepl("strata\\(", full_labs)]
  tt <- if (length(keep_labs)) terms(reformulate(keep_labs)) else terms(~ 1)
  attr(tt, ".Environment") <- baseenv()
  
  # Baseline survival at mean covariates, per stratum, on the grid. survfit
  # orders strata by factor level, so we key baselines by level order.
  sf <- survfit(fit)
  # Fixed clinical grid; ceiling() hides the exact maximum follow-up time.
  if (is.null(grid)) {
    top  <- ceiling(max(sf$time))
    if (!is.null(landmarks)) top <- max(top, ceiling(max(landmarks)))
    grid <- seq(0, top, by = step)
  }
  # Force landmark evaluation times onto the grid: a survival step function only
  # changes at event times, so an off-grid landmark would not be reproduced exactly.
  if (!is.null(landmarks)) grid <- sort(unique(c(grid, landmarks[landmarks >= 0])))
  if (is.null(sf$strata)) {
    baseline <- list(.step_surv(sf$time, sf$surv, grid))
    slevels  <- character(0)
  } else {
    grp      <- rep(seq_along(sf$strata), sf$strata)
    baseline <- lapply(split(seq_along(sf$time), grp), function(ii)
      .step_surv(sf$time[ii], sf$surv[ii], grid))
    slevels  <- fit$xlevels[[paste0("strata(", svars, ")")]]
    if (length(slevels) != length(baseline))
      stop("strata count mismatch between xlevels and survfit")
  }
  
  list(
    type        = "slim_coxph",
    coef        = coef(fit),
    means       = fit$means,
    xlevels     = fit$xlevels[intersect(names(fit$xlevels), all.vars(tt))],
    terms       = tt,
    strata_vars = svars,
    strata_levels = slevels,
    grid        = grid,
    baseline    = baseline,        # list, indexed by strata level order (or [[1]])
    n_train     = NA_integer_      # deliberately not stored from $y
  )
}

# Rebuild a CLEAN terms/formula environment on a fit from a training df, binding
# the Surv response under whatever symbol the formula actually uses (e.g.
# surv.object.train, surv.object.train.late, surv.object.train.shh). Returns the
# fit with its terms and formula environments reset and call$data cleared, so
# survfit() can reconstruct the baseline regardless of the original session.
# Used only at conversion time; the df itself is never shipped.
.rebuild_fit_env <- function(fit, data) {
  for (v in names(fit$xlevels)) {              # coerce factors to the model's levels
    nm <- sub("^strata\\((.*)\\)$", "\\1", v)
    if (!is.null(data[[nm]]))
      data[[nm]] <- factor(as.character(data[[nm]]), levels = fit$xlevels[[v]])
  }
  e <- list2env(as.list(data), parent = baseenv())
  if (all(c("OS_Time","OS_Status") %in% names(data))) {
    sv   <- survival::Surv(data$OS_Time, data$OS_Status)
    resp <- all.vars(attr(terms(fit), "variables")[[2]])[1]   # actual response symbol
    e[[resp]] <- sv
  }
  e$strata <- survival::strata
  e$Surv   <- survival::Surv
  environment(fit$terms) <- e
  if (!is.null(fit$formula)) environment(fit$formula) <- e
  if (!is.null(fit$call))    fit$call$data <- NULL             # force frame from env
  fit
}

# Names of variables wrapped in strata() in the model formula.
.strata_vars <- function(fit) {
  tt  <- terms(fit)
  sp  <- attr(tt, "specials")$strata
  if (is.null(sp)) return(character(0))
  vars <- attr(tt, "variables")
  vapply(sp, function(i) {
    e <- vars[[i + 1]]            # +1: first element of `variables` is `list`
    all.vars(e)[1]
  }, character(1))
}

# ── Predict survival for a new patient from a slim payload ────────────────────
# newdata: one-row data.frame (or named list) with the model covariates.
predict_survival <- function(slim, newdata, times = slim$grid) {
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
  
  # Design row, factors honoured via stored xlevels; strata cols are ignored
  # for the linear predictor because we select by names(coef) below.
  mm  <- model.matrix(slim$terms, newdata, xlev = slim$xlevels)
  cn  <- names(slim$coef)
  if (!all(cn %in% colnames(mm)))
    stop("newdata does not produce all model terms: ",
         paste(setdiff(cn, colnames(mm)), collapse = ", "))
  xrow <- mm[1, cn, drop = TRUE]
  
  lp  <- sum(slim$coef * (xrow - slim$means[cn]))
  
  # Select baseline stratum by factor-level order (robust to label formatting).
  if (length(slim$strata_vars) == 0) {
    S0_grid <- slim$baseline[[1]]
  } else {
    v <- as.character(newdata[[slim$strata_vars]])
    k <- match(v, slim$strata_levels)
    if (is.na(k))
      stop("strata value '", v, "' not in training levels: ",
           paste(slim$strata_levels, collapse = " | "))
    S0_grid <- slim$baseline[[k]]
  }
  
  # Map requested times onto the stored grid (step function again).
  S0_at <- .step_surv(slim$grid, S0_grid, times)
  surv  <- S0_at ^ exp(lp)
  data.frame(time = times, surv = surv)
}

# ── Aggregate a per-patient score distribution ────────────────────────────────
# Returns a DROP-IN data.frame with the SAME single column name as the input,
# holding n_q sorted quantiles instead of the raw per-patient scores. The
# highlight figures plot sorted scores against their rank, so a quantile vector
# reproduces the curve shape and percentile position without exposing any
# individual score. n is kept as an attribute for reference.
slim_dist <- function(d, n_q = 100L) {
  col <- names(d)[1]
  v   <- d[[1]]; v <- v[is.finite(v)]
  q   <- as.numeric(quantile(v, probs = seq(0, 1, length.out = n_q),
                             names = FALSE, type = 7))
  out <- setNames(data.frame(q), col)
  attr(out, "n_train") <- length(v)
  out
}

# Aggregate a per-patient predicted-survival column (surv5 / surv10) into a
# density curve for the survival-figure background. No per-patient rows.
slim_surv_density <- function(values, n = 512L, from = 0, to = 1) {
  v <- values[is.finite(values)]
  dn <- density(v, n = n, from = from, to = to)
  list(x = dn$x, y = dn$y, n = length(v))
}

# Percentile of a new score against the aggregated (sorted-quantile) distribution.
percentile_of <- function(slim_d, score) {
  q <- slim_d[[1]]
  vapply(score, function(s) max(c(0, which(q < s))) / length(q), numeric(1))
}
