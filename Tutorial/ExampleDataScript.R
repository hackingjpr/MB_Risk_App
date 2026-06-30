# --------------------------------------------------
# User settings
# --------------------------------------------------

idat_dir <- "~/Idats/Mix"

metagene_name <- "SHH"
# metagene_name <- "Group3/4 (Early)"
# metagene_name <- "Group3/4 (Late)"

# --------------------------------------------------
# Optional covariates
# --------------------------------------------------

# SHH
MYCN_Amplified <- 0

# Group3/4 Early
MYC_Amplified <- 0
Metastatic <- 0

# --------------------------------------------------
# Process IDATs
# --------------------------------------------------

temp.base <- get_basenames(idat_dir)

temp.processed <- process_idats(temp.base)

meta <- switch(
  metagene_name,
  "SHH"              = SHH,
  "Group3/4 (Early)" = G3_G4_sub,
  "Group3/4 (Late)"  = G3_G4_no_sub
)

test.res <- extract.metagene(
  as.character(meta[[1]]$genes),
  as.numeric(meta[[1]]$weights),
  beta2m(temp.processed$betas),
  as.numeric(meta[[2]])
)

test.res <- round(test.res, 3)

print(test.res)

figure.input <- test.res$Risk_Value
names(figure.input) <- rownames(test.res)

# --------------------------------------------------
# Risk distribution plot
# --------------------------------------------------

risk_plot <- switch(
  metagene_name,
  "SHH" = generate_figure_highlight_SHH(
    figure.input,
    1
  ),
  "Group3/4 (Early)" = generate_figure_highlight_G3_G4_sub(
    figure.input,
    1
  ),
  "Group3/4 (Late)" = generate_figure_highlight_G3_G4_no_sub(
    figure.input,
    1
  )
)

print(risk_plot)

# --------------------------------------------------
# Survival estimate
# --------------------------------------------------

if (metagene_name == "SHH") {
  
  fit <- survfit(
    .shh_fit,
    newdata = data.frame(
      fixedMG = as.numeric(figure.input),
      ConsensusMYCN = MYCN_Amplified
    )
  )
  
  surv_value <- summary(fit, time = 5)$surv
  
  cat(
    "\nEstimated 5-year survival (%)\n"
  )
  
  print(round(100 * surv_value, 1))
  
  surv_plot <- generate_survival_figure_shh(
    figure.input,
    MYCN_Amplified,
    1
  )
  
  print(surv_plot)
  
} else if (metagene_name == "Group3/4 (Early)") {
  
  surv_plot <- generate_survival_figure_G3_G4_sub(
    figure.input,
    MYC_Amplified,
    Metastatic,
    1
  )
  
  print(surv_plot)
  
  fit <- survfit(
    .g34early_fit,
    newdata = data.frame(
      fixedMG = as.numeric(figure.input),
      M._versus_M. = Metastatic,
      ConsensusMYC = factor(
        paste0("ConsensusMYC=", MYC_Amplified),
        levels = .g34early_fit$xlevels[["strata(ConsensusMYC)"]]
      )
    )
  )
  
  surv_value <- summary(fit, time = 5)$surv
  
  cat(
    "\nEstimated 5-year survival (%)\n"
  )
  
  print(round(100 * surv_value, 1))
  
} else {
  
  surv_plot <- generate_survival_figure_G3_G4_no_sub(
    figure.input,
    1
  )
  
  print(surv_plot)
  
  fit <- survfit(
    .g34late_fit,
    newdata = data.frame(
      fixedMG = as.numeric(figure.input)
    )
  )
  
  surv_value <- summary(fit, time = 10)$surv
  
  cat(
    "\nEstimated 10-year survival (%)\n"
  )
  
  print(round(100 * surv_value, 1))
}
