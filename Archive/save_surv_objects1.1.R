library(survival)

# ── SHH ─────────────────────────────────────────────────────────────────────
shh_train_fit <- readRDS("~/MB_Risk_App/mb_app/shh49.train.fit.rds")
shh_df        <- read.csv("~/MB_Risk_App/mb_app/df.pheno.mb.combined.shh.csv")

# Strip to needed columns only
shh_df <- shh_df[, c("OS_Time", "OS_Status", "fixedMG", "ConsensusMYCN")]

# Rebuild objects the stored call needs
surv.object.train.shh    <- Surv(shh_df$OS_Time, shh_df$OS_Status)
df.pheno.mb.combined.shh <- shh_df

# Refit using stripped data — identical coefficients, clean references
shh_fit_clean <- update(shh_train_fit, data = shh_df)
shh_fit_clean$call <- call("coxph", formula = shh_fit_clean$formula)

# Precompute training cohort survival
shh_surv <- summary(
  survfit(shh_fit_clean,
          newdata = data.frame(fixedMG       = shh_df$fixedMG,
                               ConsensusMYCN = shh_df$ConsensusMYCN)),
  time = 5
)
shh_df$surv5 <- as.numeric(shh_surv$surv)

saveRDS(shh_df,        "~/MB_Risk_App/mb_app/precomp_shh_df.rds")
saveRDS(shh_fit_clean, "~/MB_Risk_App/mb_app/precomp_shh_fit.rds")
message("SHH saved.")


#g34 late

library(survival)

g34late_fit     <- readRDS("~/MB_Risk_App/mb_app/g34late49.train.fit.rds")
g34late_df_full <- read.csv("~/MB_Risk_App/mb_app/df.pheno.mb.combined.grp3.grp4.late.csv")

g34late_df_clean <- g34late_df_full[complete.cases(g34late_df_full[, c("OS_Time", "OS_Status", "fixedMG")]), ]

surv.object.train.late               <- Surv(g34late_df_clean$OS_Time, g34late_df_clean$OS_Status)
df.pheno.mb.combined.grp3.grp4.late <- g34late_df_clean

g34late_fit_clean <- update(g34late_fit, data = g34late_df_clean)
g34late_fit_clean$call <- call("coxph", formula = g34late_fit_clean$formula)

# Now check
message("fit n: ", g34late_fit_clean$n)
message("df rows: ", nrow(g34late_df_clean))
message("fixedMG length: ", length(g34late_df_clean$fixedMG))
message("residuals length: ", length(g34late_fit_clean$residuals))

env <- environment(g34late_fit_clean$terms)
env$surv.object.train.late <- Surv(g34late_df_clean$OS_Time, g34late_df_clean$OS_Status)
env$fixedMG                <- g34late_df_clean$fixedMG

g34late_surv <- summary(
  survfit(g34late_fit_clean,
          data.frame(fixedMG = g34late_df_clean$fixedMG)),
  time = 10
)

g34late_df_save <- data.frame(
  fixedMG  = g34late_df_clean$fixedMG,
  OS_Time  = g34late_df_clean$OS_Time,
  OS_Status = g34late_df_clean$OS_Status,
  surv10   = as.numeric(g34late_surv$surv)
)

saveRDS(g34late_df_save,   "~/MB_Risk_App/mb_app/precomp_g34late_df.rds")
saveRDS(g34late_fit_clean, "~/MB_Risk_App/mb_app/precomp_g34late_fit.rds")
message("G3G4_no_sub saved.")


# ── G3G4 subgroup ───────────────────────────────────────────────────────────
# ── G3G4 subgroup ───────────────────────────────────────────────────────────
g34early_fit <- readRDS("~/MB_Risk_App/mb_app/g34early19.train.fit.rds")
g34early_df  <- read.csv("~/MB_Risk_App/mb_app/df.pheno.mb.combined.grp3.grp4.csv")

message("G3G4 early call: ")
print(g34early_fit$call)

# Strip to needed columns only
g34early_df <- g34early_df[complete.cases(g34early_df[, c("OS_Time", "OS_Status", "fixedMG", "M._versus_M.", "ConsensusMYC")]), ]

# Clean the factor columns
g34early_df$M._versus_M.[is.na(g34early_df$M._versus_M.)] <- 0
g34early_df$ConsensusMYC <- gsub("^ConsensusMYC=", "", g34early_df$ConsensusMYC)
g34early_df$ConsensusMYC[is.na(g34early_df$ConsensusMYC)] <- "0"
g34early_df$ConsensusMYC <- factor(
  paste0("ConsensusMYC=", g34early_df$ConsensusMYC),
  levels = g34early_fit$xlevels[["strata(ConsensusMYC)"]]
)

# Rebuild objects the stored call references by name
surv.object.train               <- Surv(g34early_df$OS_Time, g34early_df$OS_Status)
df.pheno.mb.combined.grp3.grp4 <- g34early_df

# Refit using clean data
g34early_fit_clean <- update(g34early_fit, data = g34early_df)
g34early_fit_clean$call <- call("coxph", formula = g34early_fit_clean$formula)

# Inject into terms environment
env <- environment(g34early_fit_clean$terms)
env$surv.object.train               <- surv.object.train
env$fixedMG                         <- g34early_df$fixedMG
env$M._versus_M.                    <- g34early_df$M._versus_M.
env$ConsensusMYC                    <- g34early_df$ConsensusMYC

# Verify survfit works
message("Testing survfit...")
print(survfit(g34early_fit_clean))

# Precompute training cohort survival
g34early_newdata <- data.frame(
  fixedMG      = g34early_df$fixedMG,
  M._versus_M. = g34early_df$M._versus_M.,
  ConsensusMYC = g34early_df$ConsensusMYC
)

g34early_surv <- summary(
  survfit(g34early_fit_clean, newdata = g34early_newdata),
  time = 5
)

# Save with OS_Time and OS_Status included for reconstruction on load
g34early_df_save <- data.frame(
  fixedMG      = g34early_df$fixedMG,
  M._versus_M. = g34early_df$M._versus_M.,
  ConsensusMYC = as.character(g34early_df$ConsensusMYC),
  OS_Time      = g34early_df$OS_Time,
  OS_Status    = g34early_df$OS_Status,
  surv5        = as.numeric(g34early_surv$surv)
)

saveRDS(g34early_df_save,   "~/MB_Risk_App/mb_app/precomp_g34early_df.rds")
saveRDS(g34early_fit_clean, "~/MB_Risk_App/mb_app/precomp_g34early_fit.rds")
message("G3G4_sub saved.")
message("All done.")
