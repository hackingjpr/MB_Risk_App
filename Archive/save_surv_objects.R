library(survival)
library(foreach)

# ── SHH ─────────────────────────────────────────────────────────────────────
shh_train_fit <- readRDS("~/MB_Risk_App/mb_app/shh49.train.fit.rds")
shh_df        <- read.csv("~/MB_Risk_App/mb_app/df.pheno.mb.combined.shh.csv")

shh_surv <- summary(
  survfit(shh_train_fit,
          data.frame(fixedMG       = shh_df$fixedMG,
                     ConsensusMYCN = shh_df$ConsensusMYCN)),
  time = 5
)
shh_df$surv5 <- as.numeric(shh_surv$surv)

saveRDS(shh_df,        "~/MB_Risk_App/mb_app/precomp_shh_df.rds")
saveRDS(shh_train_fit, "~/MB_Risk_App/mb_app/precomp_shh_fit.rds")
message("SHH saved.")



shh_train_fit <- readRDS("~/MB_Risk_App/mb_app/shh49.train.fit.rds")
shh_df_stripped <- read.csv("~/MB_Risk_App/mb_app/df.pheno.mb.combined.shh.stripped.csv")

surv.object.train.shh <- Surv(shh_df_stripped$OS_Time, shh_df_stripped$OS_Status)

# This will tell you exactly which columns are missing
test_fit <- update(shh_train_fit, data = shh_df_stripped)





# ── G3G4 no subgroup ────────────────────────────────────────────────────────
g34late_fit <- readRDS("~/MB_Risk_App/mb_app/g34late49.train.fit.rds")
g34late_df  <- read.csv("~/MB_Risk_App/mb_app/df.pheno.mb.combined.grp3.grp4.late.csv")

g34late_surv <- summary(
  survfit(g34late_fit,
          data.frame(fixedMG = g34late_df$fixedMG)),
  time = 10
)
g34late_df$surv10 <- as.numeric(g34late_surv$surv)

saveRDS(g34late_df,  "~/MB_Risk_App/mb_app/precomp_g34late_df.rds")
saveRDS(g34late_fit, "~/MB_Risk_App/mb_app/precomp_g34late_fit.rds")
message("G3G4_no_sub saved.")

# ── G3G4 subgroup ───────────────────────────────────────────────────────────
g34early_fit <- readRDS("~/MB_Risk_App/mb_app/g34early19.train.fit.rds")
g34early_df  <- read.csv("~/MB_Risk_App/mb_app/df.pheno.mb.combined.grp3.grp4.csv")

g34early_df$M._versus_M.[is.na(g34early_df$M._versus_M.)] <- 0
g34early_df$ConsensusMYC <- gsub("^ConsensusMYC=", "", g34early_df$ConsensusMYC)
g34early_df$ConsensusMYC[is.na(g34early_df$ConsensusMYC)] <- "0"
g34early_df$ConsensusMYC <- factor(
  paste0("ConsensusMYC=", g34early_df$ConsensusMYC),
  levels = g34early_fit$xlevels[["strata(ConsensusMYC)"]]
)

g34early_newdata <- data.frame(
  fixedMG      = g34early_df$fixedMG,
  M._versus_M. = g34early_df$M._versus_M.,
  ConsensusMYC = g34early_df$ConsensusMYC
)

g34early_surv <- summary(
  survfit(g34early_fit, newdata = g34early_newdata),
  time = 5
)
g34early_df$surv5 <- as.numeric(g34early_surv$surv)

saveRDS(g34early_df,  "~/MB_Risk_App/mb_app/precomp_g34early_df.rds")
saveRDS(g34early_fit, "~/MB_Risk_App/mb_app/precomp_g34early_fit.rds")
message("G3G4_sub saved.")

message("All done.")
