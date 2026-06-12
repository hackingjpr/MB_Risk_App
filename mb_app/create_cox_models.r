library('survival')
df.pheno.mb.combined.grp3.grp4 <- read.csv(file = "./df.pheno.mb.combined.grp3.grp4.csv")
surv.object.train <- Surv(df.pheno.mb.combined.grp3.grp4$OS_Time, df.pheno.mb.combined.grp3.grp4$OS_Status)
coxph(surv.object.train~fixedMG+M._versus_M.+strata(ConsensusMYC), data = df.pheno.mb.combined.grp3.grp4) -> train.fit
train.fit
saveRDS(train.fit, file = "./g34early19.train.fit.rds")

df.pheno.mb.combined.grp3.grp4.late <- read.csv(file = "./df.pheno.mb.combined.grp3.grp4.late.csv")
surv.object.train.late <- Surv(df.pheno.mb.combined.grp3.grp4.late$OS_Time, df.pheno.mb.combined.grp3.grp4.late$OS_Status)
coxph(surv.object.train.late~fixedMG, data = df.pheno.mb.combined.grp3.grp4.late) -> train.fit.late
train.fit.late
saveRDS(train.fit.late, file = "./g34late49.train.fit.rds")

df.pheno.mb.combined.shh <- read.csv(file = "./df.pheno.mb.combined.shh.csv")
surv.object.train.shh <- Surv(df.pheno.mb.combined.shh$OS_Time, df.pheno.mb.combined.shh$OS_Status)
coxph(surv.object.train.shh~fixedMG+ConsensusMYCN, data = df.pheno.mb.combined.shh) -> train.fit.shh
train.fit.shh
saveRDS(train.fit.shh, file = "./shh49.train.fit.rds")
