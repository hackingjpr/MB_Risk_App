###########
## Create dist and model rds objects from larger data frames
#Load in rds files  
SHHn49 <- readRDS(file.path("~/MB_Risk_App/mb_app/SHH_2_train_no_subgroup/n49.extract.re.rds"))
distSHHn49 <- read.csv(file.path("~/MB_Risk_App/mb_app/df.pheno.mb.combined.shh.csv"))
#Create dist and model rds objects from larger data frames
SHHn49.model <- list(SHHn49[[8]], SHHn49[[9]])
names(SHHn49.model) <- c("cpgs", "scaling")
SHHn49.dist <- data.frame(distSHHn49$Index.ID, distSHHn49$fixedMG)
colnames(SHHn49.dist) <- c("sample.name", "SHH")
SHHn49.dist <- SHHn49.dist[order(SHHn49.dist$SHH),]

#save objects to rds files
saveRDS(SHHn49.model, file.path("~/MB_Risk_App/mb_app/SHHn49.model.rds"))
saveRDS(SHHn49.dist, file.path("~/MB_Risk_App/mb_app/SHHn49.dist.rds"))



G3_G4_no_subgroup <- readRDS(file.path("~/MB_Risk_App/mb_app/G3_G4_2_train_no_subgroup.late/n49.extract.re.rds"))
distG3_G4_no_subgroup <- read.csv(file.path("~/MB_Risk_App/mb_app/df.pheno.mb.combined.grp3.grp4.late.csv"))
#create dist and model rds objects from larger data frames
G3_G4_no_subgroup.model <- list(G3_G4_no_subgroup[[8]], G3_G4_no_subgroup[[9]])
names(G3_G4_no_subgroup.model) <- c("cpgs", "scaling")
G3_G4_no_subgroup.dist <- data.frame(distG3_G4_no_subgroup$Index.ID, distG3_G4_no_subgroup$fixedMG)
colnames(G3_G4_no_subgroup.dist) <- c("sample.name", "G3_G4_no_sub")
G3_G4_no_subgroup.dist <- G3_G4_no_subgroup.dist[order(G3_G4_no_subgroup.dist$G3_G4_no_sub),]

#save objects to rds files
saveRDS(G3_G4_no_subgroup.model, file.path("~/MB_Risk_App/mb_app/G3_G4_no_subgroup.model.rds"))
saveRDS(G3_G4_no_subgroup.dist, file.path("~/MB_Risk_App/mb_app/G3_G4_no_subgroup.dist.rds"))


G3_G4_subgroup <- readRDS(file.path("~/MB_Risk_App/mb_app/G3_G4_2_train_subgroup/n19.extract.re.rds"))
distG3_G4_subgroup <- read.csv(file.path("~/MB_Risk_App/mb_app/df.pheno.mb.combined.grp3.grp4.csv"))
#Create dist and model rds objects from larger data frames
G3_G4_subgroup.model <- list(G3_G4_subgroup[[8]], G3_G4_subgroup[[9]])
names(G3_G4_subgroup.model) <- c("cpgs", "scaling")
G3_G4_subgroup.dist <- data.frame(distG3_G4_subgroup$Index.ID, distG3_G4_subgroup$fixedMG)
colnames(G3_G4_subgroup.dist) <- c("sample.name", "G3_G4_sub")
G3_G4_subgroup.dist <- G3_G4_subgroup.dist[order(G3_G4_subgroup.dist$G3_G4_sub),]
#save objects to rds files
saveRDS(G3_G4_subgroup.model, file.path("~/MB_Risk_App/mb_app/G3_G4_subgroup.model.rds"))
saveRDS(G3_G4_subgroup.dist, file.path("~/MB_Risk_App/mb_app/G3_G4_subgroup.dist.rds")) 
