### Etude des facteurs influençant l'incrément de biomasse des espèces d'arbres de Mikembo ###
###### ******************************************************************************** ######

# Useful packages and functions:
# ******************************
library(usdm)
library(vegan)
library(glmulti)

# Input of the datasets:
# **********************
data <- read.table("data_BAI_individual.txt", h = T, sep = "\t")
colnames(data)

dataEnv <- read.table("data_ENV_BAI_individual.txt", h = T, sep = "\t", row.names = 1)

# Environmental variable selection:
# *********************************
# VIF check:
vif(dataEnv)
vifcor(dataEnv, th = 0.7)
env <- dataEnv[, -c(17, 9, 19, 11, 13, 10, 18, 6, 4, 15, 20, 2, 14)]
# VIF check, after removing the collinear variables:
vif(env)

cor <- round(cor(dataEnv), 3)
write.table(cor, file = "expl-var_correlations.txt", sep = "\t")

#pca <- rda(env, scale = TRUE)
#cleanplot.pca(pca)

# Lists of results:
# *****************
# To retain the 100 best models for each species and the whole community:
list.modsel <- vector("list", 20)
names(list.modsel) <- c("AlbAnt", "BraSpi", "BraWan", 'ComCol', "DalBoe", "DipCon", 
                        "HexMon", "JulGlo", "JulPan", "MarMac", "MonKat", "ParCur", "PerAng",
                        "PhiKat", "PseMap", "PteAng", "PteTin", "StrInn", "UapNit", "ALL")
# To save the corresponding model averaging results:
list.modav <- list.modsel

# Input of the species x trait data:
trait <- read.table("sp_x_traits.txt", h = T, sep = "\t", row.names = 1)

# We thin the 'trait' dataframe to the species of interest only:
match <- c()
for (i in 1:nrow(trait)) {
  if (is.na(match(rownames(trait)[i], names(list.modsel))) == FALSE) match <- c(match, i)
}
ft <- trait[match, ]   # All traits of the 19 species (for the whole community analysis only)

# Loop for all the selected species and the whole community:
# **********************************************************
depart <- Sys.time()

for (i in 1:length(list.modsel)) {
  
  if (i != length(list.modsel)) {
    sub <- subset(data, species == levels(species)[i])
  } else sub <- data
  
  # Generation of the environmental dataset for the species of 'sub':
  matenv <- matrix(0, nrow = nrow(sub), ncol = ncol(env))
  colnames(matenv) <- colnames(env)
  for (j in 1:nrow(sub)) {
    row <- match(sub$sites[j], rownames(env))
    matenv[j, ] <- as.numeric(env[row, ])
  }
  # We bind the DBH_2010 and square of the DBH_2010 to the environmental predictors:
  subenv <- cbind(sub[, 6], matenv)
  colnames(subenv)[1] <- colnames(sub)[6]
  
  # If we analyse the whole community, then we add to subenv the individual values of the FT of
  # the species of each individual:
  if (i == length(list.modsel)) {
    matTrait <- matrix(0, nrow = nrow(sub), ncol = ncol(ft))
    colnames(matTrait) <- colnames(ft)
    for (j in 1:nrow(sub)) {
      row <- match(sub$species[j], rownames(ft))
      matTrait[j, ] <- as.numeric(ft[row, ])
    }
    subenv <- cbind(subenv, matTrait)
  }

  resp <- as.vector(sub[, 5])   # The response variable
  lm.test <- lm(resp ~ ., data = as.data.frame(subenv))
  # We limit the number of predictors to 12:
  # We repeat the genetic algorithm 20 times and use a consensus of the replicates to improve
  # convergence (Calcagno et al. 2010):
  nb.replic <- 20
  replic.list <- vector("list", nb.replic)
  for (h in 1:nb.replic) {
    set.seed(h)
    replic.list[[h]] <- glmulti(y = lm.test, level = 1, minsize = 0, maxsize = 12, 
                                method = "g", crit = "aicc", confsetsize = 100, popsize = 100, 
                                sexrate = 0.1, imm = 0.3, mutrate = 10^-3, deltaM = 0.05, 
                                deltaB = 0.05, conseq = 5, includeobjects = TRUE, 
                                plotty = TRUE, report = TRUE)
  }
  consen <- consensus(replic.list, confsetsize = 100)
  list.modsel[[i]] <- consen
}

stop <- Sys.time()
(diff <- stop - depart)

# Model averaging:
# ****************
modav.list <- vector("list", length(list.modsel))
for (i in 1:length(list.modsel)) modav.list[[i]] <- coef(list.modsel[[i]])

# Synthetic tables (one for the coefficient signs, and another for the importance values):
# ****************************************************************************************
summary.mat.sign <- matrix(nrow = length(modav.list), ncol = ncol(env)+2+ncol(matTrait))
colnames(summary.mat.sign) <- c("(Intercept)", "DBH_2010", colnames(env), colnames(matTrait))
rownames(summary.mat.sign) <- c(levels(data$species), "COMMUNITY")

summary.mat.imp <- summary.mat.sign

for (i in 1:nrow(summary.mat.sign)) {
  sel <- which(modav.list[[i]][, 4] >= 0.8)
  sign <- sapply(modav.list[[i]][, 1][sel], function(x) ifelse(x > 0, "+", "-"))
  col.ID <- sapply(names(sign), function(x) match(x, colnames(summary.mat.sign)))
  summary.mat.sign[i, col.ID] <- sign
  summary.mat.imp[i, col.ID] <- round(modav.list[[i]][, 4][sel], 3)
}

write.table(summary.mat.sign, file = "results_modav_sign.txt", sep = "\t")
write.table(summary.mat.imp, file = "results_modav_impVal.txt", sep = "\t")
