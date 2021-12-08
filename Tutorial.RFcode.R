# Random forests of classification trees for nonlinear, interacting predictors

### Prep

# Load packages
library(randomForestSRC)
library(tidyverse)

# Read the site module data
phyto.site.mods <- read.csv("bin.aug.DIRT.site.mods.csv", header = T)

# Read the species module data
phyto.sp.mods <- read.csv("bin.aug.DIRT.sp.mods.csv", header = T)

# Read the environmental data
phyto.env <- read.csv("Loewen_et_al._2020._Ecol_App._Environmental_data_-_final.csv", header = T)
phyto.env.aug <- subset(phyto.env, Month == "Aug")
rownames(phyto.env.aug) <- phyto.env.aug$Lake.Name

# Read phytoplankton trait data
phyto.traits <- read.csv("lno11694-sup-0002-supinfo02.csv", header = T)
phyto.traits$SPECIES.NAME <- gsub(" ", ".", phyto.traits$SPECIES.NAME)
phyto.traits$SPECIES.NAME <- gsub("\\.\\.", ".", phyto.traits$SPECIES.NAME)
rownames(phyto.traits) <- phyto.traits$SPECIES.NAME

# Merge site module and environmental data
phyto.site.mods <- merge(phyto.site.mods, phyto.env.aug[, 6:46], by = 0)
rownames(phyto.site.mods) <- phyto.site.mods$Row.names
phyto.site.mods <- phyto.site.mods[, -1]
phyto.site.mods$Mods <- as.factor(phyto.site.mods$Mods)

# Merge species module and trait data
phyto.sp.mods <- merge(phyto.sp.mods, phyto.traits, by = 0)
rownames(phyto.sp.mods) <- phyto.sp.mods$Row.names
phyto.sp.mods <- phyto.sp.mods[, -1]
phyto.sp.mods$Mods <- as.factor(phyto.sp.mods$Mods)
traits <- phyto.sp.mods[, c("GALD", "N.FIXING", "BUOYANCY", "SILICA_USING", "FLAGELLA", "HETEROTROPHY", "PHYCOBILINS",
                            "CHL_C", "CHL_B", "CHAIN.COLONY", "MUCILAGINOUS_SHEATH", "SPINES", "SUSPECTED_TOXIN_PRODUCING")]

### Correlations

# Load packages
library(corrplot)

# Obtain simple correlation matrix from environmental candidate variables
phyto.site.mods.cor <- cor(phyto.site.mods[, 11:51])

# Obtain a quick plot using corrplot with order of variables based on first principal component
phyto.site.mods.cor.plot <- corrplot(phyto.site.mods.cor, method = 'circle', 
                                     type = "upper", order = 'FPC', diag = FALSE)

### Variable selection

# Load packages
library(VSURF)
library(parallel)

# Create formula with all possible variables
site.names <- colnames(phyto.site.mods[, 11:51])
phyto.site.mods.all.fmla <- as.formula(paste("Mods ~ ", paste(site.names, collapse = "+")))

# Set seed for reproducibility
set.seed(100, kind = "L'Ecuyer-CMRG")

# Use VSURF to perform feature selection (in parallel) based on permutation-based variable importance
phyto.site.mods.vsurf <- VSURF(phyto.site.mods.all.fmla,
                               data = phyto.site.mods,
                               ntree = 10000,
                               mtry = sqrt(41),
                               nfor.thres = 500,
                               nmin = 1,
                               nfor.interp = 250,
                               nsd = 1,
                               nfor.pred = 250,
                               nmj = 1,
                               RFimplem = "randomForest",
                               parallel = TRUE,
                               ncores = parallel::detectCores() - 1,
                               clusterType = "PSOCK",
                               verbose = TRUE)

# Report VSURF results
summary(phyto.site.mods.vsurf)

# Plot VSURF results
plot(phyto.site.mods.vsurf)

# Create vector of variables selected after threshold step
phyto.site.mods.vsurf.thres <- phyto.site.mods.vsurf$varselect.thres
phyto.site.mods.vsurf.thres.names <- names[phyto.site.mods.vsurf.thres]

# Create formula with variables selected after threshold strep
phyto.site.mods.vsurf.thres.fmla <- as.formula(paste("Mods ~ ", paste(phyto.site.mods.vsurf.thres.names, collapse = "+")))

# Create vector of variables selected after interpretation step
phyto.site.mods.vsurf.interp <- phyto.site.mods.vsurf$varselect.interp
phyto.site.mods.vsurf.interp.names <- names[phyto.site.mods.vsurf.interp]
phyto.site.mods.vsurf.interp.names

# Create formula with variables selected after interpretation strep
phyto.site.mods.vsurf.interp.fmla <- as.formula(paste("Mods ~ ", paste(phyto.site.mods.vsurf.interp.names, collapse = "+")))

# Create vector of variables selected after prediction step
phyto.site.mods.vsurf.pred <- phyto.site.mods.vsurf$varselect.pred
phyto.site.mods.vsurf.pred.names <- names[phyto.site.mods.vsurf.pred]
phyto.site.mods.vsurf.pred.names

# Create formula with variables selected after prediction strep
phyto.site.mods.vsurf.pred.fmla <- as.formula(paste("Mods ~ ", paste(phyto.site.mods.vsurf.pred.names, collapse = "+")))

### Minimal depths of maximal subtrees 

# Use randomForestSRC to perform minimal depth variable selection
phyto.site.mods.md <- var.select(phyto.site.mods.all.fmla,
                                 data = phyto.site.mods,
                                 ntree = 10000,
                                 method = "md",
                                 splitrule = "gini",
                                 nsplit = 0,
                                 importance = "permute",
                                 na.action = "na.impute",
                                 seed = 100)

# Create data frame of maximal subtree information
phyto.site.mods.md.order <- as.data.frame(phyto.site.mods.md[["md.obj"]][["order"]])
phyto.site.mods.md.order$Variable <- row.names(phyto.site.mods.md.order)

# Find minimal depth threshold
phyto.site.mods.md$md.obj$threshold

# Plot minimal depths and selection threshold
phyto.site.mods.md.order %>%
  ggplot(aes(x = V1, y = reorder(Variable, -V1))) +
  geom_vline(xintercept = phyto.site.mods.md[["md.obj"]][["threshold"]], linetype = "dashed", 
             color = "grey", size = 0.5) +
  labs(y = "Variables", x = "Minimal depth") +
  geom_point() +
  theme_test()

# Create vector of variables selected by minimal depth
phyto.site.mods.md.topvars <- phyto.site.mods.md$topvars

# Create formula with variables selected by minimal depth
phyto.site.mods.md.topvars.fmla <- as.formula(paste("Mods ~ ", paste(phyto.site.mods.md.topvars, collapse = "+")))

### Grow random forest 

# Generate random forest using variables selected by minimal depth
site.mod.forest <- rfsrc(phyto.site.mods.md.topvars.fmla,
                         data = phyto.site.mods,
                         ntree = 10000, #defaults to 500; determines number of independent decisions trees to generate
                         splitrule = "gini", #defaults to Gini index
                         nsplit = 0, #0 random splits for deterministic splitting; defaults to 10 for each candidate splitting variable to increase speed
                         importance = "permute", #variable importance calculated using the permute method for imbalanced settings
                         na.action = "na.impute", #na.impute permits imputing missing data
                         seed = 100, #set seed for random number generator
                         block.size = 1) #defaults to 10; determines how often cumulative error rate is estimated and block size for variable importance

                         #mtry = NULL, #defaults to sqrt(p); determines number of variables randomly selected as candidates for splitting a node
                         #nodesize = NULL, #defaults to 1 for classification; determines minimum size of terminal node
                         #nodedepth = NULL, #default to no maximum depth; determines depth to which a tree should be grown 
                         #block.size = 10, #defaults to 10; determines how often cumulative error rate is estimated and block size for variable importance
                         #ensemble = "all", #defaults to all, returning both out-of-bag and inbag ensembles
                         #bootstrap = "by.root", #defaults to by.root for bootstrapping
                         #samptype = "swor", #defaults to swor for sampling without replacement (out-of-sample)
                         #samp = NULL, #bootstrap specification if specified by user 
                         #membership = FALSE, #defaults to FALSE, specifying terminal node membership and inbag information not returned
                         #sampsize = if (samptype == "swor") function(x){x * .632} else function(x){x}, #defaults to 0.632 times the sample size used for training
                         #nimpute = 1, #defaults to 1 iteration of missing data algorithm
                         #proximity = FALSE, #defaults to FALSE, specifying proximity of cases not returned
                         #distance = FALSE, #defaults to FALSE, specifying distance between cases not returned
                         #forest.wt = FALSE, #defaults to FALSE, specifying forest weight matrix not returned
                         #xvar.wt = NULL, #defaults to uniform weights; determines probability of selecting a variable for splitting
                         #split.wt = NULL, #defaults to uniform weights; large values encourage nodes to split on a specific variable
                         #case.wt = NULL, #defaults to uniform weights; observations with larger weights we be selected with higher probability in bootstrap samples
                         #forest = TRUE, #defaults to TRUE, specifying that the forest object should be returned
                         #var.used = FALSE, #defaults to FALSE, specifying that the number of times a split occurred on a variable is not returned
                         #split.depth = FALSE, #defaults to FALSE, specifying that the minimal depth of each variable is not returned
                         #do.trace = FALSE, #defaults to FALSE, specifying number of seconds between updates not returned
                         #statistics = FALSE) #defaults to FALSE, specifying split statistics not returned

print(site.mod.forest)

# Create data frame of error rate as a function of number of trees
site.mod.forest.err.rate <- as.data.frame(site.mod.forest$err.rate)
site.mod.forest.err.rate$Tree <- rownames(site.mod.forest.err.rate)
site.mod.forest.err.rate$Tree <- as.numeric(site.mod.forest.err.rate$Tree)

# Pivot data from wide to long
site.mod.forest.err.rate.long <- pivot_longer(site.mod.forest.err.rate, !Tree, names_to = "Mod", values_to = "Error.rate")

# Plot error rate as a function of number of trees
library(ggplot2)

site.mod.forest.err.rate.long %>%
  filter(Mod == "all") %>%
  ggplot(aes(x = Tree, y = Error.rate, colour = factor(Mod))) +
  geom_path() +
  theme_bw() +
  theme(legend.position='none') +
  facet_grid(cols = vars(Mod))

site.mod.forest.err.rate.long %>%
  filter(Mod != "all") %>%
  ggplot(aes(x = Tree, y = Error.rate, colour = factor(Mod))) +
  geom_path() +
  theme_bw() +
  theme(legend.position='none') +
  facet_grid(rows = vars(Mod))

# Select a single tree from the random forest
plot(get.tree.rfsrc(site.mod.forest, tree.id = 8000, class.type = "prob"))

### Variable dependence

# Create data frame of variable importance estimates
site.mod.forest.vimp <- as.data.frame(site.mod.forest$importance)
site.mod.forest.vimp$Variable <- row.names(site.mod.forest.vimp)

# Order variables (and relevel) by unconditional importance
site.mod.forest.vimp$Variable <-  as.factor(site.mod.forest.vimp$Variable)
site.mod.forest.vimp <- site.mod.forest.vimp[order(site.mod.forest.vimp$all), ]
order <- as.character(site.mod.forest.vimp$Variable)
site.mod.forest.vimp$Variable <- fct_relevel(site.mod.forest.vimp$Variable, order)

# Pivot data from wide to long
site.mod.forest.vimp.long <- pivot_longer(site.mod.forest.vimp, !Variable, names_to = "Mod", values_to = "Vimp")

# Plot variable importance based on permutation (noising up the data)
library(ggplot2)

# Overall
site.mod.forest.vimp.long %>%
  filter(Mod == "all") %>%
  ggplot(aes(x = Vimp, y = Variable, fill = Mod)) +
  geom_col() +
  labs(x = "Importance", y = "Variable", fill = "Module") +
  theme_bw() + 
  theme(legend.position='none') +
  facet_grid(cols = vars(Mod))

# And for individual classes
site.mod.forest.vimp.long %>%
  filter(Mod != "all") %>%
  ggplot(aes(x = Vimp, y = Variable, fill = Mod)) +
  geom_col() +
  labs(x = "Importance", y = "Variable", fill = "Module") +
  theme_bw() + 
  theme(legend.position='none') +
  facet_grid(cols = vars(Mod))

#site.mod.forest.vimp.long %>%
#  filter(Mod != "all") %>%
#  ggplot(aes(x = Vimp, y = reorder(Variable, Vimp), fill = Mod)) +
#  geom_col() +
#  labs(x = "Importance", y = "Variable", fill = "Module") +
#  theme_bw() + 
#  theme(legend.position='none')

# Generate marginal plots showing the total importance of each predictor for module 02
plot.variable.rfsrc(site.mod.forest, partial = FALSE, target = "Mod02")

# Generate partial plots showing the unique importance of each predictor for module 02
plot.variable.rfsrc(site.mod.forest, partial = TRUE, target = "Mod02", npts = nrow(phyto.site.mods))

### Finding variable interactions

# Finding interacting variables based on joint permutation-based variable importance
site.mod.forest.int.vimp <- find.interaction(site.mod.forest, method = "vimp", importance = "permute", 
                                             na.action = "na.impute", nrep = 20, seed = 100)

site.mod.forest.int.vimp.df <- as.data.frame(site.mod.forest.int.vimp)

head(site.mod.forest.int.vimp.df %>% arrange(-abs(Difference)), 5)

# Finding interacting variables using conditional minimal depths
site.mod.forest.int.mt <- find.interaction(site.mod.forest, method = "maxsubtree", 
                                           na.action = "na.impute", seed = 100)

as.data.frame(site.mod.forest.int.mt) %>% arrange()

### Coplots

#library(ggRandomForests)
#surface.water.temperature.pts <- quantile_pts(site.mod.forest$xvar$Surface.water.temperature, groups=6, intervals=TRUE)
#swt.grp <- cut(site.mod.forest$xvar$Surface.water.temperature, breaks = surface.water.temperature.pts)
#sum.buoyancy.frequency.pts <- quantile_pts(site.mod.forest$xvar$Sum.buoyancy.frequency, groups=6, intervals=TRUE)
#sbf.grp <- cut(site.mod.forest$xvar$Sum.buoyancy.frequency, breaks = sum.buoyancy.frequency.pts)

# Load library
library(classInt)

# Find 6 classes of sum buoyancy frequency with approximately equal size (number of observations)
sum.buoyancy.frequency.classes <- classIntervals(site.mod.forest$xvar$Sum.buoyancy.frequency,  
                                                 n = 6, style = "quantile", include.lowest = T)

# Create sum buoyancy frequency interval classes 
sbf.grp <- cut(site.mod.forest$xvar$Sum.buoyancy.frequency, breaks = sum.buoyancy.frequency.classes$brks, 
               include.lowest = T, ordered_result = TRUE)

# Obtain predicted out-of-bag probabilities for module 02
probmod02 <- as.data.frame(site.mod.forest$predicted.oob[, 2])
probmod02 <- mutate(probmod02, as.data.frame(site.mod.forest$xvar))
colnames(probmod02)[1] <- "Probability.Mod02"

# Add sum buoyancy frequency interval classes to predicted probability data
probmod02$sbf.grp <- sbf.grp

# Rename group classes
levels(probmod02$sbf.grp) <- paste("Buoyancy ", levels(probmod02$sbf.grp), sep="")

# Generate conditioning plots showing marginal predicted probablity of module 02 as a function of surface water temperature in different classes of buoyancy frequency
ggplot(probmod02, aes(x = Surface.water.temperature, y = Probability.Mod02, color = sbf.grp)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Surface water temperature °C", y = "Probability Mod02") +
  theme_bw() + 
  theme(legend.position='none') +
  facet_wrap(vars(sbf.grp))

#########

### species

# Obtain simple correlation matrix from environmental candidate variables
phyto.sp.mods.cor <- cor(traits)

# Obtain a quick plot using corrplot with order of variables based on first principal component
phyto.sp.mods.cor.plot <- corrplot(phyto.sp.mods.cor, method = 'circle', 
                                   type = "upper", order = 'FPC', diag = FALSE)

# Create formula with all possible variables
sp.names <- colnames(traits)
phyto.sp.mods.all.fmla <- as.formula(paste("Mods ~ ", paste(sp.names, collapse = "+")))

# Use randomForestSRC to perform minimal depth variable selection
phyto.sp.mods.md <- var.select(phyto.sp.mods.all.fmla,
                               data = phyto.sp.mods,
                               ntree = 10000,
                               method = "md",
                               splitrule = "gini",
                               nsplit = 0,
                               importance = "permute",
                               na.action = "na.impute",
                               seed = 100)

# Create data frame of maximal subtree information
phyto.sp.mods.md.order <- as.data.frame(phyto.sp.mods.md[["md.obj"]][["order"]])
phyto.sp.mods.md.order$Variable <- row.names(phyto.sp.mods.md.order)

# Find minimal depth threshold
phyto.sp.mods.md$md.obj$threshold

# Plot minimal depths and selection threshold
phyto.sp.mods.md.order %>%
  ggplot(aes(x = V1, y = reorder(Variable, -V1))) +
  geom_vline(xintercept = phyto.sp.mods.md[["md.obj"]][["threshold"]], linetype = "dashed", 
             color = "grey", size = 0.5) +
  labs(y = "Variables", x = "Minimal depth") +
  geom_point() +
  theme_test()

# Create vector of variables selected by minimal depth
phyto.sp.mods.md.topvars <- phyto.sp.mods.md$topvars

# Create formula with variables selected by minimal depth
phyto.sp.mods.md.topvars.fmla <- as.formula(paste("Mods ~ ", paste(phyto.sp.mods.md.topvars, collapse = "+")))

### Grow random forest 

# Generate random forest using variables selected by minimal depth
sp.mod.forest <- rfsrc(phyto.sp.mods.md.topvars.fmla,
                       data = phyto.sp.mods,
                       ntree = 10000,
                       splitrule = "gini",
                       nsplit = 0,
                       importance = "permute",
                       na.action = "na.impute",
                       seed = 100,
                       block.size = 1)

print(sp.mod.forest)

# Create data frame of error rate as a function of number of trees
sp.mod.forest.err.rate <- as.data.frame(sp.mod.forest$err.rate)
sp.mod.forest.err.rate$Tree <- rownames(sp.mod.forest.err.rate)
sp.mod.forest.err.rate$Tree <- as.numeric(sp.mod.forest.err.rate$Tree)

# Pivot data from wide to long
sp.mod.forest.err.rate.long <- pivot_longer(sp.mod.forest.err.rate, !Tree, names_to = "Mod", values_to = "Error.rate")

# Plot error rate as a function of number of trees
library(ggplot2)

sp.mod.forest.err.rate.long %>%
  filter(Mod == "all") %>%
  ggplot(aes(x = Tree, y = Error.rate, colour = factor(Mod))) +
  geom_path() +
  theme_bw() +
  theme(legend.position='none') +
  facet_grid(cols = vars(Mod))

sp.mod.forest.err.rate.long %>%
  filter(Mod != "all") %>%
  ggplot(aes(x = Tree, y = Error.rate, colour = factor(Mod))) +
  geom_path() +
  theme_bw() +
  theme(legend.position='none') +
  facet_grid(rows = vars(Mod))

# Select a single tree from the random forest
plot(get.tree.rfsrc(sp.mod.forest, tree.id = 5000, class.type = "prob"))

### Variable dependence

# Create data frame of variable importance estimates
sp.mod.forest.vimp <- as.data.frame(sp.mod.forest$importance)
sp.mod.forest.vimp$Variable <- row.names(sp.mod.forest.vimp)

# Order variables (and relevel) by unconditional importance
sp.mod.forest.vimp$Variable <-  as.factor(sp.mod.forest.vimp$Variable)
sp.mod.forest.vimp <- sp.mod.forest.vimp[order(sp.mod.forest.vimp$all), ]
order <- as.character(sp.mod.forest.vimp$Variable)
sp.mod.forest.vimp$Variable <- fct_relevel(sp.mod.forest.vimp$Variable, order)

# Pivot data from wide to long
sp.mod.forest.vimp.long <- pivot_longer(sp.mod.forest.vimp, !Variable, names_to = "Mod", values_to = "Vimp")

# Plot variable importance based on permutation (noising up the data)
library(ggplot2)

# Overall
sp.mod.forest.vimp.long %>%
  filter(Mod == "all") %>%
  ggplot(aes(x = Vimp, y = Variable, fill = Mod)) +
  geom_col() +
  labs(x = "Importance", y = "Variable", fill = "Module") +
  theme_bw() + 
  theme(legend.position='none') +
  facet_grid(cols = vars(Mod))

# And for individual classes
sp.mod.forest.vimp.long %>%
  filter(Mod != "all") %>%
  ggplot(aes(x = Vimp, y = Variable, fill = Mod)) +
  geom_col() +
  labs(x = "Importance", y = "Variable", fill = "Module") +
  theme_bw() + 
  theme(legend.position='none') +
  facet_grid(cols = vars(Mod))

# Generate marginal plots showing the total importance of each predictor for module 06
plot.variable.rfsrc(sp.mod.forest, partial = FALSE, target = "Mod06")

# Generate partial plots showing the unique importance of each predictor for module 06
plot.variable.rfsrc(sp.mod.forest, partial = TRUE, target = "Mod06", npts = nrow(phyto.sp.mods))

### Finding variable interactions

# Finding interacting variables based on joint permutation-based variable importance
sp.mod.forest.int.vimp <- find.interaction(sp.mod.forest, method = "vimp", importance = "permute",
                                           na.action = "na.impute", nrep = 20, seed = 100)

sp.mod.forest.int.vimp.df <- as.data.frame(sp.mod.forest.int.vimp)

sp.mod.forest.int.vimp.df %>% arrange(-abs(Difference))

# Finding interacting variables using conditional minimal depths
sp.mod.forest.int.mt <- find.interaction(sp.mod.forest, method = "maxsubtree",
                                         na.action = "na.impute", seed = 100)

as.data.frame(sp.mod.forest.int.mt) %>% arrange()

# Obtain predicted out-of-bag probabilities for module 06
probmod06 <- as.data.frame(sp.mod.forest$predicted.oob[, 2])
probmod06 <- mutate(probmod06, as.data.frame(sp.mod.forest$xvar))
colnames(probmod06)[1] <- "Probability.Mod06"

# Rename group classes
probmod06["SUSPECTED_TOXIN_PRODUCING"][probmod06["SUSPECTED_TOXIN_PRODUCING"] == "1"] <- "Toxin producers"
probmod06["SUSPECTED_TOXIN_PRODUCING"][probmod06["SUSPECTED_TOXIN_PRODUCING"] == "0"] <- "Not toxin producers"

# Generate conditioning plots showing marginal predicted probablity of module 06 as a function of greatest axial dimension in suspected toxin producers and non-producers
ggplot(probmod06, aes(x = GALD, y = Probability.Mod06, color = SUSPECTED_TOXIN_PRODUCING)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Greatest axial linear dimension", y = "Probability mod06") +
  theme_bw() + 
  theme(legend.position='none') +
  facet_wrap(vars(SUSPECTED_TOXIN_PRODUCING))

########



# load library
library(ggRandomForests)

# Use ggRandomForests package to obtain partial predicted probabilities for each subgroup
partial_coplot_pbc <- gg_partial_coplot(site.mod.forest, xvar = "Surface.water.temperature",
                                          groups = sbf.grp,
                                          show.plots = FALSE)

plot(partial_coplot_pbc)

# Generate conditioning plots showing partial predicted probablity of module 02 as a function of surface water temperature in different classes of buoyancy frequency
ggplot(partial_coplot_pbc, aes(x = Surface.water.temperature, y = yhat, color = group)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Surface water temperature °C", y = "Probability Mod02") +
  theme_bw() + 
  theme(legend.position='none') +
  facet_wrap(vars(group))




test <- partial(site.mod.forest, partial.xvar = "Surface.water.temperature",
                partial.values = site.mod.forest$xvar$Surface.water.temperature)

pdta1 <- get.partial.plot.data(test, target = "Surface.water.temperature")

plot(gg_interaction(site.mod.forest.int.mt), xvar=phyto.site.mods.md.topvars, panel=TRUE)




#




     

write.csv(one.mv.obj.vimp.Corrected.Richness, file = "one.mv.obj.vimp.Corrected.Richness.csv")
