lapply(c("vegan",
  "ggplot2",
  "data.table",
  "ecolFudge",
  "betapart",
  "modEvA",
  "seqinr",
  "Peptides",
  "ggparl",
  "ggridges",
  "patchwork",
  "viridis",
  "MASS"), library, character.only = T)

dat <- fread("envData.csv")

# create new sample code to match those in amino acid variant file
dat[, seqCode := gsub("-", "x", gsub("-b", "", sample))]

pmoa <- fread("pmoA_AAV.txt")
colnames(pmoa)[1] <- "AAV"

pmoa <- transDT(pmoa, transCol = "AAV", rowID = "seqCode")

aavOccs <- specnumber(pmoa[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
pmoaFilt <- pmoa[, .SD, .SDcols = c("seqCode", names(aavOccs[aavOccs > 1]))]

# get library sizes of filtered AAV table
pmoaFilt[, libSize := rowSums(pmoaFilt[, !"seqCode"])]

# merge with env metadata
pmoaDat <- merge(dat, pmoaFilt, by = "seqCode")

aavCols <- grep("AAV", colnames(pmoaDat), value = T)

# rarefy AAV table to even depth
pmoaDat[, (aavCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = aavCols]

# remove any empty AAVs
emptyAAVs <- pmoaDat[, colSums(.SD) == 0, .SDcols = aavCols]
pmoaDat[, (names(emptyAAVs[emptyAAVs == T])) := NULL]

# update remaining list of aavs
nsAavs <- aavCols[!emptyAAVs]

# create beta diversity matrix
pmoaBetaMat <- beta.pair(ifelse(pmoaDat[, .SD, .SDcols = nsAavs] != 0, 1, 0))

nmdsRes <- metaMDS(pmoaBetaMat[[3]], autotransform = F, trymax = 200)

scoreCols = c("x", "y")
pmoaDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

pmoaBfiEuclid <- vegdist(pmoaDat$bfi, "euclid")

pmoaAavBfiDat <- data.table(dissim = as.vector(pmoaBetaMat[[3]]),
  bfiSim = as.vector(pmoaBfiEuclid), type = "AAV")

# calculate pmoA AAV richness
pmoaDat[, richness := specnumber(.SD), .SDcols = nsAavs]

# test for differences in hydrophobicity
test <- read.fasta(
  "sequences/pmoA/mergedSeqs/allErrorCorrected/framebotOut/AAVs.fasta",
  seqtype = "AA", as.string = T)

testSeq <- data.table(AAV = names(test),
  seq = unlist(sapply(test, getSequence, as.string = T)))

testSeq[, ":="(hydro = hydrophobicity(seq),
  aacharge = charge(seq))]

# retain only sequences that feature in final AAV table
testSeq <- testSeq[AAV %in% nsAavs]

# sort hydrophobicity vals according to AAV table
testSeq <- testSeq[order(match(AAV, nsAavs))]

# calculated mean hydrophobicity weighted by abundance of each individual amino acid variant.
pmoaDat[, ":="(
  meanHydrophob = unlist(lapply(1:nrow(pmoaDat), function(s)
    weighted.mean(x = testSeq$hydro,
      w = pmoaDat[s, .SD, .SDcols = nsAavs]))),
  meanCharge = unlist(lapply(1:nrow(pmoaDat), function(s)
    weighted.mean(x = testSeq$aacharge,
      w = pmoaDat[s, .SD, .SDcols = nsAavs]))))]

pmoaHydro_bfi <- ggplot(pmoaDat,
    aes(x = bfi, y = meanHydrophob, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(pmoA)~hydrophobicity),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

pmoaCharge_bfi <- ggplot(pmoaDat,
    aes(x = bfi, y = meanCharge, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(pmoA)~net~charge),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

pmoaHydroLm <- lm(meanHydrophob ~ bfi, data = pmoaDat)
pmoaChargeLm <- lm(meanCharge ~ bfi, data = pmoaDat)

### otu analysis
otus <- fread("pmoAOtuTab.txt")
colnames(otus)[1] <- "OTU"

pmoaOtus <- transDT(otus, transCol = "OTU", rowID = "seqCode")

otuOccs <- specnumber(pmoaOtus[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
pmoaFiltOtus <- pmoaOtus[, .SD,
  .SDcols = c("seqCode", names(otuOccs[otuOccs > 1]))]

# get library sizes of filtered AAV table
pmoaFiltOtus[, libSize := rowSums(pmoaFiltOtus[, !"seqCode"])]

# merge with env metadata
pmoaOtuDat <- merge(dat, pmoaFiltOtus, by = "seqCode")

otuCols <- grep("OTU", colnames(pmoaOtuDat), value = T)

# rarefy AAV table to even depth
pmoaOtuDat[, (otuCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = otuCols]

# remove any empty OTUs
emptyOTUs <- pmoaOtuDat[, colSums(.SD) == 0, .SDcols = otuCols]
pmoaOtuDat[, (names(emptyOTUs[emptyOTUs == T])) := NULL]

# remove empty OTUs from column list
otuCols <- otuCols[! emptyOTUs]

# create beta diversity matrix
pmoaBetaOtuMat <- beta.pair(
  ifelse(pmoaOtuDat[, .SD, .SDcols = otuCols] != 0, 1, 0))

nmdsRes <- metaMDS(betaOtuMat[[3]], autotransform = F, trymax = 200)

scoreCols <- c("x", "y")
pmoaOtuDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

pmoaBfiOtuEuclid <- vegdist(pmoaOtuDat$bfi, "euclid")

pmoaOtuBfiDat <- data.table(dissim = as.vector(pmoaBetaOtuMat[[3]]),
  bfiSim = as.vector(pmoaBfiOtuEuclid), type = "OTU")

# pmoa negative exponential decay model
pmoaAavDecay <- decay.model(pmoaBetaMat[[3]], pmoaBfiEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)
pmoaOtuDecay <- decay.model(pmoaBetaOtuMat[[3]], pmoaBfiOtuEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)

# function to calculate prediction curve for range of data
calcDecayFit <- function(x){
  xCoord <- seq(min(x$data[, 1]), max(x$data[, 1]),  0.005)
  pred <- 1 - (1 - x$a.intercept) * exp(-x$b.slope * xCoord)
  decayFit <- data.table(x = xCoord, y = pred)
  return(decayFit)
}

# combine the two dissim data frames, then calculate fits for both decay models
pmoaDists <- rbindlist(list(pmoaAavBfiDat, pmoaOtuBfiDat))
pmoaAavFit <- calcDecayFit(pmoaAavDecay)
pmoaOtuFit <- calcDecayFit(pmoaOtuDecay)

decayLabel <- function(x){
  paste0("list(italic(R^2) ==", round(x$pseudo.r.squared, 2), ", italic('P')",
  ifelse(x$p.value >= 0.05, "==round(x$p.value, 2)",
    ifelse(x$p.value < 0.05 & x$p.value >= 0.01, "< 0.05", "< 0.01")), ")")
}

# bootstrap coefficients
pmoaAavBoot <- boot.coefs.decay(pmoaAavDecay, 1000)
pmoaOtuBoot <- boot.coefs.decay(pmoaOtuDecay, 1000)

combineBoots <- function(gene, aavFit, otuFit){
  combnBoots <- data.table(
    gene = gene,
    type = rep(c("AAV", "OTU"), each = nrow(aavFit$boot.coefs)),
    coefs = c(aavFit$boot.coefs[, 2], otuFit$boot.coefs[, 2]))
  return(combnBoots)
}

# combine all bootstrap estimates
pmoaBoot <- combineBoots(gene = "pmoA", pmoaAavBoot,
  pmoaOtuBoot)

# generate 2 cols from viridis for otu vs aav comparisons
darkCol <- "#DE7A22"
lightCol <- "#20948B"

# calculate pmoA OTU richness
pmoaOtuDat[, richness := specnumber(.SD), .SDcols = otuCols]

pmoaRichness <- rbindlist(list(AAV = pmoaDat[, .(bfi, geol, richness, month)],
  OTU = pmoaOtuDat[, .(bfi, geol, richness, month)]), idcol = "type")

pmoaAavRichness <- glm.nb(richness ~ bfi, data = pmoaRichness[type == "AAV"])
pmoaAavMonthRichness <- glm.nb(richness ~ month + bfi,
  data = pmoaRichness[type == "AAV"])
pmoaOtuRichness <- glm.nb(richness ~ bfi, data = pmoaRichness[type == "OTU"])
pmoaOtuMonthRichness <- glm.nb(richness ~ month + bfi,
  data = pmoaRichness[type == "OTU"])

lapply(list(
  pmoaAavRichness, pmoaAavMonthRichness, pmoaOtuRichness, pmoaOtuMonthRichness),
  summary)
lapply(list(
  pmoaAavRichness, pmoaAavMonthRichness, pmoaOtuRichness, pmoaOtuMonthRichness),
  Dsquared, adjust = T)

lapply(list(
  pmoaAavRichness, pmoaAavMonthRichness, pmoaOtuRichness, pmoaOtuMonthRichness), function(x) exp(coef(x)))
############################ MCRA analysis #####################################

mcra <- fread("mcrA_AAV.txt")
colnames(mcra)[1] <- "AAV"

mcra <- transDT(mcra, transCol = "AAV", rowID = "seqCode")

aavOccs <- specnumber(mcra[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
mcraFilt <- mcra[, .SD, .SDcols = c("seqCode", names(aavOccs[aavOccs > 1]))]

# get library sizes of filtered AAV table
mcraFilt[, libSize := rowSums(mcraFilt[, !"seqCode"])]

# merge with env metadata
mcraDat <- merge(dat, mcraFilt, by = "seqCode")

aavCols <- grep("AAV", colnames(mcraDat), value = T)

# rarefy AAV table to even depth
mcraDat[, (aavCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = aavCols]

# remove any empty AAVs
emptyAAVs <- mcraDat[, colSums(.SD) == 0, .SDcols = aavCols]
mcraDat[, (names(emptyAAVs[emptyAAVs == T])) := NULL]

# update remaining list of aavs
nsAavs <- aavCols[!emptyAAVs]

# create beta diversity matrix
mcraBetaMat <- beta.pair(ifelse(mcraDat[, .SD, .SDcols = nsAavs] != 0, 1, 0))

nmdsRes <- metaMDS(mcraBetaMat[[3]], autotransform = F, trymax = 200)

scoreCols = c("x", "y")
mcraDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

mcraBfiEuclid <- vegdist(mcraDat$bfi, "euclid")

mcraAavBfiDat <- data.table(dissim = as.vector(mcraBetaMat[[3]]),
  bfiSim = as.vector(mcraBfiEuclid), type = "AAV")

# calculate mcra AAV richness
mcraDat[, richness := specnumber(.SD), .SDcols = nsAavs]

# test for differences in hydrophobicity
test <- read.fasta(
  "sequences/mcrA/mergedSeqs/allErrorCorrected/framebotOut/AAVs.fasta",
  seqtype = "AA", as.string = T)

testSeq <- data.table(AAV = names(test),
  seq = unlist(sapply(test, getSequence, as.string = T)))

testSeq[, ":="(hydro = hydrophobicity(seq),
  aacharge = charge(seq))]

# retain only sequences that feature in final AAV table
testSeq <- testSeq[AAV %in% nsAavs]

# sort hydrophobicity vals according to AAV table
testSeq <- testSeq[order(match(AAV, nsAavs))]

# calculated mean hydrophobicity weighted by abundance of each individual amino acid variant.
mcraDat[, ":="(
  meanHydrophob = unlist(lapply(1:nrow(mcraDat), function(s)
    weighted.mean(x = testSeq$hydro,
      w = mcraDat[s, .SD, .SDcols = nsAavs]))),
  meanCharge = unlist(lapply(1:nrow(mcraDat), function(s)
    weighted.mean(x = testSeq$aacharge,
      w = mcraDat[s, .SD, .SDcols = nsAavs]))))]

mcraHydro_bfi <- ggplot(mcraDat,
    aes(x = bfi, y = meanHydrophob, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(mcrA)~hydrophobicity),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

mcraCharge_bfi <- ggplot(mcraDat,
    aes(x = bfi, y = meanCharge, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(mcrA)~net~charge),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

mcraHydroLm <- lm(meanHydrophob ~ bfi, data = mcraDat)
mcraChargeLm <- lm(meanCharge ~ bfi, data = mcraDat)
# compare to dna OTUS
otus <- fread("mcrAOtuTab.txt")
colnames(otus)[1] <- "OTU"

mcraOtus <- transDT(otus, transCol = "OTU", rowID = "seqCode")

otuOccs <- specnumber(mcraOtus[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
mcraFiltOtus <- mcraOtus[, .SD,
  .SDcols = c("seqCode", names(otuOccs[otuOccs > 1]))]

# get library sizes of filtered AAV table
mcraFiltOtus[, libSize := rowSums(mcraFiltOtus[, !"seqCode"])]

# merge with env metadata
mcraOtuDat <- merge(dat, mcraFiltOtus, by = "seqCode")

otuCols <- grep("OTU", colnames(mcraOtuDat), value = T)

# rarefy AAV table to even depth
mcraOtuDat[, (otuCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = otuCols]

# remove any empty Otus
emptyOTUs <- mcraOtuDat[, colSums(.SD) == 0, .SDcols = otuCols]
mcraOtuDat[, (names(emptyOTUs[emptyOTUs == T])) := NULL]

otuCols <- otuCols[! emptyOTUs]

# create beta diversity matrix
mcraBetaOtuMat <- beta.pair(
  ifelse(mcraOtuDat[, .SD, .SDcols = otuCols] != 0, 1, 0))

nmdsRes <- metaMDS(mcraBetaOtuMat[[3]], autotransform = F, trymax = 200)

scoreCols <- c("x", "y")
mcraOtuDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

mcraBfiOtuEuclid <- vegdist(mcraOtuDat$bfi, "euclid")

mcraOtuBfiDat <- data.table(dissim = as.vector(mcraBetaOtuMat[[3]]),
  bfiSim = as.vector(mcraBfiOtuEuclid), type = "OTU")

# mcra negative exponential decay model
mcraAavDecay <- decay.model(mcraBetaMat[[3]], mcraBfiEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)
mcraOtuDecay <- decay.model(mcraBetaOtuMat[[3]], mcraBfiOtuEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)

# combine the two dissim data frames, then calculate fits for both decay models
mcraDists <- rbindlist(list(mcraAavBfiDat, mcraOtuBfiDat))
mcraAavFit <- calcDecayFit(mcraAavDecay)
mcraOtuFit <- calcDecayFit(mcraOtuDecay)

# bootstrap coefficients and combine
mcraAavBoot <- boot.coefs.decay(mcraAavDecay, 1000)
mcraOtuBoot <- boot.coefs.decay(mcraOtuDecay, 1000)

mcraBoot <- combineBoots(gene = "mcrA", mcraAavBoot, mcraOtuBoot)

mcraOtuDat[, richness := specnumber(.SD), .SDcols = otuCols]

mcraRichness <- rbindlist(list(AAV = mcraDat[, .(bfi, geol, richness, month)],
  OTU = mcraOtuDat[, .(bfi, geol, richness, month)]), idcol = "type")

mcraAavRichness <- glm.nb(richness ~ bfi, data = mcraRichness[type == "AAV"])
mcraAavMonthRichness <- glm.nb(richness ~ month + bfi,
  data = mcraRichness[type == "AAV"])
mcraOtuRichness <- glm.nb(richness ~ bfi, data = mcraRichness[type == "OTU"])
mcraOtuMonthRichness <- glm.nb(richness ~ month + bfi,
  data = mcraRichness[type == "OTU"])
lapply(list(
  mcraAavRichness, mcraAavMonthRichness, mcraOtuRichness, mcraOtuMonthRichness),
  summary)
lapply(list(
  mcraAavRichness, mcraAavMonthRichness, mcraOtuRichness, mcraOtuMonthRichness),
  Dsquared, adjust = T)
lapply(list(
  mcraAavRichness, mcraAavMonthRichness, mcraOtuRichness, mcraOtuMonthRichness),
  function(x) exp(coef(x)))

################################## nirS analysis ###############################
nirs <- fread("nirS_AAV.txt")
colnames(nirs)[1] <- "AAV"

nirs <- transDT(nirs, transCol = "AAV", rowID = "seqCode")

aavOccs <- specnumber(nirs[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
nirsFilt <- nirs[, .SD, .SDcols = c("seqCode", names(aavOccs[aavOccs > 1]))]

# get library sizes of filtered AAV table
nirsFilt[, libSize := rowSums(nirsFilt[, !"seqCode"])]

# merge with env metadata
nirsDat <- merge(dat, nirsFilt, by = "seqCode")

aavCols <- grep("AAV", colnames(nirsDat), value = T)

# rarefy AAV table to even depth
nirsDat[, (aavCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = aavCols]

# remove any empty AAVs
emptyAAVs <- nirsDat[, colSums(.SD) == 0, .SDcols = aavCols]
nirsDat[, (names(emptyAAVs[emptyAAVs == T])) := NULL]

# update remaining list of aavs
nsAavs <- aavCols[!emptyAAVs]

# create beta diversity matrix
nirsBetaMat <- beta.pair(ifelse(nirsDat[, .SD, .SDcols = nsAavs] != 0, 1, 0))

nmdsRes <- metaMDS(nirsBetaMat[[3]], autotransform = F, trymax = 200)

scoreCols = c("x", "y")
nirsDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

nirsBfiEuclid <- vegdist(nirsDat$bfi, "euclid")

nirsAavBfiDat <- data.table(dissim = as.vector(nirsBetaMat[[3]]),
  bfiSim = as.vector(nirsBfiEuclid), type = "AAV")

nirsDat[, richness := specnumber(.SD), .SDcols = nsAavs]

# test for differences in hydrophobicity
test <- read.fasta(
  "sequences/nirS/mergedSeqs/allErrorCorrected/framebotOut/AAVs.fasta",
  seqtype = "AA", as.string = T)

testSeq <- data.table(AAV = names(test),
  seq = unlist(sapply(test, getSequence, as.string = T)))

testSeq[, ":="(hydro = hydrophobicity(seq),
  aacharge = charge(seq))]

# retain only sequences that feature in final AAV table
testSeq <- testSeq[AAV %in% nsAavs]

# sort hydrophobicity vals according to AAV table
testSeq <- testSeq[order(match(AAV, nsAavs))]

# calculated mean hydrophobicity weighted by abundance of each individual amino acid variant.
nirsDat[, ":="(
  meanHydrophob = unlist(lapply(1:nrow(nirsDat), function(s)
    weighted.mean(x = testSeq$hydro,
      w = nirsDat[s, .SD, .SDcols = nsAavs]))),
  meanCharge = unlist(lapply(1:nrow(nirsDat), function(s)
    weighted.mean(x = testSeq$aacharge,
      w = nirsDat[s, .SD, .SDcols = nsAavs]))))]

nirsHydro_bfi <- ggplot(nirsDat,
    aes(x = bfi, y = meanHydrophob, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(nirS)~hydrophobicity),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

nirsCharge_bfi <- ggplot(nirsDat,
    aes(x = bfi, y = meanCharge, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(nirS)~net~charge),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

nirsHydroLm <- lm(meanHydrophob ~ bfi, data = nirsDat)
nirsChargeLm <- lm(meanCharge ~ bfi, data = nirsDat)

ggsave("../figures/nirS_hydrophobicity.pdf", nirsHydro_bfi, height = 6,
  width = 7, device = "pdf")

otus <- fread("nirSOtuTab.txt")
colnames(otus)[1] <- "OTU"

nirsOtus <- transDT(otus, transCol = "OTU", rowID = "seqCode")

otuOccs <- specnumber(nirsOtus[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
nirsFiltOtus <- nirsOtus[, .SD,
  .SDcols = c("seqCode", names(otuOccs[otuOccs > 1]))]

# get library sizes of filtered AAV table
nirsFiltOtus[, libSize := rowSums(nirsFiltOtus[, !"seqCode"])]

# merge with env metadata
nirsOtuDat <- merge(dat, nirsFiltOtus, by = "seqCode")

otuCols <- grep("OTU", colnames(nirsOtuDat), value = T)

# rarefy AAV table to even depth
nirsOtuDat[, (otuCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = otuCols]

# remove any empty OTUs
emptyOTUs <- nirsOtuDat[, colSums(.SD) == 0, .SDcols = otuCols]
nirsOtuDat[, (names(emptyOTUs[emptyOTUs == T])) := NULL]

# remove empty OTUs from column list
otuCols <- otuCols[! emptyOTUs]

# create beta diversity matrix
nirsBetaOtuMat <- beta.pair(
  ifelse(nirsOtuDat[, .SD, .SDcols = otuCols] != 0, 1, 0))

nmdsRes <- metaMDS(nirsBetaOtuMat[[3]], autotransform = F, trymax = 200)

scoreCols <- c("x", "y")
nirsOtuDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

nirsBfiOtuEuclid <- vegdist(nirsOtuDat$bfi, "euclid")

nirsOtuBfiDat <- data.table(dissim = as.vector(nirsBetaOtuMat[[3]]),
  bfiSim = as.vector(nirsBfiOtuEuclid), type = "OTU")

# mcra negative exponential decay model
nirsAavDecay <- decay.model(nirsBetaMat[[3]], nirsBfiEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)
nirsOtuDecay <- decay.model(nirsBetaOtuMat[[3]], nirsBfiOtuEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)

# combine the two dissim data frames, then calculate fits for both decay models
nirsDists <- rbindlist(list(nirsAavBfiDat, nirsOtuBfiDat))
nirsAavFit <- calcDecayFit(nirsAavDecay)
nirsOtuFit <- calcDecayFit(nirsOtuDecay)

# calculate bootstraps and combine
nirsAavBoot <- boot.coefs.decay(nirsAavDecay, 1000)
nirsOtuBoot <- boot.coefs.decay(nirsOtuDecay, 1000)

nirsBoot <- combineBoots(gene = "nirS", nirsAavBoot, nirsOtuBoot)

nirsOtuDat[, richness := specnumber(.SD), .SDcols = otuCols]

nirsRichness <- rbindlist(list(AAV = nirsDat[, .(bfi, geol, richness, month)],
  OTU = nirsOtuDat[, .(bfi, geol, richness, month)]), idcol = "type")

nirsAavRichness <- glm.nb(richness ~ bfi, data = nirsRichness[type == "AAV"])
nirsOtuRichness <- glm.nb(richness ~ bfi, data = nirsRichness[type == "OTU"])
nirsAavMonthRichness <- glm.nb(richness ~ month + bfi,
  data = nirsRichness[type == "AAV"])
nirsOtuMonthRichness <- glm.nb(richness ~ month + bfi,
  data = nirsRichness[type == "OTU"])
lapply(list(
  nirsAavRichness, nirsAavMonthRichness, nirsOtuRichness, nirsOtuMonthRichness),
  summary)
lapply(list(
  nirsAavRichness, nirsAavMonthRichness, nirsOtuRichness, nirsOtuMonthRichness),
  Dsquared, adjust = T)
lapply(list(
  nirsAavRichness, nirsAavMonthRichness, nirsOtuRichness, nirsOtuMonthRichness),
  function(x) exp(coef(x)))

################################## AOB Analysis ################################
aob <- fread("AOB_AAV.txt")
colnames(aob)[1] <- "AAV"

aob <- transDT(aob, transCol = "AAV", rowID = "seqCode")

aavOccs <- specnumber(aob[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
aobFilt <- aob[, .SD, .SDcols = c("seqCode", names(aavOccs[aavOccs > 1]))]

# get library sizes of filtered AAV table
aobFilt[, libSize := rowSums(aobFilt[, !"seqCode"])]

# merge with env metadata
aobDat <- merge(dat, aobFilt, by = "seqCode")

aavCols <- grep("AAV", colnames(aobDat), value = T)

# rarefy AAV table to even depth
aobDat[, (aavCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = aavCols]

# remove any empty AAVs
emptyAAVs <- aobDat[, colSums(.SD) == 0, .SDcols = aavCols]
aobDat[, (names(emptyAAVs[emptyAAVs == T])) := NULL]

# update remaining list of aavs
nsAavs <- aavCols[!emptyAAVs]

# create beta diversity matrix
aobBetaMat <- beta.pair(ifelse(aobDat[, .SD, .SDcols = nsAavs] != 0, 1, 0))

nmdsRes <- metaMDS(aobBetaMat[[3]], autotransform = F, trymax = 200)

scoreCols = c("x", "y")
aobDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

aobBfiEuclid <- vegdist(aobDat$bfi, "euclid")

aobAavBfiDat <- data.table(dissim = as.vector(aobBetaMat[[3]]),
  bfiSim = as.vector(aobBfiEuclid), type = "AAV")

aobDat[, richness := specnumber(.SD), .SDcols = nsAavs]

# test for differences in hydrophobicity
test <- read.fasta(
  "sequences/AOB/mergedSeqs/allErrorCorrected/framebotOut/AAVs.fasta",
  seqtype = "AA", as.string = T)

testSeq <- data.table(AAV = names(test),
  seq = unlist(sapply(test, getSequence, as.string = T)))

testSeq[, ":="(hydro = hydrophobicity(seq),
  aacharge = charge(seq))]

# retain only sequences that feature in final AAV table
testSeq <- testSeq[AAV %in% nsAavs]

# sort hydrophobicity vals according to AAV table
testSeq <- testSeq[order(match(AAV, nsAavs))]

# calculated mean hydrophobicity weighted by abundance of each individual amino acid variant.
aobDat[, ":="(
  meanHydrophob = unlist(lapply(1:nrow(aobDat), function(s)
    weighted.mean(x = testSeq$hydro,
      w = aobDat[s, .SD, .SDcols = nsAavs]))),
  meanCharge = unlist(lapply(1:nrow(aobDat), function(s)
    weighted.mean(x = testSeq$aacharge,
      w = aobDat[s, .SD, .SDcols = nsAavs]))))]

aobHydroLm <- lm(meanHydrophob ~ bfi, data = aobDat)
aobChargeLm <- lm(meanCharge ~ bfi, data = aobDat)

aobHydro_bfi <- ggplot(aobDat, aes(x = bfi, y = meanHydrophob)) +
  geom_point(aes(col = geol, shape = month), size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(aob)~hydrophobicity),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

aobCharge_bfi <- ggplot(aobDat, aes(x = bfi, y = meanCharge)) +
  geom_point(aes(col = geol, shape = month), size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(aob)~net~charge),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../figures/aob_hydrophobicity.pdf", aobHydro_bfi, height = 6,
  width = 7, device = "pdf")

# AOB OTU analysis
otus <- fread("AOB_amoAOtuTab.txt")
colnames(otus)[1] <- "OTU"

aobOtus <- transDT(otus, transCol = "OTU", rowID = "seqCode")

otuOccs <- specnumber(aobOtus[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
aobFiltOtus <- aobOtus[, .SD,
  .SDcols = c("seqCode", names(otuOccs[otuOccs > 1]))]

# get library sizes of filtered AAV table
aobFiltOtus[, libSize := rowSums(aobFiltOtus[, !"seqCode"])]

# merge with env metadata
aobOtuDat <- merge(dat, aobFiltOtus, by = "seqCode")

otuCols <- grep("OTU", colnames(aobOtuDat), value = T)

# rarefy AAV table to even depth
aobOtuDat[, (otuCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = otuCols]

# remove any empty OTUs
emptyOTUs <- aobOtuDat[, colSums(.SD) == 0, .SDcols = otuCols]
aobOtuDat[, (names(emptyOTUs[emptyOTUs == T])) := NULL]

# remove empty OTUs from column list
otuCols <- otuCols[! emptyOTUs]

# create beta diversity matrix
aobBetaOtuMat <- beta.pair(
  ifelse(aobOtuDat[, .SD, .SDcols = otuCols] != 0, 1, 0))

nmdsRes <- metaMDS(aobBetaOtuMat[[3]], autotransform = F, trymax = 200)

scoreCols <- c("x", "y")
aobOtuDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

aob_OTU_nmds <- ggplot(aobOtuDat,
    aes(x = x, y = y, col = bfi, shape = geol)) +
  geom_point(size = 4) +
  scale_color_viridis() +
  labs(x = "NMDS 1", y = "NMDS 2", col = "Base flow index", shape = "Geology",
    title = "97% OTUs") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 18),
    panel.grid = element_blank())

aobBfiOtuEuclid <- vegdist(aobOtuDat$bfi, "euclid")

aobOtuBfiDat <- data.table(dissim = as.vector(aobBetaOtuMat[[3]]),
  bfiSim = as.vector(aobBfiOtuEuclid), type = "OTU")

# mcra negative exponential decay model
aobAavDecay <- decay.model(aobBetaMat[[3]], aobBfiEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)
aobOtuDecay <- decay.model(aobBetaOtuMat[[3]], aobBfiOtuEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)

# combine the two dissim data frames, then calculate fits for both decay models
aobDists <- rbindlist(list(aobAavBfiDat, aobOtuBfiDat))
aobAavFit <- calcDecayFit(aobAavDecay)
aobOtuFit <- calcDecayFit(aobOtuDecay)

# bootstrap coefficients and combine
aobAavBoot <- boot.coefs.decay(aobAavDecay, 1000)
aobOtuBoot <- boot.coefs.decay(aobOtuDecay, 1000)

aobBoot <- combineBoots(gene = "AOB_amoA", aobAavBoot, aobOtuBoot)

aobOtuDat[, richness := specnumber(.SD), .SDcols = otuCols]

aobRichness <- rbindlist(list(AAV = aobDat[, .(bfi, geol, richness, month)],
  OTU = aobOtuDat[, .(bfi, geol, richness, month)]), idcol = "type")

aobAavRichness <- glm.nb(richness ~ bfi, data = aobRichness[type == "AAV"])
aobOtuRichness <- glm.nb(richness ~ bfi, data = aobRichness[type == "OTU"])
aobAavMonthRichness <- glm.nb(richness ~ month + bfi,
  data = aobRichness[type == "AAV"])
aobOtuMonthRichness <- glm.nb(richness ~ month + bfi,
  data = aobRichness[type == "OTU"])
lapply(list(
  aobAavRichness, aobAavMonthRichness, aobOtuRichness, aobOtuMonthRichness),
  summary)
lapply(list(
  aobAavRichness, aobAavMonthRichness, aobOtuRichness, aobOtuMonthRichness), Dsquared, adjust = T)
lapply(list(
  aobAavRichness, aobAavMonthRichness, aobOtuRichness, aobOtuMonthRichness), function(x) exp(coef(x)))

################################## AOA Analysis ################################
aoa <- fread("AOA_AAV.txt")
colnames(aoa)[1] <- "AAV"

aoa <- transDT(aoa, transCol = "AAV", rowID = "seqCode")

aavOccs <- specnumber(aoa[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
aoaFilt <- aoa[, .SD, .SDcols = c("seqCode", names(aavOccs[aavOccs > 1]))]

# get library sizes of filtered AAV table
aoaFilt[, libSize := rowSums(aoaFilt[, !"seqCode"])]

# merge with env metadata
aoaDat <- merge(dat, aoaFilt, by = "seqCode")

aavCols <- grep("AAV", colnames(aoaDat), value = T)

# rarefy AAV table to even depth
aoaDat[, (aavCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = aavCols]

# remove any empty AAVs
emptyAAVs <- aoaDat[, colSums(.SD) == 0, .SDcols = aavCols]
aoaDat[, (names(emptyAAVs[emptyAAVs == T])) := NULL]

# update remaining list of aavs
nsAavs <- aavCols[!emptyAAVs]

# create beta diversity matrix
aoaBetaMat <- beta.pair(ifelse(aoaDat[, .SD, .SDcols = nsAavs] != 0, 1, 0))

nmdsRes <- metaMDS(aoaBetaMat[[3]], autotransform = F, trymax = 200)

scoreCols = c("x", "y")
aoaDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

aoaBfiEuclid <- vegdist(aoaDat$bfi, "euclid")

aoaAavBfiDat <- data.table(dissim = as.vector(aoaBetaMat[[3]]),
  bfiSim = as.vector(aoaBfiEuclid), type = "AAV")

aoaDat[, richness := specnumber(.SD), .SDcols = nsAavs]

# test for differences in hydrophobicity
test <- read.fasta(
  "sequences/AOA/qualTrimmedSeqs/allErrorCorrected/framebotOut/AAVs.fasta",
  seqtype = "AA", as.string = T)

testSeq <- data.table(AAV = names(test),
  seq = unlist(sapply(test, getSequence, as.string = T)))

testSeq[, ":="(hydro = hydrophobicity(seq),
  aacharge = charge(seq))]

# retain only sequences that feature in final AAV table
testSeq <- testSeq[AAV %in% nsAavs]

# sort hydrophobicity vals according to AAV table
testSeq <- testSeq[order(match(AAV, nsAavs))]

# calculated mean hydrophobicity weighted by abundance of each individual amino acid variant.
aoaDat[, ":="(meanHydrophob = unlist(lapply(1:nrow(aoaDat), function(s)
  weighted.mean(x = testSeq$hydro,
    w = aoaDat[s, .SD, .SDcols = nsAavs]))),
  meanCharge = unlist(lapply(1:nrow(aoaDat), function(s)
    weighted.mean(x = testSeq$aacharge,
      w = aoaDat[s, .SD, .SDcols = nsAavs]))))]

aoaHydroLm <- lm(meanHydrophob ~ bfi, data = aoaDat)
aoaChargeLm <- lm(meanCharge ~ bfi, data = aoaDat)

aoaHydro_bfi <- ggplot(aoaDat,
    aes(x = bfi, y = meanHydrophob, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(aoa)~hydrophobicity),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

aoaCharge_bfi <- ggplot(aoaDat,
    aes(x = bfi, y = meanCharge, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(aoa)~net~charge),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../figures/aoa_hydrophobicity.pdf", aoaHydro_bfi, height = 6,
  width = 7, device = "pdf")

# aoa OTU analysis
otus <- fread("AOA_amoAOtuTab.txt")
colnames(otus)[1] <- "OTU"

aoaOtus <- transDT(otus, transCol = "OTU", rowID = "seqCode")

otuOccs <- specnumber(aoaOtus[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
aoaFiltOtus <- aoaOtus[, .SD,
  .SDcols = c("seqCode", names(otuOccs[otuOccs > 1]))]

# get library sizes of filtered AAV table
aoaFiltOtus[, libSize := rowSums(aoaFiltOtus[, !"seqCode"])]

# merge with env metadata
aoaOtuDat <- merge(dat, aoaFiltOtus, by = "seqCode")

otuCols <- grep("OTU", colnames(aoaOtuDat), value = T)

# rarefy AAV table to even depth
aoaOtuDat[, (otuCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = otuCols]

# remove any empty OTUs
emptyOTUs <- aoaOtuDat[, colSums(.SD) == 0, .SDcols = otuCols]
aoaOtuDat[, (names(emptyOTUs[emptyOTUs == T])) := NULL]

# remove empty OTUs from column list
otuCols <- otuCols[! emptyOTUs]

# create beta diversity matrix
aoaBetaOtuMat <- beta.pair(
  ifelse(aoaOtuDat[, .SD, .SDcols = otuCols] != 0, 1, 0))

nmdsRes <- metaMDS(aoaBetaOtuMat[[3]], autotransform = F, trymax = 200)

scoreCols <- c("x", "y")
aoaOtuDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

aoaBfiOtuEuclid <- vegdist(aoaOtuDat$bfi, "euclid")

aoaOtuBfiDat <- data.table(dissim = as.vector(aoaBetaOtuMat[[3]]),
  bfiSim = as.vector(aoaBfiOtuEuclid), type = "OTU")

# mcra negative exponential decay model
aoaAavDecay <- decay.model(aoaBetaMat[[3]], aoaBfiEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)
aoaOtuDecay <- decay.model(aoaBetaOtuMat[[3]], aoaBfiOtuEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)

# combine the two dissim data frames, then calculate fits for both decay models
aoaDists <- rbindlist(list(aoaAavBfiDat, aoaOtuBfiDat))
aoaAavFit <- calcDecayFit(aoaAavDecay)
aoaOtuFit <- calcDecayFit(aoaOtuDecay)

# bootstrap coefficients and combine
aoaAavBoot <- boot.coefs.decay(aoaAavDecay, 1000)
aoaOtuBoot <- boot.coefs.decay(aoaOtuDecay, 1000)

aoaBoot <- combineBoots(gene = "aoa_amoA", aoaAavBoot, aoaOtuBoot)

aoaOtuDat[, richness := specnumber(.SD), .SDcols = otuCols]

aoaRichness <- rbindlist(list(AAV = aoaDat[, .(bfi, geol, richness, month)],
  OTU = aoaOtuDat[, .(bfi, geol, richness, month)]), idcol = "type")

aoaAavRichness <- glm.nb(richness ~ bfi, data = aoaRichness[type == "AAV"])
aoaOtuRichness <- glm.nb(richness ~ bfi, data = aoaRichness[type == "OTU"])
aoaAavMonthRichness <- glm.nb(richness ~ month + bfi,
  data = aoaRichness[type == "AAV"])
aoaOtuMonthRichness <- glm.nb(richness ~ month + bfi,
  data = aoaRichness[type == "OTU"])
lapply(list(
  aoaAavRichness, aoaAavMonthRichness, aoaOtuRichness, aoaOtuMonthRichness),
  summary)
lapply(list(
  aoaAavRichness, aoaAavMonthRichness, aoaOtuRichness, aoaOtuMonthRichness),
  Dsquared, adjust = T)
lapply(list(
  aoaAavRichness, aoaAavMonthRichness, aoaOtuRichness, aoaOtuMonthRichness),
  function(x) exp(coef(x)))
################################ anammox hzo analysis ##########################

hzo <- fread("hzo_AAV.txt")
colnames(hzo)[1] <- "AAV"

hzo <- transDT(hzo, transCol = "AAV", rowID = "seqCode")

aavOccs <- specnumber(hzo[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
hzoFilt <- hzo[, .SD, .SDcols = c("seqCode", names(aavOccs[aavOccs > 1]))]

# get library sizes of filtered AAV table
hzoFilt[, libSize := rowSums(hzoFilt[, !"seqCode"])]

# merge with env metadata
hzoDat <- merge(dat, hzoFilt, by = "seqCode")

aavCols <- grep("AAV", colnames(hzoDat), value = T)

# rarefy AAV table to even depth
hzoDat[, (aavCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = aavCols]

# remove any empty AAVs
emptyAAVs <- hzoDat[, colSums(.SD) == 0, .SDcols = aavCols]
hzoDat[, (names(emptyAAVs[emptyAAVs == T])) := NULL]

# update remaining list of aavs
nsAavs <- aavCols[!emptyAAVs]

# create beta diversity matrix
hzoBetaMat <- beta.pair(ifelse(hzoDat[, .SD, .SDcols = nsAavs] != 0, 1, 0))

nmdsRes <- metaMDS(hzoBetaMat[[3]], autotransform = F, trymax = 200)

scoreCols = c("x", "y")
hzoDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

hzoBfiEuclid <- vegdist(hzoDat$bfi, "euclid")

hzoAavBfiDat <- data.table(dissim = as.vector(hzoBetaMat[[3]]),
  bfiSim = as.vector(hzoBfiEuclid))

# think about using decay.model from betapart package
# test association between bfi and hzo AAV richness.
hzoDat[, richness := specnumber(.SD), .SDcols = nsAavs]

# test for differences in hydrophobicity
test <- read.fasta(
  "sequences/hzo/mergedSeqs/allErrorCorrected/framebotOut/AAVs.fasta",
  seqtype = "AA", as.string = T)

testSeq <- data.table(AAV = names(test),
  seq = unlist(sapply(test, getSequence, as.string = T)))

testSeq[, ":="(hydro = hydrophobicity(seq),
  aacharge = charge(seq))]

# retain only sequences that feature in final AAV table
testSeq <- testSeq[AAV %in% nsAavs]

# sort hydrophobicity vals according to AAV table
testSeq <- testSeq[order(match(AAV, nsAavs))]

# calculated mean hydrophobicity weighted by abundance of each individual
# amino acid variant.
hzoDat[, ":="(meanHydrophob = unlist(lapply(1:nrow(hzoDat), function(s)
  weighted.mean(x = testSeq$hydro,
    w = hzoDat[s, .SD, .SDcols = nsAavs]))),
  meanCharge = unlist(lapply(1:nrow(hzoDat), function(s)
  weighted.mean(x = testSeq$aacharge,
    w = hzoDat[s, .SD, .SDcols = nsAavs]))))]

hzoHydroLm <- lm(meanHydrophob ~ bfi, data = hzoDat)
hzoChargeLm <- lm(meanCharge ~ bfi, data = hzoDat)

hzoHydro_bfi <- ggplot(hzoDat,
    aes(x = bfi, y = meanHydrophob, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(hzo)~hydrophobicity),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

hzoCharge_bfi <- ggplot(hzoDat,
    aes(x = bfi, y = meanCharge, col = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(x = "Base flow index",
    y = expression(Weighted~mean~italic(hzo)~net~charge),
    col = "Geology") +
  scale_color_viridis(discrete = T) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

# hzo OTU analysis
otus <- fread("hzoOtuTab.txt")
colnames(otus)[1] <- "OTU"

hzoOtus <- transDT(otus, transCol = "OTU", rowID = "seqCode")

otuOccs <- specnumber(hzoOtus[, !"seqCode"], MARGIN = 2)

# remove AAVs that occur in only one sample
hzoFiltOtus <- hzoOtus[, .SD,
  .SDcols = c("seqCode", names(otuOccs[otuOccs > 1]))]

# get library sizes of filtered AAV table
hzoFiltOtus[, libSize := rowSums(hzoFiltOtus[, !"seqCode"])]

# merge with env metadata
hzoOtuDat <- merge(dat, hzoFiltOtus, by = "seqCode")

otuCols <- grep("OTU", colnames(hzoOtuDat), value = T)

# rarefy AAV table to even depth
hzoOtuDat[, (otuCols) := as.data.table(rrarefy(.SD, sample = min(libSize))),
  .SDcols = otuCols]

# remove any empty OTUs
emptyOTUs <- hzoOtuDat[, colSums(.SD) == 0, .SDcols = otuCols]
hzoOtuDat[, (names(emptyOTUs[emptyOTUs == T])) := NULL]

# remove empty OTUs from column list
otuCols <- otuCols[! emptyOTUs]

# create beta diversity matrix
hzoBetaOtuMat <- beta.pair(
  ifelse(hzoOtuDat[, .SD, .SDcols = otuCols] != 0, 1, 0))

nmdsRes <- metaMDS(hzoBetaOtuMat[[3]], autotransform = F, trymax = 200)

scoreCols <- c("x", "y")
hzoOtuDat[, (scoreCols) := as.data.table(scores(nmdsRes))]

hzoBfiOtuEuclid <- vegdist(hzoOtuDat$bfi, "euclid")

hzoOtuBfiDat <- data.table(dissim = as.vector(hzoBetaOtuMat[[3]]),
  bfiSim = as.vector(hzoBfiOtuEuclid))

# hzo negative exponential decay model
hzoAavDecay <- decay.model(betaMat[[3]], hzoBfiEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)
hzoOtuDecay <- decay.model(betaOtuMat[[3]], hzoBfiOtuEuclid,
  model.type = "exponential", y.type = "dissimilarities", perm = 1000)

# give dissim data tables a type column
hzoOtuBfiDat[, type := "OTU"]
hzoAavBfiDat[, type := "AAV"]

# combine the two dissim data frames, then calculate fits for both decay models
hzoDists <- rbindlist(list(hzoAavBfiDat, hzoOtuBfiDat))
hzoAavFit <- calcDecayFit(hzoAavDecay)
hzoOtuFit <- calcDecayFit(hzoOtuDecay)

# bootstrap coefficients and combine
hzoAavBoot <- boot.coefs.decay(hzoAavDecay, 1000)
hzoOtuBoot <- boot.coefs.decay(hzoOtuDecay, 1000)
hzoBoot <- combineBoots(gene = "hzo", hzoAavBoot, hzoOtuBoot)

# combine all dist data
distData <- rbindlist(list(
  aoa = aoaDists,
  aob = aobDists,
  hzo = hzoDists,
  nirs = nirsDists,
  mcra = mcraDists,
  pmoa = pmoaDists),
  idcol = "gene")

distData[, gene := factor(gene,
  levels = unique(gene),
  labels = c(expression(AOA~italic(amoA)),
    expression(AOB~italic(amoA)),
    expression(italic(hzo)),
    expression(italic(nirS)),
    expression(italic(mcrA)),
    expression(italic(pmoA))))]

# combine all prediction data
aavDistPredicts <- rbindlist(list(
  aoa = aoaAavFit,
  aob = aobAavFit,
  hzo = hzoAavFit,
  nirs = nirsAavFit,
  mcra = mcraAavFit,
  pmoa = pmoaAavFit),
  idcol = "gene")

aavDistPredicts[, gene := factor(gene,
  levels = unique(gene),
  labels = c(expression(AOA~italic(amoA)),
    expression(AOB~italic(amoA)),
    expression(italic(hzo)),
    expression(italic(nirS)),
    expression(italic(mcrA)),
    expression(italic(pmoA))))]

otuDistPredicts <- rbindlist(list(
  aoa = aoaOtuFit,
  aob = aobOtuFit,
  hzo = hzoOtuFit,
  nirs = nirsOtuFit,
  mcra = mcraOtuFit,
  pmoa = pmoaOtuFit),
  idcol = "gene")

otuDistPredicts[, gene := factor(gene,
  levels = unique(gene),
  labels = c(expression(AOA~italic(amoA)),
    expression(AOB~italic(amoA)),
    expression(italic(hzo)),
    expression(italic(nirS)),
    expression(italic(mcrA)),
    expression(italic(pmoA))))]

# data table to add annotation text
statLabs <- data.table(gene = rep(unique(distData$gene), each = 2),
  type = rep(c("AAV", "OTU"), times = 6),
  lab = unlist(lapply(list(aoaAavDecay,
    aoaOtuDecay,
    aobAavDecay,
    aobOtuDecay,
    hzoAavDecay,
    hzoOtuDecay,
    nirsAavDecay,
    nirsOtuDecay,
    mcraAavDecay,
    mcraOtuDecay,
    pmoaAavDecay,
    pmoaOtuDecay), decayLabel)))

statLabs[, gene := factor(gene,
  levels = unique(gene),
  labels = c(expression(AOA~italic(amoA)),
    expression(AOB~italic(amoA)),
    expression(italic(hzo)),
    expression(italic(nirS)),
    expression(italic(mcrA)),
    expression(italic(pmoA))))]

# manually edit one stat label
statLabs$lab[1] <- "list(italic(R^2) =='0.40', italic('P')< 0.01)"

fullPanel <- ggplot(distData, aes(x = bfiSim, y = dissim, col = type)) +
  geom_point(size = 1.5, alpha = 0.3) +
  scale_colour_manual(values = c(darkCol, lightCol), name = "") +
  labs(x = expression(paste(Delta, "BFI")),
    y = expression(paste("S", "\u00F8", "rensen dissimilarity"))) +
  facet_wrap(~ gene, ncol = 2, labeller = label_parsed) +
  geom_line(data = aavDistPredicts, aes(x = x, y = y), col = darkCol,
    size = 1.2) +
  geom_line(data = otuDistPredicts, aes(x = x, y = y), col = lightCol,
    size = 1.2) +
  geom_text(data = statLabs[type == "AAV"], aes(x = Inf, y = -Inf, label = lab),
    parse = T, hjust = 1.1, vjust = -0.7, col = darkCol, fontface = "bold") +
  geom_text(data = statLabs[type == "OTU"], aes(x = Inf, y = -Inf, label = lab),
    parse = T, hjust = 1.1, vjust = 0, col = lightCol, fontface = "bold") +
  theme_bw() +
  coord_cartesian(ylim = c(0.0, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 0.7))) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../figures/full_dist_panel.pdf", fullPanel, height = 10, width = 8,
  device = cairo_pdf)

# hzo richness analysis
hzoOtuDat[, richness := specnumber(.SD), .SDcols = otuCols]

hzoRichness <- rbindlist(list(AAV = hzoDat[, .(bfi, geol, richness, month)],
  OTU = hzoOtuDat[, .(bfi, geol, richness, month)]), idcol = "type")

hzoAavRichness <- glm.nb(richness ~ bfi, data = hzoRichness[type == "AAV"])
hzoMonthAavRichness <- glm.nb(richness ~ month + bfi,
  data = hzoRichness[type == "AAV"])
hzoOtuRichness <- glm.nb(richness ~ bfi, data = hzoRichness[type == "OTU"])
hzoMonthOtuRichness <- glm.nb(richness ~ month + bfi,
  data = hzoRichness[type == "OTU"])
lapply(list(
  hzoAavRichness, hzoMonthAavRichness, hzoOtuRichness, hzoMonthOtuRichness),
  summary)
lapply(list(
  hzoAavRichness, hzoMonthAavRichness, hzoOtuRichness, hzoMonthOtuRichness),
  Dsquared, adjust = T)
lapply(list(
  hzoAavRichness, hzoMonthAavRichness, hzoOtuRichness, hzoMonthOtuRichness),
  function(x) exp(coef(x)))

# assemble facet plots for AAV and OTU richness for each gene
richList <- list(aoa = aoaRichness,
  aob = aobRichness,
  hzo = hzoRichness,
  nirS = nirsRichness,
  mcrA = mcraRichness,
  pmoA = pmoaRichness)

#list of all aav and otu richness models
aavMods <- list(
    aoaAavRichness,
    aobAavRichness,
    hzoAavRichness,
    nirsAavRichness,
    mcraAavRichness,
    pmoaAavRichness
)

otuMods <- list(
    aoaOtuRichness,
    aobOtuRichness,
    hzoOtuRichness,
    nirsOtuRichness,
    mcraOtuRichness,
    pmoaOtuRichness
)

richPreds <- data.table(
  bfi = seq(min(dat$bfi), max(dat$bfi), 0.01))

aavPredicts <- lapply(aavMods, function(x)
  predict(x, newdata = richPreds, type = "link", se.fit = T))
otuPredicts <- lapply(otuMods, function(x)
  predict(x, newdata = richPreds, type = "link", se.fit = T))

test <- lapply(1:6, function(x)
  data.table(
    bfi = rep(richPreds$bfi, times = 2),
    type = rep(c("AAV", "OTU"), each = length(richPreds$bfi)),
    prediction = c(aavPredicts[[x]]$fit, otuPredicts[[x]]$fit),
    se = c(aavPredicts[[x]]$se.fit, otuPredicts[[x]]$se.fit),
    signif = c(
      rep(ifelse(
        coef(summary(aavMods[[x]]))[, 4][2] < 0.05, 1, 2),
          times = length(richPreds$bfi)),
        rep(ifelse(coef(summary(otuMods[[x]]))[, 4][2] < 0.05, 1, 2),
          times = length(richPreds$bfi)))))

# bind all predictions into one df
allPreds <- rbindlist(test, idcol = T)

# calculate prediction 95% conf intervals and
# back transform richness predictions
allPreds[, ":="(
  uppCI = exp(prediction + (1.96 * se)),
  lowCI = exp(prediction - (1.96 * se)),
  richness = exp(prediction),
  geol = NA
)]

allPreds[, signif := as.factor(signif)]

# add geom_line and geom_ribbon
richnessPlots <- list(length = 6)
for(i in 1:6){
  richnessPlots[[i]] <- ggplot(data = richList[[i]],
      aes(x = bfi, y = richness, fill = geol)) +
    geom_point(size = 3, shape = 21, alpha = 0.7) +
    geom_ribbon(data = allPreds[.id == i, ],
      aes(ymin = lowCI, ymax = uppCI),
      fill = "grey", alpha = 0.4) +
    geom_line(data = allPreds[.id == i, ],
      aes(x = bfi, y = richness, linetype = signif)) +
    scale_fill_manual(values = c("white", "grey", "darkseagreen3")) +
    scale_linetype_manual(name = "", values = c(1, 2), drop = F,
      labels = c(expression(italic(P)<0.05), expression(italic(P)>=0.05))) +
    facet_wrap(~ type, scales = "free_y") +
    labs(x = "Base flow index", y = "Richness", fill = "") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 14),
      panel.grid = element_blank(),
      strip.text.x = element_text(size = 14))
    }

# arrange signif ones into panel and write up stats. include ns ones in SI
richnessPlots[[1]] <- richnessPlots[[1]] +
  labs(title = "A") +
  theme(axis.title.x = element_blank(),
    plot.title = element_text(size = 18, hjust = -0.15, vjust = 1.2,
      margin = margin(b = -20)))
richnessPlots[[2]] <- richnessPlots[[2]] +
  labs(title = "B") +
  theme(axis.title.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 18, hjust = -0.15, vjust = 1.2,
      margin = margin(b = -20)))
richnessPlots[[3]] <- richnessPlots[[3]] +
  labs(title = "C") +
  theme(legend.position = "none",
    plot.title = element_text(size = 18, hjust = -0.15, vjust = 1.2,
      margin = margin(b = -20)))

mainRichPanel <- richnessPlots[[1]] + richnessPlots[[2]] + richnessPlots[[3]] +
  plot_layout(ncol = 1)

mainRichPanel <- mainRichPanel + plot_layout(guides = "collect")

ggsave("../figures/main_richness_panel.pdf", mainRichPanel, height = 8,
  width = 7, device = "pdf")

richnessPlots[[4]] <- richnessPlots[[4]] +
  labs(title = "A") +
  theme(axis.title.x = element_blank(),
    plot.title = element_text(size = 18, hjust = -0.15, vjust = 1.2,
      margin = margin(b = -20)))
richnessPlots[[5]] <- richnessPlots[[5]] +
  labs(title = "B") +
  theme(axis.title.x = element_blank(),
    plot.title = element_text(size = 18, hjust = -0.15, vjust = 1.2,
      margin = margin(b = -20)))
richnessPlots[[6]] <- richnessPlots[[6]] +
  labs(title = "C") +
  theme(plot.title = element_text(size = 18, hjust = -0.15, vjust = 1.2,
      margin = margin(b = -20)))

siRichPanel <- richnessPlots[[4]] + richnessPlots[[5]] + richnessPlots[[6]] +
  plot_layout(ncol = 1)

siRichPanel <- siRichPanel + plot_layout(guides = "collect")

ggsave("../figures/si_richness_panel.pdf", siRichPanel, height = 8, width = 7,
  device = "pdf")

# create ridge plot of all bootstrapped coefficients
bootCoefs <- rbindlist(
  list(aoaBoot, aobBoot, hzoBoot, nirsBoot, mcraBoot, pmoaBoot))

bootCoefs[, gene := factor(gene,
  levels = c("pmoA", "mcrA", "nirS", "hzo", "AOB_amoA", "aoa_amoA"),
  labels = c(expression(italic(pmoA)),
    expression(italic(mcrA)),
    expression(italic(nirS)),
    expression(italic(hzo)),
    expression(AOB~italic(amoA)),
    expression(AOA~italic(amoA))))]

coefPlot <- ggplot(bootCoefs,
    aes(x = coefs, y = gene, fill = type, col = type)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9, rel_min_height = 0.005)+
  theme_bw() +
  labs(x = "Boostrapped coefficient estimate", y = "", fill = "", col = "") +
  scale_fill_manual(values = c(darkCol, lightCol)) +
  scale_color_manual(values = c(darkCol, lightCol)) +
  scale_y_discrete(expand = c(0.01, 0),
    labels = parse(text = levels(bootCoefs$gene))) +
  expand_limits(y = 7) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())

ggsave("../figures/bootstrapped_coefficients.pdf", coefPlot, height = 4,
  width = 6, device = "pdf")

# merge data for aob and aoa hydro plot
hydroData <- rbindlist(list(
  AOA = aoaDat[, .(bfi, geol, month, meanHydrophob, meanCharge)],
  AOB = aobDat[, .(bfi, geol, month, meanHydrophob, meanCharge)]),
  idcol = T)

hydroData <- melt(hydroData, id.vars = c(".id", "bfi", "geol", "month"))

hydroPredData <- data.table(bfi = seq(min(dat$bfi), max(dat$bfi), 0.005))
hydroPreds <- lapply(list(
  aoaHydroLm, aoaChargeLm, aobHydroLm, aobChargeLm), function(x)
  predict(x, newdata = hydroPredData, se.fit = T))

allHydroPreds <- data.table(bfi = rep(hydroPredData$bfi, times = 4),
    .id = rep(c("AOA", "AOB"), each = 2 * nrow(hydroPredData)),
    value = c(
      hydroPreds[[1]]$fit, hydroPreds[[2]]$fit, hydroPreds[[3]]$fit,
        hydroPreds[[4]]$fit),
    variable = rep(c(rep("meanHydrophob", times = nrow(hydroPredData)),
      rep("meanCharge", times = nrow(hydroPredData))), times = 2),
    se = c(hydroPreds[[1]]$se.fit, hydroPreds[[2]]$se.fit, hydroPreds[[3]]$se.fit, hydroPreds[[4]]$se.fit))

allHydroPreds[, ":="(uppCI = value + (1.96 * se),
  lowCI = value - (1.96 * se),
  geol = NA,
  month = NA)]

# need to fix geol legend
hydroPlot <- ggplot(hydroData[variable == "meanHydrophob"],
    aes(x = bfi, y = value, fill = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.7) +
  facet_wrap(~ .id, scales = "free_y") +
  scale_shape_manual("Sample month", values = c(21, 22)) +
  scale_fill_manual("Geology", values = c("white", "grey", "darkseagreen3")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  geom_ribbon(data = allHydroPreds[variable == "meanHydrophob"],
    aes(x = bfi, y = value, ymin = lowCI, ymax = uppCI),
    alpha = 0.4, fill = "grey", colour = NA) +
  geom_line(data = allHydroPreds[variable == "meanHydrophob"],
    aes(x = bfi, y = value),
    linetype = 1) +
  labs(x = "Base flow index", y = "Average hydrophobicity", title = "A") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size = 20),
    panel.grid = element_blank())

chargePlot <- ggplot(hydroData[variable == "meanCharge"],
    aes(x = bfi, y = value, fill = geol, shape = month)) +
  geom_point(size = 4, alpha = 0.7) +
  facet_wrap(~ .id, scales = "free_y") +
  scale_shape_manual("Sample month", values = c(21, 22)) +
  scale_fill_manual("Geology", values = c("white", "grey", "darkseagreen3")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  geom_ribbon(data = allHydroPreds[variable == "meanCharge"],
    aes(x = bfi, y = value, ymin = lowCI, ymax = uppCI),
    alpha = 0.4, fill = "grey", colour = NA) +
  geom_line(data = allHydroPreds[variable == "meanCharge"],
    aes(x = bfi, y = value),
    linetype = 1) +
  labs(x = "Base flow index", y = "Average net charge", title = "B") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size = 20),
    panel.grid = element_blank())

plotLegend <- cowplot::get_legend(hydroPlot)

hydroPlot <- hydroPlot + theme(legend.position = "none")
chargePlot <- chargePlot + theme(legend.position = "none")

aoAminoPanel <-((hydroPlot / chargePlot) | plotLegend) +
  plot_layout(widths = c(1, 0.2))

ggsave("../figures/aoa_aob_hyrdo.pdf", aoAminoPanel, height = 8, width = 9,
  device = "pdf")

# aoa:aob ratio
dat[, totalAo := round(aoa) + round(aob)]
dat[, aoaProp := round(aoa)/totalAo]
aoaRatio <- glm(aoaProp ~ bfi, weights = totalAo, data = dat,
  family = binomial)
summary(aoaRatio)

dat[, hzoTotalAo := round(totalAo) + round(hzo)]
dat[, hzoProp := round(hzo)/hzoTotalAo]
hzoRatio <- glm(hzoProp ~ bfi, weights = hzoTotalAo, data = dat,
  family = binomial)
summary(hzoRatio)
Dsquared(hzoRatio)

aoaPreds <- predict(aoaRatio, newdata = richPreds, type = "link", se.fit = T)
hzoPreds <- predict(hzoRatio, newdata = richPreds, type = "link", se.fit = T)

aoaPreds <- data.table(bfi = richPreds$bfi, linkPred = aoaPreds$fit,
  linkSe = aoaPreds$se.fit)
hzoPreds <- data.table(bfi = richPreds$bfi, linkPred = hzoPreds$fit,
  linkSe = hzoPreds$se.fit)

invLogit <- function(x){
    exp(x)/(1+exp(x))
}

aoaPreds[, ":="(aoaProp = invLogit(linkPred),
  uppCI = invLogit(linkPred + (1.96 * linkSe)),
  lowCI = invLogit(linkPred - (1.96 * linkSe)),
  geol = NA)]
hzoPreds[, ":="(hzoProp = invLogit(linkPred),
  uppCI = invLogit(linkPred + (1.96 * linkSe)),
  lowCI = invLogit(linkPred - (1.96 * linkSe)),
  geol = NA)]

# prediction intervals are too close to see
aoaProb <- ggplot(dat, aes(x = bfi, y = aoaProp, fill = geol)) +
  geom_point(size = 3, alpha = 0.7, shape = 21) +
  geom_line(data = aoaPreds, aes(x = bfi, y = aoaProp), size = 1.1) +  scale_fill_manual("Geology", values = c("white", "grey", "darkseagreen3")) +
  labs(x = "Base flow index", y = "P(AOA)", title = "A") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

hzoProb <- ggplot(dat, aes(x = bfi, y = hzoProp, fill = geol)) +
  geom_point(size = 3, alpha = 0.7, shape = 21) +
  geom_line(data = hzoPreds, aes(x = bfi, y = hzoProp), size = 1.1) +  scale_fill_manual("Geology", values = c("white", "grey", "darkseagreen3")) +
  labs(x = "Base flow index", y = "P(anammox)", title = "B") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

aoPanel <- aoaProb + hzoProb + plot_layout(ncol = 2, guides = "collect")

ggsave("../figures/ammox_panel.pdf", aoPanel, height = 4.5, width = 8,
  device = "pdf")

# nmds panel
aavNmds <- lapply(list(AOA = aoaDat,
    AOB = aobDat,
    hzo = hzoDat,
    nirS = nirsDat,
    mcrA = mcraDat,
    pmoA = pmoaDat),
  function(x) x[, .(x, y, month, bfi)])

otuNmds <- lapply(list(AOA = aoaOtuDat,
    AOB = aobOtuDat,
    hzo = hzoOtuDat,
    nirS = nirsOtuDat,
    mcrA = mcraOtuDat,
    pmoA = pmoaOtuDat),
  function(x) x[, .(x, y, month, bfi)])

aavNmds <- rbindlist(aavNmds, idcol = "gene")
otuNmds <- rbindlist(otuNmds, idcol = "gene")

allNmds <- rbindlist(list(AAV = aavNmds, OTU = otuNmds), idcol = "type")

allNmds[, gene := factor(gene,
  levels = unique(gene),
  labels = c(expression(AOA~italic(amoA)),
    expression(AOB~italic(amoA)),
    expression(italic(hzo)),
    expression(italic(nirS)),
    expression(italic(mcrA)),
    expression(italic(pmoA))))]

nmdsPanel <- ggplot(allNmds, aes(x = x, y = y, col = bfi, shape = month)) +
  geom_point(size = 3, alpha = 0.7) +
  facet_grid(gene ~ type, labeller = label_parsed, scales = "free_y") +
  labs(x = "NMDS 1", y = "NMDS 2", col = "BFI", shape = "Sample month") +
  theme_bw() +
  scale_color_viridis() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../figures/nmds_panel.pdf", nmdsPanel, height = 12, width = 9,
  device = "pdf")
