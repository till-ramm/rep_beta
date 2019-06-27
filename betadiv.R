###reptile beta diversity

library(betapart)
library(picante)
library(ade4)
library(phytools)
library(geosphere)
library(vegan)
library(ggplot2)

setwd("D:/Dropbox/Africa2")

tree <- read.tree("tree.txt")
bodysize <- read.table("traits_bodysize.txt", header = T, row.names = 1)
table <- read.table("species_info.txt", header = T, row.names = 1)
loc.info <- read.table("loc_info_new.txt", header = TRUE)


hist((bodysize$bodysize_g_log10))
hist((decompo_beta$betadiv$SES_PhyloSor_turn))

p <- ggplot(bodysize, aes(x=bodysize$bodysize_g_log10)) +theme_minimal()+geom_histogram(aes(y=..density..), colour="black", fill="pink")+geom_density(color="darkblue", fill="lightblue", alpha=0.4, size=2)
p


plot(size)


#bodysize bins
size <- cut(bodysize$bodysize_g_log10, 3, include.lowest=TRUE, labels=c("small","medium", "large"))

table(size)

table <- cbind(table, size)

##different tables, different biomes


biomes <- c("all_b", "DXS", "TGSS", "TMBF", "MGS")
sets <- c("no_s", "lizards", "snakes", "small", "medium", "large")



#wb <- loadWorkbook("D:/Dropbox/Africa2/results.xlsx", create = TRUE)

## species filtering

for(b in biomes){
  if(b == "all_b") locality.list <- loc.info
  if(b != "all_b") locality.list <- loc.info[loc.info$Biome=="TMBF",]
  for(s in sets){
    if(s == "no_s") species.list <- table[,c(5:96)]
    if(s =="snakes") species.list <- table[table$clade=="s",c(5:96)]
    if(s =="lizards") species.list <- table[table$clade=="l",c(5:96)]
    if(s == "small") species.list <- table[table$size=="small",c(5:96)]
    if(s == "medium") species.list <- table[table$size=="medium",c(5:96)]
    if(s == "large") species.list <- table[table$size=="large",c(5:96)]
    
    ###locality table
    
    locality.list2 <- rownames(locality.list)
    
    locality.list2 <- t(locality.list2)
    
    submatrix <- t(species.list[colnames(species.list) %in% locality.list2])
    
    submatrix <- submatrix[rowSums(submatrix)!=0,]
    submatrix <- submatrix[,colSums(submatrix)!=0]
    
    rownames(submatrix)
    
    #spatial distances
    
    loc.info.variables <- loc.info[rownames(loc.info) %in% rownames(submatrix),]
    
    lonlat <- as.matrix(cbind(loc.info.variables$Lon, loc.info.variables$Lat))
    spatial <- dist(lonlat, method="euclidean")
    
    
    #environmental distances, current + LGM
    
    
    
    pcacurrent <-  prcomp(loc.info.variables[5:23], scale = T, center = T)
    biplot(pcacurrent)
    summary(pcacurrent)
    
    pcapast <-  prcomp(loc.info.variables[24:42], scale = T, center = T)
    summary(pcapast)
    
    currentclimate <- (dist(cbind(pcacurrent$x[,1:4]), method = "manhattan"))
    summary(currentclimate)
    
    pastclimate <- (dist(cbind(pcapast$x[,1:4]), method = "manhattan"))
    summary(pastclimate)
    
    ###tree
    
    tips_to_drop <- tree$tip.label[!(tree$tip.label %in% colnames(submatrix))]
    sub_tree <- drop.tip(tree, tip = tips_to_drop)


    ###calculate SES.PBD
    
    decompo_beta <- beta.pd.decompo(submatrix, sub_tree,type="both",output.dist=T, random=1000)
    
    
    save(decompo_beta, file=paste("decompo_beta", b, s, ".Rdata", sep=" ")) 
    
    
    #createSheet(wb, name = paste("decompo_beta", b, s, sep=" "))
    
    #writeWorksheet(wb, decompo_beta$betadiv$PhyloSor, sheet = paste("decompo_beta", b, s, sep=" "), startRow = 1, startCol = 1)
    
    #saveWorkbook(wb)
        
    #betapart measures
    
    ##tax beta
    tax.core <- betapart.core(submatrix)
    tax.dist.core <- beta.pair(tax.core, index.family = "sor")
    tax.dist.multi <- beta.multi(tax.core, index.family = "sor")
    
    save(tax.dist.core, file=paste("tax.dist.core", b, s, ".Rdata", sep=" ")) 
    save(tax.dist.multi, file=paste("tax.dist.multi", b, s, ".Rdata", sep=" ")) 
    
    
    #phylo beta
    phylo.core <- phylo.betapart.core(submatrix, sub_tree)
    phylo.dist.core <- phylo.beta.pair(phylo.core, sub_tree, index.family = "sor")
    phylo.dist.multi <- phylo.beta.multi(phylo.core, tree, index.family = "sor")
    
    save(phylo.dist.core, file=paste("phylo.dist.core", b, s, ".Rdata", sep=" ")) 
    save(phylo.dist.multi, file=paste("phylo.dist.multi", b, s, ".Rdata", sep=" ")) 
    
    
    #createSheet(wb, name = paste("tax.beta.sim", b, s, sep=" "))
    #createSheet(wb, name = paste("tax.beta.sne", b, s, sep=" "))
    #createSheet(wb, name = paste("tax.beta.sor", b, s, sep=" "))
    
    #createSheet(wb, name = paste("phylo.beta.sim", b, s, sep=" "))
    #createSheet(wb, name = paste("phylo.beta.sne", b, s, sep=" "))
    #createSheet(wb, name = paste("phylo.beta.sor", b, s, sep=" "))
    
    #createSheet(wb, name = paste("tax.beta.multi", b, s, sep=" "))
    #writeWorksheet(wb, tax.dist.multi, sheet = paste("tax.beta.multi", b, s, sep=" "), startRow = 1, startCol = 1)
    
    #createSheet(wb, name = paste("phylo.beta.multi", b, s, sep=" "))
    #writeWorksheet(wb, tax.dist.multi, sheet = paste("phylo.beta.multi", b, s, sep=" "), startRow = 1, startCol = 1)
    
    #writeWorksheet(wb, tax.dist.core$beta.sim, sheet = paste("tax.beta.sim", b, s, sep=" "), startRow = 1, startCol = 1)
    #writeWorksheet(wb, tax.dist.core$beta.sne, sheet = paste("tax.beta.sne", b, s, sep=" "), startRow = 1, startCol = 1)
    #writeWorksheet(wb, tax.dist.core$beta.sor, sheet = paste("tax.beta.sor", b, s, sep=" "), startRow = 1, startCol = 1)
    
    #writeWorksheet(wb, phylo.dist.core$beta.sim, sheet = paste("phylo.beta.sim", b, s, sep=" "), startRow = 1, startCol = 1)
    #writeWorksheet(wb, phylo.dist.core$beta.sne, sheet = paste("phylo.beta.sne", b, s, sep=" "), startRow = 1, startCol = 1)
    #writeWorksheet(wb, phylo.dist.core$beta.sor, sheet = paste("phylo.beta.sor", b, s, sep=" "), startRow = 1, startCol = 1)
    
    
    #saveWorkbook(wb)

  }
}

##neuer loop? taxonomic / phylo, decay_spatial, decay_current, decay_lgm

#distance decay model, groups vs spatial, current climate, historical climate

biomes <- c("all_b")
sets <- c("snakes", "small", "medium", "large")
predictor <-c("spatial", "current", "past")

wb <- loadWorkbook("D:/Dropbox/Africa2/results_decay2.xlsx", create = TRUE)


## species filtering

for(b in biomes){
  if(b == "all_b") locality.list <- loc.info
  if(b != "all_b") locality.list <- loc.info[loc.info$Biome==b,]
  for(s in sets){
    if(s == "no_s") species.list <- table[,c(5:96)]
    if(s =="snakes") species.list <- table[table$clade=="s",c(5:96)]
    if(s =="lizards") species.list <- table[table$clade=="l",c(5:96)]
    if(s == "small") species.list <- table[table$size=="small",c(5:96)]
    if(s == "medium") species.list <- table[table$size=="medium",c(5:96)]
    if(s == "large") species.list <- table[table$size=="large",c(5:96)]
    
    ###locality table
    
    locality.list2 <- rownames(locality.list)
    
    locality.list2 <- t(locality.list2)
    
    submatrix <- t(species.list[colnames(species.list) %in% locality.list2])
    
    submatrix <- submatrix[rowSums(submatrix)!=0,]
    submatrix <- submatrix[,colSums(submatrix)!=0]
    
    rownames(submatrix)
    
    #spatial distances
    
    loc.info.variables <- loc.info[rownames(loc.info) %in% rownames(submatrix),]
    
    lonlat <- as.matrix(cbind(loc.info.variables$Lon, loc.info.variables$Lat))
    spatial <- dist(lonlat, method="euclidean")
    
    
    #environmental distances, current + LGM
    
    
    
    pcacurrent <-  prcomp(loc.info.variables[5:23], scale = T, center = T)
    biplot(pcacurrent)
    summary(pcacurrent)
    
    pcapast <-  prcomp(loc.info.variables[24:42], scale = T, center = T)
    summary(pcapast)
    
    currentclimate <- (dist(cbind(pcacurrent$x[,1:4]), method = "manhattan"))
    summary(currentclimate)
    
    pastclimate <- (dist(cbind(pcapast$x[,1:4]), method = "manhattan"))
    summary(pastclimate)
    
    for(p in predictor){
      if(p == "spatial") predictor <- spatial
      if(p == "current") predictor <- currentclimate
      if(p == "past") predictor <- pastclimate


      fname <- paste("D:/Dropbox/Africa2/decompo_beta", b, s, ".Rdata", sep=" ")

      load(fname)

decay_model_pow <- decay.model(decompo_beta$betadiv$PhyloSor, predictor, perm = 1000, model.type="power", y.type="dissimilarities")
decay_model_exp <- decay.model(decompo_beta$betadiv$PhyloSor, predictor, perm = 1000, model.type="exponential", y.type="dissimilarities")

bootcoefs_pow <- boot.coefs.decay(decay_model_pow, 1000)
bootcoefs_exp <- boot.coefs.decay(decay_model_exp, 1000)

createSheet(wb, name = paste("m_pow", b, s, p, sep=" "))
writeWorksheet(wb, decay_model_pow$model$aic, sheet = paste("m_pow", b, s, p, sep=" "), startRow = 1, startCol = 1)
writeWorksheet(wb, decay_model_pow$pseudo.r.squared, sheet = paste("m_pow", b, s, p, sep=" "), startRow = 1, startCol = 3)
writeWorksheet(wb, decay_model_pow$pseudo.p.value, sheet = paste("m_pow", b, s, p, sep=" "), startRow = 1, startCol = 6)

createSheet(wb, name = paste("b_pow", b, s, p, sep=" "))
writeWorksheet(wb, bootcoefs_pow$boot.coefs, sheet = paste("b_pow", b, s, p, sep=" "), startRow = 1, startCol = 1)
writeWorksheet(wb, bootcoefs_pow$original.coefs, sheet = paste("b_pow", b, s, p, sep=" "), startRow = 1, startCol = 5)
writeWorksheet(wb, bootcoefs_pow$mean.boot, sheet = paste("b_pow", b, s, p, sep=" "), startRow = 1, startCol = 8)
writeWorksheet(wb, bootcoefs_pow$sd.boot, sheet = paste("b_pow", b, s, p, sep=" "), startRow = 1, startCol = 11)

createSheet(wb, name = paste("m_exp", b, s, p, sep=" "))
writeWorksheet(wb, decay_model_exp$model$aic, sheet = paste("m_exp", b, s, p, sep=" "), startRow = 1, startCol = 1)
writeWorksheet(wb, decay_model_exp$pseudo.r.squared, sheet = paste("m_exp", b, s, p, sep=" "), startRow = 1, startCol = 3)
writeWorksheet(wb, decay_model_exp$pseudo.p.value, sheet = paste("m_exp", b, s, p, sep=" "), startRow = 1, startCol = 6)

createSheet(wb, name = paste("b_exp", b, s, p, sep=" "))
writeWorksheet(wb, bootcoefs_exp$boot.coefs, sheet = paste("b_exp", b, s, p, sep=" "), startRow = 1, startCol = 1)
writeWorksheet(wb, bootcoefs_exp$original.coefs, sheet = paste("b_exp", b, s, p, sep=" "), startRow = 1, startCol = 5)
writeWorksheet(wb, bootcoefs_exp$mean.boot, sheet = paste("b_exp", b, s, p, sep=" "), startRow = 1, startCol = 8)
writeWorksheet(wb, bootcoefs_exp$sd.boot, sheet = paste("b_exp", b, s, p, sep=" "), startRow = 1, startCol = 11)

saveWorkbook(wb)

}}}


