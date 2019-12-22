library(GO.db)
getAllBPChildren <- function(goids){
  ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names = FALSE))
  ans <- ans[!is.na(ans)]
}
getAllCCChildren <- function(goids){
  ans <- unique(unlist(mget(goids, GOCCCHILDREN), use.names = FALSE))
  ans <- ans[!is.na(ans)]
}
getAllMFChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOMFCHILDREN), use.names = FALSE))
  ans <- ans[!is.na(ans)]
}
listGO <-
  list(`BP` = list(getAllBPChildren("GO:0008150")),
       `CC` = list(getAllCCChildren("GO:0005575")),
       `MF` = list(getAllMFChildren("GO:0003674")))
for(i in c(1:9)){
  listGO[["BP"]][[i+1]] <-
    getAllBPChildren(listGO[["BP"]][[i]])
  listGO[["CC"]][[i+1]] <-
    getAllCCChildren(listGO[["CC"]][[i]])
  listGO[["MF"]][[i+1]] <-
    getAllMFChildren(listGO[["MF"]][[i]])
}
nLvel <-
  unlist(
    lapply(listGO, function(GO){
      which.max(unlist(lapply(GO, length)))
    }))

library(org.Hs.eg.db)
geneGO <-
  lapply(names(listGO), function(nGO){
    gLv <-
      mget(intersect(listGO[[nGO]][[nLvel[nGO]]], keys(org.Hs.egGO2EG)),
           org.Hs.egGO2EG)
    d <-
      stack(gLv)
    return(unstack(d[, c(2, 1)]))
  })
names(geneGO) <-
  names(listGO)
saveRDS(geneGO, file = "Data/GeneGOLvl_Hs.RDS")

library(org.Mm.eg.db)
geneGO <-
  lapply(names(listGO), function(nGO){
    gLv <-
      mget(intersect(listGO[[nGO]][[nLvel[nGO]]], keys(org.Mm.egGO2EG)),
           org.Mm.egGO2EG)
    d <-
      stack(gLv)
    return(unstack(d[, c(2, 1)]))
  })
names(geneGO) <-
  names(listGO)
saveRDS(geneGO, file = "Data/GeneGOLvl_Mm.RDS")

library(org.Ce.eg.db)
geneGO <-
  lapply(names(listGO), function(nGO){
    gLv <-
      mget(intersect(listGO[[nGO]][[nLvel[nGO]]], keys(org.Ce.egGO2EG)),
           org.Ce.egGO2EG)
    d <-
      stack(gLv)
    return(unstack(d[, c(2, 1)]))
  })
names(geneGO) <-
  names(listGO)
saveRDS(geneGO, file = "Data/GeneGOLvl_Ce.RDS")

library(org.Dm.eg.db)
geneGO <-
  lapply(names(listGO), function(nGO){
    gLv <-
      mget(intersect(listGO[[nGO]][[nLvel[nGO]]], keys(org.Dm.egGO2EG)),
           org.Dm.egGO2EG)
    d <-
      stack(gLv)
    return(unstack(d[, c(2, 1)]))
  })
names(geneGO) <-
  names(listGO)
saveRDS(geneGO, file = "Data/GeneGOLvl_Dm.RDS")

library(org.Sc.sgd.db)
geneGO <-
  lapply(names(listGO), function(nGO){
    gLv <-
      mget(intersect(listGO[[nGO]][[nLvel[nGO]]], keys(org.Sc.sgdGO2ORF)),
           org.Sc.sgdGO2ORF)
    d <-
      stack(gLv)
    return(unstack(d[, c(2, 1)]))
  })
names(geneGO) <-
  names(listGO)
saveRDS(geneGO, file = "Data/GeneGOLvl_Sc.RDS")
