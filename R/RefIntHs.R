refCpx <-
  readRDS("Data/CORUM_03.09.2018_Hs_Complex_n_2357.RDS")
geneGO <-
  readRDS("Data/GeneGOLvl_Hs.RDS")
ppiB <-
  readRDS("Data/BioGrid_PPI_Hs_3.5.179_n_508899.RDS")
allGene <-
  unique(unlist(refCpx))
allGene <-
  allGene[allGene != "None"]
allS <-
  combn(allGene, 2)
allS <-
  apply(allS, 2, sort)
allPPI <-
  paste(allS[1, ], allS[2, ], sep = "~")
tpPPI <-
  unique(unlist(lapply(refCpx, function(prt){
    prt <- unique(prt)
    prt <- prt[prt != "None"]
    if(length(prt) > 1){
      s <- combn(prt, 2, FUN = sort)
      return(paste(s[1, ], s[2, ], sep = "~"))  
    }else{
      return(NA)
    }
  })))
tpPPI <- tpPPI[!is.na(tpPPI)]

rnPPI <-
  setdiff(allPPI, tpPPI)
shrPPI <-
  list()
for(goCat in c("BP", "CC", "MF")){
  goDat <-
    geneGO[[goCat]]
  goSub <-
    goDat[names(goDat) %in% allGene]
  goOv <-
    lapply(goSub, function(goI){
      lapply(goSub, function(goT){
        return(((length(intersect(goI, goT)))^2)/((length(unique(goI)))*(length(unique(goT)))))
      })
    })
  goInt <-
    lapply(goOv, function(geneOV){
      g <-
        unlist(geneOV)
      return(names(g[g!= 0]))
    })
  dInt <-
    apply(stack(goInt), 1, sort)
  shrPPI[[goCat]] <-
    paste(dInt[1, ], dInt[2, ], sep = "~")
}
tnPPI <-
  setdiff(rnPPI, unlist(shrPPI))
tnPPI <-
  setdiff(tnPPI, ppiB$PPI)
(length(tnPPI))/(length(tpPPI))
intersect(tnPPI, tpPPI)
saveRDS(list(`TP` = tpPPI, `TN` = tnPPI),
        file = "Data/RefIntHs.RDS")