version <- # 03.09.2018 Corum 3.0 current release.
  "03.09.2018"
url <-
  "http://mips.helmholtz-muenchen.de/corum/download/"
fn <-
  c("all" = 
         "allComplexes.txt.zip",
    "core" = 
         "coreComplexes.txt.zip")
i <- 2
download.file(paste0(url, fn[i]), 
              destfile = file.path("Raw", fn[i]))
unzip(file.path("Raw", fn[i]), 
      exdir = "Raw")
rawFile <-
  read.delim(gsub(".zip", "", file.path("Raw", fn[i])), 
             stringsAsFactors = FALSE)
speList <-
  c("Hs" = "Human", "Mm" = "Mouse")
for(i in seq_len(length(speList))){
  species <-
    speList[i]
  rawDat <-
    rawFile[rawFile[, 3] == species, ] # The third column: Organism
  rawCpx <-
    strsplit(rawDat[, 7], ";")
  names(rawCpx) <-
    rawDat[, 2]
  outCpx <-
    rawCpx[unlist(lapply(rawCpx, length)) > 1]
  saveRDS(outCpx,
          file = paste0("Data/CORUM_", version, 
                        "_", names(species), 
                        "_Complex_n_",
                        length(outCpx),
                        ".RDS"))
}



ComplexToKeep <-
  unlist(lapply(strsplit(Dat$`subunits(UniProt IDs)`, ";"), length)) > 1
Complex <-
  Dat[ComplexToKeep, c("ComplexName", "subunits(Entrez IDs)" )]
colnames(Complex) <-
  c("complexName", "subunitID")
listCpx <- 
  strsplit(Complex$subunitID, ";")
names(listCpx) <-
  Complex$complexName
ppiCpx <-
  lapply(listCpx, function(cpx){
    m <- 
      apply(combn(cpx, 2), 2, sort)
    s <- 
      paste(m[1, ], m[2, ], sep = "~")
    return(s)
  })
names(ppiCpx) <-
  Complex$complexName
datTP <-
  aggregate(data = stack(ppiCpx),
            `ind` ~ `values`,
            FUN = paste,
            collapse = ";")
colnames(datTP) <-
  c("PPI", "ComplexNames")  

prtSpace <-
  unique(unlist(listCpx))
mSpace <-
  combn(prtSpace, 2)
s <-
  apply(mSpace, 2, sort)
ppiSpace <-
  paste(s[1, ], s[2, ], sep = "~")
ppiTN <-
  setdiff(ppiSpace, datTP$PPI)
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
genGO <-
  select(org.Hs.eg.db, 
         keys = prtSpace, 
         keytype = c("ENTREZID"),
        columns=c("GO", "ONTOLOGY"))
genGO <-
  genGO[genGO$EVIDENCE != "IEA", ]

g <-
  unique(genGO[, c("ENTREZID", "GO")])
g$isGO <- 1
datG <-
  tidyr::spread(g, `ENTREZID`, `isGO`)
datG <-
  datG[!is.na(datG$GO), ]
mG <-
  as.matrix(datG[, -1])
rownames(mG) <-
  datG$GO
mG[is.na(mG)] <- 0
tMG <-
  t(mG)
library(philentropy)
mJaccard <-
  distance(tMG, method = "jaccard")
mDice <-
  distance(tMG, method = "dice")
# mJaccard[1:5, 1:5]
# mJaccard[lower.tri(mJaccard, diag = FALSE)][1:5]
genePair <- 
  t(combn(rownames(tMG), 2))
colnames(genePair) <-
  c("EntrezA", "EntrezB")
disGO <-
  data.frame(genePair,
             `Jaccard` = mJaccard[lower.tri(mJaccard, diag = FALSE)],
             `Dice` = mDice[lower.tri(mDice, diag = FALSE)],
             stringsAsFactors = FALSE)
# Not necessary, as rownames(mG) already sort
# s <-
#   apply(disGO[, c("EntrezA", "EntrezB")], 1, sort)
disGO[, "PPI"] <-
  paste(disGO$EntrezA, disGO$EntrezB, sep = "~")

datTPdisGO <-
  left_join(datTP, disGO, by = "PPI")

datTNdisGO <-
  disGO[disGO$PPI %in% ppiTN, ]
boxplot(list(`TP` = datTPdisGO$Dice,
             `TN` = datTNdisGO$Dice),outline = FALSE)
hist(datTNdisGO$Jaccard)
hist(datTPdisGO$Jaccard)
saveRDS(list(`Complex` = listCpx,
             `datTPdisGO` = datTPdisGO,
             `datTNdisGO` = datTNdisGO),
        file = "Hs/Refset_coreCORUM_Entrez.RDS")
