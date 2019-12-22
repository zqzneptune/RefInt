version <-
  "3.5.179" # Update monthly
fn <-
  paste0("BIOGRID-ORGANISM_", version, ".zip")
url <-
  paste0("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-", version, "/BIOGRID-ORGANISM-", version, ".tab2.zip")
download.file(url, 
              destfile = file.path("Raw", fn))
unzip(file.path("Raw", fn), 
      exdir = "/tmp/BIOGRID")
fnList <-
  c("Hs" = "Homo_sapiens",
    "Mm" = "Mus_musculus",
    "Ce" = "Caenorhabditis_elegans",
    "Dm" = "Drosophila_melanogaster",
    "Sc" = "Saccharomyces_cerevisiae_S288c")
fn <- fnList[1]
for(i in seq_len(length(fnList))){
  fn <-
    fnList[i]
  B <- 
    read.delim(paste0("/tmp/BIOGRID/BIOGRID-ORGANISM-",
                      fn,
                      "-", version, ".tab2.txt"), stringsAsFactors = FALSE)
  ppiB <-
    B[B[, 13] == "physical", ]
  ppiB[, 2] <-
    as.character(ppiB[, 2])
  ppiB[, 3] <-
    as.character(ppiB[, 3])
  s <-
    apply(ppiB[, c(2, 3)], 1, sort)
  ppiB[, "PPI"] <-
    paste(s[1, ], s[2, ], sep = "~")
  saveRDS(ppiB, 
          file = paste0("Data/BioGrid_PPI_", 
                        names(fn),
                        "_",
                        version,
                        "_n_",
                        nrow(ppiB),
                        ".RDS"))  
}
