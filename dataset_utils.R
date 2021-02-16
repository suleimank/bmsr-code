download.and.prepare.datasets <- function() {

  # Download GDSC and CCLE data
  ## set1.name <- "GDSC_2020(v1-8.2)"
  ## set2.name <- "GDSC_2020(v2-8.2)"  
  set1.name <- "GDSC_2020(v2-8.2)"  
  set2.name <- "CCLE_2015"
  cat(paste0("Downloading ", set1.name, " data. This will take a while\n"))
  set1 <- downloadPSet(set1.name, saveDir = ".")

  cat(paste0("Downloading ", set2.name, " data. This will take a while\n"))
  set2 <- downloadPSet(set2.name, saveDir = ".")

  # Restrict to breast cancer cell lines  
  cInfo <- cellInfo(set1)
  ## flag <- cInfo$Cancer.Type...matching.TCGA.label. %in% c("LAML")
  ## flag <- cInfo$Cancer.Type...matching.TCGA.label. %in% c("LUAD", "LUSC")
  flag <- cInfo$Cancer.Type...matching.TCGA.label. %in% c("BRCA")
  set1_cell_lines = cellNames(set1)[flag]

  cInfo <- cellInfo(set2)
  ## flag <- cInfo$Cancer.Type...matching.TCGA.label. %in% c("LAML")  
  ## flag <- cInfo$Hist.Subtype1 %in% c("acute_myeloid_leukaemia")
  ## flag <- cInfo$tissueid == "Lung"
  flag <- cInfo$tissueid == "Breast"
  set2_cell_lines = cellNames(set2)[flag]

  commonGenes <- intersect(fNames(set1, "rna"),
                           fNames(set2, "rna"))


  common <- intersectPSet(list('set2'=set2, 'set1'=set1),
                          intersectOn=c("drugs"),
                          strictIntersect=FALSE)
  
  set1.resp <-
    summarizeSensitivityProfiles(common$set1,
                                 cell.lines = set1_cell_lines,    
                                 sensitivity.measure='aac_recomputed',
                                 summary.stat="median",
                                 verbose=FALSE)
  
  set2.resp <-
    summarizeSensitivityProfiles(common$set2,
                                 cell.lines = set2_cell_lines,
                                 sensitivity.measure='aac_recomputed',
                                 summary.stat="median",
                                 verbose=FALSE)
  			       
  set1.expression <-
    summarizeMolecularProfiles(common$set1,
                               set1_cell_lines,
                               mDataType="rna",
                               features=commonGenes,
                               verbose=FALSE)
  
  set2.expression <-
    summarizeMolecularProfiles(common$set2,
                               set2_cell_lines,
                               mDataType="rna",
                               features=commonGenes,
                               verbose=FALSE)

  datasets <- list("set1" = list("expr" = assay(set1.expression, "exprs"),
                                 "response" = set1.resp),
                   "set2" = list("expr" = assay(set2.expression, "exprs"),
                                 "response" = set2.resp))

  return(datasets)
}

