#' @export
go_analysis <- function(
    entrez_genes,
    entrez_universe,
    annotation_package,
    ontologies = c("BP", "CC", "MF"),
    min_size = 3,
    mdlinks = TRUE,
    p_value = 0.01
    ) {

  ontologies <- match.arg(ontologies, several.ok = TRUE)
  ghgp <- new(
    "GOHyperGParams",
    geneIds = entrez_genes,
    annotation = annotation_package,
    universeGeneIds = entrez_universe,
    conditional = TRUE,
    ontology = ontologies[[1]],
    minSizeCutoff = min_size,
    pvalueCutoff = p_value
  )
  go_res <- lapply(ontologies, 
    function(ontology) {
      ontology(ghgp) <- ontology
      hyperGTest(ghgp)
    }
  )
  out <- do.call(rbind, lapply(go_res, process_res))
  out[order(out[["Pvalue"]]), ]
}



process_go_res <- function(res, mdlinks = TRUE) {
  tab <- summary(res)
  ontology <- gsub("GO(\\w+)ID", "\\1", colnames(tab)[[1]])
  colnames(tab)[[1]] <- "GOID"
  if (mdlinks) {
    tab[["GOID"]] <- mdlink(
      tab[["GOID"]],
      sprintf("http://amigo.geneontology.org/amigo/term/%s", tab[["GOID"]])
    )
  }
  tab[["Ontology"]] <- rep(ontology, nrow(tab))
  tab
}

mdlink <- function(text, link) {
  sprintf("[%s](%s)", text, link)
}
