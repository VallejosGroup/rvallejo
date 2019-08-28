#' Run GOstats functional enrichment analysis for a list of genes.
#' 
#' Runs the conditional GOstats algorithm for a set of genes in a given 
#' universe. This function runs the algorithm for multiple ontologies, and returns 
#' a table of merged results.
#' 
#' @param genes The genes of interest; ie the balls drawn from the urn in the
#' classical hypergeometric analogy.
#' @param universe The background set of genes. This corresponds to the bag
#' of balls in the classical hypergeometric analogy.
#' @param annotation_package An organism annotation package used by the GOstats
#' functions internally to annotate genes with GO terms and to traverse the GO 
#' hierarchy. eg, for humans, org.Hs.eg.db would be suitable and would require
#' Entrez IDs for the \code{genes} and \code{universe} arguments.
#' @param ontologies The ontologies to test. Options are "BP" 
#' (biological process), "CC" (cellular component) and "MF" 
#' (molecular function).
#' @param min_size The minimum size of GO terms to include. It is advised to set
#' this to at least 2; otherwise you will observe quite a few terms which are 
#' comprised of only a single gene.
#' @param mdlinks Should the GO IDs be hyperlinked using Markdown in the output
#' table?
#' @param p_value p-value cutoff. This is important to set a priori, as it 
#' affects the conditional algorithm.
#' @export
go_analysis <- function(
    genes,
    universe,
    annotation_package,
    ontologies = c("BP", "CC", "MF"),
    min_size = 3,
    mdlinks = TRUE,
    p_value = 0.01
    ) {

  ontologies <- match.arg(ontologies, several.ok = TRUE)
  if (!require(annotation_package, character.only = TRUE)) {
    stop(
      paste0(
        "You must first install the annotation_package.
        Try BiocManager::install(\"", annotation_package, "\")")
    )
  }
  ghgp <- new(
    "GOHyperGParams",
    geneIds = genes,
    annotation = annotation_package,
    universeGeneIds = universe,
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
  out <- do.call(rbind, lapply(go_res, process_go_res))
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
