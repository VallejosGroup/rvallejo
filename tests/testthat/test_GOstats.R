context("GOstats")

test_that("GOstats analysis works", {
  ## Random ENTREZ IDs (first n from the 10x PBMC v3 data).
  genes <- c("645520",
    "79501",
    "100996442",
    "729759",
    "81399",
    "400728",
    "79854",
    "284593",
    "100130417",
    "118424"
  )
  expect_error(
    expect_warning(
      go_analysis(
        genes = genes[1:5],
        universe = genes,
        annotation_package = "foobar"
      ),
      "there is no package called 'foobar'"
    ),
    "You must first install the annotation_package.")
  res <- go_analysis(
    genes = genes[1:5],
    universe = genes,
    ontologies = "BP",
    p_value = 0.5,
    annotation_package = "org.Hs.eg.db"
  )
  expect_is(res, "data.frame")
  expect_equal(res[[1, "Term"]], "sensory perception of chemical stimulus")

})
