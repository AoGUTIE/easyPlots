#' Perform Gene Ontology (GO) Enrichment Analysis
#'
#' This function performs a Gene Ontology (GO) enrichment analysis on a set of differentially expressed genes (DEGs),
#' using the `clusterProfiler` package. It can optionally translate the resulting GO term descriptions and/or filter
#' the results to include only those associated with a specified gene list.
#'
#' @param data A data frame containing at least the columns: `log2FoldChange`, `pvalue`, `SYMBOL` or `ENSEMBL`, `geneID`, `ONTOLOGY`, and `Description`.
#' @param codeType A character string specifying the gene identifier type to use in the enrichment analysis: `"SYMBOL"` or `"ENSEMBL"`.
#' @param DEG A character string indicating which DEGs to analyze: `"up"` for upregulated, `"down"` for downregulated.
#' @param Output_Language A character string with the language code (e.g., `"es"`) to which GO term descriptions should be translated.
#' @param Input_Language A character string with the language code (e.g., `"en"`) of the input GO term descriptions.
#' @param Translate Logical (`TRUE` or `FALSE`); whether to translate the GO term descriptions from `Input_Language` to `Output_Language`.
#' @param FC Numeric; the absolute log2 fold-change threshold used to define upregulated and downregulated genes.
#' @param PV Numeric; the p-value threshold used to define significantly differentially expressed genes.
#' @param maxSize Integer; maximum number of GO terms to retain for each ontology category (BP, MF, CC).
#' @param Filtering Logical (`TRUE` or `FALSE`); whether to filter the GO results to include only terms associated with genes in `GeneList`.
#' @param GeneList A character vector of gene identifiers to filter the GO results against (only used if `Filtering = TRUE`).
#'
#' @return A named list of three data frames:
#' \describe{
#'   \item{BiologicalProcess}{Data frame of enriched Biological Processes (BP).}
#'   \item{MolecularFunctions}{Data frame of enriched Molecular Functions (MF).}
#'   \item{CellularComponents}{Data frame of enriched Cellular Components (CC).}
#' }
#'
#' @details
#' The function first filters the input `data` by the specified `FC` and `PV` thresholds to select significant DEGs.
#' Depending on the `DEG` parameter, it runs the GO enrichment analysis using the upregulated or downregulated gene set.
#' Optionally, it can translate GO term descriptions and/or filter the terms based on a user-provided `GeneList`.
#' The output consists of the top `maxSize` GO terms (sorted by `Count`) for each of the three ontology categories.
#'
#' @examples
#' \dontrun{
#' result <- ontologyAnalysis(
#'   data = my_data,
#'   codeType = "SYMBOL",
#'   DEG = "up",
#'   Output_Language = "es",
#'   Input_Language = "en",
#'   Translate = TRUE,
#'   FC = 1,
#'   PV = 0.05,
#'   maxSize = 10,
#'   Filtering = TRUE,
#'   GeneList = c("TP53", "BRCA1", "EGFR")
#' )
#'
#' result$BiologicalProcess
#' }
#'
#' @export
ontologyAnalysis <- function(data,
                             codeType,
                             DEG,
                             Output_Language,
                             Input_Language,
                             Translate = TRUE,
                             FC,
                             PV,
                             maxSize,
                             Filtering = TRUE,
                             GeneList) {

  # Filtrar genes downregulados
  downregulated <- data[data$log2FoldChange < -FC & data$pvalue < PV, ]
  upregulated   <- data[data$log2FoldChange >  FC & data$pvalue < PV, ]

  if (codeType == "SYMBOL") {
    upGenes   <- upregulated$SYMBOL
    downGenes <- downregulated$SYMBOL
  } else if (codeType == "ENSEMBL") {
    upGenes   <- upregulated$ENSEMBL
    downGenes <- downregulated$ENSEMBL
  }

  if (DEG == "up") {
    onto <- clusterProfiler::enrichGO(
      gene            = upGenes,
      OrgDb           = "org.Hs.eg.db",
      keyType         = codeType,
      pAdjustMethod   = "BH",
      ont             = "ALL",
      pvalueCutoff    = 0.05,
      minGSSize       = 5,
      maxGSSize       = 800
    )
    onto <- clusterProfiler::simplify(onto)
    onto <- as.data.frame(onto)
    ontoSort <- onto[order(onto$Count, decreasing = TRUE), ]

  } else if (DEG == "down") {
    onto <- clusterProfiler::enrichGO(
      gene            = downGenes,
      OrgDb           = "org.Hs.eg.db",
      keyType         = codeType,
      pAdjustMethod   = "BH",
      ont             = "ALL",
      pvalueCutoff    = 0.05,
      minGSSize       = 5,
      maxGSSize       = 800
    )
    onto <- clusterProfiler::simplify(onto)
    onto <- as.data.frame(onto)
    ontoSort <- onto[order(onto$Count, decreasing = TRUE), ]
  }

  # Filter if necessary
  if (Filtering) {
    ontoSort <- dplyr::filter(
      ontoSort,
      vapply(strsplit(ontoSort$geneID, "/"),
             function(genes) any(genes %in% GeneList),
             logical(1))
    )
  }

  # Translate if necessary
  if (Translate) {
    ontoSort$Description <- polyglotr::google_translate(
      ontoSort$Description,
      target_language = Output_Language,
      source_language = Input_Language
    )
  }

  # Analyses
  BP <- ontoSort[ontoSort$ONTOLOGY == "BP", ]
  BP <- utils::head(BP[order(BP$Count, decreasing = TRUE), ], maxSize)
  BP <- BP[order(BP$Count, decreasing = FALSE), ]
  BP$Position <- stringr::str_wrap(BP$Description, width = 25)
  unique_descriptions <- unique(BP$Position)
  BP$Position <- factor(BP$Position, levels = unique_descriptions)
  BP$title <- "BP"

  MF <- ontoSort[ontoSort$ONTOLOGY == "MF", ]
  MF <- utils::head(MF[order(MF$Count, decreasing = TRUE), ], maxSize)
  MF <- MF[order(MF$Count, decreasing = FALSE), ]
  MF$Position <- stringr::str_wrap(MF$Description, width = 25)
  unique_descriptions <- unique(MF$Position)
  MF$Position <- factor(MF$Position, levels = unique_descriptions)
  MF$title <- "MF"

  CC <- ontoSort[ontoSort$ONTOLOGY == "CC", ]
  CC <- utils::head(CC[order(CC$Count, decreasing = TRUE), ], maxSize)
  CC <- CC[order(CC$Count, decreasing = FALSE), ]
  CC$Position <- stringr::str_wrap(CC$Description, width = 25)
  unique_descriptions <- unique(CC$Position)
  CC$Position <- factor(CC$Position, levels = unique_descriptions)
  CC$title <- "CC"

  ontoResult <- list(
    BiologicalProcess    = BP,
    MolecularFunctions   = MF,
    CellularComponents   = CC
  )

  return(ontoResult)
}
