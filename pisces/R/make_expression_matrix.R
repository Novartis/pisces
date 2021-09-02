#!/usr/bin/env Rscript

getScriptPath <- function() {
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, pe = TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if (length(script.dir) == 0) 
    stop("can't determine script dir: please call the script with Rscript")
  if (length(script.dir) > 1) 
    stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
script.dir <- getScriptPath()
wd.old <- getwd()
setwd(script.dir)
renv::activate() 
setwd(wd.old)
suppressMessages(library(tximport))
suppressMessages(library(docopt))
args <- docopt("PISCES summary expression matrix and differential expression

Usage: summarize-expression [options] [--exclude-genes=GENE]... [DIR DIR...]

Options:
  -n NAME, --name NAME                                          Output file base name [default: expression_matrix]
  -i IDX, --salmon-index SALMON_INDEX                           PISCES index to aggregate (default: defined in config.json)
  -t TYPE, --sample-type TYPE                                   Sample type (as specified in config.json) [default: human]
  -m META, --metadata METADATA_DIR                              CSV file describing contrast variables and sample names
  -r VAR, --group-by VAR                                        Column name describing variable to group samples for normalization (e.g. cell line or timepoint)
  -b VAR, --norm-by VAR                                         Column name of the main variable used for within-group normalization (e.g. treatment)
  -c FACTOR, --control-factor FACTOR                            Name of factor in '--norm-by' column used for within-group normalization (e.g. DMSO)
  -d PATSY, --deseq-formula PATSY                               `patsy` notation to be passed to DESeq2 e.g: ~ treatment + treatment:timepoint (derived from columns in `--metadata`)
  -f CSV, --contrasts CSV                                       CSV file defining contrasts in three column format (term, treatment, control)
  -s BIOTYPE, --scale-tpm BIOTYPE                               TMM normalize using genes belonging to this ENSEMBL `biotype` [default: protein_coding]
  -e TPM_CUTTOFF, --quartile-expression TPM_CUTTOFF             Exclude genes from TMM normalization and DESeq2 analysis that have expression lower than TPM_CUTTOFF in less than 25% of samples [default: 1]
  -l BP, --min-transcript-length BP                             Exclude transcripts that have a length less than --min-transcript-length [default: 0]
  -x GENE, --exclude-genes GENE                                 List of genes to exclude from TMM normalization
  -o PREFIX, --only-transcripts PREFIX                          Only summarize transcripts starting with PREFIX e.g: ENST
  --allow-missing                                               Allow summary of data when some samples are missing. Missing samples will have NA.
  --no-rescale                                                  Do not re-scale TPM of --scale-tpm BIOTYPE to 1e6
  --isoform-only                                                Write only isoform data, not gene level data
  --debug                                                       Print debugging information
  --config CONFIG                                               Internal use only

Arguments:
  DIR      Directories containing `pisces run` analysis results")

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(reshape))
suppressMessages(library(magrittr))
suppressMessages(library(edgeR))
suppressMessages(library(jsonlite))

if (args[["--debug"]]) {
  cat(str(args))
  options(error = quote(dump.frames("pisces_fault", TRUE)))
}
NormalizeSamples <- function(abundances, metadata, args) {
  groupingVar <- args[["--group-by"]]
  normVar <- args[["--norm-by"]]
  controlGroup <- args[["--control-factor"]]
  if (is.data.frame(metadata)) {
    if (is.character(groupingVar) & is.character(normVar)) {
      # --group-by --norm-by
      for (group in unique(metadata[[groupingVar]])) {
        controls <- as.character(metadata$SampleID[which(metadata[[normVar]] == 
          controlGroup & metadata[[groupingVar]] == group)])
        experimentals <- as.character(metadata$SampleID[which(metadata[[normVar]] != 
          controlGroup & metadata[[groupingVar]] == group)])
        message(paste("Group:", group))
        message("Controls:", paste(controls, collapse = ", "))
        message("Treated:", paste(experimentals, collapse = ", "))
        if (length(controls) == 0) {
          message(paste("No controls for group", group, "so median normalizing to treated samples."))
          if (length(experimentals) == 1) {
          abundances[, c(controls, experimentals)] <- log2(abundances[, 
            c(controls, experimentals)] + 1) - log2(sapply(abundances[, 
            c(experimentals)] + 1, median))
          } else {
          abundances[, c(controls, experimentals)] <- log2(abundances[, 
            c(controls, experimentals)] + 1) - log2(apply(abundances[, 
            c(experimentals)] + 1, 1, median))
          }
        } else if (length(controls) == 1) {
          abundances[, c(controls, experimentals)] <- log2(abundances[, c(controls, 
          experimentals)] + 1) - log2(sapply(abundances[, c(controls)] + 
          1, median))
        } else {
          abundances[, c(controls, experimentals)] <- log2(abundances[, c(controls, 
          experimentals)] + 1) - log2(apply(abundances[, c(controls)] + 
          1, 1, median))
        }
      }
    } else if (is.character(groupingVar) & !is.character(normVar)) {
      # --group-by
      for (group in unique(metadata[[groupingVar]])) {
        controls <- as.character(metadata$SampleID[which(metadata[[groupingVar]] == 
          group)])
        message(paste("Group:", group))
        message("Members:", paste(controls, collapse = ", "))
        message("Normalizing to median of within-group samples.")
        if (length(controls) == 1) {
          abundances[, c(controls)] <- log2(abundances[, c(controls)] + 1) - 
          log2(sapply(abundances[, c(controls)] + 1, median))
        } else {
          abundances[, c(controls)] <- log2(abundances[, c(controls)] + 1) - 
          log2(apply(abundances[, c(controls)] + 1, 1, median))
        }
      }
    } else if (!is.character(groupingVar) & is.character(normVar)) {
      # --norm-by
      controls <- as.character(metadata$SampleID[which(metadata[[normVar]] == 
        controlGroup)])
      experimentals <- as.character(metadata$SampleID[which(metadata[[normVar]] != 
        controlGroup)])
      message("Controls:", paste(controls, collapse = ", "))
      message("Treated:", paste(experimentals, collapse = ", "))
      if (length(controls) == 1) {
        abundances[, c(controls, experimentals)] <- log2(abundances[, c(controls, 
          experimentals)] + 1) - log2(sapply(abundances[, c(controls)] + 
          1, median))
      } else {
        abundances[, c(controls, experimentals)] <- log2(abundances[, c(controls, 
          experimentals)] + 1) - log2(apply(abundances[, c(controls)] + 1, 
          1, median))
      }
    } else {
      # normalize to median of all samples
      message("Normalizing to median of all samples.")
      abundances <- log2(abundances + 1) - log2(apply(abundances + 1, 1, median))
    }
  } else {
    # no metadata file
    message("Normalizing to median of all samples.")
    abundances <- log2(abundances + 1) - log2(apply(abundances + 1, 1, median))
  }
  return(abundances)
}

rescaleTPM <- function(df, columns = colnames(df), rows = 1:nrow(df), precision = 3) {
    sums.before.norm <- apply(df[rows, columns], 2, sum)
    tpm.scale.factors <- 1e+06/sums.before.norm
    scaled_df <- df
    scaled_df[rows, columns] <- data.frame(mapply("*", scaled_df[rows, 
        columns], tpm.scale.factors))
    scaled_df[, columns] <- round(scaled_df[, columns], precision)
    return(scaled_df)
    }
    
Summarize <- function(txi, tx2gene, annotation, args, metadata, species) {
    message(paste("Summarizing", species, "transcripts..."))
    
    colnames(tx2gene) <- c("transcript_id", "gene_id")
    tx.idx <- which(rownames(txi$abundance) %in% tx2gene$transcript_id)
    intron.idx <- which(stringr::str_detect(rownames(txi$abundance), "intronic_.+"))
    intergene.idx <- which(stringr::str_detect(rownames(txi$abundance), "intergene_.+"))
    
    message("Writing isoform level tables...")
    as.data.frame(txi$abundance[tx.idx,]) %>% 
        rescaleTPM() %>%
        tibble::rownames_to_column("transcript_id") %>% 
        left_join(tx2gene) %>%
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".isoforms.TPM.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
    
    as.data.frame(txi$counts[tx.idx,]) %>% 
        tibble::rownames_to_column("transcript_id") %>% 
        left_join(tx2gene) %>%
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".isoforms.counts.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
        
    as.data.frame(txi$length[tx.idx,]) %>% 
        tibble::rownames_to_column("transcript_id") %>% 
        left_join(tx2gene) %>%
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".isoforms.length.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
    
    message("Writing intron level tables...")
    as.data.frame(txi$abundance[intron.idx,]) %>% 
        rescaleTPM() %>%
        tibble::rownames_to_column("gene_id") %>% 
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".introns.TPM.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
    
    as.data.frame(txi$counts[intron.idx,]) %>% 
        tibble::rownames_to_column("gene_id") %>% 
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".introns.counts.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
        
    as.data.frame(txi$length[intron.idx,]) %>% 
        tibble::rownames_to_column("gene_id") %>% 
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".introns.length.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
        
    message("Writing intergene level tables...")
    as.data.frame(txi$abundance[intergene.idx,]) %>% 
        rescaleTPM() %>%
        tibble::rownames_to_column("gene_id") %>% 
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".intergenes.TPM.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
    
    as.data.frame(txi$counts[intergene.idx,]) %>% 
        tibble::rownames_to_column("gene_id") %>% 
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".intergenes.counts.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
        
    as.data.frame(txi$length[intergene.idx,]) %>% 
        tibble::rownames_to_column("gene_id") %>% 
        as.data.frame() %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".intergenes.length.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
        
    message("Building transcript QC measure table...")
    pct_intronic <- colSums(txi$counts[intron.idx,]) / colSums(txi$counts[c(tx.idx, intron.idx, intergene.idx),]) * 100
    pct_intergenic <- colSums(txi$counts[intergene.idx,]) / colSums(txi$counts[c(tx.idx, intron.idx, intergene.idx),]) * 100
    pct_genic <- colSums(txi$counts[tx.idx,]) / colSums(txi$counts[c(tx.idx, intron.idx, intergene.idx),]) * 100
    
    rbind(pct_genic, pct_intronic, pct_intergenic) -> qc_table
    rownames(qc_table) <- c("exonic", "intronic", "intergenic")
    qc_table %>%
        as.data.frame() %>%
        tibble::rownames_to_column("compartment") %>%
        format(nsmall = 3, scientific = F, trim = T) %>%
        write.table(paste0(args[["--name"]], ".", species, ".transcript.qc.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
        
  if (!args[["--isoform-only"]]) {
    txi.gene <- summarizeToGene(txi, tx2gene = tx2gene, ignoreTxVersion = F)
            
    message("Writing annotation data matrix...")
    write.table(as.data.frame(annotation), paste0(args[["--name"]], ".", species, 
      ".annotation.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
      
    message("Making gene-level length matrix...")  
    as.data.frame(txi.gene$abundance) %>% tibble::rownames_to_column("gene_id") %>% 
      left_join(annotation) %>% as.data.frame %>%
      write.table(format(., nsmall = 3, scientific = F, trim = T), paste0(args[["--name"]], 
      ".", species, ".length.txt"), quote = FALSE, sep = "\t", row.names = F, 
      col.names = T)
    
    message("Making gene-level TPM matrix...")
    raw_df <- as.data.frame(txi.gene$abundance) %>% tibble::rownames_to_column("gene_id") %>% 
      left_join(annotation) %>% as.data.frame
   
    if (!args[["--no-rescale"]]) {
      message(paste("Scaling TPM of", args[["--scale-tpm"]], "to 1e6..."))
        to.rescale <- which(raw_df[["biotype"]] == args[["--scale-tpm"]])
        abundance.cols <- colnames(raw_df)[which(colnames(raw_df) %in% colnames(as.data.frame(txi.gene$abundance)))]
        scaled_df <- rescaleTPM(raw_df, abundance.cols, to.rescale)
    }
    
    write.table(format(scaled_df, nsmall = 3, scientific = F, trim = T), paste0(args[["--name"]], 
      ".", species, ".TPM.txt"), quote = FALSE, sep = "\t", row.names = F, 
      col.names = T)
    normed_df <- scaled_df
    normed_df[, abundance.cols] <- NormalizeSamples(normed_df[, abundance.cols], 
      metadata, args)
    normed_df[, abundance.cols] <- round(normed_df[, abundance.cols], 3)
    message("Making log2 fold change matrix...")
    write.table(format(normed_df, nsmall = 3, scientific = F, trim = T), paste0(args[["--name"]], 
      ".", species, ".log2fc.txt"), quote = FALSE, sep = "\t", row.names = F, 
      col.names = T)

    
    
    message(paste("Calculating TMM scaling factors using", args[["--scale-tpm"]], 
      "genes."))
    ribo.genes <- which(grepl("^RP[SL]", scaled_df[, "symbol"], ignore.case = T))
    


    first_quartile <- function(x) { 
        # > first_quartile(c(0,0,0,0,0,0,0,1,10,10))
        # [1] 1
        y <- quantile(x, c(0.25, 0.5, 0.75), type=1) 
        return(y[[3]])
        }
    
    message(paste("Excluding", length(ribo.genes), "ribosomal genes from TMM scaling factor calulation."))
    failing.quartile.filter <- which(apply(scaled_df[, abundance.cols], 1, first_quartile) < 
      as.numeric(args[["quartile-expression"]]))
    message(paste(length(failing.quartile.filter), "genes failing --quartile-expression cutoff"))
    if (is.character(args[["--exclude-genes"]])) {
      exclusion.list <- which(!(scaled_df[, "symbol"] %in% args[["--exclude-genes"]]))
      message(paste("Excluding", args[["--exclude-genes"]], "from TMM scaling factor calulation."))
      message(paste(length(exclusion.list), "genes in --exclude-genes list"))
    } else {
      exclusion.list <- c()
    }
    
    if (args[["--debug"]]) {
      message("ribo.genes:")
      print(head(scaled_df[ribo.genes, "symbol"]))
      message("median.filter:")
      print(head(scaled_df[failing.quartile.filter, "symbol"]))
      message("exclusion.list:")
      print(head(scaled_df[exclusion.list, "symbol"]))
    }
    tmm.genes <- setdiff(to.rescale, ribo.genes)
    tmm.genes <- setdiff(tmm.genes, failing.quartile.filter)
    tmm.genes <- setdiff(tmm.genes, exclusion.list)
    
    
    sums.before.norm <- apply(raw_df[tmm.genes, abundance.cols], 2, sum)
    tpm.scale.factors <- 1e+06/sums.before.norm
    scaled_df <- raw_df
    if (!args[["--no-rescale"]]) {
      scaled_df[tmm.genes, abundance.cols] <- data.frame(mapply("*", scaled_df[tmm.genes, 
        abundance.cols], tpm.scale.factors))
    }
    message(paste("Using", length(tmm.genes), "genes for TMM factor calculation"))
    tmm.factors <- calcNormFactors(scaled_df[tmm.genes, abundance.cols])
    tmm_df <- scaled_df
    tmm_df[to.rescale, abundance.cols] <- data.frame(mapply("/", tmm_df[to.rescale, 
      abundance.cols], tmm.factors))
    tmm_df[, abundance.cols] <- round(tmm_df[, abundance.cols], 3)
    write.table(format(tmm_df, nsmall = 3, scientific = F, trim = T), paste0(args[["--name"]], 
      ".", species, ".TPM.TMM-scaled.txt"), quote = FALSE, sep = "\t", row.names = F, 
      col.names = T)
    normed_tmm_df <- tmm_df
    normed_tmm_df[, abundance.cols] <- NormalizeSamples(normed_tmm_df[, abundance.cols], 
      metadata, args)
    normed_tmm_df[, abundance.cols] <- round(normed_tmm_df[, abundance.cols], 
      3)
    write.table(format(normed_tmm_df, nsmall = 3, scientific = F, trim = T), 
      paste0(args[["--name"]], ".", species, ".log2fc.TMM-scaled.txt"), quote = FALSE, 
      sep = "\t", row.names = F, col.names = T)
    
    message("Writing counts matrix...")
    as.data.frame(txi.gene$counts) %>% 
        tibble::rownames_to_column("gene_id") %>% 
        inner_join(annotation) %>%
        write.table(paste0(args[["--name"]], ".", species, ".counts.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
        
    message("Writing effective length matrix...")
    as.data.frame(txi.gene$length) %>% 
        tibble::rownames_to_column("gene_id") %>% 
        inner_join(annotation) %>%
        write.table(paste0(args[["--name"]], ".", species, ".length.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
    
    if (is.character(args[["--deseq-formula"]])) {
      message(paste("DEseq formula:", args[["--deseq-formula"]]))
      formula <- as.formula(args[["--deseq-formula"]])
      terms <- rownames(attr(terms(formula), "factors"))
      } else {
        formula <- NULL
        terms <- NULL
        }
      
    if (is.character(args[["--contrasts"]])) {
        contrasts <- read.csv(args[["--contrasts"]], header = F)
        colnames(contrasts) <- c("Factor", "Experimental", "Control")
        if (is.null(formula)) {
            terms <- levels(contrasts$Factor)
            }
    } else {
        contrasts <- NULL
    }
    
    if (!is.null(formula) | !is.null(contrasts)) {
      tidy.table.file <- paste0(args[["--name"]], ".", species, ".deseq.tidy.txt")
      if (file.exists(tidy.table.file)) {
        file.remove(tidy.table.file)
      }
      write.table(data.frame(Contrast = character(), gene_id = character(), 
        log10p = double(), log2fc = double(), baseMean = double(), log2fcSE = double(), stat = double()), 
        tidy.table.file, quote = FALSE, sep = "\t", row.names = F, col.names = T, 
        append = T)
      message("Making DESeq object...")
      suppressMessages(library(DESeq2))
      contrast.metadata <- as.data.frame(metadata[, c(terms, "SampleID")])
      
      print(contrast.metadata)
      message(paste("Contrasts specified in ", args[["--contrasts"]], ":", 
        sep = ""))
      if (is.data.frame(contrasts)) {
        print(contrasts)
        contrasts <- dplyr::group_by(contrasts, Factor)
        contrasts <- dplyr::group_split(contrasts)
        lapply(contrasts, deseq_analysis, txi.gene, contrast.metadata, formula, failing.quartile.filter, tidy.table.file)
      }
      }
  }
  
  
}

deseq_analysis <- function(contrasts, txi.gene, contrast.metadata, formula, failing.quartile.filter, tidy.table.file) {
  contrasts <- as.data.frame(contrasts)
  if (!is.null(formula)) {
    deseq.dataset <- DESeqDataSetFromTximport(txi.gene, contrast.metadata, formula)
    } else {
    deseq.dataset <- DESeqDataSetFromTximport(txi.gene, contrast.metadata, as.formula(paste0("~", contrasts$Factor[1])))}
  message(paste("Filtering", length(failing.quartile.filter), "genes failing --quartile-expression cutoff from DESeq2 dataset."))
  deseq.dataset <- deseq.dataset[-failing.quartile.filter,]
  message("Running DESeq2...")
  
  deseq.dataset <- DESeq(deseq.dataset)
  
  message("Contrast terms:")
  results_names <- resultsNames(deseq.dataset)
  print(results_names)
  
  message("Generating log2fc and p-values for contrasts:")
  
  if (!is.null(contrasts)) {
    for (i in 1:nrow(contrasts)) {
      term <- trimws(as.character(contrasts[i, 1]))
      exp_factor <- trimws(as.character(contrasts[i, 2]))
      con_factor <- trimws(as.character(contrasts[i, 3]))
      write.contrasts(exp_factor, con_factor, term, deseq.dataset, tidy.table.file)
    }
  } else {
    for (term in terms) {
      variable <- metadata[, term]
      con_factor <- levels(variable)[1]
      for (exp_factor in levels(variable)[-1]) term <- trimws(as.character(contrasts[i, 
      1]))
      exp_factor <- trimws(as.character(contrasts[i, 2]))
      con_factor <- trimws(as.character(contrasts[i, 3]))
      write.contrasts(exp_factor, con_factor, term, deseq.dataset, tidy.table.file)
    }
  }
}

write.contrasts <- function(exp_factor, con_factor, term, deseq.dataset, tidy.table.file) {
    message(paste(exp_factor, con_factor, sep = " vs "))
    deseq_df <- as.data.frame(results(deseq.dataset, contrast = c(term, 
      exp_factor, con_factor)))
    transformed_results <- list(log10p = log10(deseq_df$padj + 1e-299), 
      log2fc = deseq_df$log2FoldChange, baseMean = deseq_df$baseMean, 
      log2fcSE = deseq_df$lfcSE, stat = deseq_df$stat)
    padj_name <- paste(c("log10p", paste(exp_factor, con_factor, sep = " vs ")), 
      collapse = "_")
    lfc_name <- paste(c("log2fc", paste(exp_factor, con_factor, sep = " vs ")), 
      collapse = "_")
    transformed_results <- setNames(transformed_results, c(padj_name, 
      lfc_name, "baseMean", "log2FCSE"))
    write.table(format(data.frame(Contrast = paste(exp_factor, con_factor, 
      sep = " vs "), gene_id = rownames(deseq_df), log10p = transformed_results[[1]], 
      log2fc = transformed_results[[2]], baseMean = transformed_results[[3]], 
      log2fcSE = transformed_results[[4]], stat = transformed_results[[5]]), nsmall = 3, scientific = F, 
      trim = T), tidy.table.file, quote = FALSE, sep = "\t", row.names = F, 
      col.names = F, append = T)
  }

## Start executing analysis
logfile <- file(paste0(args[["--name"]], ".summary.log"), open = "at")
sink(file = logfile, split = T, type = "output")

config <- fromJSON(args[["--config"]], flatten = TRUE)
if (is.character(args[["--salmon-index"]])) {
  this.index <- args[["--salmon-index"]]
} else {
  this.index <- config[[args[["--sample-type"]]]][["default"]]
}
this.config <- config[[args[["--sample-type"]]]][["index"]][[this.index]]
this.index_dir <- config[[args[["--sample-type"]]]][["index_dir"]]

if (args[["--debug"]]) {
  cat("Index:")
  cat(str(this.index))
  cat("Config:")
  cat(str(this.config))
}

salmon.files <- file.path(args$DIR, "salmon", this.index, "quant.sf")
if (is.character(args[["--metadata"]])) {
  message("Reading metadata file...")
  metadata <- read.csv(args$metadata, header = TRUE, check.names = FALSE)
  if ("QuantFilePath" %in% colnames(metadata)) {
    message("Found 'QuantFilePath' column in metadata file. Using these paths to salmon quant files.")
    salmon.files <- file.path(metadata$QuantFilePath)
  }
  if ("Directory" %in% colnames(metadata)) {
    message("Found 'Directory' column in metadata file. Using these paths.")
    salmon.files <- file.path(metadata$Directory, "salmon", this.index, "quant.sf")
    if (args[["--debug"]]) {
      print(paste(c("salmon", this.index, "quant.sf"), sep = "/"))
      print(as.data.frame(salmon.files))
    }
  }
  message("Using sample names from SampleID column.")
  stopifnot(length(salmon.files) == length(metadata$SampleID))
  names(salmon.files) <- metadata$SampleID
} else {
  message("Using directories as sample names.")
  names(salmon.files) <- args$DIR
  metadata <- FALSE
}
message("Using samples:")
message(cbind(as.data.frame(salmon.files), file.exists(salmon.files)))

if (!all(file.exists(salmon.files))) {
    message("Some files are missing:")
    cat(unname(salmon.files[which(!file.exists(salmon.files))]), sep = "\n", file = stderr())
    if (args[["--allow-missing"]]) {
        message("Proceeding without missing samples.")
        missing.files <- salmon.files[which(!file.exists(salmon.files))]
        salmon.files <- salmon.files[which(file.exists(salmon.files))]
    }
}

if ("readr" %in% rownames(installed.packages()) == TRUE) {
  read_fn <- readr::read_tsv
} else {
  read_fn <- read.delim
}
if (is.character(args[["--only-transcripts"]])) {
  read_fn <- function(file, ...) {
    salmon.transcripts <- readr::read_tsv(file, ...)
    filtered.transcripts <- dplyr::filter(salmon.transcripts, grepl(args[["--only-transcripts"]], 
      Name))
    message(paste(length(filtered.transcripts$Name), "transcripts in", file))
    return(filtered.transcripts)
  }
}

txi.salmon.isoform <- tximport(salmon.files, type = "salmon", txOut = TRUE, importer = read_fn)


min_length <- as.numeric(args[["--min-transcript-length"]])
min_idx <- which(apply(txi.salmon.isoform$length, 1, max) > min_length)
message(paste("Removing", length(which(apply(txi.salmon.isoform$length, 1, max) <= 
  min_length)), "transcripts less than", min_length, "bp."))
txi.salmon.isoform$abundance <- txi.salmon.isoform$abundance[min_idx, ]
txi.salmon.isoform$counts <- txi.salmon.isoform$counts[min_idx, ]
txi.salmon.isoform$length <- txi.salmon.isoform$length[min_idx, ]

for (assembly in this.config[["fastas"]]) {
  message("Reading gene and transcript annotations.")
  # Split and normalize abundances based on the assembly containing genes
  transcripts2genes <- read.delim(paste(this.index_dir, args[["--sample-type"]], this.index, paste0(basename(assembly), 
    "_transcripts_to_genes.txt"), sep = "/"), stringsAsFactors = F, header = F)
  annotation <- read.delim(paste(this.index_dir, args[["--sample-type"]], this.index, paste0(basename(assembly), 
    "_gene_annotation.txt"), sep = "/"), stringsAsFactors = F, header = F)
  annotation[] <- lapply(annotation, as.character)
  transcripts2genes[] <- lapply(transcripts2genes, as.character)
  colnames(annotation) <- c("gene_id", "biotype", "symbol", "chrom", "start", "stop", 
    "length", "frac_masked")
  Summarize(txi.salmon.isoform, transcripts2genes, annotation, args, metadata, 
    species = tools::file_path_sans_ext(basename(assembly)))
}

message("Done.")
sink()
