library("devtools")
library("svglite")
library("openxlsx")
library("edgeR")
devtools::load_all(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Documents/CleanUpRNAseq)")
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\DNA conatmination\data\Li-RNA-seq)")

## metadata
full_metadata <-  read.delim("Li.metadata.txt", 
                             header = TRUE, as.is = TRUE)
## making a EnsDb using the GTF v110 for GRCh38
ensdb_sqplite <- make_ensdb (gtf = "docs/Homo_sapiens.GRCh38.110.gtf",
                             outfile = "docs/GRCh38.V110.ensdb.sqlite",
                             organism_latin_name = "Homo_sapiens",
                             genome_version = "GRCh38",
                             Ensembl_release_version = 110)
        
## generate SAF for genomic features:
saf_list <- get_feature_saf(ensdb_sqlite = ensdb_sqplite,
                            bamfile = "results/003.samtools.sort.out/HRR589587.srt.bam",
                            mitochondrial_genome = "MT")

## read summarization
counts_summary <- summarize_reads(metadata = "docs/Li_metadata.txt",
                                  isPairedEnd = TRUE,
                                  strandSpecific = 0,
                                  saf_list = saf_list,
                                  gtf = "docs/Homo_sapiens.GRCh38.110.gtf",
                                  threads = 8)
  
save(saf_list, counts_summary, file = "01112024.Li.RNASeq.featureCounts.RData")


## split by library method: polyA(+)-RNA selection (pA) and rRNA depletion (RZ)
pA_counts_summary <- list(gene = list(counts = counts_summary$gene$counts[, 1:18],
                                      stat = counts_summary$gene$stat[, c(1:19)]),
                          exon = list(counts = counts_summary$exon$counts[, 1:18],
                                      stat = counts_summary$exon$stat[, c(1:19)]),
                          intergenic_region = list(counts = counts_summary$intergenic_region$counts[, 1:18],
                                                   annotation = counts_summary$intergenic_region$annotation,
                                                   stat = counts_summary$intergenic_region$stat[, c(1:19)]),
                          intronic_region = list(counts = counts_summary$intronic_region$counts[, 1:18],
                                                 stat = counts_summary$intronic_region$stat[, c(1:19)]),
                          rRNA = list(counts = counts_summary$rRNA$counts[, 1:18],
                                      stat = counts_summary$rRNA$stat[, c(1:19)]),
                          mitochondrion = list(counts = counts_summary$mitochondrion$counts[, 1:18, drop = FALSE],
                                               stat = counts_summary$mitochondrion$stat[, c(1:19)]),
                          gtf = list(counts = counts_summary$gtf$counts[, 1:18],
                                     counts_junction = counts_summary$gtf$counts_junction[, 1:26],
                                     stat = counts_summary$gtf$stat[, c(1:19)]))


RZ_counts_summary <- list(gene = list(counts = counts_summary$gene$counts[, 19:36],
                                     stat = counts_summary$gene$stat[, c(1, 20:37)]),
                         exon = list(counts = counts_summary$exon$counts[, 19:36],
                                     stat = counts_summary$exon$stat[, c(1, 20:37)]),
                         intergenic_region = list(counts = counts_summary$intergenic_region$counts[, 19:36],
                                                  annotation = counts_summary$intergenic_region$annotation,
                                                  stat = counts_summary$intergenic_region$stat[, c(1,20:37)]),
                         intronic_region = list(counts = counts_summary$intronic_region$counts[, 19:36],
                                                stat = counts_summary$intronic_region$stat[, c(1,20:37)]),
                         rRNA = list(counts = counts_summary$rRNA$counts[, 19:36],
                                     stat = counts_summary$rRNA$stat[, c(1,20:37)]),
                         mitochondrion = list(counts = counts_summary$mitochondrion$counts[1, 19:36, drop = FALSE],
                                              stat = counts_summary$mitochondrion$stat[, c(1, 20:37)]),
                         gtf = list(counts = counts_summary$gtf$counts[, 19:36],
                                    counts_junction = counts_summary$gtf$counts_junction[, c(1:8, 27:44)],
                                    stat = counts_summary$gtf$stat[, c(1, 20:37)]))


## salmon summary
salmon <- salmon_res(metadata="docs/Li.metadata.txt",
                     ensdb_sqlite = ensdb_sqplite)
saveRDS(salmon, file = "12182023_Li_Salmon.quant.summary.RDS")

salmon <- readRDS("12182023_Li_Salmon.quant.summary.RDS")

## split by library method: polyA(+)-RNA selection (pA) and rRNA depletion (RZ)
pA_salmon <- lapply(salmon[1:3], function(.x){
  .x <- .x[, 1:18]
})

RZ_salmon <- lapply(salmon[1:3], function(.x){
  .x <- .x[, 19:36]
})


## modify polyA(+) selection metadata
modify_pA_metadata <- function(full_metadata)
{
  metadata <- full_metadata[1:18, ]
  simple_names <- paste(c(rep(LETTERS[1:6], each =3)),
                                     1:3, sep = "")
  metadata$sample_name <- factor(metadata$sample_name,
                                 levels = unique(metadata$sample_name),
                                 labels = simple_names)
  
  metadata$group <- factor(metadata$group,
                           levels = c("pA_noDNase",
                                      "pA_0pct_gDNA",
                                      "pA_0.01pct_gDNA",
                                      "pA_0.1pct_gDNA",
                                      "pA_1pct_gDNA",
                                      "pA_10pct_gDNA"),
                           labels = c("DNase I_free","0%",
                                      "0.01%", "0.1%",
                                      "1%", "10%"))
  list(metadata = metadata, simple_names = simple_names)
}

modify_RZ_metadata <- function(full_metadata)
{
  metadata <- full_metadata[19:36, ]
  simple_names <- paste(c(rep(LETTERS[1:6], each =3)),
                                     4:6, sep = "")
  metadata$sample_name <- factor(metadata$sample_name,
                                 levels = unique(metadata$sample_name),
                                 labels = simple_names)
  
  metadata$group <- factor(metadata$group,
                           levels = c("RZ_noDNase",
                                      "RZ_0pct_gDNA",
                                      "RZ_0.01pct_gDNA",
                                      "RZ_0.1pct_gDNA",
                                      "RZ_1pct_gDNA",
                                      "RZ_10pct_gDNA"),
                           labels = c("DNase I_free", "0%",
                                      "0.01%", "0.1%",
                                      "1%", "10%"))
  list(metadata = metadata, simple_names = simple_names)
}

rename_samples <- function(counts_summary, simple_names)
{
  colnames(counts_summary$gene$counts) <- simple_names
  colnames(counts_summary$gene$stat)[2:19] <- simple_names
  
  colnames(counts_summary$exon$counts) <- simple_names
  colnames(counts_summary$exon$stat)[2:19] <- simple_names
  
  colnames(counts_summary$intergenic_region$counts) <- simple_names
  colnames(counts_summary$intergenic_region$stat)[2:19] <- simple_names
  
  colnames(counts_summary$intronic_region$counts) <- simple_names
  colnames(counts_summary$intronic_region$stat)[2:19] <- simple_names
  
  colnames(counts_summary$rRNA$counts) <- simple_names
  colnames(counts_summary$rRNA$stat)[2:19] <- simple_names
  
  colnames(counts_summary$mitochondrion$counts) <- simple_names
  colnames(counts_summary$mitochondrion$stat)[2:19] <- simple_names
  
  colnames(counts_summary$gtf$counts) <- simple_names
  colnames(counts_summary$gtf$counts_junction)[9:26] <- simple_names
  colnames(counts_summary$gtf$stat)[2:19] <- simple_names
  counts_summary
}

## Detecting gDNA contamination
library_types <- c("polyA", "RiboZero")

for (library_type in library_types) 
{
  ## Assignment statistics using GTF as annotation
  if (library_type == "polyA")
  {
    metadata_simplenames <- modify_pA_metadata(full_metadata)
    metadaata <- metadata_simplenames$metadata
    counts_summary <- pA_counts_summary
    counts_summary <- rename_samples(counts_summary, 
                                        metadata_simplenames$simple_names)
    salmon <- pA_salmon
    salmon <- lapply(salmon[1:3], function(.x){
      .x <- .x[, 1:18]
      colnames(.x) <- simple_names
      .x
    })
  } else {
    metadata_simplename <- modify_RZ_metadata(full_metadata)
    metadata <- metadata_simplenames$metadata
    counts_summary <- RZ_counts_summary
    counts_summary <- rename_samples(counts_summary, 
                                        metadata_simplenames$simple_names)
    salmon <- RZ_salmon
    salmon <- lapply(salmon[1:3], function(.x){
      .x <- .x[, 1:18]
      colnames(.x) <- simple_names
      .x
    })
  }
  out <- library_type
  if(!dir.exists(out))
  {
    dir.create(out)
  }
  mapping_stat <- counts_summary$gtf$stat
  colnames(mapping_stat)[2:ncol(mapping_stat)] <- as.character(metadata$sample_name)
  p1 <- check_read_assignment_stat(assignment_stat = mapping_stat)
  svglite(filename = file.path(out, "Fig1.Read assignment status.svg"),
          height = 3, width = 7.)
  p1 & theme(axis.text.x = element_text(size = 5))
  dev.off()
  
  p2 <- check_read_distribution(featurecounts_list = counts_summary,
                               metadata = metadata)
  svglite(filename = file.path(out, "Fig1.Read distribution across genomic features.svg"),
          height = 6, width = 8)
  p2$p & theme(axis.text.x = element_text(size = 5))
  dev.off()
  
  p3_list <- check_expr_distribution(counts = salmon$counts,
                                     normalization = "qsmooth",
                                     metadata = metadata)
  svglite(filename = file.path(out, "Fig3.Salmon.expression level distribution.svg"),
          height = 8, width = 5)
  wrap_plots(p3_list, nrow = 3, ncol = 1) & theme(axis.text.x = element_text(size = 5))
  dev.off()
  
  p3_list <- check_expr_distribution(counts = counts_summary$gtf$counts,
                                     normalization = "qsmooth",
                                     metadata = metadata)
  svglite(filename = file.path(out, "Fig3.featureCounts.expression level distribution.svg"),
          height = 8, width = 5)
  
  wrap_plots(p3_list, nrow = 3, ncol = 1) & theme(axis.text.x = element_text(size = 5))
  dev.off()
  
  p7 <- check_expressed_gene_percentage (metadata = metadata,
                                         counts = salmon$counts,
                                         abundance = salmon$abundance,
                                         min_cpm = 1,
                                         min_tpm = 1)
  svglite(filename = file.path(out, "Fig4.Salmon.Percent of expressed genes.svg"),
          height = 3, width = 4.5)
  wrap_plots(p7) & theme(axis.text.x = element_text(size = 5))
  dev.off()
  
  svglite(filename = file.path(out, "Fig5.Salmon.sample.correlation.svg"),
          height = 8, width = 8)
  check_sample_correlation(counts = salmon$counts[, seq(2, ncol(salmon$counts), by =3)])
  dev.off()
  
  ## DESeq2 exploratory analysis before correction
  p9_list <- exploratory_analysis(counts = salmon$counts,
                                  metadata = metadata)
  svglite(filename = file.path(out, "Fig6.Salmon.PCA plot.svg"),
          height = 5, width = 6.5)
  options(ggrepel.max.overlaps = Inf)
  p9_list$pca
  dev.off()
  
  svglite(filename = file.path(out, "Fig6.Salmon.distance.heatmap.svg"),
          height = 4, width = 5)
  p9_list$heatmap
  dev.off()
}  


## Differential expression analysis before and after correction
library("DESeq2")
library("edgeR")
library("limma")
library("pheatmap")
library("RColorBrewer")
library("dplyr")
library("ggfortify")
library("genefilter")
library("gplots")
library("doParallel")
library("readxl")
library("openxlsx")
library("EnhancedVolcano")
library("msigdbr")
library("ReactomePA")
library("clusterProfiler")
library("magrittr")
library("biomaRt")
library("quantro")

plotPCA <- function(pca, metadata = data.frame(group = NULL, sample_name = NULL))
{
  pc12 <- as.data.frame(pca$x[, 1:2])
  colnames(pc12) <- c("PC1", "PC2")
  pc12 <- cbind(pc12, metadata)
  pc12_var <- round(pca$sdev[1:2]^2/(sum(pca$sdev^2)) * 100, digits = 2)
  pc12$group <- factor(pc12$group, levels = unique(pc12$group))
  pc12$sample_name <- factor(pc12$sample_name, levels = unique(pc12$sample_name))
  p_pca <- ggplot(pc12, aes(x = PC1, y = PC2,
                            color = group,
                            label = sample_name)) +
    geom_text_repel(size = 2.5, show.legend = FALSE)+
    geom_point() + xlab(paste0("PC1 (", pc12_var[1], "%)")) +
    guides(color=guide_legend(title="Group"))+
    ylab(paste0("PC2 (", pc12_var[2], "%)")) + ggtitle("PCA score plot") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.title=element_text(size = 10),
          legend.text =element_text(size = 8))
  p_pca
}


get.adj <- function(expr.data=NULL, design=NULL,  meta.data=NULL)
{
  ## adjusted expression using voom()
  v <- voom(expr.data, design=design)
  fit <- lmFit(v, design=design)
  group_num <- nlevels(as.factor(meta.data$group))
  if (ncol(design) > group_num)
  {
    col.covariate <- (group_num + 1):ncol(design)
    adj.exp <- v$E -fit$coefficients[, col.covariate] %*%
      t(design[, col.covariate])
  } else {
    adj.exp <- v$E
  }
  adj.exp
}

## Heatmap showing sample distances
distancePlot <- function(sampleDists, sampleNames)
{
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- sampleNames
  colnames(sampleDistMatrix) <- sampleNames
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           fontsize = 6,
           col = colors)
  sampleDistMatrix
}


thresholds <- 3
library_types <- c("polyA", "RiboZero")
correction_methods <- c("global", "gc", "IRrate", "none")
for (library_type in library_types){
  for (threshold in thresholds){
    for (correction_method in correction_methods)
    {

      experiment_name <- paste0("0425", library_type, ".",
                                correction_method,
                                ".cpm.cutoff.", threshold)

      if (!dir.exists(experiment_name))
      {
        dir.create(experiment_name, recursive = TRUE)
      }

      ## Assignment statistics using GTF as annotation
      if (library_type == "polyA")
      {
        metadata_simplenames <- modify_pA_metadata(full_metadata)
        metadaata <- metadata_simplenames$metadata
        counts_summary <- pA_counts_summary
        counts_summary <- rename_samples(counts_summary, 
                                         metadata_simplenames$simple_names)
        salmon <- pA_salmon
        salmon <- lapply(salmon[1:3], function(.x){
          .x <- .x[, 1:18]
          colnames(.x) <- simple_names
          .x
        })
      } else {
        metadata_simplename <- modify_RZ_metadata(full_metadata)
        metadata <- metadata_simplenames$metadata
        counts_summary <- RZ_counts_summary
        counts_summary <- rename_samples(counts_summary, 
                                         metadata_simplenames$simple_names)
        salmon <- RZ_salmon
        salmon <- lapply(salmon[1:3], function(.x){
          .x <- .x[, 1:18]
          colnames(.x) <- simple_names
          .x
        })
      }
    
      intergenic_rate <- check_read_distribution(featurecounts_list = counts_summary,
                                                metadata = metadata)$IR_rate
      metadata <- merge(metadata, intergenic_rate, 
                        by.x = "sample_name",
                        b.y = "row.names",
                        all.x = TRUE,
                        sort = FALSE)
      metadata$intergenic_rate <- scale(intergenic_rate, 
                                        center = TRUE, scale = FALSE)


      ### correction
      if (correction_method == "global")
      {
        ## global correction
        counts <- global_correction(intergenic_featureCounts_res =
                                      counts_summary$intergenic_region,
                                    salmon_res = salmon,
                                    lambda = 1)
      } else if (correction_method == "gc") {
        gene_gc <- readRDS("GRCh38.gene.exon.collapsed.GC.content.RDS")
        intergenic_gc <- readRDS("GRCh38.intergenic.GC.content.RDS")
        counts <- gc_bias_correction(salmon_res = salmon,
                                     gene_gc = gene_gc,
                                     intergenic_counts =
                                       counts_summary$intergenic_region$counts,
                                     intergenic_gc = intergenic_gc,
                                     plot = FALSE)
      } else {
        counts <- round(salmon$counts)
      }


       cpms = cpm(counts)
      # keep = rowSums(cpms > 0.5) >= threshold
       keep = rowSums(cpms > 1) >= threshold
      #keep <- filterByExpr(counts, group = metadata$group)

      counts = counts[keep, ]  ## 15689 genes
      #rm("cpms", "keep")
      counts <- counts[sapply(as.data.frame(counts), sd )!= 0, ]
      all(colnames(counts) == metadata$sample_name)

      if (correction_method == "global" ||
          correction_method == "gc" ||
          correction_method == "none")
      {
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                      colData = metadata,
                                      design = ~0 + group)
        model <- model.matrix(~ 0 + group, data = metadata)
      }  else {
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                      colData = metadata,
                                      design = ~0 + group + intergenic_rate)
        model <- model.matrix(~ 0 + group + intergenic_rate, data=metadata)
      }

      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds, fitType="parametric",
                                 maxit=1000)
      ## exploratory analysis
      vsd <- vst(dds, blind = TRUE)
      sampleDists <- dist(t(assay(vsd)))


      out <- experiment_name
      pdf(file.path(out,"Fig.2.Heatmap showing raw sample distances.pdf"),
          width = 3.5, height = 3)
      sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                       sampleNames = vsd$sample_name)
      dev.off()

      ## PCA plot
      pc <- prcomp(t(assay(vsd)), scale = T, center = T, retx = TRUE)
      pdf(file = file.path(out,"Fig 3. PCA.plot.pdf"),
          width = 4.2, height = 3)
      print(plotPCA(pc, metadata))
      dev.off()

      normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))

        adj.exp <- get.adj(expr.data=normCounts,
                           design=model, meta.data=metadata)

        write.table(adj.exp, file.path(out,
                                       paste0("Table 2.", correction_method,
                                              "adjusted.gene.expression.txt")),
                    sep = "\t", quote = FALSE, row.names = TRUE)

        ## Heatmap showing sample distances after adjusting for hidden variation
        sampleDists <- dist(t(adj.exp))
        pdf(file.path(out, paste0("Fig 2.2.Heatmap showing", correction_method,
                                  "adjusted sample distances.pdf")),
            width = 3.5, height = 3)
        sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                         sampleNames = paste0(vsd$sample_name))
        dev.off()

        ## PCA using adjusted value
        pca <- prcomp(t(adj.exp), scale = T, center = T, retx = TRUE)
        pdf(file = file.path(out,"Fig 3.2.PCA.plot.adjusted.samples.pdf"),
            width = 4.2, height = 3)
        print(plotPCA(pca, metadata))
        dev.off()


        pdf(file.path(out,"Figure 4.exploratory.data.analysis.after.adjustment.pdf"),
            width = 5, height = 4)
        op <- par(mar = c(7.5, 4, 2, 2))
        matboxplot(adj.exp,
                   groupFactor = metadata$group,
                   main = "after filtering",
                   cex = 0.3, ylab = "CPM")
        par(op)

        matdensity(adj.exp,
                   groupFactor = metadata$group,
                   ylab= "Density",
                   xlab = "CPM",
                   xlim= c(-10, 15),
                   main = "after filtering")

        dev.off()

      dds <- nbinomWaldTest(dds, modelMatrix = NULL,
                            betaPrior=FALSE,
                            maxit = 50000,
                            useOptim = TRUE,
                            quiet = FALSE,
                            useT = FALSE, useQR = TRUE)

      output_DESeq <- function(i, .name, dds, contrast.matrix,
                               threshold, shrink)
      {
        res <- results(dds, alpha = 0.01,
                       contrast=contrast.matrix[i,],
                       lfcThreshold=log2(threshold),
                       format = "DataFrame",
                       altHypothesis="greaterAbs")
        if (shrink)
        {
          res <- lfcShrink(dds,
                           contrast = contrast.matrix[i,],
                           res = res,
                           format = "DataFrame",
                           lfcThreshold = log2(threshold),
                           type = "ashr")
        }

        pdf(file.path(out,paste0("Figure 8.", .name, ".DESeq.MA.plot.pdf")),
            width = 7, height = 6)
        DESeq2::plotMA(res, alpha = 0.05,
                       main = rownames(contrast.matrix)[i])

        dev.off()

        geneid_name <- read.delim("GRCh38.V110.geneid_name_biotype.txt",
                                  header = FALSE, as.is = TRUE)
        colnames(geneid_name) <- c("GeneID", "Symbol", "Biotype")
        res <- merge(data.frame(res), geneid_name,
                     by.x = "row.names", by.y = "GeneID",
                     all.x = TRUE, sort = FALSE)
        rownames(res) <- res[,1]

        res <- res[, -1]
        res <- with(res, res[order(padj),])
      }

      ## raw FLC
        contrast.matrix <- matrix(c(-1, 1, 0, 0, 0, 0,
                                    -1, 0, 1, 0, 0, 0,
                                    -1, 0, 0, 1, 0, 0,
                                    -1, 0, 0, 0, 1, 0,
                                    -1, 0, 0, 0, 0, 1),
                                  nrow = 5, byrow = TRUE)
        rownames(contrast.matrix) <- c("noDNase",
                                       "0.01pct",
                                       "0.1pct",
                                       "1pct",
                                       "10pct")

      if (correction_method == "IRrate"){
        contrast.matrix <- cbind(contrast.matrix, 0)
      }

      DESeq_out <- mapply(output_DESeq, 1:nrow(contrast.matrix),
                          rownames(contrast.matrix),
                          MoreArgs = list(dds = dds,
                                          contrast.matrix = contrast.matrix,
                                          threshold = 1, shrink = FALSE), SIMPLIFY = FALSE)

      write.xlsx(DESeq_out, file = file.path(out,
                                             "Table 3.2.DESeq.differential expressed genes.xlsx"),
                 asTable = FALSE,
                 overwrite = TRUE,
                 rowNames = TRUE,
                 sheetName = rownames(contrast.matrix))

      pdf(file.path(out, paste("Figure 5.DESeq.volcanoplot.pdf")),
          width = 6, height = 8)
      for (i in 1:nrow(contrast.matrix))
      {

        p <- EnhancedVolcano(toptable = DESeq_out[[i]],
                             lab = DESeq_out[[i]]$Symbol,
                             x = 'log2FoldChange',
                             y = 'padj',
                             legendLabSize = 5,
                             legendPosition = "top",
                             title = rownames(contrast.matrix)[i],
                             ylab=  bquote(~-Log[10]~italic(adj.P)),
                             pCutoff = 0.05,
                             FCcutoff = 1,
                             legendLabels = c("NS",
                                              expression(Log[2]~FC),
                                              "adj.p-value",
                                              expression(adj.p-value~and~log[2]~FC)),
                             legendIconSize = 4,
                             drawConnectors = TRUE,
                             widthConnectors = 0.5,
                             pointSize = 0.5,
                             labSize = 2)
        print(p)
      }
      dev.off()
    }
  }
}




### density plots to compare correction effect.
library(reshape2)
library(ggplot2)
colors_20 <- c('#3cb44b', '#ffe119', '#4363d8', '#f58231',
               '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
               '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
               '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
               '#000000')

adjusted_expression <- dir(".", "adjusted.gene.expression.txt$", recursive = TRUE)
adjusted_expression <- adjusted_expression[grepl("0425", adjusted_expression)]
# split the adjusted expression: polyA(+) selection and rRNA depletion

for (library_type in c("polyA", "RZ"))
{
  if (library_type == "polyA")
  {
    expression <- adjusted_expression[1:4]
    metadata <- full_metadata[1:18, ]
    simple_names <- paste(c(rep(LETTERS[1:6], each =3)),
                                       1:3, sep = "")
    metadata$sample_name <- factor(metadata$sample_name,
                                   levels = unique(metadata$sample_name),
                                   labels = simple_names)
  } else {
    expression <- adjusted_expression[5:8]
    metadata <- full_metadata[19:36, ]
    simple_names <- paste(c(rep(LETTERS[1:6], each =3)),
                                       4:6, sep = "")
     metadata$sample_name <- factor(metadata$sample_name,
                                   levels = unique(metadata$sample_name),
                                   labels = simple_names)
  }
  names(expression) <- c("GC", "Global", "IR_rate", "None")
  
  expr_long <- mapply(function(.x, .name){
    dat <- read.delim(.x, header = TRUE, as.is = TRUE)
    dat <- melt(as.matrix(dat),
                value.name = "expression")
    colnames(dat)[1:2] <- c("GeneID", "Sample")
    dat$correction <- .name
    dat
  }, expression, names(expression), SIMPLIFY = FALSE)
  
  polyA_long <- do.call(rbind, expr_long)

  
  expr_long <- merge(expr_long, metadata[, c("group", "sample_name")],
                      by.x = "Sample", by.y = "sample_name", all.x = TRUE)
  
  if (library_type == "polyA")
  {
    expr_long$group <- factor(expr_long$group, levels = c("pA_0pct_gDNA",
                                                            "pA_noDNase",
                                                            "pA_0.01pct_gDNA",
                                                            "pA_0.1pct_gDNA",
                                                            "pA_1pct_gDNA",
                                                            "pA_10pct_gDNA"),
                               labels = c("0%", "DNase I_free",
                                          "0.01%", "0.1%",
                                          "1%", "10%"))
    expr_long$correction <- factor(expr_long$correction,
                                    levels = c("None", "IR_rate",
                                               "GC", "Global"),
                                    labels = c("None", "IR rate",
                                               "GC%", "Global"))
  } else {
    expr_long$group <- factor(expr_long$group, levels = c("RZ_0pct_gDNA",
                                                                  "RZ_noDNase",
                                                                  "RZ_0.01pct_gDNA",
                                                                  "RZ_0.1pct_gDNA",
                                                                  "RZ_1pct_gDNA",
                                                                  "RZ_10pct_gDNA"),
                                  labels = c("0%", "DNase I_free",
                                             "0.01%", "0.1%",
                                             "1%", "10%"))
    expr_long$correction <- factor(expr_long$correction,
                                       levels = c("None", "IR_rate",
                                                  "GC", "Global"),
                                       labels = c("None", "IR rate",
                                                  "GC%", "Global"))
  }

  
  p1 <- ggplot(expr_long, aes(x = expression, color = Sample)) +
    geom_density( linewidth = 0.5) + facet_wrap(~correction, ncol = 1)
  
  pdf(file = paste0("Fig 001.Density plots showing ", library_type,
                    " library gene expression.pdf"),
      height = 6, width = 3.5)
  p1 <- p1 + xlim(c(-10, 15)) + guides(color=guide_legend(title="Sample", ncol =1)) +
    theme_bw() + xlab(expression(log[2](CPM)))+ ylab("Density")+
    scale_color_manual(values= rep(colors_20[1:6], each =3)) +
    theme(legend.key.size = unit(0.3, 'cm'),
          title = element_text(size = 8, face = "bold"),
          strip.text = element_text(size = 8, face = "bold"))
  print(p1)
  dev.off()
}



### number of DEGs by correction methods
library(readxl)
degs_files <- dir(".", "Table 3.2.DESeq.differential expressed genes.xlsx$", 
                  recursive = TRUE)
degs_files <- degs_files[grepl("0425", degs_files)]

for (library_type in c("polyA", "ribozero"))
{
  if (library_type == "polyA")
  {
    deg_f <- degs_files[1:4]
  } else {
    deg_f <- degs_files[5:8]
  }
  names(deg_f) <- c("GC", "Global", "IR_rate", "None")
  num_deg <- mapply(function(.x, .z){
    dat <- lapply(excel_sheets(.x), function(.y) {
      degs <- read_excel(.x, sheet = .y)
      padj0.05_n <- sum(!is.na(degs$padj) & degs$padj <= 0.05)
      fc1.5_padj0.05_n <- sum(!is.na(degs$padj) & degs$padj <= 0.05 &
                                !is.na(degs$log2FoldChange) &
                                abs(degs$log2FoldChange) >= log2(1.5))
      fc2_padj0.05_n <- sum(!is.na(degs$padj) & degs$padj <= 0.05 &
                              !is.na(degs$log2FoldChange) &
                              abs(degs$log2FoldChange) >= log2(2))
      c(padj0.05_n = padj0.05_n,
        fc1.5_padj0.05_n = fc1.5_padj0.05_n,
        fc2_padj0.05_n = fc2_padj0.05_n)
    })
    dat <- as.data.frame(do.call("rbind", dat))
    rownames(dat) <- excel_sheets(.x)
    dat$correction <- .z
    dat
  }, deg_f, names(deg_f), SIMPLIFY = FALSE)
  
  
  num_deg <- do.call("rbind", num_deg)
  num_deg$contamination_level <- gsub("[^.]+\\.(.+)", "\\1", rownames(num_deg))
  num_deg$contamination_level[num_deg$contamination_level == "noDNAase"] <- "noDNase"
  
  num_deg <- melt(num_deg, id.vars = c("correction", "contamination_level"),
                  valu.name = "DEGs")
  num_deg$contamination_level <- factor(num_deg$contamination_level,
                                        levels =c("noDNase",
                                                  "0.01pct",
                                                  "0.1pct",
                                                  "1pct",
                                                  "10pct"),
                                        labels = c("DNase I_free",
                                                   "0.01%", "0.1%",
                                                   "1%", "10%"))
  
  num_deg$variable <- factor(num_deg$variable,
                             levels = unique(num_deg$variable),
                             labels = c("Padj \U2264 0.05",
                                        "Padj \U2264 0.05\n|log2(FC)| \U2265 0.585",
                                        "Padj \U2264 0.05\n|log2(FC)| \U2265 1"))
  num_deg$correction <- factor(num_deg$correction,
                               levels = c("None", "IR_rate",
                                          "GC", "Global"),
                               labels = c("None", "IR rate",
                                          "GC%", "Global"))
  
  p1 <- ggplot(data = num_deg,
               aes(x = contamination_level, y = value,
                   fill = contamination_level)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.5) +
    geom_text(aes(label = value, y = value*1.10),
              angle = 0, size = 2,
              hjust = .5) +
    guides(fill = guide_legend(title="Group", nrow =1)) +
    xlab("gDNA levels")+ ylab("#DEGs") +
    facet_grid(correction ~variable, scales = "free_y") +
    theme_bw() +  theme(legend.position = "bottom",
                        strip.text = element_text(size = 8, face = "bold"),
                        legend.key.size = unit(0.3, "cm"),
                        text = element_text(size = 8),
                        title = element_text(size = 8, face = "bold"),
                        axis.title.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank())
  pdf(file = paste0("Fig 003.", library_type,
                    "Number of DEGs by correction methods and cutoffs.pdf"),
      width = 6, height = 4)
  print(p1)
  dev.off()
}


##### volcano plots
library(EnhancedVolcano)
library(readxl)
library(ggrepel)
library(patchwork)
degs_files <- dir(".", "Table 3.2.DESeq.differential expressed genes.xlsx$",
                    recursive = TRUE)
degs_files <- degs_files[grepl("0425", degs_files)]
for (library_type in c("polyA", "ribozero"))
{
  if (library_type == "polyA")
  {
    deg_f <- degs_files[1:4][c(4,3, 1, 2)]
  } else {
    deg_f <- degs_files[5:8][c(4,3, 1, 2)]
  }
  
  names(deg_f) <- c("None", "IR%","GC%", "Global")
  gDNA_levels <- c("DNase I_free", "0.01%", "0.1%", "1%",
                   "10%")
  num_deg <- mapply(function(.x, .z){

    degs <- mapply(function(.y, .level) {

      degs <- read_excel(.x, sheet = .y)
      
      degs$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
      degs$diffexpressed[degs$log2FoldChange >= 1 &
                           !is.na(degs$padj) &
                           degs$padj <= 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      degs$diffexpressed[degs$log2FoldChange <= -1 &
                           !is.na(degs$padj) &
                           degs$padj <= 0.05] <- "DOWN"
      degs$diffexpressed <- factor(degs$diffexpressed, levels = c("DOWN",
                                                                  "NO",
                                                                  "UP"))
      
      degs$delabel <- NA
      degs$delabel[degs$diffexpressed != "NO"] <-
        degs$Symbol[degs$diffexpressed != "NO"]
      degs$levels <- .level
      degs
    }, excel_sheets(.x), gDNA_levels, SIMPLIFY = FALSE)
    
    degs <- do.call("rbind",  degs)
    degs$levels <- factor(degs$levels, levels =
                            c("DNase I_free", "0.01%", "0.1%", "1%",
                              "10%"))
    
    colours <- if (all(c("NO", "DOWN", "UP") %in% degs$diffexpressed)) {
      c("blue", "grey", "red")
    } else if (all(c("NO", "DOWN") %in% degs$diffexpressed)) {
      c("blue", "grey")
    } else if (all(c("NO", "UP") %in% degs$diffexpressed)) {
      c("grey", "red")
    } else{
      "grey"
    }
    
    p <- ggplot(data=degs, aes(x=log2FoldChange, y=-log10(padj),
                               col=diffexpressed, label=delabel)) +
      geom_point(size = 1) +
      theme_bw() +
      geom_text_repel(size = 1.5, force = 2, max.overlaps =20,
                      show.legend  = F) +
      scale_color_manual(values= colours) +
      geom_vline(xintercept=c(-log2(1.5), log2(1.5)),
                 col="black", linetype = "dashed", linewidth = 0.5,) +
      geom_hline(yintercept=-log10(0.05), linewidth = 0.5, 
                 col="black", linetype = "dashed") +
      xlab(expression(log[2](FoldChange))) +
      ylab(expression(-log[10](adjP))) + ggtitle(.z) +
      guides(color = guide_legend(title = "Diff. expressed")) +
      facet_wrap(~levels, ncol = 1, scales = "free_y") +
      theme(strip.text = element_text(size = 6, face = "bold"),
            plot.title = element_text(size = 8, hjust = 0.5, 
                                      face = "bold"),
            axis.text = element_text(size = 4),
            legend.text=element_text(size = 5),
            axis.title = element_text(size = 5, face = "bold"),
            legend.key.size = unit(0.4, "cm"))
    p
  }, deg_f, names(deg_f), SIMPLIFY = FALSE)
  
  pdf(paste0("Fig005.0426.",library_type, ".DEG.volcano.plots.pdf"),
      width = 12, height = 8)
  print(wrap_plots(num_deg,
                   ncol =4, 
                   guides = "collect",
                   tag_levels = "A") & theme(legend.position = 'bottom'))
  dev.off()
  
}
