library("devtools")
library(svglite)
pkg <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Documents\CleanUpRNAseq)"
devtools::load_all(pkg)

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\DNA conatmination\data\Michelle.RNA-seq\03082024)")

## making a EnsDb using the GTF v110 for GRCh38
ensdb_sqlite <- make_ensdb(gtf = "Homo_sapiens.GRCh38.110.gtf",
                       outfile = "GRCh38.V110.ensdb.sqlite",
                       organism_latin_name = "Homo_sapiens",
                       genome_version = "GRCh38",
                       Ensembl_release_version = 110)
#
# ## generate SAF for genomic features:
saf_list <- get_feature_saf(ensdb_sqlite = "GRCh38.V110.ensdb.sqlite",
      bamfile = "results/003.samtools.sort.out/K084CD7PCD1P_S17_L001.srt.bam",
                            mitochondrial_genome = c("MT", "chrM"),
                            chloroplast_genome = c("chrPltd", "Pltd"))

## Summarize read counts
counts_summary <- summarize_reads(
                            metadata = "Michelle.CD1A.RNAseq.txt",  
                            isPairedEnd = TRUE,
                            strandSpecific = 0,
                            saf_list = saf_list,
                            gtf = "docs/Homo_sapiens.GRCh38.110.gtf",
                            threads = 8)
salmon <- salmon_res(metadata = metadata,
                     ensdb_sqlite = ensdb_sqlite,)
metadata <- read.delim("CD1A.RNAseq.metadata.txt", header = TRUE, as.is = TRUE)
out <- "./full_analysis"

if (!dir.exists(out)){
  dir.create(out)
}

saveRDS(counts_summary, file = file.path(out,
                                         "Kevin_CD1A_14_full.counts.RDS"))
saveRDS(salmon, file.path(out, "Kevin_14sample.salmon.count.RDS"))
## using 14 samples data
counts_summary <- readRDS(file.path(out, "Kevin_CD1A_14_full.counts.RDS"))

mapping_stat <- counts_summary$gtf$stat

p1 <- check_read_assignment_stat(assignment_stat = mapping_stat)

svglite(filename = file.path(out, "Fig1. Read assignment status.svg"),
        height = 3, width = 4.5)
p1 & theme(axis.text.x = element_text(size = 6))
dev.off()


## read distribution across different genomic regions
p3 <- check_read_distribution(featurecounts_list = counts_summary,
                             metadata = metadata)

svglite(filename = file.path(out, "Fig2. Read distribution.svg"),
        height = 6, width = 7)
p3$p & theme(axis.text.x = element_text(size = 6))
dev.off()

## metadata reformatting
metadata$sample_name <- factor(metadata$sample_name,
                               levels = metadata$sample_name)
metadata$group <- gsub("CD1A\\(-\\)", "CD1A_N", metadata$group)
metadata$group <- gsub("CD1A\\(\\+\\)", "CD1A_P", metadata$group)
metadata$group <- factor(metadata$group, levels = unique(metadata$group))
metadata$batch <- factor(metadata$batch)

## before correction, QC
expression_qc <- function(counts = NULL,
                          abundance = NULL,
                             metadata = NULL,
                             count_type = c("featureCounts",
                                            "SalmonQuant",
                                            "correctedSalmon"))
{
  out2 <- file.path(out, count_type)
  if (!dir.exists(out2))
  {
    dir.create(out2)
  }

  count_type <- match.arg(count_type,
                          choices = c("featureCounts",
                                      "SalmonQuant",
                                      "correctedSalmon"),
                          several.ok = FALSE)

  p5_list <- check_expr_distribution(counts = counts,
                                     metadata = metadata)
  svglite(filename = file.path(out2, paste0("Fig3.", count_type,
                                           ".Expression.level.distribution.svg")),
          height = 8, width = 6)
  print(wrap_plots(p5_list, nrow = 3, ncol = 1) &
    theme(axis.text.x = element_text(size = 6)))
  dev.off()


  p7 <- check_expressed_gene_percentage (metadata = metadata,
                                         counts = counts,
                                         abundance = abundance,
                                         min_cpm = 1,
                                         min_tpm = 1)

  svglite(filename = file.path(out2, paste0("Fig4.", count_type,
                                           "Percent of expressed genes.svg")),
          height = 3, width = 4)
  print(wrap_plots(p7) & theme(axis.text.x = element_text(size = 6)))
  dev.off()


  #svglite(filename = file.path(out2, paste0("Fig5.", count_type,
  #                                         ".sample correlation.svg")),
  tiff(filename = file.path(out2, paste0("Fig5.", count_type,
                                           ".sample correlation.tiff")),
  units="in", width=10.5, height=10.5, res=600)
  check_sample_correlation(counts = counts)
  dev.off()


  ## DESeq2 exploratory analysis before correction
  p9_list <- exploratory_analysis(counts = counts,
                                  metadata = metadata)
  svglite(filename = file.path(out2, paste0("Fig6.", count_type,
                                           ".PCA plot.svg")),
          height = 5, width = 6.5)
  options(ggrepel.max.overlaps = Inf)
  print(p9_list$pca)
  dev.off()

  svglite(filename = file.path(out2, paste0("Fig6.", count_type,
                                           ".distance heatmap.svg")),
          height = 5, width = 5.2)
  print(p9_list$heatmap)
  dev.off()
}

## featureCounts-based QC
stopifnot(all(colnames(counts_summary$gtf$counts) == metadata$sample_name))
expression_qc(counts = counts_summary$gtf$counts,
              metadata = metadata,
              count_type = "featureCounts")

## Salmon quant-based QC
salmon <- readRDS(file.path(out, "Kevin_14sample.salmon.count.RDS"))
stopifnot(all(colnames(salmon$counts) == metadata$sample_name))
expression_qc(counts = salmon$counts,
              abundance = salmon$abundance,
              metadata = metadata,
              count_type = "SalmonQuant")

## global correction
corrected_counts <- global_correction(intergenic_featureCounts_res =
                                        counts_summary$intergenic_region,
                                      salmon_res = salmon,
                                      lambda = 1.0)

expression_qc(counts = corrected_counts,
              abundance = NULL,
              metadata = metadata,
              count_type = "correctedSalmon")


#### GC-aware correction
gene_gc <- readRDS("../GRCh38.gene.exon.collapsed.GC.content.RDS")
colnames(gene_gc)[1] <- "gc_content"
saveRDS(gene_gc, "../GRCh38.gene.exon.collapsed.GC.content.RDS")
IR_gc <- readRDS("../GRCh38.intergenic.GC.content.RDS")

GCcorrected_counts <-
    gc_bias_correction(salmon_res = salmon,
                       gene_gc = gene_gc,
                       intergenic_counts =
                         counts_summary$intergenic_region$counts,
                       intergenic_gc = IR_gc,
                       plot = FALSE)


### check differential genes
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

counts_summary <- readRDS(file.path(out, "Kevin_CD1A_14_full.counts.RDS"))
metadata$sample_name <- factor(metadata$sample_name,
                               levels = metadata$sample_name)
metadata$group <- gsub("CD1A\\(-\\)", "CD1A_N", metadata$group)
metadata$group <- gsub("CD1A\\(\\+\\)", "CD1A_P", metadata$group)
metadata$group <- factor(metadata$group, levels = unique(metadata$group))
metadata$batch <- factor(metadata$batch)

stopifnot(all(colnames(counts_summary$gtf$counts) == metadata$sample_name))

salmon <- readRDS(file.path(out, "Kevin_14sample.salmon.count.RDS"))
stopifnot(all(colnames(salmon$counts) == metadata$sample_name))
corrected_counts <- global_correction(intergenic_featureCounts_res =
                                        counts_summary$intergenic_region,
                                      salmon_res = salmon,
                                      lambda = 1.0)
### GC-aware correction
gene_gc <- readRDS("../GRCh38.gene.exon.collapsed.GC.content.RDS")
colnames(gene_gc)[1] <- "gc_content"
saveRDS(gene_gc, "../GRCh38.gene.exon.collapsed.GC.content.RDS")
IR_gc <- readRDS("../GRCh38.intergenic.GC.content.RDS")

GCcorrected_counts <-
  gc_bias_correction(salmon_res = salmon,
                     gene_gc = gene_gc,
                     intergenic_counts =
                       counts_summary$intergenic_region$counts,
                     intergenic_gc = IR_gc,
                     plot = FALSE)

plotPCA <- function(pca, metadata = data.frame(group = NULL, sample_name = NULL))
{
  pc12 <- as.data.frame(pca$x[, 1:2])
  colnames(pc12) <- c("PC1", "PC2")
  pc12 <- cbind(pc12, metadata)
  pc12_var <- round(pca$sdev[1:2]^2/(sum(pca$sdev^2)) * 100, digits = 2)
  pc12$group <- factor(pc12$group, levels = unique(pc12$group))
  pc12$sample_name <- factor(pc12$sample_name,
                             levels = unique(pc12$sample_name))
  p_pca <- ggplot(pc12, aes(x = PC1, y = PC2,
                            color = group,
                            label = sample_name)) +
    geom_text_repel(size = 2.5)+
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

IR_rate <- check_read_distribution(featurecounts_list =
                                          counts_summary,
                                  metadata = metadata)$IR_rate

metadata$intergenic_rate <- merge(metadata, IR_rate,
                                  by.x = "sample_name", by.y = "row.names",
                                  all.x = TRUE,
                                  sort = FALSE)

## scale for GLM stable
if (sd(metadata$intergenic_rate) > 5) {
  metadata$intergenic_rate <- 
    scale(intergenic_rate, center = TRUE, scale = FALSE)
}

thresholds  <- 2
correction_methods <- c("global", "gc", "IRrate", "none")

  for (threshold in thresholds){
    for (correction_method in correction_methods)
    {
      experiment_name <- paste0(use_sva, ".",
                                correction_method,
                                ".cpm.", threshold)

      if (!dir.exists(experiment_name))
      {
        dir.create(experiment_name, recursive = TRUE)
      }

      ### correction
      if (correction_method == "global")
      {
        ## global correction
        counts <- corrected_counts
      } else if (correction_method == "gc") {
        counts <- GCcorrected_counts
      } else {
        counts <- round(salmon$counts)
      }

      cpms = cpm(counts)
      # keep = rowSums(cpms > 1) >= threshold
      # keep = rowSums(cpms > 0.5) >= threshold
      keep = rowSums(cpms > median(10/(colSums(counts)/10^6))) >= threshold

      counts = counts[keep, ]
      rm("cpms", "keep")
      counts <- counts[sapply(as.data.frame(counts), sd )!= 0, ]
      all(colnames(counts) == metadata$sample_name)
      if (correction_method == "global" ||
          correction_method == "gc" ||
          correction_method == "none")
      {
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                      colData = metadata,
                                      design = ~0 + group + batch)
      }  else {
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                      colData = metadata,
                                      design = ~0 + group + batch + intergenic_rate)
      }
      
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds, fitType="parametric",
                                 maxit=1000)
      ## exploratory analysis
      vsd <- vst(dds, blind = TRUE)
      sampleDists <- dist(t(assay(vsd)))
      
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
      
      out <- experiment_name
      
      pdf(file.path(out,"Fig.2.Heatmap showing raw sample distances.pdf"),
          width = 6.5, height = 5)
      sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                       sampleNames = vsd$sample_name)
      dev.off()
      
      ## PCA plot
      pc <- prcomp(t(assay(vsd)), scale = T, center = T, retx = TRUE)
      pdf(file = file.path(out,"Fig 3. PCA.plot.pdf"), width = 6.5, height = 5)
      print(plotPCA(pc, metadata))
      dev.off()
      
      normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
      
      ## adjust log(cpm)
      get.adj <- function(expr.data=NULL, 
                          design=NULL,
                          meta.data=NULL)
      {
        ## adjusted expression using voom()
        v <- voom(expr.data, design=design)
        fit <- lmFit(v, design=design)
        group_num <- nlevels(as.factor(meta.data$group))
        if (ncol(design) > group_num)
        {
          col.covariate <- (group_num + 1):ncol(design)
          adj.exp <- v$E - fit$coefficients[, col.covariate] %*%
            t(design[, col.covariate])
        } else {
          adj.exp <- v$E
        }
        adj.exp
      }
      
      if (correction_method == "global" ||
          correction_method == "gc" ||
          correction_method == "none")
      {
        model <- model.matrix(~ 0 + group + batch, data=metadata)
        adj.exp <- get.adj(expr.data=normCounts,
                           design=model, meta.data = metadata)
      } else {
        model <- model.matrix(~ 0 + group + batch + intergenic_rate,
                              data=metadata)
        adj.exp <- get.adj(expr.data=normCounts,
                           design=model, meta.data = metadata)
      }
      
      
      ## adjust expression for hidden variations for EDA plots
      write.table(adj.exp, file.path(out, "Table 2.batch-adjusted.gene.expression.txt"),
                  sep = "\t", quote = FALSE, row.names = TRUE)
      
      ## Heatmap showing sample distances after adjusting for hidden variation
      sampleDists <- dist(t(adj.exp))
      pdf(file.path(out, "Fig 2.2.Heatmap showing batch-adjusted sample distances.pdf"),
          width = 6, height = 5)
      sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                       sampleNames = paste0(vsd$sample_name))
      dev.off()
      
      ## PCA using adjusted value
      pca <- prcomp(t(adj.exp), scale = T, center = T, retx = TRUE)
      pdf(file = file.path(out,"Fig 3.2.PCA.plot.batch.sva.adjusted.samples.pdf"),
          width = 6.5, height = 5)
      print(plotPCA(pca, metadata))
      dev.off()
      
      pdf(file.path(out,"Figure 4.exploratory.data.analysis.before.qsmooth.pdf"),
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
                 main = "after filtering")
      #legend('topright', legend = levels(metadata$group),
      #       col = 1:length(unique(metadata$group)),lty = 1)
      
      dev.off()
      
      if(correction_method == "global" ||
         correction_method == "gc" ||
         correction_method == "none")
      {
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                      colData = metadata,
                                      design = ~ 0 + group + batch)
      } else {  ## IR_rate
        dds <- 
          DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                 colData = metadata,
                                 design = ~ 0 + group + batch + intergenic_rate)
      }
      
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds, fitType="parametric",
                                 maxit=1000)
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
      
      contrast.matrix <- matrix(c(-1, 1, 0, 0, 0),
                                nrow = 1, byrow = TRUE)
      rownames(contrast.matrix) <- "CD1A_N-CD1AP"
      
      if (correction_method == "IRrate"){
        contrast.matrix <- cbind(contrast.matrix, 0)
      }
      
      DESeq_out <- mapply(output_DESeq, 1:nrow(contrast.matrix),
                          rownames(contrast.matrix),
                          MoreArgs = list(dds = dds,
                                          contrast.matrix = contrast.matrix,
                                          threshold = 1, shrink = FALSE), SIMPLIFY = FALSE)
      
      write.xlsx(DESeq_out, file =
                   file.path(out, "Table 3.2.DESeq.differential expressed genes.xlsx"),
                 asTable = FALSE,
                 overwrite = TRUE,
                 rowNames = TRUE,
                 sheetName = rownames(contrast.matrix))
      
      #### volcano plots
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
                             legendIconSize = 5,
                             arrowheads = FALSE,
                             drawConnectors = TRUE,
                             widthConnectors = 0.5,
                             pointSize = 0.8,
                             #xlim = ifelse(i ==1, c(-5, 18), c(-6, 18)),
                             labSize = 2)
        print(p)
      }
      dev.off()
    }
  }



### density plots to compare correction effect.
library(reshape2)
library(ggplot2)
adjusted_expression <- dir(".", "adjusted.gene.expression.txt$", 
                           recursive = TRUE)

adjusted_expression <- adjusted_expression[grepl("^noSVA", 
                                                 adjusted_expression)]
names(adjusted_expression) <- c("Batch(+);GC%",
                                "Batch(+);Global",
                                "Batch(+);IR_rate",
                                "Batch(+);None",
                                "Batch(-);None")

adjusted_expression_long <- mapply(function(.x, .name){
  
  #.x = adjusted_expression[1]
  dat <- read.delim(.x, header = TRUE, as.is = TRUE, check.names = FALSE)
  dat <- melt(as.matrix(dat),
              value.name = "expression")
  colnames(dat)[1:2] <- c("GeneID", "Sample")
  dat$correction <- .name
  dat
}, adjusted_expression, names(adjusted_expression), SIMPLIFY = FALSE)

adjusted_expression_long <- do.call(rbind, adjusted_expression_long)

adjusted_expression_long <- merge(adjusted_expression_long,
                                  metadata[, c("group", "sample_name")],
                                  by.x = "Sample", 
                                  by.y = "sample_name", 
                                  all.x = TRUE)
adjusted_expression_long$group <- factor(adjusted_expression_long$group,
                                         levels = c("CD1A_N",
                                                    "CD1A_P"),
                                         labels = c("CD1A(-)", "CD1A(+)"))
adjusted_expression_long$correction <- factor(adjusted_expression_long$correction,
                                              levels = c("Batch(-);None",
                                                         "Batch(+);None",
                                                         "Batch(+);IR_rate",
                                                         "Batch(+);GC%",
                                                         "Batch(+);Global"),
                                              
                                              labels = c("Batch(-);None",
                                                         "Batch(+);None",
                                                         "Batch(+);IR%",
                                                         "Batch(+);GC%",
                                                         "Batch(+);Global"))


p1 <- ggplot(adjusted_expression_long, aes(x = expression, color = Sample)) +
  geom_density( linewidth = 1) + facet_wrap(~correction, ncol = 1)

pdf(file = paste0("Fig 001.Density plots showing SMART-seq4", 
                  " library gene expression-",
                  sva, ".pdf"),
    height = 6, width = 5.5)
p1 <- p1 + xlim(c(-10, 15)) + guides(color=guide_legend(title="Sample")) +
  theme_bw() + xlab(expression(log[2](CPM)))+ ylab("Density")+
  theme(legend.key.size = unit(0.3, 'cm'),
        strip.text = element_text(size = 8, face = "bold"),
        text = element_text(size = 8),
        title = element_text(size = 8, face = "bold"),)
print(p1)
dev.off()


### number of DEGs by correction methods
degs_files <- dir(".", "Table 3.2.DESeq.differential expressed genes.xlsx$", 
                  recursive = TRUE)

degs_files <- degs_files[grepl("^noSVA|clean", degs_files)]
names(degs_files) <- c("Batch(+);10samples_1",
                       "Batch(+);IR%;12samples",
                       "Batch(+);IR%;10samples_2",
                           "Batch(+);GC%",
                           "Batch(+);Global",
                           "Batch(+);IR_rate",
                           "Batch(+);None")


  library(readxl)
  num_deg <- mapply(function(.x, .z){
    dat <- lapply(excel_sheets(.x), function(.y) {
      degs <- read_excel(.x, sheet = .y)
      padj0.05_up <- sum(!is.na(degs$padj) & degs$padj <= 0.05 &
                           degs$log2FoldChange > 0)
      padj0.05_dn <- sum(!is.na(degs$padj) & degs$padj <= 0.05 & 
                           degs$log2FoldChange < 0)
      fc1.5_padj0.05_up <- sum(!is.na(degs$padj) & degs$padj <= 0.05 &
                                !is.na(degs$log2FoldChange) &
                                degs$log2FoldChange >= log2(1.5))

      fc1.5_padj0.05_dn <- sum(!is.na(degs$padj) & degs$padj <= 0.05 &
                                 !is.na(degs$log2FoldChange) &
                                 degs$log2FoldChange <= -log2(1.5))

      fc2_padj0.05_up <- sum(!is.na(degs$padj) & degs$padj <= 0.05 &
                              !is.na(degs$log2FoldChange) &
                              degs$log2FoldChange >= log2(2))
      fc2_padj0.05_dn <- sum(!is.na(degs$padj) & degs$padj <= 0.05 &
                               !is.na(degs$log2FoldChange) &
                               degs$log2FoldChange <= -log2(2))

      c(padj0.05_up = padj0.05_up,
        padj0.05_dn = padj0.05_dn,
        fc1.5_padj0.05_up = fc1.5_padj0.05_up,
        fc1.5_padj0.05_dn = fc1.5_padj0.05_dn,
        fc2_padj0.05_up = fc2_padj0.05_up,
        fc2_padj0.05_dn = fc2_padj0.05_dn)
    })
    dat <- as.data.frame(do.call("rbind", dat))
    rownames(dat) <- excel_sheets(.x)
    dat$correction <- .z
    dat
  }, degs_files, names(degs_files), SIMPLIFY = FALSE)


  num_deg <- do.call("rbind", num_deg)
  num_deg <- melt(num_deg, id.vars = c("correction"),
                  value.name = "DEGs")

  num_deg$DEGs <- ifelse(grepl("_dn", num_deg$variable), -num_deg$DEGs, num_deg$DEGs)
  num_deg$cutoff <- gsub("_(up|dn)$", "", num_deg$variable)
  num_deg$cutoff  <- factor(num_deg$cutoff ,
                             levels = unique(num_deg$cutoff),
                             labels = c("Padj \U2264 0.05",
                                        "Padj \U2264 0.05\n|log2(FC)| \U2265 0.585",
                                        "Padj \U2264 0.05\n|log2(FC)| \U2265 1"))

num_deg$correction <- factor(num_deg$correction,
                                 levels = c("Batch(+);IR%;12samples",
                                            "Batch(+);IR%;10samples_2",
                                            "Batch(+);10samples_1",
                                            "Batch(+);IR_rate",
                                            "Batch(+);GC%",
                                            "Batch(+);Global",
                                            "Batch(+);None"),

                                 labels = c("Batch(+);IR%;12samples",
                                            "Batch(+);IR%;10samples_2",
                                            "Batch(+);10samples_1",
                                            "Batch(+);IR%",
                                            "Batch(+);GC%",
                                            "Batch(+);Global",
                                            "Batch(+);None"))

  #num_deg2$value[num_deg2$value == 0] <- 1
  p1 <- ggplot(data = num_deg,
               aes(x = cutoff, y = DEGs,
                   fill = correction)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.6) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_text(aes(label = ifelse(DEGs > 0, DEGs, -DEGs),
                  y = ifelse(DEGs > 0, DEGs + 10, DEGs -200)),
              position = position_dodge(0.6),
              angle = 90, size = 2.5,
              hjust = 0) +
    guides(fill = guide_legend(title="Correction")) +
    xlab("Threshold")+ ylab("#DEGs") +
    expand_limits(y = 0) + scale_y_continuous(breaks =c(-1000, -500, 0, 500, 1000),
                                            labels =c(1000, 500, 0, 500, 1000),
                                            limits = c(-1200, 1500)) +
    theme_bw() +  theme(legend.key.size = unit(0.25, "cm"),
                        text = element_text(size = 8),
                        title = element_text(size = 8,
                                             face = "bold"))
  svg(filename = paste0("Fig 003.Number of DEGs by correction methods and cutoffs",
                    ".svg"),
      width = 5, height = 3)
  print(p1)
  dev.off()



##### volcano plots
library(EnhancedVolcano)
library(readxl)
library(ggrepel)
library(patchwork)
  degs_files <- dir(".", "Table 3.2.DESeq.differential expressed genes.xlsx$",
                    recursive = TRUE)
  
  degs_files <- degs_files[grepl("^noSVA|clean", degs_files)]
  names(degs_files) <- c("Batch(+);10samples",
                         "Batch(+);GC%",
                         "Batch(+);Global",
                         "Batch(+);IR_rate",
                         "Batch(+);None")
  
  degs <- mapply(function(.y, .name) {
    #.y = excel_sheets(.x)[1]
    degs <- read_excel(.y, sheet = 1)
    
    degs$diffexpressed <- "NO"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
    degs$diffexpressed[degs$log2FoldChange >= 1 &
                         !is.na(degs$padj) &
                         degs$padj <= 0.05] <- "UP"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    degs$diffexpressed[degs$log2FoldChange <= -1 &
                         !is.na(degs$padj) &
                         degs$padj <= 0.05] <- "DOWN"
    degs$diffexpressed <- factor(degs$diffexpressed,
                                 levels = c("DOWN",
                                            "NO",
                                            "UP"))
    
    
    degs$delabel <- NA
    degs$delabel[degs$diffexpressed != "NO"] <-
      degs$Symbol[degs$diffexpressed != "NO"]
    degs$correction <- .name
    degs
  }, degs_files, names(degs_files), SIMPLIFY = FALSE)
  
  degs <- do.call("rbind",  degs)
  
  
  degs$correction <- factor(degs$correction,
                            
                            levels = c("Batch(+);10samples",
                                       "Batch(+);IR_rate",
                                       "Batch(+);None",
                                       "Batch(+);GC%",
                                       "Batch(+);Global"),
                            
                            labels = c("Batch(+);10samples",
                                       "Batch(+);IR%",
                                       "Batch(+);None",
                                       "Batch(+);GC%",
                                       "Batch(+);Global"))
  
  
  colours <- if (all(c("NO", "DOWN", "UP") %in% degs$diffexpressed)) {
    c("blue", "grey", "red")
  } else if (all(c("NO", "DOWN") %in% degs$diffexpressed)) {
    c("blue", "grey")
  } else if (all(c("NO", "UP") %in% degs$diffexpressed)) {
    c("grey", "red")
  } else{
    "grey"
  }
  
  p <- ggplot(data=degs, aes(x = log2FoldChange, y = -log10(padj),
                             color = diffexpressed, label = delabel)) +
    geom_point(size = 1) +
    theme_bw() +
    geom_text_repel(size = 2, force = 2, max.overlaps =20,
                    show.legend  = F) +
    scale_color_manual(values= colours) +
    geom_vline(xintercept=c(-log2(1.5), log2(1.5)),
               col="black", linetype = "dashed") +
    geom_hline(yintercept=-log10(0.05),
               col="black", linetype = "dashed") +
    xlab(expression(log[2](FoldChange))) +
    ylab(expression(-log[10](adjP))) +
    guides(color = guide_legend(title = "Diff. expressed")) +
    facet_wrap(~correction, nrow = 3, scales = "free") +
    theme(strip.text = element_text(size = 8, face = "bold"),
          legend.position="top",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 8),
          legend.text=element_text(size = 8),
          axis.title = element_text(size = 10, face = "bold"),
          legend.key.size = unit(0.4, "cm"))
  pdf(paste0("Fig005.Volcano plots showing DEGs-", sva,".pdf"),
      width = 8, height = if (sva=="SVA") {11} else {10})
  print(p)
  dev.off()
  


### Sankey plot showing DEGs
if (!requireNamespace("ggalluvial")) {
  BiocManager::install("ggalluvial")
}
library("ggalluvial")

## DEGs from different methods
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\DNA conatmination\data\Michelle.RNA-seq\03082024)")

degs_files <- dir(".", "Table 3.2.DESeq.differential expressed genes.xlsx",
                   recursive = TRUE, full.names = TRUE)


degs_files <- degs_files[grepl("noSVA|clean", degs_files)]
names(degs_files) <- 
  c("Batch(+);10samples_1",
    "Batch(+);IR%;12samples",
    "Batch(+);IR%;10samples_2",
    "Batch(+);GC%",
    "Batch(+);Global",
    "Batch(+);IR_rate",
    "Batch(+);None")

library(readxl)

dat <- mapply(function(.x, .name) {
    degs <- read_excel(.x, sheet = 1)
    degs$status <- ifelse(!is.na(degs$padj) & degs$padj <= 0.05 &
                            !is.na(degs$log2FoldChange) &
                            degs$log2FoldChange >= log2(2),
                          "UP", ifelse(!is.na(degs$padj) & degs$padj <= 0.05 &
                                         !is.na(degs$log2FoldChange) &
                                         degs$log2FoldChange <= -log2(2), 
                                       "DOWN", "NO"))
    degs$correction <- .name
    degs <- degs[, c(1, 10,11)]
    colnames(degs)[1] <- "GeneID"
    degs
}, degs_files, names(degs_files), SIMPLIFY = FALSE)
all_genes <- unique(unlist(sapply(dat, "[[", 1)))
dat <- lapply(dat, function(.x){
       tmp_df <- data.frame(GeneID = all_genes[!all_genes %in% .x$GeneID],
                            status = "NO",
                            correction = .x$correction[1])
       .x <- rbind(.x, tmp_df)
       .x
})

all_nonDE <- intersect(dat[[1]]$GeneID[dat[[1]]$status == "NO"],
             intersect(dat[[2]]$GeneID[dat[[2]]$status == "NO"],
             intersect(dat[[3]]$GeneID[dat[[3]]$status == "NO"],
             intersect(dat[[4]]$GeneID[dat[[4]]$status == "NO"],
             intersect(dat[[5]]$GeneID[dat[[5]]$status == "NO"],
             intersect(dat[[6]]$GeneID[dat[[6]]$status == "NO"],
                       dat[[7]]$GeneID[dat[[7]]$status == "NO"]))))))

## 500 of 28664/30014 gene non-DE for all nonSVA analysis
## 500 of 28966/30014 gene non-DE for all SVA analysis
exluded_nonDE <- all_nonDE[500:length(all_nonDE)]

dat_cmb <- do.call('rbind', dat)
dat_cmb$status <- factor(dat_cmb$status, levels = c("UP",
                                                    "NO",
                                                    "DOWN"))

## 10 samples
dat_cmb0 <- dat_cmb
samples  <- c(10, 12, 14)[3]
for (s in samples)
{
  if (s == 10)
  {
    dat_cmb <- dat_cmb0[dat_cmb0$correction != "Batch(+);IR%;12samples", ]
    dat_cmb$correction <- factor(dat_cmb$correction, levels = c(

      "Batch(+);IR%;10samples_2",
      "Batch(+);10samples_1",
      "Batch(+);IR_rate",
      "Batch(+);GC%",
      "Batch(+);Global",
      "Batch(+);None"),
      labels = c(
        "Batch(+);IR%;10samples_2",
        "Batch(+);10samples_1",
        "Batch(+);IR%",
        "Batch(+);GC%",
        "Batch(+);Global",
        "Batch(+);None"))
  } else if (s == 12) {
    dat_cmb <-  dat_cmb0[dat_cmb0$correction != "Batch(+);IR%;10samples_2", ]
    dat_cmb$correction <- factor(dat_cmb$correction, levels = c(
      "Batch(+);IR%;12samples",
      "Batch(+);10samples_1",
      "Batch(+);IR_rate",
      "Batch(+);GC%",
      "Batch(+);Global",
      "Batch(+);None"),
      labels = c(
        "Batch(+);IR%;12samples",
        "Batch(+);10samples_1",
        "Batch(+);IR%",
        "Batch(+);GC%",
        "Batch(+);Global",
        "Batch(+);None"))
  } else {
    dat_cmb <- dat_cmb0
    dat_cmb$correction <- factor(dat_cmb$correction,
                                 levels = c(
      "Batch(+);IR%;12samples",
      "Batch(+);IR%;10samples_2",
      "Batch(+);10samples_1",
      "Batch(+);IR_rate",
      "Batch(+);GC%",
      "Batch(+);Global",
      "Batch(+);None"),
      labels = c(
        "Batch(+);IR%;12samples",
        "Batch(+);IR%;10samples_2",
        "Batch(+);10samples_1",
        "Batch(+);IR%",
        "Batch(+);GC%",
        "Batch(+);Global",
        "Batch(+);None"))
  }
  dat_cmb$freq <- 1
  dat_cmb <- dat_cmb[!dat_cmb$GeneID%in%exluded_nonDE, ]

  p <- ggplot(dat_cmb,
               aes(x = correction,
                   stratum = status,
                   alluvium = GeneID,
                   y = freq,
                   fill = status,
                   label = status)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_discrete(type = c("#f58231", "#aaffc3", "#42d4f4")) +
    geom_flow(alpha = 0.6) +
    geom_stratum(alpha = 1) + xlab("Correction Method") +
    ylab("Count") +
    geom_text(stat = "stratum", size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
          plot.margin = margin(0.1, 0.8, 0.1, 0.1, unit ="cm"),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10, face = "bold"))

  #pdf("Fig 010.noSVA.DEGs.Sankey.plot.pdf", height = 3, width = 4)
  pdf(paste0("Fig 011.DEGs.Sankey.plot.", s, "samples.pdf"),
      height = 3, width = 7.5)
  print(p)
  dev.off()
}

## gene rank files
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\DNA conatmination\data\Michelle.RNA-seq\03082024)")

library(readxl)
out <- "Rnk.files"
if (!dir.exists(out))
{
  dir.create(out)
}
degs_files <- dir(".", "Table 3.2.DESeq.differential expressed genes.xlsx",
                  recursive = TRUE, full.names = TRUE)

degs_files <- degs_files[grepl("noSVA", degs_files)]
names(degs_files) <- c("Batch(+);GC%",
                       "Batch(+);Global",
                       "Batch(+);IR_rate",
                       "Batch(+);None")
dat <- mapply(function(.x, .name) {
  degs <- read_excel(.x, sheet = 1)
  rnk <- with(degs, degs[order(degs$log2FoldChange, decreasing = TRUE), ])
  colnames(rnk)[1] <- "GeneID"
  rnk <- rnk[, c("Symbol", "log2FoldChange")]
  write.table(rnk, file = file.path(out,paste0(sva, ".", .name, ".rnk")),
              col.names = FALSE, row.names =FALSE,
              sep = "\t", quote = FALSE)
}, degs_files, names(degs_files), SIMPLIFY = FALSE)

  

### DEG analysis without the pairs with heavy contamination
genes <- read.delim("IR_rate.genes.for.DE.analysis.txt", header = FALSE)[, 1]

counts_summary <- readRDS("Kevin_CD1A_14_full.counts.RDS")
metadata <- counts_summary$metadata
metadata$sample_name <- factor(metadata$sample_name,
                               levels = metadata$sample_name)
metadata$group <- gsub("CD1A\\(-\\)", "CD1A_N", metadata$group)
metadata$group <- gsub("CD1A\\(\\+\\)", "CD1A_P", metadata$group)
metadata$group <- factor(metadata$group, levels = unique(metadata$group))
metadata$batch <- factor(metadata$batch)

stopifnot(all(colnames(counts_summary$gtf$counts) == metadata$sample_name))

salmon <- readRDS("Kevin_14sample.salmon.count.RDS")
stopifnot(all(colnames(salmon$counts) == metadata$sample_name))

metadata <- metadata[metadata$batch != 2, ]
# metadata <- metadata[metadata$batch != 3, ]

metadata$batch <- factor(metadata$batch)

counts <- round(salmon$counts)
counts <- counts[, metadata$sample_name]
#threshold <- 2
stopifnot(all(colnames(counts) == metadata$sample_name))

# cpms = cpm(counts)
# keep = rowSums(cpms > 1) >= threshold
# keep = rowSums(cpms > 0.5) >= threshold
# keep = rowSums(cpms > median(10/(colSums(counts)/10^6))) >= threshold

counts = counts[genes, ]

#counts <- counts[sapply(t(as.data.frame(counts)), sd )!= 0, ]
all(colnames(counts) == metadata$sample_name)


dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                colData = metadata,
                                design = ~0 + group + batch)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType="parametric",
                           maxit=1000)

out <- "10.cleansamples.out"
if (!dir.exists(out)){
  dir.create(out)
}

## exploratory analysis
vsd <- vst(dds, blind = TRUE)
sampleDists <- dist(t(assay(vsd)))

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

pdf(file.path(out,"Fig.2.Heatmap showing raw sample distances.pdf"),
    width = 6.5, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = vsd$sample_name)
dev.off()

## PCA plot
vsd_filter <- t(assay(vsd))
vsd_filter <- vsd_filter[, unname(sapply(as.data.frame(vsd_filter), sd))!= 0]
pc <- prcomp(vsd_filter, scale = T, center = T, retx = TRUE)
pdf(file = file.path(out,"Fig 3. PCA.plot.pdf"), width = 6.5, height = 5)
print(plotPCA(pc, metadata))
dev.off()

normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))

# using SVA to adjust for hidden variations

## adjust log(cpm)
get.adj <- function(expr.data=NULL, design=NULL, meta.data=NULL)
{
  ## adjusted expression using voom()
  v <- voom(expr.data, design=design)
  fit <- lmFit(v, design=design)
  group_num <- nlevels(as.factor(meta.data$group))
  if (ncol(design) > group_num)
  {
    col.covariate <- (group_num + 1):ncol(design)
    adj.exp <- v$E - fit$coefficients[, col.covariate] %*%
      t(design[, col.covariate])
  } else {
    adj.exp <- v$E
  }
  adj.exp
}

model <- model.matrix(~ 0 + group + batch, data=metadata)
adj.exp <- get.adj(expr.data=normCounts,
                   design=model, meta.data = metadata)



## adjust expression for hidden variations for EDA plots
write.table(adj.exp, file.path(out, "Table 2.batch-adjusted.gene.expression.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

## Heatmap showing sample distances after adjusting for hidden variation
sampleDists <- dist(t(adj.exp))
pdf(file.path(out, "Fig 2.2.Heatmap showing batch-adjusted sample distances.pdf"),
    width = 6, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = paste0(vsd$sample_name))
dev.off()

## PCA using adjusted value
pca <- prcomp(t(adj.exp), scale = T, center = T, retx = TRUE)
pdf(file = file.path(out,"Fig 3.2.PCA.plot.batch.sva.adjusted.samples.pdf"),
    width = 6.5, height = 5)
print(plotPCA(pca, metadata))
dev.off()

pdf(file.path(out,"Figure 4.exploratory.data.analysis.before.qsmooth.pdf"),
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
           main = "after filtering")
#legend('topright', legend = levels(metadata$group),
#       col = 1:length(unique(metadata$group)),lty = 1)

dev.off()


dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                  colData = metadata,
                                  design = ~ 0 + group + batch)


dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType="parametric",
                           maxit=1000)
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

contrast.matrix <- matrix(c(-1, 1, 0, 0),
                          nrow = 1, byrow = TRUE)
rownames(contrast.matrix) <- "CD1A_P-CD1A_N"

DESeq_out <- mapply(output_DESeq, 1:nrow(contrast.matrix),
                    rownames(contrast.matrix),
                    MoreArgs = list(dds = dds,
                                    contrast.matrix = contrast.matrix,
                                    threshold = 1, shrink = FALSE), SIMPLIFY = FALSE)

write.xlsx(DESeq_out, file =
             file.path(out, "Table 3.2.DESeq.differential expressed genes.xlsx"),
           asTable = FALSE,
           overwrite = TRUE,
           rowNames = TRUE,
           sheetName = rownames(contrast.matrix))

#### volcano plots
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
                       legendIconSize = 5,
                       arrowheads = FALSE,
                       drawConnectors = TRUE,
                       widthConnectors = 0.5,
                       pointSize = 0.8,
                       #xlim = ifelse(i ==1, c(-5, 18), c(-6, 18)),
                       labSize = 2)
  print(p)
}
dev.off()




## Hallmark GSEA
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\DNA conatmination\data\Michelle.RNA-seq\03082024\Rnk.files\Hallmark.gsea.out)")

gsea_f <- dir(".", ".txt$")
names(gsea_f) <- gsub("nonSVA.(.+?)_h.all.(.+?).txt", "\\1_\\2", gsea_f)

gsea <- mapply(function(.x, .name){
  dat <- read.delim(.x, header = TRUE, as.is = TRUE)
  dat$correction <- gsub("(.+)_.+", "\\1", .name)
  dat$direction <- gsub(".+_(.+)", "\\1", .name)
  dat <- dat[, c(1, 4, 6, 8, 13, 14)]
},
       gsea_f, names(gsea_f), SIMPLIFY = FALSE)

sig.genesets_neg <- unique(do.call("c", sapply(gsea[seq(1,10, by =2)], function(.x){
  .x$NAME[.x[,4] <= 0.25]
})))

sig.neg_all <- lapply(gsea[seq(1,10, by =2)],
                           function(.x){
                             dat <-.x[.x$NAME %in% sig.genesets_neg, ]
                           })
sig.neg_all <- do.call("rbind", sig.neg_all)


sig.genesets_pos <- unique(do.call("c", sapply(gsea[seq(2,10, by =2)], function(.x){
  .x$NAME[.x[,4] <= 0.25]
})))

sig.pos_all <- lapply(gsea[seq(2,10, by =2)],
                      function(.x){
                        dat <- .x[.x$NAME %in% sig.genesets_pos, ]
                      })
sig.pos_all <- do.call("rbind", sig.pos_all)
sig.all <- rbind(sig.pos_all, sig.neg_all)
sig.all$NAME <- gsub("_", " ", gsub("HALLMARK_", "", sig.all$NAME))
sig.all$direction <- ifelse(sig.all$direction == "pos", "Postively enriched",
                            "Negatively enriched")
sig.all$correction <- factor(sig.all$correction,
                             levels = c("Batch(+);10samples",
                                        "Batch(+);IR_rate",
                                        "Batch(+);None",
                                        "Batch(+);GC%",
                                        "Batch(+);Global"),
                             labels = c("Batch(+);10samples",
                                        "Batch(+);IR%",
                                        "Batch(+);None",
                                        "Batch(+);GC%",
                                        "Batch(+);Global"))

library(ggplot2)
library(scales)

p <- ggplot(sig.all, aes(x = NES, y = NAME, color = -log10(FDR.q.val + 10^(-5)))) +
  geom_point()+ geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dotted") +
  facet_grid(direction~correction,
                          scales = "free_y",
                          space = "free") +
  scale_color_gradient2(low = muted("red"), mid = "white",
                        midpoint = 2.5,  high = muted("blue"),
                      name = expression(-log[10](FDR))) +
   ylab(NULL) + xlab("Normalized enrichment score")+
  theme(text = element_text(size = 6, face = "bold"),
        strip.text = element_text(size = 6, face = "bold"),
        legend.key.width = unit(0.3,"cm"))
pdf("Figure 007.GSEA.Hallmark.pdf", height = 3.3, width = 7)
print(p)
dev.off()



### 9 cleansample one contaminated
library("sva")
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
library(quantro)

counts_summary <- readRDS("Kevin_CD1A_14_full.counts.RDS")
metadata <- counts_summary$metadata
metadata$sample_name <- factor(metadata$sample_name,
                               levels = metadata$sample_name)
metadata$group <- gsub("CD1A\\(-\\)", "CD1A_N", metadata$group)
metadata$group <- gsub("CD1A\\(\\+\\)", "CD1A_P", metadata$group)
metadata$group <- factor(metadata$group, levels = unique(metadata$group))
metadata$batch <- factor(metadata$batch)

stopifnot(all(colnames(counts_summary$gtf$counts) == metadata$sample_name))

salmon <- readRDS("Kevin_14sample.salmon.count.RDS")
stopifnot(all(colnames(salmon$counts) == metadata$sample_name))

## 9 clean 1 dirty
# excluded_samples <-  c("CD1A(-)_m2_1",
#                        "CD1A(+)_m2_1",
#                        "CD1A(-)_m3_2",
#                        "CD1A(+)_m3_2")

## 11 clean 1 dirty-1
#excluded_samples <- c("CD1A(-)_m2_1",
#                      "CD1A(+)_m2_1")
##  11 clean 1 dirty-2
excluded_samples <- c("CD1A(-)_m3_2",
                      "CD1A(+)_m3_2")


metadata <- metadata[!metadata$sample_name %in% excluded_samples, ]

metadata$batch <- factor(metadata$batch)

salmon$counts <- round(salmon$counts)
salmon$counts <- salmon$counts[, metadata$sample_name]

intergenic_stats <- check_read_distribution(featurecounts_list =
                                          counts_summary,
                                           metadata = metadata)$IR_rate
intergenic_stats <- assign_stats[!rownames(assign_stats) %in% excluded_samples, ,
                                 drop = FALSE]

intergenic_stats <- intergenic_stats[as.character(metadata$sample_name), ,
                                     drop =FALSE]

stopifnot(all(metadata$sample_name == rownames(intergenic_stats)))
metadata$intergenic_rate <- intergenic_stats$assigned_percent

## scale for GLM stable
if (sd(metadata$intergenic_rate) > 5) {
  metadata$intergenic_rate <- scale(metadata$intergenic_rate,
                                    center = TRUE, scale = FALSE)
}
## SVA, Batch, but no intergenic rate
thresholds  <- 2
correction_methods <- c("IRrate")
use_svas <- c("noSVA")
for (use_sva in use_svas){
  for (threshold in thresholds){
    for (correction_method in correction_methods)
    {
      # experiment_name <- paste0("9clean1dirty", ".", use_sva, ".",
      #                           correction_method,
      #                           ".cpm.", threshold)

      experiment_name <- paste0("10clean2dirty", ".", use_sva, ".",
                                 correction_method,
                                 ".cpm.", threshold)

      # experiment_name <- paste0("9clean1dirty-2", ".", use_sva, ".",
      #                           correction_method,
      #                           ".cpm.", threshold)



      if (!dir.exists(experiment_name))
      {
        dir.create(experiment_name, recursive = TRUE)
      }

      ### correction
      if (correction_method == "global")
      {
        ## global correction
        counts <- corrected_counts
      } else if (correction_method == "gc") {
        counts <- GCcorrected_counts
      } else {
        counts <- round(salmon$counts)
      }

      # cpms = cpm(counts)
      # keep = rowSums(cpms > 1) >= threshold
      # keep = rowSums(cpms > 0.5) >= threshold
      # keep = rowSums(cpms > median(10/(colSums(counts)/10^6))) >= threshold

      counts = counts[genes, ]
      #rm("cpms", "keep")
      #counts <- counts[sapply(as.data.frame(counts), sd )!= 0, ]
      all(colnames(counts) == metadata$sample_name)

      if (correction_method == "global" ||
          correction_method == "gc" ||
          correction_method == "none")
      {
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                      colData = metadata,
                                      design = ~0 + group + batch)
      }  else {
        dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                      colData = metadata,
                                      design = ~0 + group + batch + intergenic_rate)
      }

      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds, fitType="parametric",
                                 maxit=1000)
      ## exploratory analysis
      vsd <- vst(dds, blind = TRUE)
      sampleDists <- dist(t(assay(vsd)))

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

      out <- experiment_name

      pdf(file.path(out,"Fig.2.Heatmap showing raw sample distances.pdf"),
          width = 6.5, height = 5)
      sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                       sampleNames = vsd$sample_name)
      dev.off()

      ## PCA plot
      vsd_filter <- t(assay(vsd))
      vsd_filter <- vsd_filter[, unname(sapply(as.data.frame(vsd_filter), sd))!= 0]
      pc <- prcomp(vsd_filter, scale = T, center = T, retx = TRUE)
      pdf(file = file.path(out,"Fig 3. PCA.plot.pdf"), width = 6.5, height = 5)
      print(plotPCA(pc, metadata))
      dev.off()

      normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))

      # using SVA to adjust for hidden variations
      if (use_sva == "SVA")
      {
        get.sva <- function (expr.data=NULL, meta.data=NULL,
                             mod = NULL, mod0 = NULL)
        {
          num.sva <-svaseq(as.matrix(expr.data), mod, mod0)$n.sv
          if (num.sva >=1)
          {
            sv <- svaseq(as.matrix(expr.data), mod, mod0, n.sv=num.sva)$sv
            colnames(sv)<- paste0("sv", seq_len(num.sva))
            meta.data.sva <-cbind(meta.data, sv, num.sva = num.sva)
          } else {
            meta.data.sva <- cbind(meta.data, num.sva = 0)
          }
          meta.data.sva
        }
      }


      ## adjust log(cpm)
      get.adj <- function(expr.data=NULL, design=NULL, meta.data=NULL)
      {
        ## adjusted expression using voom()
        v <- voom(expr.data, design=design)
        fit <- lmFit(v, design=design)
        group_num <- nlevels(as.factor(meta.data$group))
        if (ncol(design) > group_num)
        {
          col.covariate <- (group_num + 1):ncol(design)
          adj.exp <- v$E - fit$coefficients[, col.covariate] %*%
            t(design[, col.covariate])
        } else {
          adj.exp <- v$E
        }
        adj.exp
      }

      if (correction_method == "global" ||
          correction_method == "gc" ||
          correction_method == "none")
      {

        ## get sva: 4 hidden variables
        if (use_sva == "SVA")
        {
          meta.sva <- get.sva(expr.data = normCounts,
                              meta.data = metadata,
                              mod = model.matrix(~0 + group + batch,
                                                 data = metadata),
                              mod0 = model.matrix(~ batch, data = metadata))
          if (meta.sva$num.sva[1] > 0)
          {
            model <- model.matrix(as.formula(paste0("~ 0 + group + batch +",
                                                    paste(paste0("sv",
                                                                 seq_len(meta.sva$num.sva[1])),
                                                          collapse  = "+"))),
                                  data=meta.sva)
          } else {
            model <- model.matrix(~ 0 + group + batch, data=meta.sva)
          }
          adj.exp <- get.adj(expr.data=normCounts,
                             design=model, meta.data = meta.sva)
        } else {
          model <- model.matrix(~ 0 + group + batch, data=metadata)
          adj.exp <- get.adj(expr.data=normCounts,
                             design=model, meta.data = metadata)
        }
      } else {  # IR_rate
        ## get sva: 4 hidden variables
        if (use_sva == "SVA")
        {
          meta.sva <- get.sva(expr.data = normCounts,
                              meta.data = metadata,
                              mod = model.matrix(~0 + group + batch + intergenic_rate,
                                                 data = metadata),
                              mod0 = model.matrix(~ batch + intergenic_rate, data = metadata))
          if (meta.sva$num.sva[1] > 0) {
            model <- model.matrix(as.formula(paste0("~ 0 + group + batch + intergenic_rate +",
                                                    paste(paste0("sv", 1:meta.sva$num.sva[1]),
                                                          collapse  = "+"))),
                                  data=meta.sva)
          } else {
            model <- model.matrix(~0 + group + batch + intergenic_rate,
                                  data = meta.sva)
          }
          adj.exp <- get.adj(expr.data=normCounts,
                             design=model, meta.data=meta.sva)
        } else {
          model <- model.matrix(~ 0 + group + batch + intergenic_rate,
                                data=metadata)
          adj.exp <- get.adj(expr.data=normCounts,
                             design=model, meta.data = metadata)
        }
      }


      ## adjust expression for hidden variations for EDA plots
      write.table(adj.exp, file.path(out, "Table 2.batch-adjusted.gene.expression.txt"),
                  sep = "\t", quote = FALSE, row.names = TRUE)

      ## Heatmap showing sample distances after adjusting for hidden variation
      sampleDists <- dist(t(adj.exp))
      pdf(file.path(out, "Fig 2.2.Heatmap showing batch-adjusted sample distances.pdf"),
          width = 6, height = 5)
      sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                       sampleNames = paste0(vsd$sample_name))
      dev.off()

      ## PCA using adjusted value
      pca <- prcomp(t(adj.exp), scale = T, center = T, retx = TRUE)
      pdf(file = file.path(out,"Fig 3.2.PCA.plot.batch.sva.adjusted.samples.pdf"),
          width = 6.5, height = 5)
      print(plotPCA(pca, metadata))
      dev.off()

      pdf(file.path(out,"Figure 4.exploratory.data.analysis.before.qsmooth.pdf"),
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
                 main = "after filtering")
      #legend('topright', legend = levels(metadata$group),
      #       col = 1:length(unique(metadata$group)),lty = 1)

      dev.off()

      if(correction_method == "global" ||
         correction_method == "gc" ||
         correction_method == "none")
      {
        if (use_sva == "SVA" && meta.sva$num.sva[1] > 0)
        {
          dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                        colData = meta.sva,
                                        design = as.formula(paste0("~ 0 + group + batch +",
                                                                   paste(paste0("sv", 1:meta.sva$num.sva[1]),
                                                                         collapse  = "+"))))
        } else {
          dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                        colData = metadata,
                                        design = ~ 0 + group + batch)
        }
      } else {  ## IR_rate
        if (use_sva == "SVA" && meta.sva$num.sva[1] > 0)
        {
          dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                        colData = meta.sva,
                                        design = as.formula(paste0("~ 0 + group + batch + intergenic_rate +",
                                                                   paste(paste0("sv", 1:meta.sva$num.sva[1]),
                                                                         collapse  = "+"))))
        } else {
          dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                        colData = metadata,
                                        design = ~ 0 + group + batch + intergenic_rate)
        }
      }

      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds, fitType="parametric",
                                 maxit=1000)
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


      ## 10 samples
      contrast.matrix <- matrix(c(-1, 1, 0, 0),
                                nrow = 1, byrow = TRUE)

      ## 12 samples
      #contrast.matrix <- matrix(c(-1, 1, 0, 0, 0),
      #                          nrow = 1, byrow = TRUE)
      rownames(contrast.matrix) <- "CD1A_N-CD1AP"

      if (correction_method == "IRrate"){
        contrast.matrix <- cbind(contrast.matrix, 0)
      }

      if (use_sva == "SVA" && meta.sva$num.sva[1] >= 1)
      {
        for (i in seq_len(meta.sva$num.sva[1]))
        {
          contrast.matrix <- cbind(contrast.matrix, rep(0, nrow(contrast.matrix)))
        }
      }

      DESeq_out <- mapply(output_DESeq, 1:nrow(contrast.matrix),
                          rownames(contrast.matrix),
                          MoreArgs = list(dds = dds,
                                          contrast.matrix = contrast.matrix,
                                          threshold = 1, shrink = FALSE), SIMPLIFY = FALSE)

      write.xlsx(DESeq_out, file =
                   file.path(out, "Table 3.2.DESeq.differential expressed genes.xlsx"),
                 asTable = FALSE,
                 overwrite = TRUE,
                 rowNames = TRUE,
                 sheetName = rownames(contrast.matrix))

      #### volcano plots
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
                             legendIconSize = 5,
                             arrowheads = FALSE,
                             drawConnectors = TRUE,
                             widthConnectors = 0.5,
                             pointSize = 0.8,
                             #xlim = ifelse(i ==1, c(-5, 18), c(-6, 18)),
                             labSize = 2)
        print(p)
      }
      dev.off()
    }
  }
}



#### plot library size based on Salmon results
library("ggplot2")
frag_files <- dir(".", "flenDist.txt$", recursive = TRUE, full.names = TRUE)
names(frag_files) <- gsub("./(.+?)/.+", "\\1", frag_files)

frag_dat <- mapply(function(.x, .name){
  frag <- t(read.delim(.x, header = FALSE))
  frag <- as.data.frame(cbind(frag, 0:1000))
  colnames(frag) <- c("density", "size")
  frag$sample <- .name
  frag
}, frag_files, names(frag_files), SIMPLIFY = FALSE)

frag_dat <- do.call(rbind, frag_dat)

metadata <- read.delim("../../docs/Michelle.RNAseq-salmon.quant.txt")
metadata$samples <- gsub("results/010.salmon.quant.out/(.+?)/quant.sf",
                         "\\1", metadata$salmon_quant_file)

frag_dat <- merge(frag_dat, metadata[, c(2,5)], by.x = "sample", by.y ="samples",
                  all.x = TRUE)
p <- ggplot(frag_dat, aes(x = size, y = density, color = sample_name)) +
  geom_point(size = 0.5) +
  #geom_smooth(method = "loess") +
  geom_line(linewidth = 0.5 ) + facet_wrap(~sample_name, nrow = 5)
pdf("Kevin.RNAseq.lib.size-nonsmoothed.pdf", width = 15, height = 15)
p + xlim(c(0, 1000))
dev.off()


