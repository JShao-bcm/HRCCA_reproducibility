library(edgeR)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(stringr)
library(tidyr)

#conda env: mamba activate /pub/jianms1/local/miniconda/envs/Seurat
#adata <- LoadH5Seurat("../Endothelial_cell_filter_subclass_final_rawcounts.h5seurat", meta.data = FALSE, misc = FALSE)
#metadata <- read.table("../Endothelial_cell_filter_subclass_final_rawcounts_obs.csv.gz", header = TRUE, sep = ",", row.names = 1)
#adata <- AddMetaData(object = adata, metadata = metadata)
setwd("/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_clean/annotation/Regional_DEGs/Majorclass")

obj=readRDS("/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_clean/annotation/snRNA_clean_all_rawcounts_final.rds")

# Filter to keep only cells with percent.mt less than or equal to 5
adata <- subset(obj, subset = percent.mt <= 5)

# Verify the filtering by checking the range of the new object
summary(adata$percent.mt)

# Check how many cells were removed
print(paste("Cells before:", ncol(obj)))
print(paste("Cells after:", ncol(adata)))

# 2. Filter low quality Endothelial cells
#adata <- subset(adata, subset = Endothelial_subclass %in% c('Stressed ECs', 'Endothelial cells LowQual'), invert = TRUE)
adata$Age <- as.character(adata$Age)

pca <- read.table("/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/DEG_scRNA_snRNA/DEG_edgeR/sample_PCA_for_DEG.txt", header=TRUE, sep="\t")

#celltypes <- c("Arterial", "Arteriole", "Fenestrated capillary", "Vein", "Venule") 
celltypes <- unique(adata@meta.data$snRNA_majorclass)

for(ct in celltypes){
    print(paste0("Now processing cell type: ", ct))
    adata_sub <- subset(adata, subset = snRNA_majorclass == ct) 
    
    # --- NEW: Calculate Percent Expressed per Gene per Region ---
    mac_cells <- which(adata_sub$Region == "Macular")
    per_cells <- which(adata_sub$Region == "Peripheral")
    
    # Using the 'counts' slot for accurate detection of any RNA molecules
    sc_counts <- GetAssayData(adata_sub, assay = "RNA", slot = "counts")
    
    # Calculate % of cells with counts > 0
    pct_mac <- rowSums(sc_counts[, mac_cells] > 0) / length(mac_cells) * 100
    pct_per <- rowSums(sc_counts[, per_cells] > 0) / length(per_cells) * 100
    
    # Store in a reference dataframe for later merging
    pct_df <- data.frame(
        gene = rownames(adata_sub),
        pct_macular = round(pct_mac, 2),
        pct_peripheral = round(pct_per, 2)
    )
    # ----------------------------------------------------------

    # Aggregate counts: Sum counts by Donor and Region
    pb <- AggregateExpression(adata_sub, 
                              group.by = c("Donor", "Region"), 
                              assays = "RNA", 
                              slot = "counts",
                              return.seurat = TRUE)
    
    cell_counts <- as.data.frame(table(adata_sub$Donor, adata_sub$Region))
    colnames(cell_counts) <- c("Donor", "Region", "n_cells")
    cell_counts <- cell_counts %>% filter(n_cells >= 20) %>% mutate(sampleid_pb = paste0(Donor, "_", Region))

    valid_samples <- intersect(colnames(pb), cell_counts$sampleid_pb)
    if(length(valid_samples) == 0) next
    pb <- subset(pb, cells = valid_samples)

    age_map <- adata_sub@meta.data %>% select(Donor, Age) %>% distinct()
    sex_map <- adata_sub@meta.data %>% select(Donor, Gender) %>% distinct()

    metadata_pb <- data.frame(ids = colnames(pb)) %>%
        mutate(
            Donor = str_extract(ids, ".*(?=_)"), 
            Region = str_extract(ids, "[^_]+$"), 
            Donor_snRNA = paste0(Donor, "_snRNA")
        ) %>%
        left_join(pca, by = c("Donor_snRNA" = "sampleid")) %>% 
        left_join(age_map, by = "Donor") %>% 
        left_join(sex_map, by = "Donor")

    rownames(metadata_pb) <- metadata_pb$ids
    pb <- AddMetaData(pb, metadata_pb)

    # Prepare edgeR inputs
    counts_mat <- GetAssayData(pb, assay = "RNA", slot = "counts")
    colnames(counts_mat) <- as.character(colnames(counts_mat))

    meta <- pb@meta.data %>%
        mutate(
            sampleid = rownames(.),
            donor = factor(Donor),
            region = factor(Region, levels = c("Peripheral", "Macular")),
            lib_size = colSums(counts_mat)
        )

    # Filter for paired donors
    donors_with_both <- meta %>%
        group_by(donor) %>%
        filter(any(region == "Macular") & any(region == "Peripheral")) %>%
        pull(donor) %>% unique()
    
    if(length(donors_with_both) < 2){
        message(paste("Skipping", ct, "- insufficient pairs."))
        next
    }
    
    meta_paired <- meta %>% filter(donor %in% donors_with_both)
    counts_paired <- counts_mat[, meta_paired$sampleid]
    
	if(identical(colnames(counts_paired), meta_paired$sampleid)){print("smapleid matched")}
	
    meta_paired$donor <- droplevels(meta_paired$donor)
    meta_paired$region <- factor(meta_paired$region, levels = c("Peripheral", "Macular"))

    # edgeR Analysis
    y <- DGEList(counts = counts_paired, samples = meta_paired)
    y <- normLibSizes(y)
    
    # CPM Filtering
    #keep.genes <- rowMeans(cpm(y)) >= 5
	#keep.genes <- rowSums(cpm(counts_paired) > 0.1) >= 5
	keep.genes <- rowSums(cpm(counts_paired) >= 1) >= (ncol(counts_paired) * 0.5)
    y <- y[keep.genes, , keep.lib.sizes = FALSE]
    y <- normLibSizes(y)

    log2_tmm <- cpm(y, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
    
    write.csv(as.data.frame(log2_tmm), paste0("Log2_TMM_Normalized_Expression_",ct,".csv"))
    write.csv(y$samples, paste0("Metadata_",ct,".csv"))
    
    # Design and Model
    design <- model.matrix(~ donor + region , data = y$samples) #+ scale(lib_size)
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = "regionMacular")
	
    # --- NEW: Final Table Cleanup and Percentage Merge ---
    all.degs <- topTags(qlf, n = Inf)$table
    all.degs$gene <- rownames(all.degs)
    
    # Merge the percentages we calculated at the start of the loop
    all.degs <- all.degs %>%
        left_join(pct_df, by = "gene") %>%
        mutate(Total_Pairs = length(donors_with_both)) %>%
        select(gene, everything()) # Move gene name to the first column
    
    # Final Output
    ct_formatted <- gsub(" ", "_", ct)
	# Save the fit object to a file
	saveRDS(fit, file = paste0("snRNA_majorclass_", ct_formatted, "_Paired_CPM1_Macular_vs_Periphery_PCT_fit.rds"))
	## re-Load the saved model
	#fit_loaded <- readRDS("edgeR_fit_paired_region.rds")
	## Extract the DEGs for your specific region
	## Note: You still use the 'fit' object inside glmQLFTest
	#qlf_new <- glmQLFTest(fit_loaded, coef = "regionMacular")
	## View your top genes
	#topTags(qlf_new)
    outfile <- paste0("snRNA_majorclass_", ct_formatted, "_Paired_CPM1_Macular_vs_Periphery_PCT.csv")
    write.csv(all.degs, file = outfile, row.names = FALSE, quote = FALSE)
}
