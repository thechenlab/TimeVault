library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(purrr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(biomaRt)

# Set data and result path
data_path <- '/data/qiyu/mRNAcap/data'
result_path <- '/data/qiyu/mRNAcap/result'
figure_path <- '/data/qiyu/mRNAcap/figure'

# Set input sample information
batch_list <- c('241220_4SUdata_22J5GLLT4', '250104_241205_4SUexp', '250125_4SUtrack_Nova')
batch = batch_list[1]
if (batch == batch_list[1]) {
    timepoint_list = c('No', '3hr', '6', '9', '12', '24', '36', '48')
    replicte_num = 3
} else if (batch == batch_list[2]) {
    timepoint_list = c('No', '3hr', '6', '9', '12', '24', '36', '48')
    replicte_num = 3
} else if (batch == batch_list[3]) {
    timepoint_list = c('No', '3hr', '6', '9', '12', '24', '36', '48')
    replicte_num = 3
} else {
    stop("batch is not in the batch_list")
}


# Create metadata for input samples
replicate_df <- data.frame(
    timepoint = as.factor(timepoint_list),
    batch = rep(seq(replicte_num), length(timepoint_list)),
    sample = paste0(timepoint_list, rep(1:replicte_num, each = length(timepoint_list)))
)
all_rna_path <- paste0(data_path, "/", batch, '/4_counts/all/')
new_rna_path <- paste0(data_path, "/", batch, '/4_counts/new/')
sample_list <- list.files(path = all_rna_path, pattern = "matrix_counts_genename.txt$", full.names = FALSE)
sample_list <- sapply(sample_list, function(x) strsplit(x, "_matrix_counts_genename.txt")[[1]][1]) %>% {names(.) <- NULL;.}
if (! (all(sample_list %in% replicate_df$sample) && all(replicate_df$sample %in% sample_list))) {
    stop("sample_list does not match with input timepoint_list")
}


# Quadratic model
half_life_quadratic <- function(exp_values, time_values) {
    exp_values <- as.numeric(exp_values)
    time_values <- as.numeric(time_values)
    if (length(exp_values) < 3) return(NA)

    # Fit a quadratic polynomial model: y = a * t^2 + b * t + c
    model <- lm(exp_values ~ poly(time_values, 2, raw = TRUE))
    coef <- coef(model)
    a <- coef[3]  # Quadratic term coefficient
    b <- coef[2]  # Linear term coefficient
    c <- coef[1]  # Intercept
    t_peak <- -b / (2 * a)
    y_peak <- a * t_peak^2 + b * t_peak + c
    half_y <- y_peak / 2

    time_half <- tryCatch({
        roots <- polyroot(c(c - half_y, b, a)) 
        roots_real <- Re(roots[abs(Im(roots)) < 1e-6])
        half_life <- roots_real[roots_real > t_peak]
        if (length(half_life) == 0) return(NA)
        return(half_life[1] - t_peak)
    }, error = function(e) {
        return(NA)
    })

    return(time_half)
}


# Read DEG files for new and all samples
read_deg <- function(path, x) {
    deg_path <- file.path(path, paste0(x, '_matrix_counts_genename.txt'))
    deg <- read.table(deg_path, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
    colnames(deg) <- x
    return(deg)
}
new_deg <- map(sample_list, ~ read_deg(new_rna_path, .x)) %>% bind_cols() %>% as.data.frame()   
all_deg <- map(sample_list, ~ read_deg(all_rna_path, .x)) %>% bind_cols() %>% as.data.frame()
ratio_df <- map2(new_deg, all_deg, ~ ifelse(.y == 0, 0, .x / .y)) %>% bind_cols() %>% as.data.frame()
rownames(ratio_df) <- rownames(new_deg)

if (!all(colnames(ratio_df) %in% replicate_df$sample)) {
    stop("Column names of ratio_df do not match replicate_df$sample")
}
time_values <- seq(length(timepoint_list)) - 1
sample_groups <- split(replicate_df$sample, replicate_df$timepoint)
average_ratio <- lapply(sample_groups, function(samples) {rowMeans(ratio_df[, samples, drop = FALSE], na.rm = TRUE)}) %>% as.data.frame()   
rownames(average_ratio) <- rownames(ratio_df)
colnames(average_ratio) <- names(sample_groups)
average_ratio <- average_ratio[, timepoint_list]
colnames(average_ratio) <- time_values  


# Compute half-life for each gene in average_ratio 
half_life_df <- average_ratio %>%
    filter(ncol(select(., where(is.numeric))) == length(time_values)) %>%
    mutate(half_life = pmap_dbl(select(., where(is.numeric)), ~{
        values <- c(...)
        if (length(values) == length(time_values)) {
            half_life_quadratic(values, time_values)
        } else {
            NA_real_
        }
    })) 
colnames(half_life_df) <- c(timepoint_list, "half_life")
dir.create(file.path(result_path, batch), recursive = TRUE, showWarnings = FALSE)
write.csv(half_life_df, file.path(result_path, batch, "half_life_average_df.csv"), row.names = TRUE)


# Comput half-life for each replicates
replicate_half_life_list <- lapply(unique(replicate_df$replicate), function(r) {
    replicate_name_list = unique(replicate_df[replicate_df$replicate == r, "sample"])
    half_life_df = ratio_df[, replicate_name_list] %>%
        dplyr::filter(ncol(dplyr::select(., where(is.numeric))) == length(time_values)) %>%
        dplyr::mutate(half_life = pmap_dbl(dplyr::select(., where(is.numeric)), ~{
            values <- c(...)
            if (length(values) == length(time_values)) {
                half_life_quadratic(values, time_values)
            } else {
                NA_real_
            }
        })) %>% dplyr::select(half_life)
    colnames(half_life_df) <- r
    return(half_life_df)
})

replicate_half_life <- bind_cols(replicate_half_life_list)
write.csv(replicate_half_life, file.path(result_path, batch, "replicate_half_life_df.csv"), row.names = TRUE)


# all NA genes to do enrichment
gene_list <- half_life_df %>% dplyr::filter(!is.na(half_life)) %>% rownames

# biomaRt to get entrez id
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
biomaRt_symbol_entrez <- getBM(attributes = c("external_gene_name", "entrezgene_id"), filters = "external_gene_name", values = gene_list, mart = ensembl)
biomaRt_ensembl_entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = "ensembl_gene_id", values = gene_list, mart = ensembl)
colnames(biomaRt_symbol_entrez) <- c("gene", "ENTREZID")
colnames(biomaRt_ensembl_entrez) <- c("gene", "ENTREZID")
biomaRt_entrez_df <- bind_rows(biomaRt_symbol_entrez, biomaRt_ensembl_entrez)
biomaRt_entrez_df$ENTREZID <- as.character(biomaRt_entrez_df$ENTREZID)

# bitr to get entrez id
symbol_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>% rename(gene = SYMBOL) 
ensembl_entrez <- bitr(gene_list[grepl("^ENSG", gene_list)], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>% rename(gene = ENSEMBL)
colnames(gene_entrez_df) <- c("gene", "ENTREZID")
gene_entrez_df <- bind_rows(symbol_entrez, ensembl_entrez)
gene_entrez_df$ENTREZID <- as.character(gene_entrez_df$ENTREZID)

# combine biomaRt and bitr results
entrez_df <- bind_rows(biomaRt_entrez_df, gene_entrez_df) %>% 
    distinct() %>%
    group_by(gene) %>%
    dplyr::filter(!is.na(ENTREZID) | n() == 1) %>%  
    slice(1) %>%  
    ungroup()
entrez_ids <- entrez_df %>% dplyr::filter(!is.na(ENTREZID)) %>% pull(ENTREZID)

# GO enrichment
go_result <- as.data.frame(enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE))
go_df_grouped <- go_result %>% 
    separate_rows(geneID, sep = "/") %>% 
    as.data.frame() %>%
    dplyr::select(geneID, ID, ONTOLOGY, Description) %>% 
    group_by(geneID) %>%
    summarise(
        Description = paste(unique(Description), collapse = "; "),
        GO_terms = paste(unique(ID), collapse = "; "), 
        ONTOLOGY = paste(unique(ONTOLOGY), collapse = "; "),
        .groups = "drop"
    ) 


# KEGG enrichment
kegg_result <- as.data.frame(enrichKEGG(gene = entrez_ids, organism = "hsa"))
kegg_df_grouped <- kegg_result %>% 
    separate_rows(geneID, sep = "/") %>% 
    as.data.frame() %>%
    dplyr::select(geneID, ID, Description, category, subcategory) %>% 
    group_by(geneID) %>%
    summarise(
        Description = paste(unique(Description), collapse = "; "), 
        Category = paste(unique(category), collapse = "; "), 
        Subcategory = paste(unique(subcategory), collapse = "; "), 
        KEGG_ID = paste(unique(ID), collapse = "; "), 
        .groups = "drop"
    ) %>% 
    left_join(entrez_df, ., by = c("ENTREZID" = "geneID"))

# save results
write.csv(go_df_grouped, file.path(result_path, batch, "go_annotation_df.csv"), row.names = TRUE)
write.csv(kegg_df_grouped, file.path(result_path, batch, "kegg_annotation_df.csv"), row.names = TRUE)