###########################################################
################### MtAD - Kathy ######################
###########################################################

# ================================================
# ================================================
# ====== Step 0. Read in, data preprocessing =====
# ================================================
# ================================================

library(tidyverse)
library(purrr)
library(data.table)
library(pheatmap)

# Read in files
deg <- readRDS("/Users/spectremac/Desktop/mtAD/DEG_rds/AD_nonAD.sex.deseq2.rds")
dream <- readRDS("/Users/spectremac/Desktop/mtAD/DEG_rds/continuous_variable.sex.DREAM.rds")

# ====================================
# below is the sanity check
#class(deg)
#str(deg)

#class(cont_var_sex)
#str(cont_var_sex)
#length(cont_var_sex)  
#nm <- names(cont_var_sex) 
#nm

#names(cont_var_sex)[1]
#head(cont_var_sex[[1]], 10)
# ====================================

# make it readable
dream_df <- imap_dfr(dream, ~ mutate(.x, list_name = .y), .id = "idx") %>%
  as_tibble()

#dream_df

# data cleaning - pathology (amyloid.scaled -> amyloid)
dream_df <- dream_df %>%
  mutate(path_short = sub("\\.scaled$", "", pathology))

# Unify celltype names
dream_df <- dream_df %>%
  mutate(
    cell_clean = celltype |> 
      str_replace_all("\\.", "_") |>      # let . be _
      str_replace("^InN\\b", "Inh") |>    # InN → Inh
      str_replace("^InN_", "Inh_"),       # InN_ → Inh_
    # extract main types（Ast/Exc/Inh/Mic/Oli/OPC/Vas/T_cells）
    cell_major = case_when(
      str_starts(cell_clean, "Ast")      ~ "Ast",
      str_starts(cell_clean, "Exc")      ~ "Exc",
      str_starts(cell_clean, "Inh")      ~ "Inh",
      str_starts(cell_clean, "Mic")      ~ "Mic",
      str_starts(cell_clean, "Oli")      ~ "Oli",
      str_starts(cell_clean, "OPC")      ~ "OPC",
      str_starts(cell_clean, "Vas")      ~ "Vas",
      str_starts(cell_clean, "T_cells")  ~ "T_cells",
      TRUE ~ cell_clean
    )
  )

# see cleaning result
table(dream_df$cell_major)
#      Ast      Exc      Inh      Mic      Oli      OPC    T_cells    Vas 
#    1653781  6419500 12985817  1278812  1137824   544610   255337  1670165 

dream_df %>% count(celltype, cell_clean, cell_major, sort = TRUE) %>% head(20)


# ================================================
# ================================================
# ====== Step 1. Identify mitochondria genes =====
# ================================================
# ================================================

#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")   
#}
#BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "GO.db"),
#                     update = FALSE, ask = FALSE)


suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(stringr)
})

## Mt related GO terms
go_terms <- c(
  "GO:0005739",  # mitochondrion
  "GO:0005743",  # mitochondrial inner membrane
  "GO:0005741",  # mitochondrial outer membrane
  "GO:0005759",  # mitochondrial matrix
  "GO:0005758",  # intermembrane space
  "GO:0070469",  # respiratory chain
  "GO:0005740"   # mitochondrial envelope
)

## GO ALL infer human genes
mt_genes_go <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_terms,
  keytype = "GOALL",
  columns = "SYMBOL"
) %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL) %>%
  sort()

## Collect all the gene names
get_all_genes <- function() {
  if (exists("cont_var_sex")) {
    unique(unlist(lapply(cont_var_sex, `[[`, "gene")))
  } else if (file.exists("Continuous_variable.sex.DREAM.rds")) {
    x <- readRDS("Continuous_variable.sex.DREAM.rds")
    unique(unlist(lapply(x, `[[`, "gene")))
  } else if (exists("dream_df")) {
    unique(dream_df$gene)
  } else {
    stop("error")
  }
}
all_genes <- get_all_genes()

## MT- genes
mt_from_prefix <- sort(unique(all_genes[str_detect(all_genes, "^MT[-\\.]")]))

## Merge and limit in the dataset
mt_set_all   <- sort(unique(c(mt_genes_go, mt_from_prefix)))
mt_set_in_df <- sort(intersect(mt_set_all, all_genes))

cat("# GO infer MT related genes", length(mt_genes_go), "\n")
# 1828 

cat("# MT- genes", length(mt_from_prefix), "\n")
# MT- genes 26 

cat("# Combine (filter duplicate)", length(mt_set_all), "\n")
# Combine (filter duplicate) 1841 

cat("# Truly exists", length(mt_set_in_df), "\n")
# # Truly exists 1770

#head(mt_set_in_df)

#df <- data.frame(id = mt_set_in_df, stringsAsFactors = FALSE)
#df




# ================================================
# ================================================
# = Step 1.5 Boxplot - cell type/# sig. mt genes =
# ================================================
# ================================================

#################################################
################ All Pathways ###################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    

# all path
path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
                 "cognwo_demog_slope","CR_score")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Box plot ----------
counts_box <- dream_mt_sig %>%
  dplyr::group_by(cell_major, path_short) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop")

# Show 0 sig. genes
counts_box <- counts_box %>%
  tidyr::complete(cell_major, path_short,
                  fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::filter(!is.na(cell_major), !is.na(path_short))

ggplot(counts_box,
       aes(x = cell_major, y = n_sig_mt_genes)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = path_short), width = 0.12, height = 0,
              alpha = 0.7, size = 1.8, show.legend = TRUE) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes per pathology",
       title = sprintf("Distribution across pathologies (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min),
       color = "Pathology") +
  theme_bw(base_size = 12)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# make sure sex has a consistent order (optional)
# 1) Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  # (optional) include explicit zeros for combos with no sig genes
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))

# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


#============== Newly added ===============#


library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)  

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    

# all path
path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
                 "cognwo_demog_slope","CR_score")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)


##### Heatmap  #####

counts_hp <- dream_mt_sig %>%
  dplyr::filter(!is.na(cell_major), !is.na(path_short)) %>%
  dplyr::group_by(cell_major, path_short) %>%
  dplyr::summarise(
    n_sig   = dplyr::n_distinct(gene),
    med_lfc = stats::median(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # fill in 0s for combos with no sig genes
  tidyr::complete(cell_major, path_short, fill = list(n_sig = 0, med_lfc = NA))

# Order rows/cols by total counts 
row_order <- counts_hp %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(total = sum(n_sig), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total)) %>% dplyr::pull(cell_major)

col_order <- counts_hp %>%
  dplyr::group_by(path_short) %>%
  dplyr::summarise(total = sum(n_sig), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total)) %>% dplyr::pull(path_short)

counts_hp <- counts_hp %>%
  dplyr::mutate(
    cell_major = factor(cell_major, levels = row_order),
    path_short = factor(path_short, levels = col_order)
  )


ggplot(counts_hp, aes(x = path_short, y = cell_major, fill = n_sig)) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = ifelse(n_sig > 0, n_sig, "")), size = 3) +
  scale_fill_gradientn(
    colours = c("#ffffcc", "#fed976", "#fd8d3c", "#e31a1c"),
    trans = "sqrt",
    breaks = c(0,1,5,10,25,50,100,200,400),
    name = "# sig. mt genes"
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    barheight = grid::unit(200, "pt"),
    barwidth  = grid::unit(15,  "pt")
  )) +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  labs(x = "Pathology", y = "Cell type (major)", fill = "# sig. mt genes",
       title = "Significant mitochondrial genes per pathology × cell type") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8)
  )



##### Bubble Plot (indicate direction) #####

# Same counts object: counts_hp (n_sig, med_lfc)
cap <- 0.8
counts_hp <- counts_hp %>%
  dplyr::mutate(med_lfc_cap = pmin(pmax(med_lfc, -cap), cap))

ggplot(counts_hp, aes(x = path_short, y = cell_major)) +
  geom_point(aes(size = n_sig, color = med_lfc_cap), alpha = 0.9) +
  scale_size_area(max_size = 12, breaks = c(1,5,10,25,50,100,200)) +
  scale_color_gradient2(low = "#2c7bb6", mid = "grey90", high = "#d7191c",
                        midpoint = 0, name = "median logFC") +
  labs(x = "Pathology", y = "Cell type (major)", size = "# sig. mt genes",
       title = "Counts (size) and direction (color) by pathology × cell type") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

##### Heatmap  #####
# incorprate sex
if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale")
library(ggnewscale)

ggplot() +
  ## FEMALE panel -----------------------------------------------------------
geom_tile(
  data = dplyr::filter(counts_hp_sex, sex == "female"),
  aes(path_short, cell_major, fill = n_sig),
  color = "white", size = 0.3
) +
  geom_text(
    data = dplyr::filter(counts_hp_sex, sex == "female"),
    aes(path_short, cell_major, label = ifelse(n_sig > 0, n_sig, "")),
    size = 2.7
  ) +
  scale_fill_gradientn(
    colours = c("#ffffcc", "#fed976", "#fd8d3c", "#e31a1c"),
    trans = "sqrt",
    breaks = c(0,1,5,10,25,50,100,200,400),
    name = "# sig. mt genes",
    guide = guide_colourbar(
      title.position = "top",
      barheight = grid::unit(200, "pt"),
      barwidth  = grid::unit(15,  "pt")
    )
  ) +
  ggnewscale::new_scale_fill() +
  
  ## MALE panel -------------------------------------------------------------
geom_tile(
  data = dplyr::filter(counts_hp_sex, sex == "male"),
  aes(path_short, cell_major, fill = n_sig),
  color = "white", size = 0.3
) +
  geom_text(
    data = dplyr::filter(counts_hp_sex, sex == "male"),
    aes(path_short, cell_major, label = ifelse(n_sig > 0, n_sig, "")),
    size = 2.7
  ) +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5"),
    trans = "sqrt",
    breaks = c(0,1,5,10,25,50,100,200,400),
    name = "# sig. mt genes",
    guide = guide_colourbar(
      title.position = "top",
      barheight = grid::unit(200, "pt"),
      barwidth  = grid::unit(15,  "pt")
    )
  ) +
  
  facet_wrap(~ sex, nrow = 1) +
  labs(x = "Pathology", y = "Cell type (major)", fill = "# sig. mt genes",
       title = "Sig. mt genes per pathology × cell type (by sex)") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )












#################################################
#################### gpath ######################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("gpath")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))

pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


#################################################
################### amyloid #####################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("amyloid")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



#################################################
################### nft #####################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("nft")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



#################################################
################### tangles #####################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("tangles")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


#################################################
################### tangles #####################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("tangles")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


#################################################
################### plaq_d #####################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("plaq_d")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)


  # ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


#################################################
################### plaq_n #####################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("plaq_n")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



#################################################
################ cogn_global_lv #################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("cogn_global_lv")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



#################################################
################ cognep_demog_slope #################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("cognep_demog_slope")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")




#################################################
################ cogng_demog_slope #################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("cogng_demog_slope")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")




#################################################
################ cognpo_demog_slope #################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("cognpo_demog_slope")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


#################################################
################ cognps_demog_slope #################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("cognps_demog_slope")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  # (optional) include explicit zeros for combos with no sig genes
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")



#################################################
################ cognse_demog_slope #################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("cognse_demog_slope")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  # (optional) include explicit zeros for combos with no sig genes
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


#################################################
################ cognwo_demog_slope #################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("cognwo_demog_slope")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


#################################################
################ CR_score #################
#################################################

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    
# all path
#path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
#                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
#                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
#                 "cognwo_demog_slope","CR_score")  

# gpath
path_subset <- c("CR_score")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

#table(dream_mt_sig$celltype)

# ----------  Bar plot ----------
counts_bar <- dream_mt_sig %>%
  dplyr::group_by(cell_major) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sig_mt_genes))

ggplot(counts_bar,
       aes(x = reorder(cell_major, n_sig_mt_genes), y = n_sig_mt_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_mt_genes), vjust = -0.3, size = 3) +
  labs(x = "Cell type (major)",
       y = "# significant mitochondrial genes",
       title = sprintf("Per cell type: # sig. mt genes (FDR<%.2f, |logFC|≥%.1f)", alpha, lfc_min)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  expand_limits(y = max(counts_bar$n_sig_mt_genes) * 1.08)

# ---------- Sex-driven bar/box plots ----------
# bar plot
# Build counts per cell_major × sex
counts_bar_sex <- dream_mt_sig %>%
  dplyr::group_by(cell_major, sex) %>%
  dplyr::summarise(n_sig_mt_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
  # (optional) include explicit zeros for combos with no sig genes
  tidyr::complete(cell_major, sex, fill = list(n_sig_mt_genes = 0)) %>%
  dplyr::mutate(sex = factor(tolower(sex), c("female","male")))


# a little padding for the y-axis so labels fit
pad <- max(1, ceiling(0.07 * max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE)))

ggplot(counts_bar_sex,
       aes(x = cell_major, y = n_sig_mt_genes, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig_mt_genes),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(counts_bar_sex$n_sig_mt_genes, na.rm = TRUE) + pad),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cell type (major)", y = "# sig. mt genes", fill = "Sex",
       title = "Per cell type × sex") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")





# ================================================
# == Step 2. Heatmap for gpath, sig. mt genes  ==
# ================================================

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")

alpha   <- 0.05  
lfc_min <- 0.2  

## filter mt-genes + gpath
dream_mt0 <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df,
                path_short == "gpath") %>%
  dplyr::mutate(sex = tolower(sex))

## keep sig.
dream_sig <- dream_mt0 %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)

## cell_major (gene×sex×cell_major select median logFC）
dream_mt <- dream_sig %>%
  dplyr::group_by(gene, sex, cell_major) %>%
  dplyr::summarise(logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    sex = factor(sex, c("female","male")),
    col_key = paste(cell_major, sex, sep = "|")    # sex
  )

## only keep sig. combinations, others NA
mat_mt <- dream_mt %>%
  dplyr::select(gene, col_key, logFC) %>%
  tidyr::pivot_wider(names_from = col_key, values_from = logFC) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

## drop all NAs
if (ncol(mat_mt) > 0) {
  mat_mt <- mat_mt[, colSums(!is.na(mat_mt)) > 0, drop = FALSE]
}


if (any(mat_mt == 0, na.rm = TRUE)) {
  message("Converting 0s to NA (only appropriate if 0s were placeholders).")
  mat_mt[mat_mt == 0] <- NA
}

# drop rows that are entirely NA
mat_mt <- mat_mt[rowSums(!is.na(mat_mt)) > 0, , drop = FALSE]


## filter rows, pick sd highest genes
min_nonNA <- 1   # at least 1 col has value 
topN <- 150
keep_genes <- apply(mat_mt, 1, function(x) sum(!is.na(x))) >= min_nonNA
mat_mt2 <- mat_mt[keep_genes, , drop = FALSE]

v <- apply(mat_mt2, 1, stats::var, na.rm = TRUE)
ord_all <- order(v, decreasing = TRUE, na.last = NA)  # drop rows with var = NA
n_keep  <- min(topN, length(ord_all))                 # cap at available rows
mat_mt2 <- mat_mt2[ord_all[seq_len(n_keep)], , drop = FALSE]

# drop any rows that somehow ended up all-NA (paranoia check)
mat_mt2 <- mat_mt2[rowSums(!is.na(mat_mt2)) > 0, , drop = FALSE]

## Winsorize (replace extreme outliers)
cap <- 2
mat_cap <- pmin(pmax(mat_mt2, -cap), cap)

# verify row names are real
stopifnot(!any(is.na(rownames(mat_cap))), !any(rownames(mat_cap) == "NA", na.rm = TRUE))

# ======================= Newly added: =========================
# compute distances - pairwise complete observations
# replace any remaining NA distances with a neutral value (dist=1)


row_dist <- (function(m) {
  R <- cor(t(m), use = "pairwise.complete.obs", method = "pearson")
  D <- as.dist(1 - pmin(pmax(R, -1), 1))
  D[is.na(D)] <- 1  # no overlap -> treat as uncorrelated
  D
})(mat_cap)

col_dist <- (function(m) {
  R <- cor(m, use = "pairwise.complete.obs", method = "pearson")
  D <- as.dist(1 - pmin(pmax(R, -1), 1))
  D[is.na(D)] <- 1
  D
})(mat_cap)

# If there are too few rows/cols, disable clustering on that axis
cluster_rows_ok <- nrow(mat_cap) > 1
cluster_cols_ok <- ncol(mat_cap) > 1
# ================================================================


## annotation (only CellType & Sex)
ann_col <- data.frame(
  CellType = sub("\\|.*$", "", colnames(mat_cap)),
  Sex      = sub("^.*\\|", "", colnames(mat_cap)),
  row.names = colnames(mat_cap),
  stringsAsFactors = FALSE
)
ann_col$Sex <- tolower(ann_col$Sex)

## color
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
library(RColorBrewer)

cell_lvls <- sort(unique(ann_col$CellType))
sex_lvls  <- sort(unique(ann_col$Sex))

base_cols_cell <- RColorBrewer::brewer.pal(min(12, max(3, length(cell_lvls))), "Set3")
pal_cell <- setNames(grDevices::colorRampPalette(base_cols_cell)(length(cell_lvls)), cell_lvls)
pal_sex  <- setNames(c(female = "#e41a1c", male = "#377eb8")[sex_lvls], sex_lvls)

ann_colors <- list(CellType = pal_cell, Sex = pal_sex)

lim <- max(abs(mat_cap), na.rm = TRUE)
breaks <- seq(-lim, lim, length.out = 101)
col_fun <- colorRampPalette(c("#313695","#74add1","#ffffbf","#f46d43","#a50026"))

## plot
pheatmap(
  mat_cap,
  clustering_distance_rows = if (cluster_rows_ok) row_dist else "euclidean",
  clustering_distance_cols = if (cluster_cols_ok) col_dist else "euclidean",
  cluster_rows = cluster_rows_ok,
  cluster_cols = cluster_cols_ok,
  clustering_method = "average",
  color = col_fun(length(breaks)-1),
  breaks = breaks,
  na_col = "grey90",
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = ifelse(nrow(mat_cap) > 80, 6, 8),
  border_color = NA,
  annotation_col = ann_col,
  annotation_colors = ann_colors
)

sum(is.na(rownames(mat_cap)))         # should be 0
sum(rownames(mat_cap) == "NA")        # should be 0


#=====================A======================#
# Directional decomposition: Count/ratio of up vs down
# cell_major × pathology × sex 

# parameters
alpha   <- 0.05   
lfc_min <- 0.2    

# all path
path_subset <- c("amyloid","gpath","nft","tangles","plaq_d","plaq_n",
                 "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
                 "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
                 "cognwo_demog_slope","CR_score")  

# filter mt- and sig.
dream_mt <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(
    sex = tolower(sex),
    cell_major = factor(cell_major, c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells"))
  )

dream_mt <- dream_mt %>%
  dplyr::filter(path_short %in% path_subset)

dream_mt_sig <- dream_mt %>%
  dplyr::filter(!is.na(adj.P.Val), adj.P.Val < alpha,
                !is.na(logFC), abs(logFC) >= lfc_min)


dir_counts <- dream_mt_sig %>%
  mutate(direction = ifelse(logFC > 0, "up", "down")) %>%
  group_by(cell_major, path_short, direction) %>%
  summarise(n = n_distinct(gene), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(log2_up_down = log2((up + 0.5)/(down + 0.5)))  # 稳健比值（可热图）


ggplot(dir_counts, aes(path_short, cell_major)) +
  geom_tile(aes(fill = up + down), colour = "white", size = .3) +
  geom_text(aes(label = ifelse(up+down>0, up+down, "")), size = 2.6) +
  scale_fill_viridis_c(option = "C", trans = "sqrt") +
  facet_grid(~ "Counts") +
  theme_bw(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(dir_counts, aes(path_short, cell_major)) +
  geom_tile(aes(fill = log2_up_down), colour = "white", size = .3) +
  scale_fill_gradient2(low="#2b8cbe", mid="white", high="#e34a33", midpoint=0) +
  facet_grid(~ "log2(up/down)") +
  theme_bw(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))



#=====================B======================#
# 2) Recurrent Gene Map Repetitive Significance Across Pathological Variables

library(tidytext)

min_recur <- 2   
topN <- 20       

recur_plot <- recur %>%
  dplyr::filter(n_path >= min_recur) %>%
  dplyr::group_by(cell_major) %>%
  dplyr::arrange(dplyr::desc(n_path), dplyr::desc(abs(pos - neg)), gene, .by_group = TRUE) %>%
  dplyr::slice_head(n = topN) %>%
  dplyr::ungroup()

ggplot(recur_plot,
       aes(x = n_path,
           y = reorder_within(gene, n_path, cell_major),
           fill = (pos - neg))) +
  geom_col() +
  scale_y_reordered() +
  scale_fill_gradient2(low = "#2b8cbe", mid = "grey90", high = "#e34a33",
                       midpoint = 0, name = "(+) − (−)") +
  facet_wrap(~ cell_major, scales = "free_y", ncol = 2) +
  labs(x = "# pathology with sig.", y = NULL,
       title = "Recurrent mitochondrial genes across pathologies") +
  theme_bw(11) +
  theme(axis.text.y = element_text(size = 7))   


# alternative method

min_recur <- 2
topN_by_major <- c(
  Inh = 25, Exc = 15, Oli = 20, Ast = 15, Mic = 12, OPC = 10, Vas = 10, T_cells = 10
)

recur_plot <- recur %>%
  dplyr::filter(n_path >= min_recur) %>%
  dplyr::left_join(
    tibble::tibble(cell_major = names(topN_by_major),
                   topN = as.integer(topN_by_major)),
    by = "cell_major"
  ) %>%
  dplyr::group_by(cell_major) %>%
  dplyr::arrange(dplyr::desc(n_path), dplyr::desc(abs(pos - neg)), gene, .by_group = TRUE) %>%
  dplyr::slice_head(n = dplyr::first(topN)) %>%
  dplyr::ungroup()

ggplot(recur_plot,
       aes(x = n_path,
           y = reorder_within(gene, n_path, cell_major),
           fill = (pos - neg))) +
  geom_col() +
  scale_y_reordered() +
  scale_fill_gradient2(low = "#2b8cbe", mid = "grey90", high = "#e34a33",
                       midpoint = 0, name = "(+) − (−)") +
  facet_wrap(~ cell_major, scales = "free_y", ncol = 2) +
  labs(x = "# pathology with sig.", y = NULL,
       title = "Recurrent mitochondrial genes across pathologies") +
  theme_bw(11) +
  theme(axis.text.y = element_text(size = 7))


#=====================C======================#
# 3) Statistical Test for Gender Differences by Cell Type

# two columns
sex_hit <- dream_mt_sig %>%
  mutate(sex = tolower(sex)) %>%
  group_by(cell_major, gene, sex) %>%
  summarise(sig = any(TRUE), .groups="drop") %>%
  tidyr::pivot_wider(names_from = sex, values_from = sig, values_fill = FALSE)

## Fisher test per cell type
sex_diff <- sex_hit %>%
  group_by(cell_major) %>%
  summarise(
    both = sum(female & male),
    f_only = sum(female & !male),
    m_only = sum(!female & male),
    none = sum(!female & !male),
    .groups="drop"
  ) %>%
  rowwise() %>%
  mutate(p_fisher = fisher.test(matrix(c(both, f_only, m_only, none), nrow=2))$p.value,
         logOR = suppressWarnings(log((f_only+0.5)/(m_only+0.5)))) %>%
  ungroup()


ggplot(sex_diff, aes(x = cell_major, y = -log10(p_fisher), fill = logOR)) +
  geom_col(width=.7) +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  scale_fill_gradient2(low="#2b8cbe", mid="grey90", high="#e34a33", midpoint=0,
                       name="log OR (F vs M)") +
  labs(x="Cell type (major)", y="-log10 p (Fisher)",
       title="Sex difference in #sig. mt genes (per cell type)") +
  theme_bw(11) + theme(axis.text.x = element_text(angle=45, hjust=1))


totals <- dream_df %>%
  dplyr::filter(gene %in% mt_set_in_df) %>%
  dplyr::mutate(sex = tolower(sex)) %>%
  dplyr::distinct(gene, cell_major, sex) %>%
  dplyr::count(cell_major, sex, name = "n_total")

sig_counts <- dream_mt_sig %>%                
  dplyr::mutate(sex = tolower(sex)) %>%
  dplyr::distinct(gene, cell_major, sex) %>%
  dplyr::count(cell_major, sex, name = "n_sig")

tab <- totals %>%
  dplyr::left_join(sig_counts, by = c("cell_major","sex")) %>%
  tidyr::complete(cell_major, sex, fill = list(n_total = 0, n_sig = 0)) %>%
  tidyr::pivot_wider(names_from = sex, values_from = c(n_total, n_sig), values_fill = 0) %>%
  dplyr::mutate(
    a = n_sig_female,
    b = pmax(n_total_female - n_sig_female, 0),
    c = n_sig_male,
    d = pmax(n_total_male   - n_sig_male,   0),
    a2 = a + 0.5, b2 = b + 0.5, c2 = c + 0.5, d2 = d + 0.5,
    logOR  = log((a2 * d2) / (b2 * c2)),
    pval   = purrr::pmap_dbl(list(a,b,c,d),
                             ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))$p.value),
    nlog10p = -log10(pval)
  )

ggplot(tab, aes(x = cell_major, y = nlog10p, fill = logOR)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
  scale_fill_gradient2(low="#2b8cbe", mid="grey90", high="#e34a33",
                       midpoint = 0, name = "log OR (F vs M)") +
  labs(x = "Cell type (major)", y = "-log10 p (Fisher)",
       title = "Sex difference in #sig. mt genes (per cell type)") +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#========================================================================#
#========================================================================#
#===================== adding mt_pathway scores =========================#
#========================================================================#
#========================================================================#


mt_pathway <- read.csv("/Users/spectremac/Desktop/PFC427_adata_filtered_obs.mt_pathway.csv.gz")

summary(mt_pathway)

# pseudo-bulk

library(data.table)
library(dplyr)
library(tidyr)


# mt <- readRDS("mt_pathway.rds")  
mt <- mt_pathway  

mt <- mt %>%
  mutate(
    cell_major = factor(major.celltype.0,
                        levels = c("Ast","Exc","Inh","Mic","Oli","OPC","Vas","T_cells")),
    sex = ifelse(msex == 1, "male", "female"),   
    projid = as.character(projid)
  )

path_cols <- grep(
  "(MOOTHA|REACTOME_|KEGG_|HALLMARK_OXIDATIVE_PHOSPHORYLATION|WP_|GO(BP|CC|MF)_|GALLUZZI_|WONG_|BIOCARTA_)",
  names(mt), value = TRUE
)

#qc_cols <- c("pct_counts_mt","total_counts","doublet_score","batch")

num_qc   <- intersect(c("pct_counts_mt","total_counts","doublet_score"), names(mt))
cols_keep <- c("projid","cell_major","sex", path_cols, num_qc)

mt_small <- mt %>%
  dplyr::select(all_of(cols_keep)) %>%         
  dplyr::mutate(
    sex        = tolower(sex),
    projid     = as.factor(projid),           
    cell_major = as.factor(cell_major)
  )

pb <- mt_small %>%
  dplyr::group_by(projid, cell_major, sex) %>%
  dplyr::summarise(
    dplyr::across(all_of(path_cols), ~median(.x, na.rm = TRUE)),
    dplyr::across(all_of(num_qc),    ~median(.x, na.rm = TRUE)),
    .groups = "drop"
  )


library(lme4)
library(broom.mixed)
library(purrr)


# candidates
pheno_candidates <- c(
  "gpath","amyloid","nft","tangles","plaq_d","plaq_n",
  "cogn_global_lv","cognep_demog_slope","cogng_demog_slope",
  "cognpo_demog_slope","cognps_demog_slope","cognse_demog_slope",
  "cognwo_demog_slope","CR_score" 
)


present  <- intersect(pheno_candidates, names(mt))
missing  <- setdiff(pheno_candidates, present)
if (length(missing)) message("Missing phenotype columns in mt: ", paste(missing, collapse=", "))


pheno_donor <- mt %>%
  dplyr::select(projid, dplyr::any_of(present)) %>%
  dplyr::distinct(projid, .keep_all = TRUE)


pheno_cols <- intersect(pheno_candidates, names(mt))
stopifnot(length(pheno_cols) > 0)  # stop when no one match

pheno_donor <- mt %>%
  dplyr::distinct(projid, dplyr::across(dplyr::all_of(pheno_cols)), .keep_all = TRUE) %>%
  dplyr::select(projid, dplyr::all_of(pheno_cols))

# combine to pseudo-bulk
pb2 <- pb %>% left_join(pheno_donor, by = "projid")

# for each pathway × cell_major × sex run score ~ pheno + age + batch + pct_counts_mt + (1|projid)
run_assoc <- function(df, response, predictor){
  fml <- as.formula(paste(response, "~ scale(", predictor, ") + scale(pct_counts_mt) + (1|projid)"))
  fit <- suppressWarnings(lmer(fml, data = df, REML = FALSE))
  tidy(fit) %>% filter(term == paste0("scale(", predictor, ")")) %>%
    mutate(response = response, predictor = predictor)
}

#==============================================#

suppressPackageStartupMessages({
  library(lme4)
  library(broom)
  library(broom.mixed)
  library(dplyr)
  library(purrr)
})

.ok_subset <- function(df, predictor, min_n = 30) {
  if (nrow(df) < min_n) return(FALSE)
  x <- df[[predictor]]
  if (all(is.na(x))) return(FALSE)
  if (sd(x, na.rm = TRUE) == 0) return(FALSE)
  TRUE
}

# choose lm or lmer 
run_assoc <- function(df, response, predictor){
  df <- df %>% dplyr::select(projid, !!response := all_of(response),
                             !!predictor := all_of(predictor),
                             pct_counts_mt)
  
  if (!.ok_subset(df, predictor)) return(NULL)
  
  n_id  <- dplyr::n_distinct(df$projid)
  n_obs <- nrow(df)
  has_re <- n_obs > n_id              
  
  fml_base <- as.formula(
    paste(response, "~ scale(", predictor, ") + scale(pct_counts_mt)")
  )
  
  if (has_re) {
    fml <- as.formula(paste0(deparse(fml_base), " + (1|projid)"))
    fit <- suppressWarnings(lme4::lmer(fml, data = df, REML = FALSE))
    tab <- broom.mixed::tidy(fit, effects = "fixed")
  } else {
    fit <- stats::lm(fml_base, data = df)
    tab <- broom::tidy(fit)
  }
  
  term_target <- paste0("scale(", predictor, ")")
  out <- tab %>%
    dplyr::filter(term == term_target) %>%
    dplyr::mutate(
      response  = response,
      predictor = predictor,
      n_obs = n_obs,
      n_id  = n_id,
      used_RE = has_re
    )
  out
}


# ================== run all combinations ==================
results <- list()
for (cm in levels(pb2$cell_major)){
  for (sx in c("female","male")){
    df_sub <- pb2 %>% dplyr::filter(cell_major == cm, sex == sx)
    if (nrow(df_sub) < 30) next
    
    res_cm_sx <- purrr::map_dfr(path_cols, function(pw){
      purrr::map_dfr(pheno_cols, ~ run_assoc(df_sub, response = pw, predictor = .x)) %>%
        dplyr::mutate(cell_major = cm, sex = sx, pathway = pw)
    })
    
    results[[paste(cm, sx, sep="|")]] <- res_cm_sx
  }
}

results_df <- dplyr::bind_rows(results)


library(ggplot2)
library(dplyr)
library(purrr)

res_raw <- if (exists("results_df")) results_df else dplyr::bind_rows(results)

if (!"pathway" %in% names(res_raw) && "response" %in% names(res_raw)) {
  res_raw <- dplyr::rename(res_raw, pathway = response)
}

res_path <- res_raw %>%
  dplyr::mutate(
    sex  = tolower(.data$sex),
    beta = .data$estimate,
    se   = .data$std.error,
    z    = dplyr::if_else(!is.na(.data$statistic), .data$statistic, .data$beta / .data$se),
    p    = dplyr::if_else(!is.na(.data$p.value), .data$p.value,
                          2 * pnorm(abs(.data$z), lower.tail = FALSE))
  ) %>%
  dplyr::group_by(.data$cell_major, .data$sex, .data$predictor) %>%
  dplyr::mutate(FDR = p.adjust(.data$p, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::select(.data$cell_major, .data$sex, .data$predictor, .data$pathway,
                .data$beta, .data$se, .data$p, .data$FDR)

if (!"pathway" %in% names(res_path)) {
  if ("response" %in% names(res_path)) {
    res_path <- dplyr::rename(res_path, pathway = response)
  } else {
    stop("res_path 中找不到 'pathway' 列。当前列名：",
         paste(names(res_path), collapse = ", "))
  }
}


hp <- res_path %>%
  dplyr::filter(!is.na(pathway)) %>%
  dplyr::mutate(
    sgn  = sign(beta),
    stat = -log10(pmax(FDR, 1e-300))
  ) %>%
  dplyr::group_by(cell_major, sex, predictor, pathway) %>%
  dplyr::summarise(stat = max(stat, na.rm = TRUE) * dplyr::first(sgn),
                   .groups = "drop")



########### The plot is too dense ##########


keep_pat <- regex("(OXPHOS|OXIDATIVE|TCA|CITRIC|FATTY[_ ]?ACID|RESPIRATORY|MITOPHAGY|TRANSLATION|BIOGENESIS)",
                  ignore_case = TRUE)

p <- hp %>%
  dplyr::filter(stringr::str_detect(pathway, keep_pat)) %>%
  ggplot(aes(x = predictor, y = pathway, fill = stat)) +
  geom_tile(color = "white", size = .2) +
  facet_grid(cell_major ~ sex, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e34a33", midpoint = 0,
                       name = "±log10(FDR)") +
  labs(x = "Pathology", y = "Mito pathway", title = "Pathway activity vs pathology") +
  theme_bw(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

print(p)


########################## change to this ###############################


stopifnot(all(c("cell_major","sex","predictor","pathway","beta","FDR") %in% names(res_path)))

library(dplyr); library(stringr)

hp <- res_path %>%
  mutate(
    sex  = tolower(sex),
    stat = sign(beta) * -log10(pmax(FDR, 1e-300))   # -log10(FDR)
  )


hp_strength <- hp %>%
  group_by(cell_major, sex, pathway) %>%
  summarise(max_abs = max(abs(stat), na.rm = TRUE), .groups = "drop")


class(hp_strength); names(hp_strength); head(hp_strength)
conflicted::conflict_prefer("filter","dplyr")
conflicted::conflict_prefer("slice_max","dplyr")


N <- 25

top_path <- hp_strength %>%
  dplyr::filter(is.finite(.data$max_abs) & .data$max_abs >= 1.3) %>%
  dplyr::group_by(.data$cell_major, .data$sex) %>%
  dplyr::slice_max(order_by = .data$max_abs, n = N, with_ties = FALSE) %>%
  dplyr::ungroup()


hp_plot <- hp %>%
  semi_join(top_path, by = c("cell_major","sex","pathway")) %>%   
  left_join(top_path, by = c("cell_major","sex","pathway")) %>%   
  mutate(

    pathway_label = stringr::str_replace_all(pathway, "^(REACTOME|GOBP|GO|WP)_", ""),
    pathway_label = stringr::str_replace_all(pathway_label, "_", " "),
    pathway_label = stringr::str_trunc(pathway_label, 60),
    panel         = interaction(cell_major, sex, drop = TRUE),
    pathway_ord   = tidytext::reorder_within(pathway_label, max_abs, panel)
  )

ggplot(hp_plot, aes(x = predictor, y = pathway_ord, fill = stat)) +
  geom_tile(color = "white", size = 0.2) +
  facet_grid(cell_major ~ sex, scales = "free_y", space = "free_y") +
  scale_y_reordered() +
  scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e34a33",
                       midpoint = 0, name = "±log10(FDR)") +
  labs(x = "Pathology", y = "Mito pathway", title = "Pathway activity vs pathology") +
  theme_bw(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank(),
        strip.text  = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))



#=======================================================#
#=================== Ast/Oli/Vas =======================#
#=======================================================#

# Zoom the pathway heatmap to Ast / Oli / Vas (top-N per panel)

keep_cells <- c("Ast","Oli","Vas")
N <- 25  # top pathways per (cell_major, sex)

hp3 <- res_path %>%
  mutate(sex = tolower(sex),
         stat = sign(beta) * -log10(pmax(FDR, 1e-300))) %>%
  filter(cell_major %in% keep_cells)

hp_strength3 <- hp3 %>%
  group_by(cell_major, sex, pathway) %>%
  summarise(max_abs = max(abs(stat), na.rm = TRUE), .groups="drop")

top_path3 <- hp_strength3 %>%
  filter(is.finite(max_abs), max_abs >= 1.3) %>%
  group_by(cell_major, sex) %>%
  slice_max(max_abs, n = N, with_ties = FALSE) %>%
  ungroup()

hp_plot3 <- hp3 %>%
  semi_join(top_path3, by=c("cell_major","sex","pathway")) %>%
  left_join(top_path3, by=c("cell_major","sex","pathway")) %>%
  mutate(pathway_label = pathway |>
           gsub("^(REACTOME|GOBP|GO|WP)_","",x=_) |>
           gsub("_"," ",x=_))

ggplot(hp_plot3, aes(predictor, reorder(pathway_label, max_abs), fill = stat)) +
  geom_tile(color="white", size=.2) +
  facet_grid(cell_major ~ sex, scales="free_y", space="free_y") +
  scale_fill_gradient2(low="#2b8cbe", mid="white", high="#e34a33", midpoint=0,
                       name="±log10(FDR)") +
  labs(x="Pathology", y="Mito pathway", title="Top pathways (Ast/Oli/Vas)") +
  theme_bw(11) + theme(axis.text.x=element_text(angle=45,hjust=1),
                       panel.grid=element_blank())


# Sex contrast at the pathway level (Ast / Oli / Vas)
# How different are females vs males for each pathway?

keep_cells <- c("Ast","Oli","Vas")
N <- 30   

pw_sex <- hp3 %>%
  group_by(cell_major, pathway, sex) %>%
  summarise(
    m  = mean(stat, na.rm = TRUE),             # mean across pathologies
    sd = sd(stat,   na.rm = TRUE),
    n  = sum(is.finite(stat)),                 # #pathologies with data
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = c(m, sd, n)) %>%
  mutate(
    delta = m_female - m_male,                                 # F − M
    se    = sqrt((sd_female^2 / pmax(n_female, 1)) +
                   (sd_male^2   / pmax(n_male,   1))),           # avoid /0
    ci    = 1.96 * se
  )

top_pw <- pw_sex %>%
  filter(cell_major %in% keep_cells) %>%
  group_by(cell_major) %>%
  slice_max(order_by = abs(delta), n = N, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    label = pathway |>
      gsub("^(REACTOME|GOBP|GO|WP)_", "", x = _) |>
      gsub("_", " ", x = _) |>
      stringr::str_trunc(60),
    yord = tidytext::reorder_within(label, delta, cell_major)
  )

ggplot(top_pw, aes(x = delta, y = yord, fill = delta)) +
  geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = delta - ci, xmax = delta + ci), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  facet_wrap(~ cell_major, scales = "free_y", ncol = 1) +
  scale_y_reordered() +
  scale_fill_gradient2(low = "#2b8cbe", mid = "grey90", high = "#e34a33", midpoint = 0,
                       name = "Female − Male") +
  labs(x = "Female − Male (mean ±log10(FDR) with sign)", y = NULL,
       title = "Sex contrast per pathway (Ast/Oli/Vas)") +
  theme_bw(11) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))


# Formal sex × pathology interaction at the pathway level
# This tests whether the association between a pathology and a pathway score differs by sex

library(lme4)
library(broom)
library(broom.mixed)
library(dplyr)
library(purrr)

##################### too dense ##################

covars <- c("pct_counts_mt")
cov_str <- if (length(covars)) paste("+", paste0("scale(", covars, ")", collapse=" + ")) else ""

run_interact <- function(df, response, predictor) {
  # keep only needed cols & clean
  need <- c("projid","sex", response, predictor, covars)
  df <- df %>% dplyr::select(dplyr::any_of(need)) %>%
    dplyr::mutate(sex = factor(tolower(sex), levels = c("female","male"))) %>%
    dplyr::filter(is.finite(.data[[response]]), is.finite(.data[[predictor]]))
  
  # must have both sexes
  if (nlevels(droplevels(df$sex)) < 2) {
    return(tibble(estimate = NA_real_, std.error = NA_real_, p.value = NA_real_,
                  pathway = response, predictor = predictor))
  }
  
  # projid?
  has_repl <- any(dplyr::count(df, projid)$n > 1)
  
  if (has_repl) {
    fml <- as.formula(paste0(response, " ~ sex * scale(", predictor, ") ",
                             cov_str, " + (1|projid)"))
    fit <- suppressWarnings(lme4::lmer(fml, data = df, REML = FALSE))
    tt  <- broom.mixed::tidy(fit)
  } else {
    fml <- as.formula(paste0(response, " ~ sex * scale(", predictor, ") ", cov_str))
    fit <- stats::lm(fml, data = df)
    tt  <- broom::tidy(fit)
  }
  
  term_name <- paste0("sexmale:scale(", predictor, ")")
  out <- tt %>% dplyr::filter(term == term_name)
  
  # if lm doesn’t carry p.value broom version, compute 
  if (!"p.value" %in% names(out) && all(c("estimate","std.error") %in% names(out))) {
    out$p.value <- 2*pnorm(abs(out$estimate/out$std.error), lower.tail = FALSE)
  }
  
  out %>% dplyr::transmute(estimate, std.error = std.error, p.value,
                           pathway = response, predictor = predictor)
}

# run only on Ast/Oli/Vas
res_int <- purrr::map_dfr(keep_cells, function(cm) {
  df <- pb2 %>% dplyr::filter(cell_major == cm)
  purrr::map_dfr(unique(res_path$pathway), function(pw) {
    purrr::map_dfr(unique(res_path$predictor), ~ run_interact(df, pw, .x)) %>%
      dplyr::mutate(cell_major = cm)
  })
}) %>%
  dplyr::mutate(FDR = p.adjust(p.value, method = "BH"),
                stat = sign(estimate) * -log10(pmax(FDR, 1e-300)))

ggplot(res_int, aes(predictor, pathway, fill = stat)) +
  geom_tile(color = "white", size = .2) +
  facet_grid(cell_major ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e34a33", midpoint = 0,
                       name = "Sex×Pathology ±log10(FDR)") +
  labs(x = "Pathology", y = "Mito pathway",
       title = "Sex × pathology interaction (Ast/Oli/Vas)") +
  theme_bw(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid   = element_blank())

####################################


############ Change the plot #########

library(dplyr)
library(stringr)
library(tidytext)  

N_top <- 25         
thr   <- 1.2        

nrow(res_int_top)
table(hp_strength$cell_major, hp_strength$max_abs >= 1.3)

# strength per pathway within each panel
hp_strength <- res_int %>%
  filter(is.finite(stat)) %>%
  group_by(cell_major, pathway) %>%
  summarise(max_abs = max(abs(stat), na.rm = TRUE), .groups = "drop")

top_pw_raw <- hp_strength %>%
  group_by(cell_major) %>%
  slice_max(order_by = max_abs, n = N_top, with_ties = FALSE) %>%
  ungroup()

top_pw <- top_pw_raw %>% filter(max_abs >= thr)
if (nrow(top_pw) == 0) {
  message("No pathways pass threshold; falling back to top N without threshold.")
  top_pw <- top_pw_raw
}

res_int_top <- res_int %>%
  semi_join(top_pw, by = c("cell_major","pathway")) %>%
  left_join(select(top_pw, cell_major, pathway, max_abs), by = c("cell_major","pathway")) %>%
  filter(!is.na(cell_major), !is.na(pathway)) %>%
  mutate(
    pathway_label = pathway |>
      str_replace_all("^(REACTOME|GOBP|GO|WP)_","") |>
      str_replace_all("_"," ") |>
      str_trunc(55),
    y_ord = reorder_within(pathway_label, max_abs, cell_major)
  ) %>%
  droplevels()

# order predictors by overall strength
pred_order <- res_int_top %>%
  group_by(predictor) %>%
  summarise(s = max(abs(stat), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(s)) %>% pull(predictor)
res_int_top$predictor <- factor(res_int_top$predictor, pred_order)

# ---- replot ----
ggplot(res_int_top, aes(predictor, y_ord, fill = stat)) +
  geom_tile(color = "white", size = 0.2) +
  facet_grid(cell_major ~ ., scales = "free_y", space = "free_y") +
  scale_y_reordered() +
  scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e34a33",
                       midpoint = 0, name = "Sex×Pathology ±log10(FDR)",
                       limits = c(-4, 4), oob = scales::squish) +
  labs(x = "Pathology", y = "Mito pathway",
       title = "Sex × pathology interaction — top pathways") +
  theme_bw(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid  = element_blank(),
        strip.text  = element_text(face = "bold"))
############ Change the plot #########



# 4) Link pathway hits back to gene-level drivers (from dream_df)
# For each top pathway in Ast/Oli/Vas, 
# find mt genes whose sign(logFC) agrees with the pathway direction across pathologies

# pick the top K pathways per cell and per sex from hp_plot3
K <- 5
top_pw <- hp_plot3 %>%
  group_by(cell_major, sex) %>%
  slice_max(max_abs, n=K, with_ties=FALSE) %>%
  ungroup() %>% select(cell_major, sex, pathway)

# reduce dream_df to those cell types, mt genes, and same pathologies
dream_sub <- dream_df %>%
  filter(cell_major %in% keep_cells,
         gene %in% mt_set_in_df,
         path_short %in% unique(res_path$predictor)) %>%
  mutate(sex = tolower(sex))

# summarize gene direction vs pathway
leading <- top_pw %>%
  left_join(hp3, by=c("cell_major","sex","pathway")) %>%
  select(cell_major, sex, predictor, pathway, pw_stat = stat) %>%
  left_join(dream_sub, by=c("cell_major","sex","predictor"="path_short")) %>%
  group_by(cell_major, sex, pathway, gene) %>%
  summarise(agree = mean(sign(pw_stat) == sign(logFC), na.rm=TRUE),
            n = dplyr::n(),
            .groups="drop") %>%
  filter(n >= 3) %>%  
  group_by(cell_major, sex, pathway) %>%
  slice_max(agree, n = 15, with_ties = FALSE) %>%
  ungroup()


leading %>% arrange(cell_major, sex, pathway, desc(agree)) %>% print(n=100)




# visualize #
# Dot-heatmap 


leading2 <- leading %>%
  mutate(
    sex  = factor(tolower(sex), levels = c("female","male")),
    path_label = pathway |>
      str_replace_all("^(REACTOME|GOBP|GO|WP)_", "") |>
      str_replace_all("_", " ") |>
      str_trunc(50),
    gene_label = gene
  )

path_pat <- "(MITOPHAGY|FUSION|ATP|ELECTRON TRANSPORT|RNA 3 END|RNA DEGRADATION|PROTON MOTIVE FORCE)"
leading_focus <- leading2 %>%
  filter(str_detect(path_label, regex(path_pat, ignore_case = TRUE)))

topN <- 10 
top_dot <- leading_focus %>%
  group_by(cell_major, sex, path_label) %>%
  slice_max(order_by = agree, n = topN, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(panel = interaction(cell_major, sex, path_label, drop = TRUE),
         yord  = reorder_within(gene_label, agree, panel))

ggplot(top_dot,
       aes(x = path_label, y = yord, color = agree, size = n)) +
  geom_point(alpha = 0.9) +
  scale_y_reordered() +
  facet_grid(cell_major ~ sex, scales = "free_y", space = "free_y") +
  scale_color_viridis_c(end = 0.95, name = "Agreement") +
  scale_size(range = c(1.5, 5), name = "Evidence (n)") +
  labs(x = "Pathway", y = "Top genes",
       title = "Leading gene drivers per pathway (top 10 by agreement)") +
  theme_bw(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"))



#  Lollipop per pathway 


plot_one_path <- function(cell = "Oli",
                          path = "mitophagy",   
                          topN = 15,
                          ignore_case = TRUE) {
  
  stopifnot(all(c("cell_major","sex","pathway","gene","agree","n") %in% names(leading)))
  
  dat0 <- leading %>%
    mutate(
      sex  = factor(tolower(sex), c("female","male")),
      path_clean = pathway |>
        str_replace_all("^(REACTOME|GOBP|GO|WP)_", "") |>
        str_replace_all("_", " ")
    )
  
  hits <- dat0 %>%
    filter(cell_major == cell,
           str_detect(path_clean, regex(path, ignore_case = ignore_case)))
  
  if (nrow(hits) == 0) {
    message("No rows for cell = '", cell, "' and path pattern = '", path, "'.")
    message("Examples in this cell:\n",
            paste0(dat0 %>% filter(cell_major == cell) %>% distinct(path_clean) %>%
                     pull() %>% head(10), collapse = "\n"))
    return(invisible(NULL))
  }
  
  df <- hits %>%
    group_by(sex) %>% slice_max(agree, n = topN, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(yord = reorder_within(gene, agree, sex))
  
  ggplot(df, aes(x = agree, y = yord, color = sex)) +
    geom_segment(aes(x = 0, xend = agree, yend = yord), alpha = .35) +
    geom_point(size = 2) +
    scale_y_reordered() +
    facet_wrap(~ sex, nrow = 1, scales = "free_y") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "Agreement (proportion of predictors)",
         y = NULL,
         title = paste(cell, "—", unique(hits$path_clean)[1])) +
    theme_bw(11) + theme(panel.grid = element_blank())
}

plot_one_path(cell = "Oli", path = "mitophagy")
plot_one_path(cell = "Vas", path = "proton motive force")
plot_one_path(cell = "Ast", path = "^GOBP_MITOCHONDRIAL_RNA_3_END_PROCESSING$")




