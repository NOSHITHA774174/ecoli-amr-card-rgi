library(tidyverse)
library(pheatmap)
setwd("~/amr_project")

# 1. Which result file belongs to which strain
strain_info <- tribble(
  ~id,                                    ~strain,                          ~group,
  "GCF_000005845.2_ASM584v2_genomic",     "K-12 MG1655",                    "Reference",
  "GCF_000006665.1_ASM666v1_genomic",     "EDL933",                         "Reference",
  "GCF_000007445.1_ASM744v1_genomic",     "CFT073",                         "Reference",
  "GCF_000008865.2_ASM886v2_genomic",     "Sakai",                          "Reference",
  "GCF_000013265.1_ASM1326v1_genomic",    "UTI89",                          "Reference",
  "GCF_048403405.1_ASM4840340v1_genomic", "755 (urine, Taiwan)",            "Clinical",
  "GCF_054047845.1_ASM5404784v1_genomic", "405 (urine, Houston)",           "Clinical",
  "GCF_051619495.1_ASM5161949v1_genomic", "AH19B041 (blood, Melbourne)",    "Clinical",
  "GCF_024918915.1_ASM2491891v1_genomic", "p11B (blood, Houston, ST131)",   "Clinical"
)

# 2. Read every result file, keep reliable hits (Perfect/Strict)
read_one <- function(i) {
  f <- file.path("results", paste0(strain_info$id[i], ".txt"))
  read.delim(f, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(strain = strain_info$strain[i], group = strain_info$group[i])
}
amr <- map_dfr(seq_len(nrow(strain_info)), read_one) %>%
  filter(Cut_Off %in% c("Perfect", "Strict"))

# 3. Bar chart: genes per genome, colored by group
counts <- amr %>% count(strain, group)
print(counts)
p1 <- ggplot(counts, aes(reorder(strain, -n), n, fill = group)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3) +
  scale_fill_manual(values = c(Reference = "grey60", Clinical = "firebrick")) +
  labs(x = NULL, y = "Resistance genes detected", fill = NULL,
       title = "CARD/RGI resistance genes: reference vs clinical E. coli") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("plots/genes_per_genome_v2.png", p1, width = 9, height = 6, dpi = 300)

# 4. KEY RESULT: genes in clinical isolates that NO reference strain has
ref_genes <- amr %>% filter(group == "Reference") %>% pull(Best_Hit_ARO) %>% unique()
new_genes <- amr %>%
  filter(group == "Clinical", !Best_Hit_ARO %in% ref_genes) %>%
  group_by(Best_Hit_ARO) %>%
  summarise(isolates = n_distinct(strain),
            drug_class = first(`Drug Class`),
            mechanism = first(`Resistance Mechanism`), .groups = "drop") %>%
  arrange(desc(isolates))
cat("\nGenes found only in clinical isolates:", nrow(new_genes), "\n")
print(new_genes, n = Inf, width = Inf)
write.csv(new_genes, "results/clinical_only_genes.csv", row.names = FALSE)

# 5. Your three lab drugs: how many genes target each class, per strain
classes <- c("penicillin beta-lactam", "tetracycline", "aminoglycoside")
drug_tab <- map_dfr(classes, ~ amr %>%
                      filter(str_detect(`Drug Class`, fixed(.x))) %>%
                      count(strain) %>% mutate(drug_class = .x)) %>%
  complete(strain = strain_info$strain, drug_class = classes, fill = list(n = 0)) %>%
  pivot_wider(names_from = drug_class, values_from = n)
print(drug_tab, n = Inf, width = Inf)

# 6. Heatmap of only the genes that differ between strains
mat <- amr %>% distinct(strain, Best_Hit_ARO) %>% mutate(present = 1) %>%
  pivot_wider(names_from = strain, values_from = present, values_fill = 0) %>%
  column_to_rownames("Best_Hit_ARO") %>% as.matrix()
var_mat <- mat[rowSums(mat) < ncol(mat), , drop = FALSE]
ann <- data.frame(group = strain_info$group, row.names = strain_info$strain)
png("plots/variable_genes_heatmap.png", width = 2400, height = 3000, res = 300)
pheatmap(var_mat, color = c("white", "steelblue"), legend = FALSE,
         annotation_col = ann, fontsize_row = 6, cluster_cols = FALSE,
         main = "Genes that differ between strains")
dev.off()