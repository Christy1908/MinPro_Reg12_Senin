# ====================================================================
# SCRIPT TRANSKRIPTOMIK: THE SCAR-FREE LUNG (IPF vs Normal)
# Deskripsi: Analisis Ekspresi Diferensial dan Pengayaan Fungsional
# Dataset: GSE47460 (Platform Microarray Agilent)
# ====================================================================

# 0. Penyesuaian batas waktu (timeout) koneksi untuk mengunduh dataset berukuran besar
options(timeout = 3600)

# 1. Memuat paket GEOquery untuk akuisisi data dari database NCBI GEO
library(GEOquery)

# 2. Mengunduh dataset matriks ekspresi GSE47460 beserta anotasinya
gset <- getGEO("GSE47460", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# 3. Ekstraksi matriks ekspresi gen dan metadata fenotipe sampel
ex <- exprs(gset)
metadata <- pData(gset)

# 4. Standardisasi kolom fenotipe untuk pengelompokan kondisi klinis
metadata$Group <- metadata$characteristics_ch1

# 5. Filtrasi sampel: Eksklusi pasien COPD, mempertahankan kelompok Control dan Interstitial lung disease
metadata_bersih <- metadata[metadata$Group == "disease state: Control" | metadata$Group == "disease state: Interstitial lung disease", ]

# 6. Penamaan ulang label kelompok menjadi format biner (Normal vs Fibrosis) sebagai faktor
metadata_bersih$Group <- ifelse(metadata_bersih$Group == "disease state: Control", "Normal", "Fibrosis")
metadata_bersih$Group <- as.factor(metadata_bersih$Group)

# 7. Sinkronisasi matriks ekspresi dengan metadata sampel yang telah difiltrasi
ex_bersih <- ex[, rownames(metadata_bersih)]

# ====================================================================
# ANALISIS EKSPRESI DIFERENSIAL (DIFFERENTIAL EXPRESSION ANALYSIS)
# ====================================================================
library(limma)

# 1. Pembuatan matriks desain eksperimental
design <- model.matrix(~ 0 + metadata_bersih$Group)
colnames(design) <- c("Fibrosis", "Normal")

# 2. Penentuan matriks kontras (Fibrosis vs Normal)
contrast.matrix <- makeContrasts(Fibrosis - Normal, levels = design)

# 3. Pencocokan model linier dan komputasi statistik Bayes empiris
fit <- lmFit(ex_bersih, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 4. Ekstraksi hasil gen yang terekspresi diferensial (DEGs) tanpa batasan jumlah
hasil_deg <- topTable(fit2, adjust.method="fdr", sort.by="B", number=Inf)

# ====================================================================
# ANOTASI FITUR GEN (FEATURE ANNOTATION)
# ====================================================================
kamus_gen <- fData(gset)

# Penggabungan hasil statistik dengan anotasi probe dan pengurutan berdasarkan signifikansi (P-Value)
hasil_deg_lengkap <- merge(hasil_deg, kamus_gen, by=0)
hasil_deg_lengkap <- hasil_deg_lengkap[order(hasil_deg_lengkap$P.Value), ]

# ====================================================================
# VISUALISASI 1: VOLCANO PLOT
# ====================================================================
library(ggplot2)

# Kategorisasi status regulasi gen berdasarkan ambang batas logFC dan P-Value
hasil_deg_lengkap$Status <- "Not Significant"
hasil_deg_lengkap$Status[hasil_deg_lengkap$logFC > 1 & hasil_deg_lengkap$P.Value < 0.05] <- "Upregulated (Fibrosis)"
hasil_deg_lengkap$Status[hasil_deg_lengkap$logFC < -1 & hasil_deg_lengkap$P.Value < 0.05] <- "Downregulated (Normal)"

ggplot(hasil_deg_lengkap, aes(x = logFC, y = -log10(P.Value), color = Status)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Upregulated (Fibrosis)" = "red",
                                "Downregulated (Normal)" = "blue",
                                "Not Significant" = "grey80")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Differential Expression in IPF",
       subtitle = "Dataset GSE47460 (193 Fibrosis vs 91 Normal)",
       x = "Log2 Fold Change", y = "-Log10 P-Value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme(legend.position = "bottom")

# ====================================================================
# VISUALISASI 2: HEATMAP (Z-SCORE SCALING & HIERARCHICAL CLUSTERING)
# ====================================================================
library(pheatmap)
library(RColorBrewer)

# Seleksi 50 DEGs teratas untuk kejelasan visualisasi
top50_gen <- head(hasil_deg_lengkap, 50)
mat_heatmap <- ex_bersih[top50_gen$Row.names, ]

# Penanganan duplikasi Gene Symbol untuk penamaan baris
rownames(mat_heatmap) <- make.unique(as.character(top50_gen$GENE_SYMBOL))

# Pengaturan parameter anotasi kolom dan palet warna
annotation_col <- data.frame(Kondisi = metadata_bersih$Group)
rownames(annotation_col) <- colnames(mat_heatmap)
ann_colors <- list(Kondisi = c(Normal = "#000080", Fibrosis = "#B22222"))
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

pheatmap(mat_heatmap, 
         scale = "row",              
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         color = heatmap_colors,
         show_colnames = FALSE,      
         show_rownames = TRUE,       
         fontsize_row = 6,           
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         main = "Top 50 DEGs: IPF vs Normal")

# ====================================================================
# VISUALISASI 3: ANALISIS PENGAYAAN FUNGSIONAL (GO & KEGG)
# ====================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Ekstraksi Gene Symbol dari gen yang mengalami upregulasi signifikan
gen_merah <- hasil_deg_lengkap$GENE_SYMBOL[hasil_deg_lengkap$logFC > 1 & hasil_deg_lengkap$P.Value < 0.05]
gen_merah <- as.character(gen_merah)
gen_merah <- gen_merah[!is.na(gen_merah) & gen_merah != ""]

# Konversi Gene Symbol menjadi Entrez ID untuk kompatibilitas database
gen_entrez <- bitr(gen_merah, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Analisis Gene Ontology (Biological Process)
ego <- enrichGO(gene          = gen_entrez$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

barplot(ego, showCategory=10) + 
  ggtitle("GO Biological Process: Upregulated Genes in Fibrosis") + 
  theme_minimal()

# Analisis KEGG Pathway
ekegg <- enrichKEGG(gene         = gen_entrez$ENTREZID,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)

dotplot(ekegg, showCategory=10) + 
  ggtitle("KEGG Pathways Enrichment") + 
  theme_minimal()

# Ambil 10 teratas
top10_gen <- head(hasil_deg_lengkap, 10)

# Export ke CSV
write.csv(top10_gen, "Top_10_Gen_Target_Fibrosis.csv", row.names = FALSE)

getwd()

# Simpan data lengkap (ribuan gen) untuk lampiran
write.csv(hasil_deg_lengkap, "Hasil_Analisis_Transkriptomik_Lengkap.csv", row.names = FALSE)

# Simpan data Top 10 gen untuk tabel utama makalah
write.csv(top10_gen, "Top_10_Gen_Target_Fibrosis.csv", row.names = FALSE)