# SCRIPT PARA EL ANÁLISIS GRÁFICO DE LOS RESULTADOS DE EGGNOG Y MOTUS

rm(list = ls())

# Importamos el archivo de OTUs hecho en Python y creamos una tabla con la información.
df_OTUs <- readr::read_tsv("/Users/Sonsoles/Desktop/TFG/B1T1/resultados_reales/archivos_phyloseq_eggnog/phyloseq_OTUs.tsv")
colnames(df_OTUs)[1] <- 'Index'
OTUs_table <- df_OTUs[,-1]
rownames(OTUs_table) <- df_OTUs$Index

# Importamos el archivo de taxas hecho en Python y creamos una tabla con la información.
df_taxa <- readr::read_tsv("/Users/Sonsoles/Desktop/TFG/B1T1/resultados_reales/archivos_phyloseq_eggnog/phyloseq_taxas.tsv")
colnames(df_taxa)[1] <- 'Index'
taxa_table <- df_taxa[,-1]
rownames(taxa_table) <- df_taxa$Index

# Importamos el archivo de muestras hecho en Python y creamos una tabla con la información.
df_samples <- readr::read_tsv("/Users/Sonsoles/Desktop/TFG/B1T1/resultados_reales/archivos_phyloseq_eggnog/phyloseq_samples.tsv")
colnames(df_samples)[1] <- 'Index'
samples_table <- df_samples[,-1]
rownames(samples_table) <- df_samples$Index


# Creamos una matriz de OTUs.
otumat <- data.matrix(OTUs_table)

# Creamos una matriz de taxa.
taxmat <- as.matrix(taxa_table)


# Procedemos a hacer un estudio de Phyloseq.
library("phyloseq")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
plot_bar(physeq, fill="Domain")
plot_bar(physeq, fill="Phylum")
plot_bar(physeq, fill="Class")
plot_bar(physeq, fill="Order")
plot_bar(physeq, fill="Family")
plot_bar(physeq, fill="Genus")
plot_bar(physeq, fill="Species")

# MERGE THE OTUs WITH THE SAME TAXONOMY
ps.domain = tax_glom(physeq, taxrank="Domain", NArm=FALSE)
plot_bar(ps.domain, fill="Domain")
ps.phylum = tax_glom(physeq, taxrank="Phylum", NArm=FALSE)
plot_bar(ps.phylum, fill="Phylum")
ps.class = tax_glom(physeq, taxrank="Class", NArm=FALSE)
plot_bar(ps.class, fill="Class")
ps.order = tax_glom(physeq, taxrank="Order", NArm=FALSE)
plot_bar(ps.order, fill="Order")
ps.family = tax_glom(physeq, taxrank="Family", NArm=FALSE)
plot_bar(ps.family, fill="Family")
ps.genus = tax_glom(physeq, taxrank="Genus", NArm=FALSE)
plot_bar(ps.genus, fill="Genus")
ps.species = tax_glom(physeq, taxrank="Species", NArm=FALSE)
plot_bar(ps.species, fill="Species")

# Gráficos de alfa y beta-diversity.
library("ggplot2")
plot_richness(physeq, measures=c("Observed","Chao1", "Shannon", "Simpson"))

