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
rownames(samples_table) <- colnames(OTUs_table)

# Creamos una variable con los tipos de muestra.
samples_type <- readr::read_tsv("/Users/Sonsoles/Desktop/TFG/B1T1/datos_reales/family_samples.tsv")
samples_table$Individual <- samples_type$Individual

# Creamos las matriz de OTUs y de taxa.
otumat <- data.matrix(OTUs_table)
taxmat <- as.matrix(taxa_table)

# Procedemos a hacer un estudio de Phyloseq.
library("phyloseq")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(samples_table)

sample_names(SAM) <- sample_names(OTU)

physeq = phyloseq(OTU, TAX, SAM)

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
plot_richness(physeq, x="Individual", measures=c("Observed","Chao1", "Shannon", "Simpson"))

# Gráficos de ordination
library("plyr")
wh0 = genefilter_sample(physeq, filterfun_sample(function(x) x > 5), A=0.5*nsamples(physeq))
GP1 = prune_taxa(wh0, physeq)
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))
phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)
human = get_variable(GP1, "Individual") 
sample_data(GP1)$human <- factor(human)
GP.ord <- ordinate(GP1, "NMDS", "bray")

p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
print(p1)
p1 + facet_wrap(~Phylum, 3)

p2 = plot_ordination(GP1, GP.ord, type="samples", color="Individual", shape="human") 
p2 + geom_polygon(aes(fill=Individual)) + geom_point(size=5) + ggtitle("samples")

p3 = plot_ordination(GP1, GP.ord, type="biplot", color="Individual", shape="Phylum", title="biplot")
# Some stuff to modify the automatic shape scale
GP1.shape.names = get_taxa_unique(GP1, "Class")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values=GP1.shape)

p4 = plot_ordination(GP1, GP.ord, type="split", color="Phylum", shape="human", label="Individual", title="split") 
p4
gg_color_hue <- function(n){
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
color.names <- levels(p4$data$Phylum)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols["samples"] <- "black"
  p4 + scale_color_manual(values=p4cols)
  





