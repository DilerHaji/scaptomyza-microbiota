qiime phylogeny raxml-rapid-bootstrap \
	--i-alignment duncan-aligned-filtered.qza \
	--p-seed 1234 \
	--p-rapid-bootstrap-seed 12345 \
	--p-bootstrap-replicates 100 \
	--p-n-threads 32 \
	--p-raxml-version SSE3 \
	--o-tree duncan-unrooted-tree.qza





setwd("/Users/dilerhaji/Desktop/Duncan_16S/R")
library(phyloseq)


#tree <- read_qza("../rooted-tree.qza")
#tree <- tree$data

taxonomy <- as.matrix(tax_gza_to_phyloseq("duncan-taxonomy.qza"))
table <- read_qza("../duncan-dada-table.qza")
table <- table$data
metadata <- read.csv("duncan_metadata.csv")
metadata <- sample_data(metadata)
rownames(metadata) <- metadata$sampleid
phy <- phyloseq(otu_table(table, taxa_are_rows = TRUE), tax_table(taxonomy), sample_data(metadata))
sample_data(phy)$depth_original <- as.numeric(sample_sums(phy))

## metadata stuff
[1] "sampleid"       "sample"         "conc"           "type"          
 [5] "dateCollected"  "dateExtracted"  "gDNAID"         "stage"         
 [9] "sex"            "locality"       "genus"          "species"       
[13] "depth_original"



#####################
exporting taxonomy table
#####################


phy_table <- dtable(phy)

write.csv(phy_table, "phy_table.csv")







#####################
Library size plot
#####################

# Library size plot across samples
df <- as.data.frame(sample_data(phy)) 
df$LibrarySize <- sample_sums(phy)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
lib_plot <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=type)) + geom_point() + 
	xlab("Order of library sizes") + 
	geom_hline(yintercept = 50000, col = "red") + 	
	geom_text_repel(data = df, aes(label = sampleid), size = 2)

ggsave("lib_plot.png", lib_plot, scale = 0.7)

phy <- subset_samples(phy, !sampleid %in% data.frame(df[df$LibrarySize < 50000, "sampleid"])$sampleid)






####################
Removing non-bacteria 
###################

# adding filter information to spreadsheet 

rownames(phy_table)
filtered_bacteria <- ifelse(rownames(phy_table) %in% taxa_names(subset_taxa(phy, !kingdom %in% c("Bacteria"))), yes = 1, no = 0)
filtered_mito <- ifelse(rownames(phy_table) %in% taxa_names(subset_taxa(phy, family %in% c("Mitochondria"))), yes = 1, no = 0)
filtered_chloro <- ifelse(rownames(phy_table) %in% taxa_names(subset_taxa(phy, order %in% c("Chloroplast"))), yes = 1, no = 0)
table(filtered_bacteria)
table(filtered_mito)
table(filtered_chloro)
filtered <- filtered_bacteria+filtered_mito+filtered_chloro
phy_table$filter_nonbacteria <- filtered
write.csv(phy_table, "phy_table.csv")


phy1 <- subset_taxa(phy, kingdom %in% c("Bacteria")) #186
phy1 <- subset_taxa(phy1, !family %in% c("Mitochondria"))	#191
phy1 <- subset_taxa(phy1, !order %in% c("Chloroplast"))	#74








####################
Ordination of everything after above filtering 
####################

ord <- ordinate(phy1, method = "PCoA", distance = "bray")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], species = sample_data(phy)$species, id = sample_data(phy)$sampleid)
ord1 <- ggplot(ordDF, aes(x = PC1, y = PC2, color = species)) + 
	geom_point(alpha = 0.8, size = 3) +
	geom_text_repel(data = ordDF, aes(label = id), size = 2) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Paired") +
	#stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("ord1.png", ord1, scale = 0.8)


ord <- ordinate(phy1, method = "PCoA", distance = "bray")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], stage = sample_data(phy)$stage, sex = sample_data(phy)$sex)
ord2 <- ggplot(ordDF, aes(x = PC1, y = PC2, color = stage, shape = sex)) + 
	geom_point(alpha = 0.8, size = 3) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Set1") +
	#stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("ord2.png", ord2, scale = 0.8)


ord <- ordinate(phy1, method = "PCoA", distance = "bray")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], date = sample_data(phy)$dateExtracted, conc = sample_data(phy)$conc)
ord3 <- ggplot(ordDF, aes(x = PC1, y = PC2, color = date, size = conc)) + 
	geom_point(alpha = 0.8) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Set1") +
	#stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("ord3.png", ord3, scale = 0.8)



















### Top most abundant ASVs and their taxonomic designations 


phy_bar <- phy1


taxorder <- taxa_sums(phy_bar)[order(-taxa_sums(phy_bar))]
taxname <- names(taxorder[1:20])

tax_table(phy_bar)[!rownames(tax_table(phy_bar)) %in% taxname, "genus"] <- NA
phy_bar <- relphy(phy_bar)



barall <- plot_bar2(phy_bar, fill = "genus") + facet_wrap(~stage, scales = "free") + scale_fill_brewer(palette = "Set1")

ggsave("phy_bar.png", barall, scale = 1)














### Filtering contaminants 

#### DECONTAM 
library(decontam)


## Freq based (probs better)

phy2 <- subset_samples(phy1, conc > 0)
phy2 <- subset_taxa(phy2, taxa_sums(phy) > 0)
freq <- isContaminant(phy2, method="frequency", conc="conc", threshold = 0.15)
table(freq$contaminant) #66
hist(freq$p, breaks = 100)
head(which(freq$contaminant))
plot_frequency(phy3, taxa_names(phy3)[c(which(freq$contaminant))], conc="conc")

decontam2 <- rownames(freq[which(freq$contaminant),])

tax_table(subset_taxa(phy3, taxa_names(phy3) %in% decontam2))

phy2 <- prune_taxa(!taxa_names(phy1) %in% decontam2, phy1) 







							  ## Prev based for B4
							  ## We performed a series of Decontam (R package) filtering steps. First, we used Decontam's prevalence based filtering to remove bacterial taxa that were prevalent across controls (threshold = 0.6) 
							  ## Based on the plot, increase prevelence of bacteria acroos gut samples showed a qualitative increase with prevalence across control samples, suggesting that contaminants ar widespread in gut samples. 


							  sample_data(phy2)$is.neg <- sample_data(phy2)$type == "Blank"
							  prev <- isContaminant(phy2, method="prevalence", neg="is.neg", threshold=0.1)
							  table(prev$contaminant)
							  ps.pa <- transform_sample_counts(phy2, function(abund) 1*(abund>0))
							  ps.pa.neg <- prune_samples(sample_data(ps.pa)$type == "Blank", ps.pa)
							  ps.pa.pos <- prune_samples(sample_data(ps.pa)$type == "Fly", ps.pa)
							  df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=prev$contaminant)
							  ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_jitter(alpha = 0.5) +
								xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")



							  decontam1 <- rownames(b4_prev[which(b4_prev$contaminant),])

							  phy1 <- prune_taxa(!taxa_names(phylo) %in% decontam1, phylo) #-30





#### Before and after filtering plots

phy_bar <- phy2


taxorder <- taxa_sums(phy_bar)[order(-taxa_sums(phy_bar))]
taxname <- names(taxorder[1:20])

tax_table(phy_bar)[!rownames(tax_table(phy_bar)) %in% taxname, "genus"] <- NA
phy_bar <- relphy(phy_bar)

barall <- plot_bar2(phy_bar, fill = "genus") + facet_wrap(~stage, scales = "free") + scale_fill_brewer(palette = "Set1")

ggsave("phy_bar2.png", barall, scale = 1)




#### After filtering pcoa plots

ord <- ordinate(phy2, method = "PCoA", distance = "bray")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], species = sample_data(phy)$species, id = sample_data(phy)$sampleid)
ord1 <- ggplot(ordDF, aes(x = PC1, y = PC2, color = species)) + 
	geom_point(alpha = 0.8, size = 3) +
	geom_text_repel(data = ordDF, aes(label = id), size = 2) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Paired") +
	#stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("ord1_2.png", ord1, scale = 0.8)






### Boxplots of pairwise differences : stage, species, sex


pdis <- phyloseq::distance(phy2, method="bray", type="samples")
pdis <- as.matrix(pdis)
pdis[lower.tri(pdis, diag = TRUE)] <- "diag"
sam <- sample_data(phy2)
colnames(pdis) <- rownames(pdis) <- paste(as.character(sam$sampleid), as.character(sam$species), as.character(sam$type), as.character(sam$stage), as.character(sam$sex), as.character(sam$locality), sep = "_")

ptemp <- data.frame(cbind(pdis, ID = rownames(pdis)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)

pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])
head(pdis2)
write.csv(pdis2, "pdis.csv")


dis <- pdis2[pdis2$pdis != "diag",]

dis$species_same <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 3)) == unlist(lapply(str_split(dis$id, "_"), "[", 3)), yes = 0, no = 1)
dis$species_same <- as.factor(dis$species_same)


ggplot()


as.numeric(as.factor(unlist(lapply(str_split(pdis2$ID, "_"), "[", 2))))

dis$samplid <- abs(as.numeric(as.factor(unlist(lapply(str_split(pdis2$ID, "_"), "[", 2)))) - as.numeric(as.factor(unlist(lapply(str_split(pdis2$id, "_"), "[", 2)))))
dis$type <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 4)) == unlist(lapply(str_split(dis$id, "_"), "[", 4)), yes = 0, no = 1)
dis$stage <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 5)) == unlist(lapply(str_split(dis$id, "_"), "[", 5)), yes = 0, no = 1)
dis$sex <- unlist(lapply(str_split(dis$ID, "_"), "[", 6))
dis$locality <- unlist(lapply(str_split(dis$id, "_"), "[", 6))












### Phylogenetic tree of most abundant taxa with related Genbank taxa included
# Do ASVs cluster in group suggesting that all individuals have one strain 
# Wolbachia tree 
