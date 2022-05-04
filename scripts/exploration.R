# 16S analysis

#### Loading the needed packages ####
library("phyloseq")
library("ggplot2")
library("edgeR")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")

#### Loading the biom files ####

## The object obtained by silva database
silva <- import_biom("medicago-sativa/miscellaneous/bioms/silva.biom")
# Prunning the data
silva@tax_table@.Data <- substring(silva@tax_table@.Data, 4)
colnames(silva@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


## The object obtained by pdp
rdp <- import_biom("medicago-sativa/miscellaneous/bioms/rdp.biom")
# Prunning the data
rdp@tax_table@.Data <- substring(rdp@tax_table@.Data, 4)
colnames(rdp@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

## The object obtained by greengenes
greeng <- import_biom("medicago-sativa/miscellaneous/bioms/greengenes.biom")
# Prunning the data
greeng@tax_table@.Data <- substring(greeng@tax_table@.Data, 4)
colnames(greeng@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#### Contrasting the taxonomic assignation ####
c.datab <- data.frame(row.names = c("Silva","RDP","Greeng"),
                      DataB = c("Silva","RDP","Greeng"),
                      OTUs = c(sum(sample_sums(silva)),sum(sample_sums(rdp)),
                               sum(sample_sums(greeng))),
                      NonBacteria = c(summary(silva@tax_table@.Data[,1] == "Bacteria")[2],
                                      summary(rdp@tax_table@.Data[,1] == "Bacteria")[2],
                                      summary(greeng@tax_table@.Data[,1] == "Bacteria")[2]),
                      Phyla = c(summary(unique(silva@tax_table@.Data[,2]))[1],
                                  summary(unique(rdp@tax_table@.Data[,2]))[1],
                                  summary(unique(greeng@tax_table@.Data[,2]))[1]),
                      Species= c("Na","Na","Yes"))

#Defining the order of the colors
spe.color <- palette.colors(n = 2, palette = "Dark2")
names(spe.color)<- c("Yes","Na")

#Plotting the number of OTUs
ggplot(data = c.datab, aes(y = OTUs, x = DataB, fill = Species))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values = spe.color)+
  theme_bw() + theme(text = element_text(size = 30))

#Plotting the non-Bacterial phyla
##Building matrix for heatmap
heat.fra <- data.frame(DataB = c(rep(x = "Silva", times = 5),rep(x = "RDP", times = 5),
                                 rep(x = "Greeng", times = 5)),
                       Phyla = rep(x = unique(silva@tax_table@.Data[,1]), times = 3),
                       Presence = c(1,1,1,1,1,
                                    1,0,0,0,1,
                                    1,0,0,0,1))


ggplot(data = heat.fra, mapping = aes(y= Phyla, x = DataB)) +
  geom_tile(aes(fill = Presence), colour = "grey", size = 2) +
  scale_fill_gradient(high = "#5ab4ac" ,low = "#000000" )+
  theme_bw() + theme(text = element_text(size = 30))

#### Normalizing th data ####

#### Trimming the data ####
edgeRnorm = function(physeq, ...) {
  require("edgeR")
  require("phyloseq")
  # physeq = simlist[['55000_1e-04']] z0 = simlisttmm[['55000_1e-04']] physeq
  # = simlist[['1000_0.2']] z0 = simlisttmm[['1000_0.2']] Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # See if adding a single observation, 1, everywhere (so not zeros) prevents
  # errors without needing to borrow and modify calcNormFactors (and its
  # dependent functions) It did. This fixed all problems.  Can the 1 be
  # reduced to something smaller and still work?
  x = x + 1
  # Now turn into a DGEList
  y = edgeR::DGEList(counts = x, remove.zeros = TRUE)
  # Perform edgeR-encoded normalization, using the specified method (...)
  z = edgeR::calcNormFactors(y, ...)
  # A check that we didn't divide by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  # Don't need the following additional steps, which are also built-in to some
  # of the downstream distance methods. z1 = estimateCommonDisp(z) z2 =
  # estimateTagwiseDisp(z1)
  return(z)
}
# Normalizing each of the three datasets

z.silva<- edgeRnorm(silva, method = "TMM")
silva.totu <- otu_table(z.silva[[1]], taxa_are_rows = TRUE)
n.silva <- merge_phyloseq(silva.totu, silva@tax_table)

z.rdp <- edgeRnorm(rdp, method = "TMM")
rdp.totu <- otu_table(z.rdp[[1]], taxa_are_rows = TRUE)
n.rdp<- merge_phyloseq(rdp.totu, rdp@tax_table)

z.greeng <- edgeRnorm(greeng, method = "TMM")
greeng.totu <- otu_table(z.silva[[1]], taxa_are_rows = TRUE)
n.greeng <- merge_phyloseq(greeng.totu, greeng@tax_table)

# Transforming the absolute counts into relative ones
r.silva <- transform_sample_counts(physeq = silva, function(x) x*100/sum(x))
r.rdp <- transform_sample_counts(physeq = rdp, function(x) x*100/sum(x))
r.greeng <- transform_sample_counts(physeq = greeng, function(x) x*100/sum(x))


#### Comparing the number of detected OTUs of different bacterial genera ####

# Looking for the majority taxa
maj.greeng <- filter_taxa(r.greeng, function(x) mean(x) > 1, TRUE)

# For Clavibacter
cla.silva <- subset_taxa(silva, Genus == "Clavibacter")
cla.rdp <- subset_taxa(rdp, Genus == "Clavibacter")
cla.greeng <- subset_taxa(greeng, Genus == "Clavibacter")
cla.greeng@tax_table@.Data[1,7] <- "NotIdentified"


#Dataframe of the Clavibacter OTUs
clavi <- data.frame(DataB = c("Silva","RDP","Greeng"),
                    OTUs = c(sum(sample_sums(cla.silva)),sum(sample_sums(cla.rdp)),
                             sum(sample_sums(cla.greeng))))

# Plotting
ggplot(data = clavi, aes(y = OTUs, x = DataB, fill= DataB))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(type = "",palette = "Set2", aesthetics = "fill")+
  theme_bw() + theme(text = element_text(size = 30))

#Heatmap of the Clavibacter reads found in Greengenes
plot_heatmap(cla.greeng, taxa.label = "Species") +
  theme_bw() + theme(text = element_text(size = 30))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

