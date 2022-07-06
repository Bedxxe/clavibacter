#### Loading the needed packages ####
library("phyloseq")
library("ggplot2")
library("edgeR")
library("DESeq2")
library("RColorBrewer")
library("stringr")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("ggrepel")
library("grid")
library("gridExtra")
library("ggmap")
library("cowplot")
library("ggsn")
library("scatterpie")
library("svglite")

#### Creating a palettes of colors ####
# Since we can have as much as 20 OTUs in our data. we will create a vector
# with colors that can differentiate between each other. We will create a vector 
# of 28 objects
manual.colors <- c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                 'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                 'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                 "tan1", 'wheat4', '#DDAD4B', 'chartreuse', 
                 'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                 "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8")

# Also I want to create a palettle for the host plants that we will be using
diverse.colors <- c("#1b9e77","#d95f02","#7570b3","#e7298a",
                    "#66a61e","#e6ab02","#a6761d","#666666")

#### Loading the biom files ####
#Capsicum
c.choi <- import_biom("../capsicum/choi-2020/taxonomy/biom-files/choi-2020.biom")
c.misce <- import_biom("../capsicum/miscelaneous-capsicum/taxonomy/biom-files/miscelaneous-capsicum.biom")
c.new <- import_biom("../capsicum/newberry-2020/taxonomy/biom-files/newberry-2020.biom")

#Tuberosum
t.shi <- import_biom("../tuberosum/shi-2019/taxonomy/biom-files/shi-2019.biom")

#### Adjustments the biom data ####
#Capsicum
c.choi@tax_table@.Data <- substring(c.choi@tax_table@.Data, 4)
colnames(c.choi@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(c.choi@otu_table@.Data) <- substring(colnames(c.choi@otu_table@.Data), 3)

c.misce@tax_table@.Data <- substring(c.misce@tax_table@.Data, 4)
colnames(c.misce@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(c.misce@otu_table@.Data) <- substring(colnames(c.misce@otu_table@.Data), 3)

c.new@tax_table@.Data <- substring(c.new@tax_table@.Data, 4)
colnames(c.new@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(c.new@otu_table@.Data) <- substring(colnames(c.new@otu_table@.Data), 3)

#Tuberosum
t.shi@tax_table@.Data <- substring(t.shi@tax_table@.Data, 4)
colnames(t.shi@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(t.shi@otu_table@.Data) <- substring(colnames(t.shi@otu_table@.Data), 3)
#### Loading Metadata ####

#Capsicum
#Metadata from Capsicum/Choi-2020
mraw.choi <- read.csv2(file = "../capsicum/choi-2020/metadata/SraRunTable.txt", sep = ",")
##Metadata from Capsicum/Miscelaneous 
mraw.misce <- read.csv2(file = "../capsicum/miscelaneous-capsicum/metadata/SraRunTable.txt", sep = ",")
##Metadata from Capsicum/Newberry-2020
mraw.new <- read.csv2(file = "../capsicum/newberry-2020/metadata/SraRunTable.txt", sep = ",")
mraw.new <- mraw.new[! mraw.new$Run %in% c("SRR10527395"),]

#Tuberosum
##Metadata from Tuberosum/Newberry-2020
mraw.shi <- read.csv2(file = "../tuberosum/shi-2019/metadata/SraRunTable.txt", sep = ",")
mraw.shi <- mraw.shi[! mraw.shi$Run %in% c("SRR7448326","SRR7457712"),]

## Trimming the metadata
#Capsicum
met.choi <- data.frame(row.names = mraw.choi$Run, 
                       Plant = rep("Capsicum", times = 12),
                       Sample= mraw.choi$Run, 
                       Environment = rep("Field", times = 12),
                       Isolation= c(rep("Stem",4),rep("Root",4),rep("Stem",2),rep("Root",2)) ,
                       State= c("D","H","H","H","D","D","D","H","D","D","H","H"),
                       Latitude = rep(x = 25.811389, times = 12),
                       Longitude = rep(x= 106.523333, times = 12 ))
met.misce <- data.frame(row.names = mraw.misce$Run, 
                        Plant = rep("Capsicum", times = 4),
                        Sample = mraw.misce$Run,
                        Environment = c("Field","Field","Greenhouse","Field"),
                        Isolation = c(rep("Na", times = 4)),
                        State = c(rep("Na", times = 4)),
                        Latitude = c(rep("37.470000", times = 3),mraw.misce$geographic_location_.latitude.[4]),
                        Longitude = c(rep("127.970000", times = 3),mraw.misce$geographic_location_.longitude.[4]))
met.new <- data.frame(row.names = mraw.new$Run, 
                      Plant = rep("Capsicum", times = 3),
                      Sample = mraw.new$Run,
                      Environment = mraw.new$env_local_scale,
                      Isolation = mraw.new$env_medium,
                      State = rep("Na", times = length(mraw.new$Run)),
                      Latitude = c("32.601389","32.601389","33.20654") ,
                      Longitude = c("85.353611","85.353611","87.534607"))

#Tuberosum
met.shi <- data.frame(row.names = mraw.shi$Run, 
                      Plant = rep("Tuberosum", times = 18), 
                      Sample = mraw.shi$Run,
                      Environment = rep("Field", times = 18),
                      Isolation = rep("Soil", times = 18),
                      State = rep("D", times = 18),
                      Latitude = rep("34.2487", times = 18), 
                      Longitude = rep("119.8167", times = 18))


#### Merging the metadata with the phyloseq objects ####
#Capsicum
c.choi <- merge_phyloseq(c.choi, sample_data(met.choi))
c.misce <- merge_phyloseq(c.misce,sample_data(met.misce))
c.new <- merge_phyloseq(c.new,sample_data(met.new))
##Merging all the three individual object in one that contains all the Capsicum data
capsi <- merge_phyloseq(c.choi, c.misce, c.new)

#Tuberosum
tuber <- merge_phyloseq(t.shi, sample_data(met.shi))

##Merging all 4 complete objects in one
cbact <- merge_phyloseq(capsi,tuber)
##Concentrating the metadata in one object
all.meta<- rbind(met.choi, met.misce, met.new, met.shi)



#### Normalization of the data ####
#Let's explore how the depth of the samples is distributed

# We will create a vector to allocate the colors for each of the host plants that we have
plant.colors <- diverse.colors[1:length(unique(all.meta$Plant))]
names(plant.colors) <- unique(all.meta$Plant)

#Allocating the taxonomic information into a data.frame
dprof <- data.frame(Samples = colnames(cbact@otu_table@.Data),
                    Reads = sample_sums(cbact),
                    Plant = cbact@sam_data@.Data[[2]])
#Plotting the obtained data
ggplot(data = dprof, mapping = aes(x = Samples, y = Reads))+
  geom_bar(stat = "identity", aes( fill = Plant)) +
  scale_fill_manual(values = plant.colors)+
  theme(text = element_text(size = 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_hline(yintercept = mean(sample_sums(cbact)), color = "cyan3",
             size = 1.5, alpha = 0.6)
#ggsave(filename = paste("sampleDepth-03-01","png",sep = "."),
#       scale = 3,dpi = 600, path = "figures/")

#-- Defining the normalization method --#
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
z<- edgeRnorm(cbact, method = "TMM")
#-- Merging all the objects in the new normalized phyloseq object --#
nor.cb <- merge_phyloseq(otu_table(z@.Data[[1]], taxa_are_rows = TRUE),
                         tax_table(cbact@tax_table@.Data),
                         cbact@sam_data)
rm(z)

#### Obtaining only the Clavibacter data ####
clabac <- subset_taxa(nor.cb, Genus == "Clavibacter")
#cla.percent <- transform_sample_counts(clabac, function(x) x*100/ sum(x))
# Putting the data into a data.frame and trimming it
clavi.df <- psmelt(clabac)
clavi.df$Species[clavi.df$Species==""] <- "sp."

# Transformation to relative abundance
clavi.df$ClaRelative <- sample_sums(clabac)/sample_sums(cbact)
#Plotting this data
ggplot(data = clavi.df, mapping = aes(y= ClaRelative, x = Sample, fill = Species, color = Plant))+
  geom_bar(position = "stack", stat = "identity", size=1.5)+
  scale_color_manual(values = plant.colors)+
  ylab(label = "Clavibacter-RelAbundance")+
  theme(text = element_text(size= 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Beta diversity of the Clavi reads
clav.ord <- ordinate(physeq = clabac, method = "NMDS", distance = "bray")
plot_ordination(physeq = clabac, ordination = clav.ord,
                color = "Plant") +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Plant))+
  geom_point(size=4 , alpha = 0.5)+
  scale_fill_manual(values = plant.colors)+
  scale_color_manual(values = plant.colors)



#### Preparing the data for plotting ####
#First let's try to define from where to where the map will be located
## Coercing the coordinates to numeric
all.meta$Latitude <- as.numeric(all.meta$Latitude)
all.meta$Longitude <- as.numeric(all.meta$Longitude)
# Also for clavi.df
clavi.df$Latitude <- as.numeric(clavi.df$Latitude)
clavi.df$Longitude <- as.numeric(clavi.df$Longitude)

## Extracting the max and min value for all samples with a margin to allow the visualization of the data
samp.sites <- data.frame(longitude= c(min(all.meta$Longitude)-2,max(all.meta$Longitude)+2),
                         latitude = c(min(all.meta$Latitude)-1,max(all.meta$Latitude)+1))

# We will create a list were we will allocate the information of
# each one of the samples, the Species that were identified and 
# generate a new column where we will allocate the Genus and species 
# information together
a<- list() # An empty list
b<- NULL # An empty vector to locate the different Genus+Species diversity
for (i in 1:length(unique(clavi.df$Sample))) {
  #Obtaining the data to work with inside a list
  a[[i]]<- clavi.df[clavi.df$Sample==unique(clavi.df$Sample)[i],]
  # Pasting the name of the Genus and Species in one new column
  a[[i]]$OTUName <-paste(a[[i]]$Genus, a[[i]]$Species, sep = "-")
  # We will define a vector to allocate all the different OTUs in case that not all the samples have the same ones
  b <- c(b,unique(a[[i]]$OTUName))
}
rm(i)

#We will take colors from the manual.colors vector as many unique OTUs in our data
samp.colors <- rev(manual.colors)[1:(length(unique(b)))]
# We will name this colors with the names of the unique OTUs
names(samp.colors) <- unique(b)

# We will create a set of vectors to prepare the data to do the pie graphs.
# First, I will obtain the name, Longitude,Latitude, and plant-hots name
scatter.samples <- unique(clavi.df$Sample)
scatter.Longitude <- NULL
scatter.Latitude <- NULL
scatter.Plant <- NULL
for (i in 1:length(unique(clavi.df$Sample))) {
  scatter.Longitude[i] <- a[[i]]$Longitude[1]
  scatter.Latitude[i] <- a[[i]]$Latitude[1]
  scatter.Plant[i] <- a[[i]]$Plant[1]
}
rm(i)
#Let's put the 4 vector into a data.frame
scatter.data <- data.frame(scatter.samples, scatter.Longitude,scatter.Latitude, scatter.Plant)
colnames(scatter.data) <- c("Sample","Longitude","Latitude","Plant")

# Now, we need to put all the abundances of the different OTUs in a new data.frame
# I will create an empty data.frame
OTU.scatter <- data.frame(matrix(ncol = length(unique(b)), 
                                 nrow = length(unique(clavi.df$Sample))))
colnames(OTU.scatter) <- unique(b)

#And fill it with the information
for (i in 1:length(unique(b))) {
  O.Name <- NULL #Vector to put the abundance information of each OTU
  for (j in 1:length(unique(clavi.df$Sample))) {
    O.Name <- c(O.Name,a[[j]]$ClaRelative[a[[j]]$OTUName==unique(b)[i]])
  }
  OTU.scatter[i] <- O.Name #Filling each column
}
rm(i,j,O.Name,scatter.samples,
   scatter.Plant,scatter.Longitude,scatter.Latitude)

#Finally, we will bind the two data.frames
scatter.data <- cbind(scatter.data,OTU.scatter)
rm(OTU.scatter)

#### Plotting the data on a map and pie-graphs ####
# We have all the information we need to generate the plots

#As an example, we will plot the first sample on our "a" list
ggplot(data = a[[1]], mapping = aes(x = Sample, y = ClaRelative, fill= OTUName, color = Plant))+
  geom_bar(stat = "identity", width = 1, size= 1.5) +
  theme_bw()+
  scale_fill_manual(values = samp.colors) +
  scale_color_manual(values = plant.colors) +
  coord_polar("y", start=0)+
  ggtitle(paste(unique(a[[1]]$Plant),unique(a[[1]]$Sample), sep = " "))+
  theme(plot.title = element_text(size = 30),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
  
# We will get this plots for all the samples and save them inside a new set of directories
# Creating the folders that will contain the plots
dir.create("pies/")
dir.create("pies/labeled/")
dir.create("pies/just-pie")

# First let's create the labeled pie graphs
l.pies <- list()
for (i in 1:length(unique(clavi.df$Sample))) {
  l.pies[[i]] <- ggplot(data = a[[i]], mapping = aes(x = Sample, y = Abundance, fill= OTUName,color = Plant))+
    geom_bar(stat = "identity", width = 1,size= 0.6) +
    theme_bw()+
    scale_fill_manual(values = samp.colors) +
    scale_color_manual(values = plant.colors) +
    coord_polar("y", start=0)+
    ggtitle(paste(unique(a[[i]]$Plant),unique(a[[i]]$Sample), sep = " "))+
    theme(plot.title = element_text(size = 20),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  ggsave(filename = paste(unique(clavi.df$Sample)[i],"svg",sep = "."),dpi = 900,device = "svg",
         path = "pies/labeled/", width = 30, height = 20, units = "cm")
}
rm(i)

# Now we can do de ones without labels
n.pies <- list()
for (i in 1:length(unique(clavi.df$Sample))) {
  n.pies[[i]] <- ggplot(data = a[[i]], mapping = aes(x = Sample, y = Abundance, fill= OTUName,color = Plant))+
    geom_bar(stat = "identity", width = 1,size= 0.6) +
    theme_bw()+
    scale_fill_manual(values = samp.colors) +
    scale_color_manual(values = plant.colors) +
    coord_polar("y", start=0)+
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          legend.position="none")
  ggsave(filename = paste(unique(clavi.df$Sample)[i],"svg",sep = "."),dpi = 900, device = "svg", 
         path = "pies/just-pie/", width = 20, height = 20, units = "cm")
}
rm(i)

# Let's Let's put this pie-graphs over the map

# We will use the map_data function to load the information we need.
mi.world <- map_data(map = "world")

# Then we will use the geom_polygon function from ggplot2 to plot the map
ggplot(mi.world, mapping = aes(x = long, y = lat, group = group))+
  geom_polygon(fill= "gray",colour = "black")+
  xlab("Longitude") + ylab("Latitude") + theme_bw()+ 
  coord_fixed(xlim = samp.sites$longitude,
              ylim = samp.sites$latitude,
              expand = TRUE)
ggsave(filename = paste("03-05-nakedMap","png",sep = "."),dpi = 900, 
       path = "figures/", width = 30, height = 20, units = "cm")

#Finally, we will use the scatterpie package to put the pies over this map

ggplot(mi.world, mapping = aes(x = long, y = lat, group = group))+
  geom_polygon(fill= "gray",colour = "black")+
  xlab("Longitude") + ylab("Latitude") + theme_bw()+ 
  coord_fixed(xlim = samp.sites$longitude,
              ylim = samp.sites$latitude,
              expand = TRUE) +
  geom_scatterpie(data = scatter.data, cols = unique(b), 
                  pie_scale = 0.6, size = 0.3,
                  aes(x= Longitude , y = Latitude, color = Plant))+
  scale_fill_manual(values = samp.colors) +
  scale_color_manual(aesthetics = c("color"),values = plant.colors)
ggsave(filename = paste("03-06-MapWPies","png",sep = "."),dpi = 1200, 
       path = "figures/", width = 30, height = 15, units = "cm")


