# install.packages("BiocManager")
# BiocManager::install("phyloseq", force = TRUE)
# BiocManager::install("Biostrings")
# BiocManager::install("metagenomeSeq", force = TRUE)
# install.packages("tidyverse")
# installed.packages("ggplot2")
# installed.packages("ggpubr")
# BiocManager::install("decontam", force = TRUE)
# BiocManager::install("indicspecies", force = TRUE)
# BiocManager::install("ggpubr", force = TRUE)
# BiocManager::install("DESeq2" = TRUE)
# install.packages("tidyverse")
# install.packages("dunn.test")
# install.packages("writexl")
# install.packages("pairwiseAdonis", force = TRUE)
# install.packages("remotes")
# remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# install.packages("cluster")
# install.packages("devtools")
# install.packages("BiocStyle")
# install.packages("DECIPHER")
# install.packages("phangorn")
# install.packages("knitr")


## load packages
library(knitr)
library(BiocStyle)
library(dada2)
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(ggpubr)
library(decontam)
library(metagenomeSeq)
library(indicspecies)
library(dada2); packageVersion("dada2")
library(DESeq2)
library(ggrepel)
library(pairwiseAdonis)
library(nlme)
library(dunn.test)
library(writexl)

getwd()
#change path to file directory
path <- "zr15761.rawdata.240210"
list.files(path)

# Forward and reverse fastq filenames have format: ID_SAMPLENAME_R1_001.fastq.gz and ID_SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, 
                        pattern="_R1.fastq.gz", 
                        full.names = TRUE))

fnFs

fnRs <- sort(list.files(path, 
                        pattern="_R2.fastq.gz", 
                        full.names = TRUE))

fnRs
# Extract sample names, assuming filenames have format: ID_SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

sample.names

#Plot quality
profile_fnFs <- dada2::plotQualityProfile(fnFs)
profile_fnRs <- dada2::plotQualityProfile(fnRs)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Quality trim. Trim based on quality and trim primers if still present in reads
out <- filterAndTrim(fnFs,
                     filtFs,
                     fnRs,
                     filtRs,
                     truncLen = c(320, 180),
                     trimLeft = c(16, 24),
                     maxN = 0,
                     maxEE= c(2, 2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = FALSE,
                     matchIDs = TRUE) # On Windows set multithread = FALSE

saveRDS(out, "out.rds")
out <- readRDS("out.rds")
dim(out)
out

# learning error rate
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errR, nominalQ = TRUE)

#dereplication
derepFs <- derepFastq(filtFs, verbose = FALSE)
derepRs <- derepFastq(filtRs, verbose = FALSE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
dadaFs[[1]] # show the first component 

#merge paired reads
mergers <- mergePairs(dadaFs,
                      derepFs,
                      dadaRs,
                      derepRs,
                      verbose = FALSE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# construct asv table
seqtab <- dada2::makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(dada2::getSequences(seqtab)))

#remove chimera
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab,
                                           method = "consensus",
                                           multithread = TRUE,
                                           verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#track the pipeline
getN <- function(x) {sum(dada2::getUniques(x))}

track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN),
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

track

#assign taxonomy
taxa <- dada2::assignTaxonomy(seqtab.nochim,
                              "4587955/silva_nr99_v138.1_train_set.fa.gz",
                              multithread = TRUE)

taxa <- dada2::addSpecies(taxa,
                          "4587955/silva_species_assignment_v138.1.fa.gz")

#inspect taxonomic assignment
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#construct phylogenetic tree
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = FALSE)

## using the package phangorn to build the tree using maximum likelihood
phangAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit <- pml(treeNJ, data = phangAlign)
fitGTR <- update(fit, k = 4, inv = 0.2)


fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload = TRUE)


### Handing off to phyloseq
## load meta data - make sure metadeta is well formated and sample names matches rownames of reads
samdf <- read.csv("metadata.csv", header = TRUE) # set header true if your file has an header

samples.out <- rownames(seqtab.nochim)

rownames(samdf) <- samples.out

#load trimmed sequence, tree file, taxa, and metadata into phyloseq
ps <- phyloseq(otu_table(seqtab.nochim, 
                         taxa_are_rows = FALSE),
               sample_data(samdf),
               tax_table(taxa),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

ps

# removing chloroplast or taxa not assigned at the domain level
physeq.clean <- ps

# a good general rule of thumb is samples below 1000 reads could be eliminated, although this isn't a hard rule, you can remove at 10,000 or more if you want
# Remove samples with less than 1000 reads
physeq.no.chloro <- prune_samples(sample_sums(physeq.clean) >= 1000, physeq.clean)

physeq.no.chloro

plot_richness(physeq.no.chloro, x="treatment", measures=c("Shannon", "Simpson"), color="treatment")

# Now lets look at the read distribution per sample and decide if we need to get rid of some samples because of low sequence depth
sort(sample_sums(physeq.no.chloro), decreasing = T) # read distribution

# how many total reads are we now working with?
sum(sample_sums(physeq.no.chloro))

# What is our mean and median read depth per sample? 
mean(sample_sums(physeq.no.chloro))

median(sample_sums(physeq.no.chloro))

#Lets make a histogram of the read distribution and put the median read depth
read.depths <- data.frame(sample_sums(physeq.no.chloro))

colnames(read.depths) <- "read.depth"

read.depth.plot <- ggplot(read.depths, aes(read.depth)) +
  geom_histogram(fill = "lightblue", 
                 # bins = 31,
                 color = "black") + 
  geom_vline(xintercept = median(sample_sums(physeq.no.chloro)), linetype = "dashed") + 
  theme_classic(20) + 
  theme(plot.margin = margin(t= 30, b = 30, l = 30, r = 30)) +
  xlab("Read Depth")

read.depth.plot

# Rarefaction analysis 
sam.data <- data.frame(physeq.no.chloro@sam_data)

physeq.no.chloro@otu_table

bOTU.table <- otu_table(physeq.no.chloro) %>%
  as.data.frame() %>%
  as.matrix()

raremax <- min(rowSums(bOTU.table)) ## without the "t" tranform
rare.fun <- rarecurve((bOTU.table), step = 1000, sample = raremax, tidy = T) ## without the "t" tranform

sam.data$sample_id <- row.names(sam.data)

bac.rare.curve.extract2 <- left_join(sam.data, 
                                     rare.fun, 
                                     by = c("sample_id" = "Site"))

bac.rare.curve.extract2$treatment_type <- paste0(bac.rare.curve.extract2$Treatment, "-", bac.rare.curve.extract2$Sample.Type)

bac.rare <- ggplot(bac.rare.curve.extract2,
                   aes(x = Sample,
                       y = Species,
                       group = Sample.ID,
                       color = treatment_type)) +
  # geom_point() +
  geom_line(linewidth = 1) +
  xlab("Reads") +
  ylab("Number of OTUs") +
  ggtitle("") +
  theme_classic(base_size = 20) +
  geom_vline(xintercept = median(sample_sums(physeq.no.chloro)),
             linetype = "dashed") +
  scale_color_discrete(name = "",
                       breaks = c("Basal-Fecal",
                                  "Basal-Feed",
                                  "A20-Fecal",
                                  "A20-Feed",
                                  "D2-Fecal",
                                  "D2-Feed",
                                  "RAS-Water"),
                       labels = c(
                         "Basal-Fecal",
                         "Basal-Feed",
                         "CFP-20-Fecal",
                         "CFP-20-Feed",
                         "FSC-2-Fecal",
                         "FSC-2-Feed",
                         "RAS-Water")) +
  
  theme(legend.key.spacing.y = unit(1, "cm"),
        legend.key.width = unit(2, "cm"))

# plot
bac.rare


# Normalize Sampling reads based on cumulative sum scaling (CSS normalization) for distance analysis
physeq.prop <- transform_sample_counts(physeq.no.chloro, function(x) x/sum(x))

MGS <- phyloseq_to_metagenomeSeq(physeq.no.chloro)

p <- metagenomeSeq::cumNormStatFast(MGS)

MGS <- metagenomeSeq::cumNorm(MGS, p = p)

metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample

norm.bacteria <- metagenomeSeq::MRcounts(MGS, norm = T)

norm.bacteria.OTU <- phyloseq::otu_table(norm.bacteria, taxa_are_rows = TRUE)

physeq.css <- phyloseq(otu_table(norm.bacteria.OTU, taxa_are_rows = FALSE),
                       sample_data(samdf),
                       tax_table(taxa),phy_tree(fitGTR$tree))

physeq.css

# Most abundant OTUs
most.abundant <- data.frame(sort(rowSums(physeq.no.chloro@otu_table), decreasing = TRUE))
most.abundant

# Now Lets Analyse
# ALPHA DIVERSITY
physeq.no.chloro@sam_data$shannon <- estimate_richness(physeq.no.chloro, measures = c("Shannon"))$Shannon

physeq.no.chloro@sam_data$invsimpson <- estimate_richness(physeq.no.chloro, measures = c("InvSimpson"))$InvSimpson

physeq.no.chloro@sam_data$richness <- estimate_richness(physeq.no.chloro, measures = c("Observed"))$Observed

physeq.no.chloro@sam_data$even <- physeq.no.chloro@sam_data$shannon/log(physeq.no.chloro@sam_data$richness)

sample.data.bac <- data.frame(physeq.no.chloro@sam_data)
sample.data.bac$Treatment <- as.factor(sample.data.bac$Treatment)

# Richness over treatment
sample.data.bac$treatment_type <- paste0(sample.data.bac$Treatment, "-", sample.data.bac$Sample.Type)



library(ggsci)
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
grep(pattern = "Feed|RAS",
     x = sample.data.bac$treatment_type_tank,
     value = TRUE) -> kolay

sample.data.bac_ok <- sample.data.bac[!(sample.data.bac$treatment_type_tank %in% kolay), ]


sample.data.bac_ok$Treatment <- factor(sample.data.bac_ok$Treatment,
                                       levels = c("Basal", "A20", "D2"))

richness.treatment <- ggplot(sample.data.bac_ok,
                             aes(x = Treatment, 
                                 y = richness,
                                 fill = Treatment)) + 
  geom_boxplot(outlier.color = NA,
               linewidth = 1,
               show.legend = FALSE) +
  
  scale_x_discrete(name = "Treatment",
                   labels = c("Basal",
                              "CFP-20",
                              "FSC-2")) +
  
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1300)) +
  
  geom_jitter(size = 3,
              # position = "jitter",
              width = 0.1,
              color = "#fe019a",
              show.legend = FALSE) +
  
  ylab("Richness") +
  
  # scale_fill_aaas() +
  
  scale_fill_d3(guide = "none") +
  
  # stat_compare_means(method = "anova",
  #                    size = 9,
  #                    vjust = 3) +
  
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               aes(color = "Mean"),
               size = 5,
               fill = "cyan",
               # color = "cyan",
               show.legend = TRUE) +
  
  scale_color_manual(name = "",
                     
                     values = "blue") +
  
  xlab("Treatment") +
  
  # facet_wrap(~ Treatment) +
  
  theme_bw() +
  
  theme(axis.text.x = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5),
        
        axis.text.y = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5)) +
  
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 1, 
                                                    r = 0,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(axis.title.y = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 0, 
                                                    r = 1,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "cm")) +
  
  theme(legend.position = "right") +
  
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    size = 1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  
  
  # theme(axis.line = element_line(color = "black",
  #                                linewidth = 1)) +
  
  theme(axis.ticks = element_line(color = "black",
                                  linewidth = 1)) +
  
  theme(axis.ticks.length = unit(0.25, "cm")) +
  
  theme(legend.position = c(0.95, 0.95),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box = "horizontal",
        legend.justification = c("right", "top")) +
  
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15)) 

richness.treatment

# Anova comparison
rich_trt_ok <- aov(richness ~ Treatment, 
                   data = sample.data.bac)

summary(rich_trt_ok)
TukeyHSD(rich_trt_ok)

library(agricolae)
HSD.test(rich_trt_ok, "treatment_type", console = TRUE)

shapiro.test(rich_trt_ok$residuals)

# bartlett.test(richness ~ treatment_type, 
#               data = sample.data.bac)

anova_table <- anova(rich_trt_ok)
print(anova_table)

tukey_results <- TukeyHSD(rich_trt_ok)

print(tukey_results)


# Shannon diversity by treatment
shan.treatment <- ggplot(sample.data.bac_ok,
                         aes(x = Treatment, 
                             y = shannon,
                             fill = Treatment)) + 
  geom_boxplot(outlier.color = NA,
               linewidth = 1,
               show.legend = FALSE) +
  
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 6)) +
  
  scale_x_discrete(name = "Treatment",
                   labels = c("Basal",
                              "CFP-20",
                              "FSC-2")) +
  
  geom_jitter(size = 3,
              # position = "jitter",
              width = 0.1,
              color = "#fe019a",
              show.legend = FALSE) +
  
  ylab("Shannon") +
  
  # scale_fill_aaas() +
  
  scale_fill_d3(guide = "none") +
  
  # stat_compare_means(method = "anova",
  #                    size = 9,
  #                    vjust = 3) +
  
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               aes(color = "Mean"),
               size = 5,
               fill = "cyan",
               # color = "cyan",
               show.legend = TRUE) +
  
  scale_color_manual(name = "",
                     values = "blue") +
  
  xlab("Treatment") +
  
  # facet_wrap(~ Treatment) +
  
  theme_bw() +
  
  theme(axis.text.x = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5),
        
        axis.text.y = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5)) +
  
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 1, 
                                                    r = 0,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(axis.title.y = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 0, 
                                                    r = 1,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "cm")) +
  
  theme(legend.position = "right") +
  
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    size = 1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  
  
  # theme(axis.line = element_line(color = "black",
  #                                linewidth = 1)) +
  
  theme(axis.ticks = element_line(color = "black",
                                  linewidth = 1)) +
  
  theme(axis.ticks.length = unit(0.25, "cm")) +
  
  theme(legend.position = c(0.95, 0.95),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box = "horizontal",
        legend.justification = c("right", "top")) +
  
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15)) 

shan.treatment

# Anova comparison
shan_trt <- aov(shannon ~ Treatment, 
                data = sample.data.bac)

anova_table <- anova(shan_trt)

print(anova_table)

tukey_results <- TukeyHSD(shan_trt)

# evenness diversity by treatment
even.treatment <- ggplot(sample.data.bac_ok,
                         aes(x = Treatment, 
                             y = even,
                             fill = Treatment)) + 
  geom_boxplot(outlier.color = NA,
               linewidth = 1,
               show.legend = FALSE) +
  
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  
  scale_x_discrete(name = "Treatment",
                   labels = c("Basal",
                              "CFP-20",
                              "FSC-2")) +
  
  geom_jitter(size = 3,
              # position = "jitter",
              width = 0.1,
              color = "#fe019a",
              show.legend = FALSE) +
  
  ylab("Evenness") +
  
  # scale_fill_aaas() +
  
  scale_fill_d3(guide = "none") +
  
  # stat_compare_means(method = "anova",
  #                    size = 9,
  #                    vjust = 3) +
  
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               aes(color = "Mean"),
               size = 5,
               fill = "cyan",
               # color = "cyan",
               show.legend = TRUE) +
  
  scale_color_manual(name = "",
                     values = "blue") +
  
  xlab("Treatment") +
  
  # facet_wrap(~ Treatment) +
  
  theme_bw() +
  
  theme(axis.text.x = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5),
        
        axis.text.y = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5)) +
  
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 1, 
                                                    r = 0,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(axis.title.y = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 0, 
                                                    r = 1,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "cm")) +
  
  theme(legend.position = "right") +
  
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    size = 1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  
  
  # theme(axis.line = element_line(color = "black",
  #                                linewidth = 1)) +
  
  theme(axis.ticks = element_line(color = "black",
                                  linewidth = 1)) +
  
  theme(axis.ticks.length = unit(0.25, "cm")) +
  
  theme(legend.position = c(0.95, 0.95),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box = "horizontal",
        legend.justification = c("right", "top")) +
  
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15)) 

even.treatment

#ANOVA comparison
even_trt <- aov(even ~ Treatment, 
                data = sample.data.bac)

anova_table <- anova(even_trt)

print(anova_table)

tukey_results <- TukeyHSD(even_trt)

#BETA DIVERSITY
# Principle coordinates analysis with Bray-Curtis distances
ordination.pcoa <- ordinate(physeq.css,
                            "PCoA", 
                            "bray") # calculate the resemblance and ordinate using PCoA
ordination.pcoa$vectors # positions of your points on the pcoa graph
ordination.pcoa$values #values to calculate the variance explained on each axis (dimension)

pcoa <- plot_ordination(physeq.css, 
                        ordination = ordination.pcoa,
                        axes = 1:2,
                        type = "Sample.Type", 
                        color = "treatment_type",
                        shape = "Treatment",
                        justDF = FALSE) +
  theme_bw() 
# scale_color_manual(values = cbbPalette)

pcoa + geom_point(size = 3) + 
  stat_ellipse(geom = "polygon",
               aes(fill = treatment_type),
               alpha = 0.2) +
  xlim(-1, 1) + ylim(-1, 1)


pcoa

pcoa$data -> pcoa_bray
pcoa_fecal <- pcoa_bray |> subset(Sample.Type == "Fecal")

## PLOT

pcoa_fecal$Treatment <- factor(pcoa_fecal$Treatment,
                               levels = c("Basal", "A20", "D2"),
                               labels = c("Basal", "CFP-20",
                                          "FSC-2"))

pcoa_fecal_bray <- ggplot(data = pcoa_fecal, 
                          mapping = aes(x = Axis.1,
                                        y = Axis.2,
                                        color = Treatment,
                                        shape = Treatment)) +
  
  geom_polygon(stat = "ellipse", 
               level = 0.95,
               aes(fill = Treatment),
               alpha = 0.1,
               show.legend = F) +
  
  
  
  geom_point(size = 6, 
             show.legend = F,
             alpha = 1) +
  
  # scale_color_discrete(name = "Treatment") +
  
  
  scale_color_manual(name = "Treatment",
                     values = c("red", "darkgreen", "blue")) +
  
  scale_shape_manual(values = c(17, 15, 19)) +
  
  # scale_shape_discrete(name = "Treatment") +
  
  scale_x_continuous(name = paste0("PCoA 1 (", 
                                   round(ordination.pcoa$values[1, 2] * 100, 1),
                                   "%)"),
                     limits = c(-1, 1)
  ) +
  
  scale_y_continuous(name = paste0("PCoA 2 (", 
                                   round(ordination.pcoa$values[2, 2] * 100, 1),
                                   "%)"),
                     limits = c(-1, 1)
  ) +
  
  # ggtitle(label = "PCoA beta diversity for fish fecal (Bray-Curtis)") +
  
  theme_bw() +
  
  theme(axis.ticks.length = unit(0.25, "cm")) +
  
  theme(axis.text = element_text(colour = "black",
                                 face = "bold")) +
  
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, "cm"))  +
  
  theme(axis.text.x = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5),
        
        axis.text.y = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5)) +
  
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 1, 
                                                    r = 0,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(axis.title.y = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 0, 
                                                    r = 1,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(
    # legend.title = element_text(size = 20, face = "bold", margin = margin(t = 0, 
    #                                                 r = 0,
    #                                                 b = 0.5, 
    #                                                 l = 0, 
    #                                                 "cm")),
    legend.title=element_blank(),
    legend.title.position = "left",
    legend.text = element_text(size = 20, face = "plain"),
    legend.key.spacing.x = unit(0.3, "cm"),
    legend.key.spacing.y = unit(0.3, "cm"),
    legend.key.size = unit(1, "cm")) +
  
  theme(plot.title = element_text(size = "25",
                                  face = "bold",
                                  colour = "darkblue",
                                  angle = 0,
                                  hjust = 0.5,
                                  vjust = 0,
                                  margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")))  






# Principle coordinates analysis with unWeighted-unifrac distances
physeq.css@phy_tree
unifrac <- UniFrac(physeq.css, weighted = FALSE)
ordination.unifrac.pcoa <- ordinate(physeq.css, "PCoA", distance = unifrac) # calculate the resemblance and ordinate using PCoA
ordination.unifrac.pcoa$vectors # positions of your points on the pcoa graph
ordination.unifrac.pcoa$values #values to calculate the variance explained on each axis (dimension)


# unweighted plot and centriods
pcoa.unifrac <- plot_ordination(physeq.css,
                                ordination = ordination.unifrac.pcoa,
                                type = "Sample.Type",
                                color = "treatment_type",
                                shape = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values = cbbPalette) +
  stat_ellipse(geom = "polygon",
               level = 0.95,
               aes(fill = treatment_type),
               alpha = 0.2)
pcoa.unifrac + geom_point(size = 5)

unifrac.data <- pcoa.unifrac$data



# PERMANOVA - testing for differences in centroids
#Bray-curtis centriods comprison and testing
set.seed(9999) # for reproducibility
prok.dist.bray <- phyloseq::distance(physeq.css, "bray") # create bray-curtis distance matrix

result.bray <- adonis2(prok.dist.bray~treatment_type, 
                       as(sample_data(physeq.css), "data.frame")) #Are there significant changes?

result.bray

pairwise.adonis2(prok.dist.bray~treatment_type, 
                 as(sample_data(physeq.css), "data.frame"),
                 nperm = 2000)


# unweighted UniFrac centriod comparison and testing
set.seed(9999)

unifrac <- UniFrac(physeq.css, weighted = FALSE)

result.UNWEIGHT <- adonis2(unifrac_F~treatment_type, 
                           as(sample_data(physeq.css), "data.frame")) #Are there significant changes?

pairwise.adonis2(prok.dist.bray~treatment_type, 
                 as(sample_data(physeq.css), "data.frame"),
                 nperm = 2000)


### Absolute Abundance Plot ###
top20 <- names(sort(taxa_sums(physeq.no.chloro),
                    decreasing=TRUE))[1:20]

ps.top20 <- transform_sample_counts(physeq.no.chloro,
                                    function(OTU) OTU/sum(OTU))

ps.top20 <- prune_taxa(top20, ps.top20)

bar_abundnace <- plot_bar(ps.top20,
                          x = "Treatment",
                          fill = "Genus",
                          facet_grid= ~ Sample.Type)

bar_abundnace  + theme_bw() + geom_col(color = NA)

bar_abundnace$data

### Relative Abundance plot ###

ps@otu_table -> check_otu_ps
check_otu_ps_wide <- as.matrix(check_otu_ps)
check_otu_ps_wide <- as.data.frame(check_otu_ps_wide)
check_otu_ps_wide$maso <- row.names(check_otu_ps_wide)

check_otu_ps_wide <- check_otu_ps_wide[, c(ncol(check_otu_ps_wide), 1:(ncol(check_otu_ps_wide)-1))]

row.names(check_otu_ps_wide) <- NULL

library(reshape2)
reshape2::melt(check_otu_ps_wide, id = "maso") -> check_otu_ps_long

check_otu_ps_long_clean <- check_otu_ps_long |> subset(value != 0)

ps@sam_data -> check_sam_ps

check_sam_ps$maso <- row.names(check_sam_ps)

row.names(check_sam_ps) <- NULL

check_sam_ps <- as.data.frame(check_sam_ps)

check_sam_ps <- as.matrix(check_sam_ps)

check_sam_ps <- as.data.frame(check_sam_ps)

relative_abudance_1 <- merge(x = check_otu_ps_long_clean,
                             y = check_sam_ps[, c("maso", "Sample.ID", "Treatment", "treatment_type", "Sample.Type")],
                             all.x = TRUE,
                             by = "maso")

ps@tax_table -> check_tax_ps

check_tax_ps <- as.matrix(check_tax_ps)
check_tax_ps <- as.data.frame(check_tax_ps)

# function for checking NA in dataset
sapply(X = check_tax_ps,
       FUN = function(x){  length(which(is.na(x)))  } )

as.data.frame(table(check_tax_ps$Genus)) -> genus_1
# View(genus_1)

as.data.frame(table(check_tax_ps$Phylum)) -> phylum_1
# View(phylum_1)

check_tax_ps_clean <- check_tax_ps[complete.cases(check_tax_ps$Phylum), ]

# sapply(X = check_tax_ps_clean,
#        FUN = function(x){  length(which(is.na(x)))  } )

check_tax_ps_clean$variable <- row.names(check_tax_ps_clean)
row.names(check_tax_ps_clean) <- NULL

relative_abudance_2 <- merge(x = relative_abudance_1,
                             y = check_tax_ps_clean,
                             all.x = TRUE,
                             by = "variable")

relative_abudance_3 <- relative_abudance_2[, -1]

relative_abudance_3 <- relative_abudance_3[complete.cases(relative_abudance_3$Phylum), ]

sapply(X = relative_abudance_3,
       FUN = function(x){  length(which(is.na(x)))  } )


relative_abudance_3$treatment_type <- factor(relative_abudance_3$treatment_type,
                                             levels = c("Basal-Fecal",
                                                        "A20-Fecal",
                                                        "D2-Fecal",
                                                        "Basal-Feed",
                                                        "A20-Feed",
                                                        "D2-Feed",
                                                        "RAS-Water"))

### chẻ ra

sort(table(relative_abudance_3$Phylum, useNA = "always"),
     decreasing = TRUE)[1:17] -> phylum_1

names(phylum_1) -> top_phylum

relative_abudance_3_top <- relative_abudance_3[relative_abudance_3$Phylum %in% top_phylum, ]

relative_abudance_3_other <- relative_abudance_3[!(relative_abudance_3$Phylum %in% top_phylum), ]

relative_abudance_3_other$Phylum <- "Others"

relative_abudance_3_top <- rbind(relative_abudance_3_top, relative_abudance_3_other)

relative_abudance_3_top$Treatment <- factor(relative_abudance_3_top$Treatment,
                                            levels = c("Basal", "A20", "D2", "RAS"))

sort(unique(relative_abudance_3_top$Phylum))[-14] -> sort_phylum

relative_abudance_3_top$Phylum <- factor(relative_abudance_3_top$Phylum,
                                         levels = c(sort_phylum, "Others"))

relative_abudance_fecal <- relative_abudance_3_top |> subset(Sample.Type == "Fecal")

ggplot(data = relative_abudance_fecal,
       mapping = aes(x = Treatment,
                     y = value,
                     fill = Phylum)) +
  
  geom_bar(
    position = "fill",
    stat = "identity",
    # position = "dodge",
    width = 0.8) +
  
  # ggtitle(label = "FISH FECAL") +
  
  scale_x_discrete(name = "Treatment",
                   expand = expansion(mult = 0, add = 0.5),
                   labels = c("Basal", "CFP-20", "FSC-2")) +
  
  scale_y_continuous(name = "Relative abundance",
                     expand = expansion(mult = 0, add = 0),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  
  scale_fill_manual(name = "Phylum",
                    values = c("lightblue",
                               "#baffc9",
                               "#aaaaaa", 
                               "#1b85b8",                                             
                               "lightgoldenrodyellow",
                               "navajowhite",
                               "darkseagreen3",
                               "aquamarine3",
                               "lightgoldenrod2",
                               "salmon1",
                               "lightsteelblue3",
                               "plum3",
                               "rosybrown3",
                               "lightpink4",
                               "#2bbed8",
                               "pink2",
                               "#ffad60",
                               "#386f54")) +
  
  theme_bw() +
  
  theme(axis.text.x = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5),
        
        axis.text.y = element_text(size = 20,
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5)) +
  
  theme(axis.title.x = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 1, 
                                                    r = 0,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(axis.title.y = element_text(face = "bold", 
                                    colour = "black", 
                                    size = 20,
                                    margin = margin(t = 0, 
                                                    r = 1,
                                                    b = 0, 
                                                    l = 0, 
                                                    "cm"))) +
  
  theme(plot.margin = margin(t = 1, r = 1, b = 6, l = 1, "cm")) +
  
  theme(legend.position = "right") +
  
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  
  # theme(axis.line = element_line(color = "black",
  #                                linewidth = 1)) +
  
  theme(axis.ticks = element_line(color = "black",
                                  linewidth = 1)) +
  
  theme(axis.ticks.length = unit(0.25, "cm")) +
  
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.title.position = "top",
        legend.text = element_text(size = 15),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.key.spacing.y = unit(0.3, "cm"),
        legend.key.size = unit(0.8, "cm")) +
  
  theme(legend.position = c(-0.1, -0.37),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.margin = margin(b = 0, l = 0, t = 0, r = 0, unit = "mm"),
        legend.background = element_rect(fill = "transparent",
                                         color = "transparent",
                                         linewidth = 0.5),
        legend.box.margin = margin(b = 0, l = 8, t = 0, r = 0, unit = "mm"),
        legend.justification = c("left", "bottom")) +
  
  theme(plot.title = element_text(size = "25",
                                  face = "bold",
                                  colour = "darkblue",
                                  angle = 0,
                                  hjust = 0.5,
                                  vjust = 0,
                                  margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt"))) 

#END#
