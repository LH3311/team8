# Group 8
# MICB 475 Group Project February 11, 2024
# Script adapted from "Module13_phyloseq.R" by Evelyn Sun

# Load in Packages
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### Load in data ####

metafp <- "pd_export/parkinsons_metadata.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "pd_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "pd_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "pd_export/tree.nwk"
phylotree <- read.tree(phylotreefp)


#### Format files and create phylsoleq object ####

# Substitute the ASV's for Index #'s in ASV (OTU) table and save as a matrix
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

# Substitute the Sample Names for Index #'s in metadata and save as a date.frame
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'#SampleID'
SAMP <- sample_data(samp_df)
class(SAMP)

# Convert taxon strings to a table with separate taxa rank columns, subustitute the Feature ID's for Index #'s
# and save as a matrix
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

# Merge all into a phyloseq object
pd <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Filter and rarify data ####

# Remove non-bacterial sequences
pd_filt <- subset_taxa(pd, Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Remove ASVs that have less than 5 counts total
pd_filt_nolow <- filter_taxa(pd_filt, function(x) sum(x)>5, prune = TRUE)

# Remove samples with less than 100 reads
pd_filt_nolow_reads <- prune_samples(sample_sums(pd_filt_nolow)>100, pd_filt_nolow)

# Remove samples where Glucose is na
pd_final <- subset_samples(pd_filt_nolow_reads, !is.na(Glucose) )

# Rarefy samples
rarecurve <- rarecurve(t(as.data.frame(otu_table(pd_final))), cex=0.1)

# Need to decide on Sample Size (5000) and re-run
pd_rare <- rarefy_even_depth(pd_final, rngseed = 1, sample.size = 5000) 
pd_rare

##### Saving #####
save(pd_final, file="pd_final.RData")
save(pd_rare, file="pd_rare.RData")

