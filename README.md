# june2024
install.packages("vegan")
library(vegan)
packageVersion("vegan")
setwd("D:/Marine_Iguanas_Project/MARINE_IGUANAS/JUNI2024/")
list.files()
table=read.table("rbcLOTUs_Table_june2024.txt", 
                 header = T,
                 row.names = 1,
                 check.names = FALSE)
head(table)[1:2,]
dim(table)
metadata = read.table("rbcL_Metadata_june2024.txt",
                      header = T,
                      row.names = 1,
                      check.names = FALSE)
head(metadata)
dim(metadata)

######RAREFACTION CURVES----
rarecurve(t(table),step=50,cex=1,label = TRUE ,ylab = "No.of ASVs",xlab = "No.of rbcL sequences")

#Export rarefaction plot to pdf----
dev.copy2pdf(out.type = "pdf")
dev.off()
#rename the pdf
pdf_filename = paste("rarefaction_curves_rbcL_gemapisr.pdf")
file.rename(from = "Rplot.pdf", to = pdf_filename)

#######CORRELATION BETWEEN NUMBER OF SEQS PER SAMPLE AND NUMBER OF ASVs----

library(stats)
ASVs_vs_seqs = read.table("rbcL_gemapisr_SumSeqs_vs_SumASVs.txt", header = T, check.names = FALSE)
head(ASVs_vs_seqs)
dim(ASVs_vs_seqs)

#Spearman Rank Correlation Statistics (Spearman R=)
cor(ASVs_vs_seqs$SumSeqs, ASVs_vs_seqs$SumASVs, method = "spearman")

#Get the p-value (here 'rho' is the same as Sperman R value from the previous command)
cor.test(ASVs_vs_seqs$SumSeqs, ASVs_vs_seqs$SumASVs, method = "spearman", exact = FALSE)

# If the number of ASVs and number of sequences per sample demonstrate significant correlation, then
# we have to take this into account in our analyses -> 
# either doing rarefaction of the data set or taking the number of sequences per sample as a covariate.
# [this correlation is almost always significant. But in any case, sequencing depth per sample should be comparable]
# 0.00 to 0.19 very weak correlation
# 0.20 to 0.39 weak correlation
# 0.40 to 0.69 moderate correlation
# 0.7 to 0.89 strong correlation
# higher is very strong correlation

#REMEMBER!
# we are considering tag-sequencing or amplicon data. This type of sequence data is compositional, 
#because it does not represent true absolute abundances. The total number of sequences in your dataset is arbitrary 
#(set by sequence depth), thus it is inappropriate to make conclusions about the abundance of a given species or taxa 
#(what was targeted in the sequencing effort).
#It is important to acknowledge this in both your an analysis and interpretation of the data.


### Rarefy OTU/ASV table to even sequencing depth----

# ! INSPECT THE RAREFACTION CURVES and AVS/OTU TABLE to decide the rarefaction depth

#load libs
library('base')
library('phyloseq')
library('vegan')

#check version
packageVersion("phyloseq")

###import OTU table 'as.matrix' for phyloseq package
otumat = as.matrix(as.data.frame(read.table("rbcLOTUs_Table_june2024.txt", 
                                            header = TRUE, 
                                            row.names = 1)))

###Convert imported OTU table to phyloseq item
OTU = otu_table(otumat, taxa_are_rows = TRUE)
OTU_physeq = phyloseq(OTU)

###rarefy OTU table to even sequencing depth;
###select depth in option "sample.size = "
raref_tab <- rarefy_even_depth(OTU_physeq, sample.size = min(sample_sums(OTU_physeq)),
                               rngseed = FALSE, 
                               replace = TRUE, 
                               trimOTUs = TRUE, 
                               verbose = TRUE) 

#sample.size = min(sample_sums(OTU_physeq)) --> rarefy to the depth according to sample that has lowest number of seqs

###read what this command does
#?rarefy_even_depth

### Write the rarefied OTU table to a file (tab delimited txt file)
OTU_rarefied = as(otu_table(raref_tab), "matrix")
###Convert to data.frame
OTUdf = as.data.frame(OTU_rarefied)
###Write to file
write.table(OTUdf, file = "rbcL_rarefied_OTU_table.txt", 
            sep = "\t", 
            col.names = NA,
            quote = FALSE)

raref_table=read.table("rbcL_rarefied_OTU_table.txt",
                     header = T,
                     row.names = 1,
                     check.names = FALSE)

### rarefaction plot per sample for rarefied data set #notworkin
rarecurve(t(raref_table),step=50,cex=1,label = TRUE, ylab = "No.of ASVs",xlab = "No.of rbcL sequences")

######################
### PERMANOVA test----
######################
# test the differences between 'treatments/groups' 

###read in the rarefied OTU table, if the table was rarefied
table = read.table("rbcLOTUs_Table_june2024.txt", 
                   header = TRUE,
                   row.names = 1, 
                   check.names = FALSE)
dim(table)
head(table)
###read metadata, DOUBLE-CHECK THAT SAMPLE LIST MATCH BETWEEN THE OTU TABLE AND METADATA FILE
metadata = read.table("rbcL_Metadata_june2024.txt", 
                      header = TRUE, 
                      row.names = 1, 
                      check.names = FALSE) 
dim(metadata)
head(table)

####Hellinger transformation----

#Hellinger transformation of the rarefied OTU/ASV table (sequence counts transformation prior analyses)
table_hell= decostand(t(table),method = "hellinger") 
#note that the table is transposed for Hellinger (t(table))


###export table after hellinger transformation for its use in primer
###Write to file
table_hell_trans=t(table_hell)
write.table(table_hell_trans, file = "rbcL_OTU_table_HELL.txt", 
            sep = "\t", 
            col.names = NA,
            quote = FALSE)

# Check if there is a significant effect of the treatment on the communities
# Variable to check = 'Subspecies' in metadata file
PERMANOVA = adonis2(table_hell ~ Subspecies+Island+Location, metadata, 
                    permutations = 9999, 
                    method = "bray")

#see the results
PERMANOVA

#plot permuted F-values
densityplot(permustats(PERMANOVA))

#permanova is not assumption free, so we also need to check for similar dispersion

### check the PERMANOVA assumption, assumption of homogeneity ###

###we first calculate the distances between samples
###As our data are basically abundance-based data, we have to calculate Bray-Curtis
###distances between samples instead of Euclidean distances. So:

dist.table = vegdist(table_hell, method= "bray")

###calculate multivariate dispersions----
mod <- betadisper(dist.table, metadata$Subspecies)

####perform anova
anova(mod)

#### Permutation test
permutest(mod,pairwise = TRUE, permutations=9999)

###Tukey Honest Significant Differences: NOT SUITABLE IF HOMOGENEITY IF VARIANCES IS NOT MET

#mod.HSD <- TukeyHSD(mod)
#mod.HSD

#psig=as.numeric(apply(mod.HSD$group[,2:3],1,prod)>=0)+1
#op=par(mar=c(4.2,12,4,2))
#plot((mod.HSD), las =1, col=psig)


###GAMES HOWELL POST HOC TEST----
install.packages("rstatix")
library(rstatix)

help(.games_howell_test)



## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)

## with data ellipses instead of hulls
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.9) # 90% data ellipse


## Draw a boxplot of the distances to centroid for each group
boxplot(mod)


mod3 <- betadisper(dist.table, metadata$Subspecies, type = "centroid")
mod3
permutest(mod3, permutations = 99)
anova(mod3)
plot(mod3)
boxplot(mod3)
plot(TukeyHSD(mod3),las=1)



-------------------------
anova(betadisper(dist.table, metadata$Subspecies))

###prueba
mod <- betadisper(dist.table, metadata$Subspecies)

plot(mod)

anova(mod)
permutest(mod, pairwise = TRUE, permutations = 99)

mod.HSD <- TukeyHSD(mod)
mod.HSD

psig=as.numeric(apply(mod.HSD$group[,2:3],1,prod)<=0)+1
op=par(mar=c(4.2,12,4,2))
plot((mod.HSD), las =1, col=psig)


TukeyHSD(mod, conf.level = .95)
plot (TukeyHSD(mod, conf.level = .95))
####
anova(betadisper(dist.table, metadata$Subspecies))
# if P val is larger than 0.05, then OK! 
# Otherwise the PERMANOVA difference may be due to variability in dispersion between groups

#Betadisper first calculates the average distance of group members to the group
#centroid in multivariate space (generated by a distance matrix). Then, an ANOVA is
#done to test if the dispersions (variances) of groups are different.


##################################################################
### Pairwise PERMANOVA----
##################################################################
install.packages('devtools')
library(devtools)
install.packages('pairwiseAdonis')
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") #install pairwiseAdonis
library(pairwiseAdonis)

#import OTU/ASV table
table=read.table("rbcLOTUs_Table_june2024.txt", header=T, row.names=1, check.names = FALSE)

#import variables/factors
metadata = read.table("rbcL_Metadata_june2024.txt", header=T, row.names=1, check.names = FALSE)

#Hellinger transformation of the OTU/ASV table (sequence counts transformation prior analyses)
table_hell=decostand(t(table),"hell")

PW_PERMANOVA = pairwise.adonis(table_hell, metadata$Island, perm = 9999)  #default sim.method is "bray"
PW_PERMANOVA #print results

############
### PLOTS----
############

###1. PCoA plot
library(phyloseq)
library(ggplot2)

#using Hellinger transformed table
table_hell_phy = otu_table(t(table_hell), 
                           taxa_are_rows = T)

#load metadata for phyloseq package
metadata_phy = sample_data(metadata)

#merge OTU table and metadata as phyloseq object
table_env_physeq = phyloseq(table_hell_phy, metadata_phy)

# generating and visualizing the PCoA plot with phyloseq
pcoa_plot <- ordinate(table_hell_phy, 
                      method = "PCoA", 
                      distance = "bray")

# draw the plot with ggplot2
PCoA_elipse <- plot_ordination(table_env_physeq, pcoa_plot, color = "Island") + 
  geom_point(size = 3) + 
  ggtitle("PCoA") + 
  theme_bw() 

print(PCoA_elipse)
PCoA_elipse+stat_ellipse()

#Examine more PCoA axes
# axes 1 and 3

PCoA_elipse_13 <- plot_ordination(table_env_physeq, pcoa_plot, color = "Island", axes = c(1, 3)) + 
  geom_point(size = 3) + 
  ggtitle("PCoA (Axes 1 and 3)") + 
  theme_bw()

print(PCoA_elipse_13)
PCoA_elipse_13 + stat_ellipse()

#axes 2 and 3

PCoA_elipse_23 <- plot_ordination(table_env_physeq, pcoa_plot, color = "Island", axes = c(2, 3)) + 
  geom_point(size = 3) + 
  ggtitle("PCoA (Axes 2 and 3)") + 
  theme_bw()

print(PCoA_elipse_23)
PCoA_elipse_23 + stat_ellipse()

#visualize how much variance each axis explains
#extract eigenvalues and plot the variance

# Extract eigenvalues
eigenvalues <- pcoa_plot$values$Eigenvalues
variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Create a bar plot of variance explained
variance_data <- data.frame(Axis = 1:length(variance_explained), Variance = variance_explained)
ggplot(variance_data, aes(x = Axis, y = Variance)) +
  geom_bar(stat = "identity") +
  ggtitle("Variance Explained by PCoA Axes") +
  xlab("PCoA Axis") +
  ylab("Variance Explained (%)") +
  theme_minimal()

#Enhanced Plotting with ggplot2

library(phyloseq)
library(ggplot2)

# Extract ordination scores
ordination_scores <- as.data.frame(vegan::scores(pcoa_plot))
ordination_scores$Subspecies <- as.factor(sample_data(table_env_physeq)$Subspecies)

# Plot Axes 1 and 2
ggplot(ordination_scores, aes(x = Axis.1, y = Axis.2, color = Subspecies)) +
  geom_point(size = 3) +
  stat_ellipse() +
  ggtitle("PCoA (Axes 1 and 2)") +
  theme_bw() +
  xlab(paste0("Axis 1 (", round(variance_explained[1], 2), "%)")) +
  ylab(paste0("Axis 2 (", round(variance_explained[2], 2), "%)"))

# Plot Axes 1 and 3
ggplot(ordination_scores, aes(x = Axis.1, y = Axis.3, color = Subspecies)) +
  geom_point(size = 3) +
  stat_ellipse() +
  ggtitle("PCoA (Axes 1 and 3)") +
  theme_bw() +
  xlab(paste0("Axis 1 (", round(variance_explained[1], 2), "%)")) +
  ylab(paste0("Axis 3 (", round(variance_explained[3], 2), "%)"))

# Extract ordination scores if pcoa_plot is an ordination object
ordination_scores <- as.data.frame(pcoa_plot$vectors)  # If pcoa_plot has ordination scores

# Extract sample data
sample_data_df <- as.data.frame(sample_data(table_env_physeq))

# Combine ordination scores with sample metadata
ordination_scores <- cbind(ordination_scores, sample_data_df)

library(ggplot2)

# Define a function to plot specified axes
plot_pcoa_axes <- function(ordination_scores, x_axis, y_axis) {
  ggplot(ordination_scores, aes_string(x = x_axis, y = y_axis, color = "Subspecies")) +
    geom_point(size = 3) +
    stat_ellipse() +
    ggtitle(paste("PCoA (Axes", gsub("Axis.", "", x_axis), "and", gsub("Axis.", "", y_axis), ")")) +
    theme_bw() +
    xlab(paste(x_axis)) +
    ylab(paste(y_axis))
}

# Plot Axes 1 and 2
print(plot_pcoa_axes(ordination_scores, "Axis.1", "Axis.2"))

# Plot Axes 1 and 3
print(plot_pcoa_axes(ordination_scores, "Axis.1", "Axis.3"))

# Plot Axes 2 and 3
print(plot_pcoa_axes(ordination_scores, "Axis.2", "Axis.3"))

####3D plot

library(plotly)

plot_3d <- plot_ly(ordination_scores, x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, color = ~Subspecies) %>%
  add_markers(size = 3) %>%
  layout(scene = list(xaxis = list(title = 'Axis 1'),
                      yaxis = list(title = 'Axis 2'),
                      zaxis = list(title = 'Axis 3')))

plot_3d

##########wow0
# Install packages if not already installed
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("magrittr", quietly = TRUE)) install.packages("magrittr")
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
if (!requireNamespace("phyloseq", quietly = TRUE)) install.packages("phyloseq")

# Load necessary libraries
library(dplyr)
library(magrittr)
library(plotly)
library(phyloseq)

# Assuming `table_env_physeq` is your phyloseq object
# Calculate PCoA
ordination <- ordinate(table_env_physeq, method = "PCoA", distance = "bray")

# Extract ordination scores for plotting
ordination_scores <- as.data.frame(ordination$vectors)
sample_data_df <- as.data.frame(sample_data(table_env_physeq))
ordination_scores <- cbind(ordination_scores, sample_data_df)

# Plot the PCoA in 3D
plot_3d <- plot_ly(ordination_scores, 
                   x = ~Axis.1, 
                   y = ~Axis.2, 
                   z = ~Axis.3, 
                   color = ~Subspecies) %>%
  add_markers(size = 3) %>%
  layout(scene = list(xaxis = list(title = 'Axis 1'),
                      yaxis = list(title = 'Axis 2'),
                      zaxis = list(title = 'Axis 3')))

# Print the plot
plot_3d

#####wow1

# Load ggplot2
library(ggplot2)

# Define a function to plot specified axes
plot_pcoa_axes <- function(ordination_scores, x_axis, y_axis) {
  ggplot(ordination_scores, aes_string(x = x_axis, y = y_axis, color = "Subspecies")) +
    geom_point(size = 3) +
    stat_ellipse() +
    ggtitle(paste("PCoA (Axes", gsub("Axis.", "", x_axis), "and", gsub("Axis.", "", y_axis), ")")) +
    theme_bw() +
    xlab(paste(x_axis)) +
    ylab(paste(y_axis))
}

# Plot Axes 1 and 2
print(plot_pcoa_axes(ordination_scores, "Axis.1", "Axis.2"))

# Plot Axes 1 and 3
print(plot_pcoa_axes(ordination_scores, "Axis.1", "Axis.3"))

# Plot Axes 2 and 3
print(plot_pcoa_axes(ordination_scores, "Axis.2", "Axis.3"))


#######wow2



#-------------
##### el cod q viene tiene las labels pero no corresponden a la subspecie correcta
plot_ordination(table_env_physeq, pcoa_plot, color = "Subspecies") + 
  geom_point(size = 3) + labs(col = "Subspecies") + 
  geom_text(aes(label = rownames(metadata), 
                hjust = 0.4, vjust = -0.6)) + 
  ggtitle("PCoA") + 
  theme_bw()   
#####

#edit figure with ggplot as you like (http://r-statistics.co/Complete-Ggplot2-Tutorial-Part1-With-R-Code.html)
#ggplot2 pdf manual -> https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf

### 2. NMDS plot
nmds_plot <- ordinate(table_hell_phy, 
                      method = "NMDS", 
                      distance = "bray")

#check the NMDS stress 
nmds_plot$stress

# draw the plot with ggplot2
NMDS_elipse <- plot_ordination(table_env_physeq, nmds_plot, color = "Subspecies") + 
  geom_point(size = 3)+
  ggtitle("NMDS")+
  theme_bw()

print(NMDS_elipse)
NMDS_elipse+stat_ellipse()

#####islands
NMDS_elipse <- plot_ordination(table_env_physeq, nmds_plot, color = "Island") + 
  geom_point(size = 3)+
  ggtitle("NMDS")+
  theme_bw()

print(NMDS_elipse)

NMDS_elipse+stat_ellipse()

#####################################
#####NICHE OVERLAP----
#######################################

install.packages("spaa")

library(spaa)
## to use niche overlap in the spaa package, a data matrix with each column for species 
####and eaach row for each plot? I dont know what they men with plot...
table_ASVs_subspecies=read.table("rbcL_gemapisr_ASV_vs_subspecies.txt", 
                 header = T,
                 row.names = 1,
                 check.names = FALSE)
head(table_ASVs_subspecies)[1:2,]
dim(table_ASVs_subspecies)
?niche.overlap

####table_ASVs_subspecies_trans <- t(table_ASVs_subspecies) 

niche.overlap(table_ASVs_subspecies, method="schoener")

####or use only presence/absence
table_binary <- ifelse(table_ASVs_subspecies >0,1,0)
niche.overlap(table_binary, method="schoener")

####or use table hellinger transformed
table_ASV_subs_hell <- decostand(t(table_ASVs_subspecies),"hell")
#but transpose the table back
niche.overlap(t(table_ASV_subs_hell), method="schoener")


##################################
### indicator species analyses ###
##################################
#load lib
install.packages("indicspecies")
library(indicspecies)
packageVersion("indicspecies")

#Check documentation here -> https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html

# OTU table for indicspecies analyses: OTUs are in columns, samples are in rows!
# table_hell (formed above) is in suitable format

head(table_hell)[1:2,]

#Define groups as based on the data from metadata file
#SORT METADATA FILE so that 'sample treatments' are grouped together (needed for later indicspecies analyses)
# e.g. first 28 rows (samples) are from godzilla and the following from mertensi
metadata = read.table("rbcL_gemapisr_metadata.txt",
                      header = T,
                      row.names = 1,
                      check.names = FALSE)
head(metadata)
dim(metadata)

groups = c(rep("godzilla", 33),
           rep("mertensi", 40),
           rep("hayampi", 45),
           rep("nanus", 33),
           rep("sielmanni", 26))

#read the rarefied OTU table #or not
#raref_table = read.table("rbcL_rarefied_OTU_table.txt", 
#                         header = TRUE, 
 #                        row.names = 1, 
  #                       check.names = FALSE)

table=read.table("rbcL_gemapisr_ASVstable.txt", 
                 header = T,
                 row.names = 1,
                 check.names = FALSE)
head(table)[1:2,]
dim(table)

#Sort the OTU table as based on metadata file, so match the grouping for the analyses
reorder_indexes <- match(rownames(metadata), colnames(table))
reorder_indexes

# Reorder the OTU table to have the sample names in the same order as the metadata data frame
table_ordered  <- table[ , reorder_indexes]

#Check if the sample names match between OTU table and metadata
colnames(table_ordered)
rownames(metadata)

#Write to file if needed
write.table(table_ordered, file = "rbcL_table_ordered.txt", 
            sep = "\t", 
            col.names = NA,
            quote = FALSE)


#Hellinger transformation of the table_ordered for the analyses
table_hell = decostand(t(table_ordered),"hell") 


###Indicator species analysis using multipatt
indval_results = multipatt(table_hell, 
                           groups, 
                           control = how(nperm = 9999))

#basic summary of the results
#List of ASVs/OTUs associated to each group
summary(indval_results)

#Examining the indicator value components
summary(indval_results, indvalcomp = TRUE)
# A = the probability value (0-1) that the ASV/OTU belongs to the target site group.
# Value of 1 means that this ASV/OTU in found ONLY in this group.
# This conditional probability is called the specificity.

# B = the probability of finding the species in sites belonging to the site group.
# Values of 1 means that all samples in this group contains this ASV/OTU. 
# This second conditional probability is called the fidelity or sensitivity.

#write the results to a file "indicspecies.results.txt"
capture.output(summary(indval_results, indvalcomp=TRUE), file = "indicspecies.results.txt")

#Find ubiquitous taxa, 
# if p value = NA then its ubiquitous (i.e. found everywhere)
write.table(indval_results$sign, file = "ubiquitous_taxa.txt", sep = "\t") 
