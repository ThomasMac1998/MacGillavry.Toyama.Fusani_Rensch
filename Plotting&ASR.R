# load required packages 
library(caper)
library(ape)
library(phytools)
library(ggplot2)
library(ggtree) 
library(patchwork)
library(scales)
library(dplyr)
library(scico)

# Set working directory: ### 
setwd("/Users/ThomasMac/Desktop/2_Projects/B_Current/MacGillavry&Toyama_Tails/MacGillavry&Toyama")

# Import data and tree: ### 

# This tree is from Ligon et al. (2018) 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
# Let's take a look
plot(tree, cex = 0.5)

# Read the dataframe 
data <- read.csv("Tail.data.csv")

# Quick inspection 
str(data) # We've got data for 39 species 

# Calculate sexual dimorphism in total tail length and body mass: ### 

# Calculate the sexual dimorphism index for tarsus length and body mass: ### 
#data$tarsus_dimorphism <- ifelse(data$tarsus_M > data$tarsus_F,
#                                      (data$tarsus_M / data$tarsus_F) - 1,
#                                      -1 * ((data$tarsus_F / data$tarsus_M) - 1))
# Let's look at the distribution of the sexual dimorphism index: 
hist(data$SSD_tarsus)

#data$weight_dimorphism <- ifelse(data$weight_M > data$weight_F,
#                                 (data$weight_M / data$weight_F) - 1,
#                                 -1 * ((data$weight_F / data$weight_M) - 1))
# Let's look at the distribution of the sexual dimorphism index: 
hist(data$SSD_weight) # Interesting, this looks bimodal! 

# Calculate relative tail length dimorphism: 
#data$rel_tail_length_dim <- ifelse(data$rel_tail_length_M > data$rel_tail_length_F,
#                                 (data$rel_tail_length_M / data$rel_tail_length_F) - 1,
#                                 -1 * ((data$rel_tail_length_F / data$rel_tail_length_M) - 1))
# Let's look at the distribution of the sexual dimorphism index: 
hist(data$rel_tail_dim) 

# Calculate species average tarsus length and body mass: ### 
# data$average_tarsus <- (data$tarsus_M+data$tarsus_F)/2
# data$average_weight <- (data$weight_M+data$weight_F)/2
hist(data$avg_tarsus)
hist(data$avg_weight)

# Create comparative dataset for caper 
comparative.df <- comparative.data(phy = tree, 
                                    data = data, 
                                    names.col = species, 
                                    vcv = TRUE, 
                                    na.omit = FALSE, 
                                    warn.dropped = TRUE) 
str(comparative.df) # Tree with 39 tips 

# Check to make sure if the tree is attached to the data properly: 
head(data, 4)
head(comparative.df$data, 4)
# Looks good. 

######################################################################################### 
######################################################################################### 

# Continuous trait reconstructions: ### 

# Load full dataset: 
data <- read.csv("Paradisaeidae.tail_data.csv")

# Read the tree again 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
tree <- drop.tip(tree, "Parotia_helenae")

## Specify row names so the data matches the tree
rownames(data) <- data$species

## Read the tree file
## The concatenated tree, already ultrametricized in Geneious 
print(tree, printlen = 2)

## First, let's plot the tree 
plotTree(tree,fsize=0.9,ftype="i",lwd=1)

## Convert to a vector 
xx<-setNames(data$tail_max_M,rownames(data))
xx

## Now we want to reconstruct ancestral states, and plot these results on a phylogeny. 

## Change data frame to a vector
xx<-as.matrix(xx)[,1]
xx

## Estimate ancestral states 
fit <- fastAnc(tree, xx, vars=TRUE, CI=TRUE)
fit

## We can also calculate 95% CIs 
fit$CI[1,]
range(xx)

### Plot again in ggtree ### 

# Make a dataframe with trait values at the tips
td <- data.frame(
  node = nodeid(tree, names(xx)),
  trait = xx)

# Make a dataframe with estimated trait values at the nodes
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree, d, by = 'node')

# Plot a Phenogram: 

# Create the plot with the entire tree but only display tip labels for selected species
tree <- groupClade(tree, 69)

treeplot1 <- ggtree(tree, yscale = "trait", size = 0.5, aes(color = group)) +
  theme_classic() + 
  scale_x_continuous(name = "Time (myr)", 
                     breaks = seq(0, 100, by = 5)) + 
  scale_y_continuous(name = "Male tail length (mm)", limits = c(0, 900), 
                     breaks = seq(0, 900, by = 100)) + 
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position="none",
    legend.justification = c(0, 1), 
    plot.title = element_text(size = 9, face = "bold")
  ) +
#  geom_point2(aes(subset = (node == 69)), shape = 21, size = 2, fill = 'darkgoldenrod2') + 
  scale_color_manual(values = c("black", "darkgoldenrod2"), 
                     labels = c("No", "Yes")) +  
  labs(color = "Tail wires") + 
  labs(color = "Tail wires")  

######################################################################################### 

### Reconstruction of female raw tail length values: #

# Read the tree again 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
tree <- drop.tip(tree, "Parotia_helenae")

## Convert to a vector 
xx<-setNames(data$tail_max_F,rownames(data))
xx

## Now we want to reconstruct ancestral states, and plot these results on a phylogeny. 

## Change data frame to a vector
xx<-as.matrix(xx)[,1]
xx

## Estimate ancestral states 
fit <- fastAnc(tree, xx, vars=TRUE, CI=TRUE)
fit

## We can also calculate 95% CIs 
fit$CI[1,]
range(xx)

### Plot again in ggtree ### 

# Make a dataframe with trait values at the tips
td <- data.frame(
  node = nodeid(tree, names(xx)),
  trait = xx)

# Make a dataframe with estimated trait values at the nodes
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree, d, by = 'node')

# Plot a Phenogram: 

# Create the plot with the entire tree but only display tip labels for selected species
tree <- groupClade(tree, 69)

treeplot2 <- ggtree(tree, yscale = "trait", size = 0.5, aes(color = group)) +
  theme_classic() + 
  scale_x_continuous(name = "Time (myr)", 
                     breaks = seq(0, 100, by = 12.5)) + 
  scale_y_continuous(name = "Female tail length (mm)", limits = c(0, 900), 
                     breaks = seq(0, 900, by = 100), position = "right") + 
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position="none",
    legend.justification = c(0, 1),
    plot.title = element_text(size = 9, face = "bold")  
  ) +
#  geom_point2(aes(subset = (node == 69)), shape = 21, size = 2, fill = 'darkgoldenrod2') + 
  scale_color_manual(values = c("black", "darkgoldenrod2"), 
                     labels = c("No", "Yes")) +  
  labs(color = "Tail wires") + 
  scale_x_reverse() 

###########################################################################################
###########################################################################################

### Reconstruction of male relative tail length values: #
str(data)
# Load data with relative tail length values: 
data <- read.csv("Tail.data.csv", row.names = 2)
# Read the tree again 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
# Remove NA value species (those without body mass information) from data frame and tree
tree <- drop.tip(tree, c("Parotia_helenae", 
                         "Astrapia_nigra", 
                         "Paradigalla_carunculata", 
                         "Paradisaea_apoda", 
                         "Paradisaea_decora", 
                         "Paradisaea_guilielmi"))

# Check which node we need to group for later plotting. 
plot(tree)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree))
# node 28 in this tree. 

## Convert to a vector 
xx<-setNames(data$rel_tail_length_M,rownames(data))
xx

## Now we want to reconstruct ancestral states, and plot these results on a phylogeny. 

## Change data frame to a vector
xx<-as.matrix(xx)[,1]
xx

## Estimate ancestral states 
fit <- fastAnc(tree, xx, vars=TRUE, CI=TRUE)
fit

## We can also calculate 95% CIs 
fit$CI[1,]
range(xx)

### Plot again in ggtree ### 

# Make a dataframe with trait values at the tips
td <- data.frame(
  node = nodeid(tree, names(xx)),
  trait = xx)

# Make a dataframe with estimated trait values at the nodes
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree, d, by = 'node')

# Plot a Phenogram: 

# Create the plot with the entire tree but only display tip labels for selected species
tree <- groupClade(tree, 62)

range(data$rel_tail_length_M)

treeplot3 <- ggtree(tree, yscale = "trait", size = 0.5, aes(color = group)) +
  theme_classic() + 
  scale_x_continuous(name = "Time (myr)", 
                     breaks = seq(0, 100, by = 5)) + 
  scale_y_continuous(name = "Male relative tail length", limits = c(0.75, 2.25), 
                     breaks = seq(0.75, 2.25, by = 0.25)) + 
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position="none",
    legend.justification = c(0, 1),
    plot.title = element_text(size = 9, face = "bold")  
  ) +
#  geom_point2(aes(subset = (node == 62)), shape = 21, size = 2, fill = 'darkgoldenrod2') + 
  scale_color_manual(values = c("black", "darkgoldenrod2"), 
                     labels = c("No", "Yes")) +  
  labs(color = "Tail wires") 

###########################################################################################

### Reconstruction of female relative tail length values: #

# Read the tree again 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
# Remove NA value species (those without body mass information) from data frame and tree
tree <- drop.tip(tree, c("Parotia_helenae", 
                         "Astrapia_nigra", 
                         "Paradigalla_carunculata", 
                         "Paradisaea_apoda", 
                         "Paradisaea_decora", 
                         "Paradisaea_guilielmi"))

## Convert to a vector 
xx<-setNames(data$rel_tail_length_F,rownames(data))
xx

## Now we want to reconstruct ancestral states, and plot these results on a phylogeny. 

## Change data frame to a vector
xx<-as.matrix(xx)[,1]
xx

## Estimate ancestral states 
fit <- fastAnc(tree, xx, vars=TRUE, CI=TRUE)
fit

## We can also calculate 95% CIs 
fit$CI[1,]
range(xx)

### Plot again in ggtree ### 

# Make a dataframe with trait values at the tips
td <- data.frame(
  node = nodeid(tree, names(xx)),
  trait = xx)

# Make a dataframe with estimated trait values at the nodes
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree, d, by = 'node')

# Plot a Phenogram: 

# Create the plot with the entire tree but only display tip labels for selected species
tree <- groupClade(tree, 62)

range(data$rel_tail_length_M)

treeplot4 <- ggtree(tree, yscale = "trait", size = 0.5, aes(color = group)) +
  theme_classic() + 
  scale_x_continuous(name = "Time (myr)", 
                     breaks = seq(0, 100, by = 12.5)) + 
  scale_y_continuous(name = "Female relative tail length", limits = c(0.75, 2.25), position = "right", 
                     breaks = seq(0.75, 2.25, by = 0.25)) + 
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position="none",
    legend.justification = c(0, 1),
    plot.title = element_text(size = 9, face = "bold")  
  ) +
#  geom_point2(aes(subset = (node == 62)), shape = 21, size = 2, fill = 'darkgoldenrod2') + 
  scale_color_manual(values = c("black", "darkgoldenrod2"), 
                     labels = c("No", "Yes")) +  
  labs(color = "Tail wires") + 
  scale_x_reverse()

###########################################################################################
###########################################################################################

# Stitch the plots together: ### 
p1 <- (treeplot1|plot_spacer()|treeplot2) + plot_layout(widths = c(4.5, -3, 4.5))
p2 <- (treeplot3|plot_spacer()|treeplot4) + plot_layout(widths = c(4.5, -3, 4.5))
print(p1)
print(p2)

###########################################################################################
###########################################################################################

# We can also plot the distribution of tail lengths across all birds: 

data.all <- read.csv("AllBirds_tail.data.csv")
str(data.all)

# Filter only males 
data.M <- filter(data.all, Sex == "M")

# Filter only females 
data.F <- filter(data.all, Sex == "F")

# Now for only passerines 

# Males 
Passerines.M <- filter(data.M, Order2 == "Passeriformes")

# Females 
Passerines.F <- filter(data.F, Order2 == "Passeriformes")

# Filter species with Family2 == 'Paradisaeidae'
paradisaeidae_species <- Passerines.M[Passerines.M$Family2 == 'Paradisaeidae', ]
#menuridae_species <- Passerines.M[Passerines.M$Family2 == 'Menuridae', ]
paradisaeidae_species.2 <- Passerines.F[Passerines.F$Family2 == 'Paradisaeidae', ]
#menuridae_species.2 <- Passerines.F[Passerines.F$Family2 == 'Menuridae', ]

# Create a new variable to distinguish 'Paradisaeidae' from other families
Passerines.M$Family2_grouped <- ifelse(Passerines.M$Family2 == "Paradisaeidae", 
                                       "Paradisaeidae", 
                                       "Other Families")

# For males
M <- ggplot(Passerines.M, aes(mean_Tail.Length)) +
  geom_histogram(binwidth = 10, fill = "white", colour = "grey80") +
  theme_classic() + 
  scale_y_continuous(breaks = seq(0, 900, by = 300), expand = c(0, 0)) + 
  theme(legend.position = "none") +
  labs(title = "Males (n = 5950)", 
       x = NULL, 
       y = "# species") + 
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12), 
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    plot.title = element_text(size = 10, face = "bold")
  ) +
  geom_segment(data = paradisaeidae_species, 
               aes(x = mean_Tail.Length, xend = mean_Tail.Length, y = 0, yend = 100, color = Wire.Tails)) + 
  scale_x_continuous(limits = c(0, 900), breaks = seq(0, 900, by = 100)) + 
  scale_color_manual(values = c("Y" = "darkgoldenrod2", "N" = "black"))  + 
  geom_boxplot(data = Passerines.M, 
               aes(x = mean_Tail.Length, y = 400, fill = Family2_grouped), 
               width = 100, 
               size = 0.5, 
               outlier.shape = NA, 
               position = position_dodge(width = 175), 
               alpha = 0.75) +
  scale_fill_manual(values = c("Paradisaeidae" = "darkgrey", "Other Families" = "white")) 

# Create a new variable to distinguish 'Paradisaeidae' from other families
Passerines.F$Family2_grouped <- ifelse(Passerines.F$Family2 == "Paradisaeidae", 
                                       "Paradisaeidae", 
                                       "Other Families")

# For females
F <- ggplot(Passerines.F, aes(mean_Tail.Length)) +
  geom_histogram(binwidth = 10, fill = "white", colour = "grey80") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none") +
  labs(title = "Females (n = 5618)", 
       x = "Mean tail length (mm)", 
       y = "# species") + 
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12), 
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    plot.title = element_text(size = 10, face = "bold")
  ) +
  geom_segment(data = paradisaeidae_species.2, 
               aes(x = mean_Tail.Length, xend = mean_Tail.Length, y = 0, yend = 100, color = Wire.Tails)) +  # Color by Wire.tails
  scale_x_continuous(limits = c(0, 900), breaks = seq(0, 900, by = 100)) + 
  scale_color_manual(values = c("Y" = "darkgoldenrod2", "N" = "black"))  + 
  geom_boxplot(data = Passerines.F, 
               aes(x = mean_Tail.Length, y = 400, fill = Family2_grouped), 
               width = 100, 
               size = 0.5, 
               outlier.shape = NA, 
               position = position_dodge(width = 175), 
               alpha = 0.75) +
  scale_fill_manual(values = c("Paradisaeidae" = "darkgrey", "Other Families" = "white")) 

# Plots together
(M/F) + plot_annotation(tag_levels = 'A')

######################################################################################### 
######################################################################################### 

df <- read.csv("Tail.data.csv")

# Rensch's rule for body size
p1 <- ggplot(df, aes(x = avg_weight, y = SSD_weight, shape = Clade)) +
  geom_abline(intercept = 0, slope = 0, color = "black", linetype = 3, size = 0.5) + 
  geom_point(fill = "grey90", size = 2.5, alpha = 1) +   
  scale_shape_manual(values = c("Core" = 21, "Basal" = 24)) +  
  labs(
    title = "",    
    y = "Body mass dimorphism",              
    x = "Species body mass (g)"             
  ) +
  theme_bw() + 
  theme(
    axis.title = element_text(size = 10, face = "bold"),          
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis ticks by 45 degrees
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = log10(seq(50, 450, by = 50)),  
                     labels = seq(50, 450, by = 50)) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1), 
                     labels = seq(0, 0.6, by = 0.1), 
                     limits = c(-0.1, 0.6)) + 
  geom_abline(intercept = 0.62658, slope = -0.20618, color = "black", linetype = 2, size = 0.5) + 
  geom_abline(intercept = -0.97768, slope = 0.56167, color = "black", linetype = 1, size = 0.5) 

print(p1)

# For polygamous:
#a = -1.2105356, b = 0.6707531

#For monogamous:
#a = 0.5585915, b = -0.1784704

# The ancova doesn't give you p-values for each regression line, 
# The key result is the significant difference in the slopes of both lines:
# t = 2.22, p = 0.034

######################################################################################### 
######################################################################################### 

# Rensch's rule for relative tail length: ### 

# Plot the original data and the predicted regression line
p2 <- ggplot(df, aes(x = avg_weight, y = rel_tail_dim, fill = tail_wires, shape = Clade)) +
  geom_abline(intercept = 0, slope = 0, color = "black", linetype = 3, size = 0.5) + 
  geom_point(col = "white", size = 2.5, alpha = 1) +   
  labs(
    title = "",    
    y = "Relative tail length dimorphism",              
    x = "Species body mass (g)"             
  ) +
  theme_bw() + 
  theme(
    axis.title = element_text(size = 10, face = "bold"),          
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis ticks by 45 degrees
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = log10(seq(50, 450, by = 50)),  
                     labels = seq(50, 450, by = 50)) +
  scale_y_continuous(breaks = seq(-0.25, .75, by = 0.25), 
                     labels = seq(-0.25, .75, by = 0.25), 
                     limits = c(-0.25, .75)) +
  scale_fill_manual(values = c("Y" = "darkgoldenrod2", "N" = "black")) + 
  scale_shape_manual(values = c("Core" = 21, "Basal" = 24)) + 
  # Plot means of each species (intercept only without slope): 
  geom_abline(intercept = mean(df$rel_tail_dim[df$category == "wired"]), slope = 0, color = "darkgoldenrod2", linetype = 1, size = 0.5) + 
  geom_abline(intercept = mean(df$rel_tail_dim[df$category == c("long", "short")]), slope = 0, color = "black", linetype = 1, size = 0.5) 

# New version of the plot but without monogamy shown, instead with tail 
# categories shown in the same colour scheme as Fig. 2A. 
p2 <- ggplot(df, aes(x = avg_weight, y = rel_tail_dim, fill = category, col = category)) +
  geom_abline(intercept = 0, slope = 0, color = "black", linetype = 3, size = 0.5) + 
  geom_point(shape = 21, size = 2.5, alpha = 1) +   
  labs(
    title = "",    
    y = "Relative tail length dimorphism",              
    x = "Species body mass (g)"             
  ) +
  theme_bw() + 
  theme(
    axis.title = element_text(size = 10, face = "bold"),          
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis ticks by 45 degrees
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = log10(seq(50, 450, by = 50)),  
                     labels = seq(50, 450, by = 50)) +
  scale_y_continuous(breaks = seq(-0.25, .75, by = 0.25), 
                     labels = seq(-0.25, .75, by = 0.25), 
                     limits = c(-0.25, .75)) +
  scale_fill_manual(values = c("wired" = "darkgoldenrod2", "long" = "black", "short" = "white")) + 
  scale_colour_manual(values = c("wired" = "black", "long" = "black", "short" = "grey60")) +
  scale_shape_manual(values = c("Core" = 21, "Basal" = 24)) + 
  # Plot means of each species (intercept only without slope): 
  geom_abline(intercept = mean(df$rel_tail_dim[df$category == "wired"], na.rm = TRUE), slope = 0, color = "darkgoldenrod2", linetype = 1, size = 0.5) + 
  geom_abline(intercept =  mean(df$rel_tail_dim[df$category == "long"], na.rm = TRUE), slope = 0, color = "black", linetype = 1, size = 0.5) + 
  geom_abline(intercept =  mean(df$rel_tail_dim[df$category == "short"], na.rm = TRUE), slope = 0, color = "grey60", linetype = 1, size = 0.5) 

print(p2)

# Also make a boxplot to add:
p3 <- ggplot(df, aes(x = factor(category, levels = c("short", "long", "wired")), 
                     y = rel_tail_dim, 
                     col = category)) +
  geom_boxplot(linetype = 1, size = 0.75) + 
  labs(
    title = "",    
    y = NULL,              
    x = NULL             
  ) +
  theme_minimal() + 
  theme(
    axis.title = element_text(size = 10, face = "bold"),          
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis ticks by 45 degrees
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = seq(-0.25, .75, by = 0.25), 
                     labels = seq(-0.25, .75, by = 0.25), 
                     limits = c(-0.25, .75)) +
  scale_colour_manual(values = c("wired" = "darkgoldenrod2", "long" = "black", "short" = "grey")) + 
  geom_abline(intercept = 0, slope = 0, color = "black", linetype = 3, size = 0.5) + 
  scale_x_discrete(labels = c("short" = "Short", "long" = "Long", "wired" = "Wires"))

######################################################################################### 
######################################################################################### 

# Stitch plots together 
(p1 | p2 | p3) + 
  plot_annotation(tag_levels = 'A') +  
  plot_layout(widths = c(1, 1, 0.5))

######################################################
######################################################

# ASR for different tail elaboration categories: # 

# Read the data with tail length categories: # 
X<-read.csv("categories.csv",row.names=1)
str(X)
X$category <- as.factor(X$category)
TA<-setNames(X[, 1],rownames(X))
TA <- TA[!names(TA)]
TA

# Read the tree again: # 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
tree <- drop.tip(tree, "Parotia_helenae") 
# Removed Eastern Parotia as poorly known and often considered to be subspecies of P. lawesii

# Estimate ancestral states under a ER model
fitER<-ace(TA,tree,model="ER",type="discrete")
fitER

fitER$lik.anc

# Estimate ancestral states under a ARD model
fitARD<-ace(TA,tree,model="ARD",type="discrete")
fitARD

fitARD$lik.anc

# Estimate ancestral states under a SYM model
fitSYM<-ace(TA,tree,model="SYM",type="discrete")
fitSYM

fitSYM$lik.anc

### Model comparisons 

# fit of models using AIC
AIC<-setNames(sapply(list(fitER,fitSYM,fitARD),AIC),c("ER","SYM","ARD"))
AIC
aic.w(AIC) 

# Plot the ARD reconstruction
cols<-setNames(c("black", "white", "darkgoldenrod2"),levels(TA)) # Set colours
plotTree(tree, fsize=0.45, ftype="i",lwd=1.5, offset=0.25) # Plot the tree
nodelabels(node=1:tree$Nnode+Ntip(tree), # Plot ancestral states as pies
           pie=fitER$lik.anc,piecol=cols,cex=1) 
tiplabels(pie=to.matrix(TA[tree$tip.label], # Plot tip states
                        levels(TA)),piecol=cols,cex=0.5)

######################################################
######################################################

# Maybe some plots for the graphical abstract? 



