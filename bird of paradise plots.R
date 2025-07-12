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
setwd("")

######################################################################################### 
######################################################################################### 

# Continuous trait reconstructions: ### 

# Load full dataset and drop NA species for relative male tail length:
data <- read.csv("tail_data.csv") 
# This dataset includes our full dataset for species where we could calculate
# both species average mass as well as relative tail length for both sexes. 
# Read the tree again and remove species without relative tail length values: 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
tree <- drop.tip(tree, c("Parotia_helenae", 
                         "Astrapia_nigra", 
                         "Paradigalla_carunculata", 
                         "Paradisaea_apoda", 
                         "Paradisaea_guilielmi", 
                         "Paradisaea_decora"))

## Specify row names so the data matches the tree
rownames(data) <- data$species

## Read the tree file
## The concatenated tree, already ultrametricized in Geneious 
print(tree, printlen = 2)

## First, let's plot the tree 
plotTree(tree,fsize=0.9,ftype="i",lwd=1)

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

# Plot the tree with ancestral states: 
p <- ggtree(tree, layout='rectangular', ladderize = FALSE, size=2.2) + 
  geom_tree(aes(color=trait), continuous = 'colour', size=1.8) +  
  scale_color_gradientn(colours=c("dodgerblue","white", 'firebrick'), 
                        limits = c(min(data$rel_tail_length_M), max(data$rel_tail_length_M))) + # Can remove this limits section here... 
  geom_tiplab(
    aes(label = gsub("_", " ", label), fontface = "italic"),
    color = "black",
    size = 2.5,
    hjust = -0.07
  ) +
  xlim(-5, 40) +
  theme(
    legend.position = c(0, 0.8),           # top-left inside plot (x=0, y=1)
    legend.justification = c(0, 0.8),      # align legend box to top-left corner
    legend.direction = "vertical",       # or "horizontal" if preferred
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)) 

p

range(data$rel_tail_length_M)
#########################################################################################
######################################################################################### 

# Continuous trait reconstructions: ### 

# Load full dataset and drop NA species for relative male tail length:
data <- read.csv("tail_data.csv") 
# This dataset includes our full dataset for species where we could calculate
# both species average mass as well as relative tail length for both sexes. 
# Read the tree again and remove species without relative tail length values: 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
tree <- drop.tip(tree, c("Parotia_helenae", 
                         "Astrapia_nigra", 
                         "Paradigalla_carunculata", 
                         "Paradisaea_apoda", 
                         "Paradisaea_guilielmi", 
                         "Paradisaea_decora"))

## Specify row names so the data matches the tree
rownames(data) <- data$species

## Read the tree file
## The concatenated tree, already ultrametricized in Geneious 
print(tree, printlen = 2)

## First, let's plot the tree 
plotTree(tree,fsize=0.9,ftype="i",lwd=1)

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

p2 <- ggtree(tree, layout='rectangular', ladderize = FALSE, size=2.2) + 
  geom_tree(aes(color=trait), continuous = 'colour', size=1.8) +  
  scale_color_gradientn(colours=c("dodgerblue","white", 'firebrick'), 
                        limits = c(min(data$rel_tail_length_M), max(data$rel_tail_length_M))) + # Can remove this limits section here... 
  geom_tiplab(
    aes(label = gsub("_", " ", label), fontface = "italic"),
    color = "black",
    size = 2.5,
    hjust = -0.07
  ) +
  xlim(-5, 40) +
  theme(
    legend.position = c(0, 0.8),           # top-left inside plot (x=0, y=1)
    legend.justification = c(0, 0.8),      # align legend box to top-left corner
    legend.direction = "vertical",       # or "horizontal" if preferred
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)) 

p2


######################################################################################### 

# Load full dataset and drop NA species for relative male tail length:
data <- read.csv("tail_data.csv")
# Read the tree again and remove species without relative tail length values: 
tree <- read.nexus("Ligon.et.al._UltrametricTree")
tree <- drop.tip(tree, c("Parotia_helenae", 
                         "Astrapia_nigra", 
                         "Paradigalla_carunculata", 
                         "Paradisaea_apoda", 
                         "Paradisaea_decora", 
                         "Paradisaea_guilielmi"))

## Specify row names so the data matches the tree
rownames(data) <- data$species

## Read the tree file
## The concatenated tree, already ultrametricized in Geneious 
print(tree, printlen = 2)

## First, let's plot the tree 
plotTree(tree,fsize=0.9,ftype="i",lwd=1)

## Convert to a vector 
xx<-setNames(data$rel_tail_dim,rownames(data))
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

p3 <- ggtree(tree, layout='rectangular', ladderize = FALSE, size=2.2) + 
  geom_tree(aes(color=trait), continuous = 'colour', size=1.8) +  
  scale_color_gradientn(colours=c("darkgreen", "white", 'purple')) + 
  geom_tiplab(
    aes(label = gsub("_", " ", label), fontface = "italic"),
    color = "black",
    size = 2.5,
    hjust = -0.07
  ) +
  xlim(-5, 40) +
  theme(
    legend.position = c(0, 0.8),           # top-left inside plot (x=0, y=1)
    legend.justification = c(0, 0.8),      # align legend box to top-left corner
    legend.direction = "vertical",       # or "horizontal" if preferred
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)) 

print(p3)

###########################################################################################
###########################################################################################

###########################################################################################
###########################################################################################

# Final touches and stitch the plots together: ### 
p1 <- p + ggtitle("A — Male relative tail length")
p2 <- p2 + ggtitle("B — Female relative tail length")
p3 <- p3 + ggtitle("C — Sexual dimorphism in relative tail length")

(p1|p2|p3) 

# This creates the basis of Figure 3. 

###########################################################################################
###########################################################################################

# Making Figure 2: 

data.all <- read.csv("Figure3data.csv")

# Color gradient for hex bins
colorful_gradient <- c("dodgerblue", "yellow", "red") 

fill_scale <- scale_fill_gradientn(
  colors = colorful_gradient,
  limits = c(0, 50),  # Ensures both plots share the same scale
  breaks = c(10, 20, 30, 40, 50),
  name = "Count"
)

# Filter for Paradisaeidae
paradise_only <- data.all %>% filter(Family2 == "Paradisaeidae")

# Compute convex hulls
hull_male <- paradise_only %>%
  filter(!is.na(Mass), !is.na(mean_Tail.LengthM)) %>%
  mutate(x = log10(Mass), y = log10(mean_Tail.LengthM)) %>%
  slice(chull(x, y))

hull_female <- paradise_only %>%
  filter(!is.na(Mass), !is.na(mean_Tail.LengthF)) %>%
  mutate(x = log10(Mass), y = log10(mean_Tail.LengthF)) %>%
  slice(chull(x, y))

# Plot for Males
tl1 <- ggplot(data.all, aes(x = log10(Mass), y = log10(mean_Tail.LengthM))) + 
  ylab(expression(log[10] * " Tail length (mm)")) + 
  xlab("") + 
  theme_classic() + 
  ggtitle("Males") + 
  scale_y_continuous(breaks = seq(1, 3, by = 0.5), 
                     labels = seq(1, 3, by = 0.5), 
                     limits = c(1, 3)) + 
  geom_hex(aes(fill = ..count..), bins = 50) + 
  geom_point(data = paradise_only, 
             aes(x = log10(Mass), y = log10(mean_Tail.LengthM)), 
             shape = 21, color = "black", size = .5, stroke = 1) +
  geom_polygon(data = hull_male, aes(x = x, y = y), 
               fill = "grey50", alpha = 0.25, color = "black", 
               linetype = "solid", size = 0.3) +
  fill_scale

# Plot for Females
tl2 <- ggplot(data.all, aes(x = log10(Mass), y = log10(mean_Tail.LengthF))) + 
  ylab(expression(log[10] * " Tail length (mm)")) + 
  xlab(expression(log[10] * " Species mass (g)")) + 
  theme_classic() + 
  ggtitle("Females") + 
  scale_y_continuous(breaks = seq(1, 3, by = 0.5), 
                     labels = seq(1, 3, by = 0.5), 
                     limits = c(1, 3)) + 
  geom_hex(aes(fill = ..count..), bins = 50) + 
  geom_point(data = paradise_only, 
             aes(x = log10(Mass), y = log10(mean_Tail.LengthF)), 
             shape = 21, color = "black", size = .5, stroke = 1) +
  geom_polygon(data = hull_female, aes(x = x, y = y), 
               fill = "grey50", alpha = 0.25, color = "black", 
               linetype = "solid", size = 0.3) +
  fill_scale

# Combine plots with shared legend
(tl1 / tl2) + plot_layout(guides = "collect")



######################################################################################### 
######################################################################################### 

df <- read.csv("tail_data.csv")
# Rename vars:
df$Clade[df$Clade == "Basal"] <- "monogamous"
df$Clade[df$Clade == "Core"] <- "polygynous"

# Rensch's rule for body size
p1 <- ggplot(df, aes(x = avg_weight, y = SSD_weight, fill = Clade)) +
  geom_abline(intercept = 0, slope = 0, color = "black", linetype = 3, size = 0.5) + 
  geom_point(size = 1.5, shape = 21, alpha = .7) +   
  labs(
    title = "",    
    y = "Body weight dimorphism",              
    x = "Species body weight (g)"             
  ) +
  theme_classic() + 
  theme(
    axis.title = element_text(size = 10),          
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "right", 
    legend.title = element_text(size = 7, face = "bold")
  ) +
  scale_x_continuous(breaks = log10(seq(50, 450, by = 50)),  
                     labels = seq(50, 450, by = 50)) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1), 
                     labels = seq(0, 0.6, by = 0.1), 
                     limits = c(-0.1, 0.6)) + 
  geom_abline(intercept = 0.62658, slope = -0.20618, color = "black", linetype = 2, size = 0.5) + 
  geom_abline(intercept = -0.97768, slope = 0.56167, color = "black", linetype = 1, size = 0.5) + 
  scale_fill_manual(values = c("polygynous" = "black", "monogamous" = "white")) 

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
p2 <- ggplot(df, aes(x = avg_weight, y = rel_tail_dim, shape = category)) +
  geom_abline(intercept = 0, slope = 0, color = "black", linetype = 3, size = 0.5) + 
  geom_point(size = 1.5, alpha = 0.7, fill = "black") +  # FIX: ensure outline is black
  labs(
    title = "",    
    y = "Relative tail length dimorphism",              
    x = "Species body weight (g)",
    shape = "Tail phenotype"
  ) +
  theme_classic() + 
  theme(
    axis.title = element_text(size = 10),          
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "right", 
    legend.title = element_text(size = 7, face = "bold")
  ) +
  scale_x_continuous(breaks = log10(seq(50, 450, by = 50)),  
                     labels = seq(50, 450, by = 50)) +
  scale_y_continuous(breaks = seq(-0.25, 0.75, by = 0.25), 
                     labels = seq(-0.25, 0.75, by = 0.25), 
                     limits = c(-0.25, 0.75)) +
  scale_shape_manual(values = c("long" = 24, "short" = 21)) 

######################################################################################### 
######################################################################################### 

# Stitch plots together 
(p1 / p2) + plot_annotation(tag_levels = 'A') 
# This creates Figure 5. 

######################################################
######################################################

# Make a new plot for the wire-tailed birds: 
# Used to make Figure 4: 
library(ggbeeswarm)
df<-read.csv("tail_data.csv")
df$category2 <- factor(df$category2, levels = c("short", "long", "wired"))
t1 <- ggplot(df, aes(x = category2, y = rel_tail_dim, fill = category2)) +
  geom_violin(alpha = .75) + 
  geom_beeswarm(priority='ascending', cex = 2, alpha = .8, size = 1.2, shape = 21, fill = "white") + 
  theme_minimal() + 
  scale_fill_manual(values = c("white", "grey", "goldenrod1")) + 
  theme(legend.position = "none") + 
  stat_summary(fun = mean, geom = "point", width = 0.1, color = "black") +
  labs(x = "Tail category", y = "SD in relative tail length") 

tree <- groupClade(tree, 62)
t2 <- ggtree(tree, layout='rectangular', ladderize = FALSE, size=.5) + 
  geom_tree(aes(color = as.factor(group))) + 
scale_color_manual(values=c("black", "goldenrod1")) + 
  theme(legend.position = "none") 

(t2|plot_spacer()/t1) + plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1, 1.5))






