############################################################################
##### Title: Non-Euclidean Data Analysis With Metric Statistics (HDSR) #####
##### Description: Toy Example for Fréchet regression ######################
##### Manuscript reference : Section 2 #####################################
##### Figure reference : Figures 1&2 #######################################
##### Date: 09/18/25 #######################################################
##### Author: Yidong Zhou ##################################################
############################################################################

# This script creates Figures 1 and 2, demonstrating
# Fréchet regression on a toy example with weighted networks. The example
# shows how Fréchet regression can be applied to predict network structures
# based on a covariate X.
#
# Main components:
# 1. Data setup: Create three weighted triangle networks with different edge weights
# 2. Figure 1: Display the original three networks (G1, G2, G3) with their weights
# 3. Figure 2: Show predicted networks at X = 2 and X = 5 using Fréchet regression
#
# The example uses a simple triangle network structure to illustrate the concept
# of Fréchet regression in a non-Euclidean space (the space of weighted networks).

# Load required libraries
library(ggplot2)    # For creating plots
library(latex2exp)  # For rendering LaTeX expressions in plots
# =============================================================================
# DATA SETUP FOR FIGURE 1: ORIGINAL NETWORKS
# =============================================================================

# Define the number of networks and their corresponding covariate values
n <- 3                    # Number of networks in the dataset
X <- c(1, 4, 7)          # Covariate values for the three networks

# Create data frame for plotting the three triangle networks
# Each network is a triangle with vertices at (0,0), (2,0), and (1,√3)
# The edge weights are determined by the covariate X
df <- data.frame(
  # Edge coordinates: each triangle has 3 edges, repeated for 3 networks
  x = rep(c(0, 2, 1), each = 3),                    # Starting x-coordinates of edges
  y = rep(c(0, 0, sqrt(3)), each = 3),              # Starting y-coordinates of edges
  xend = rep(c(1, 0, 2), each = 3),                 # Ending x-coordinates of edges
  yend = rep(c(sqrt(3), 0, 0), each = 3),           # Ending y-coordinates of edges
  
  # Edge weights: different for each network based on covariate X
  # Pattern: [X, 2*X, X] for each network, representing different edge weight relationships
  weight = c(X, 2 * X, X),                          # Edge weights for all three networks
  
  # Node labels for each network
  labelnode = rep(c('$v_1$', '$v_2$', '$v_3$'), each = 3),
  
  # Weight labels (converted to fractions for better readability)
  labelweight = as.character(MASS::fractions(c(X, 2 * X, X))),
  
  # Network labels for faceting
  sample = factor(rep(1:n, 3), levels = 1:n, 
                 labels = c(TeX('$G_1$'), TeX('$G_2$'), TeX('$G_3$'))),
  
  # Node label positions (slightly offset from vertices for visibility)
  xnode = rep(c(-0.1, 2.23, 1), each = 3),         # X-coordinates for node labels
  ynode = rep(c(-0.1, -0.1, sqrt(3) + 0.15), each = 3), # Y-coordinates for node labels
  
  # Weight label positions (positioned near the middle of each edge)
  xweight = rep(c(0.3, 1, 1.7), each = 3),         # X-coordinates for weight labels
  yweight = rep(c(1, -0.2, 1), each = 3)            # Y-coordinates for weight labels
)
# =============================================================================
# FIGURE 1: ORIGINAL NETWORKS G1, G2, G3
# =============================================================================

# Create Figure 1 showing the three original weighted triangle networks
# This figure demonstrates the input data for Fréchet regression
ggplot(data = df) +
  # Draw edges as line segments with thickness proportional to weight
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, linewidth = weight / 2), 
               color = '#F8766D') +  # Red color for edges
  
  # Add node labels (v1, v2, v3) using LaTeX formatting
  geom_text(aes(x = xnode, y = ynode, label = TeX(labelnode, output = 'character')), 
            size = 8, parse = TRUE) +
  
  # Add weight labels on each edge
  geom_text(aes(x = xweight, y = yweight, label = labelweight), size = 6) +
  
  # Create separate panels for each network (G1, G2, G3)
  facet_wrap(vars(sample), nrow = 1, labeller = label_parsed, strip.position = 'bottom') +
  
  # Use line width as specified (no legend for line width)
  scale_linewidth_identity(guide = 'none') +
  
  # Set axis labels (these are placeholders for the actual plot)
  labs(x = 'Number of nodes m', y = 'Run time (s)') +
  
  # Set plot limits to ensure proper display
  xlim(-0.3, 2.3) +
  ylim(-0.3, sqrt(3) + 0.2) +
  
  # Use fixed aspect ratio to maintain triangle shape
  coord_fixed(ratio = 1) +
  
  # Remove default theme elements for clean appearance
  theme_void() +
  
  # Customize text size and strip positioning
  theme(text = element_text(size=25), 
        strip.text = element_text(margin = margin(0.2,0,0,0, "cm")))

# =============================================================================
# DATA SETUP FOR FIGURE 2: PREDICTED NETWORKS
# =============================================================================

# Create data frame for plotting the predicted networks at X = 2 and X = 5
# These predictions are obtained using Fréchet regression on the original data
# The predicted weights follow the same pattern as the original data but with
# the predicted covariate values: [X_pred, 2*X_pred, X_pred]
df <- data.frame(
  # Edge coordinates: same triangle structure as Figure 1
  x = rep(c(0, 2, 1), each = 2),                    # Starting x-coordinates of edges
  y = rep(c(0, 0, sqrt(3)), each = 2),              # Starting y-coordinates of edges
  xend = rep(c(1, 0, 2), each = 2),                 # Ending x-coordinates of edges
  yend = rep(c(sqrt(3), 0, 0), each = 2),           # Ending y-coordinates of edges
  
  # Predicted edge weights: [2, 5, 4, 10, 2, 5]
  # This represents predictions at X = 2 and X = 5 following the pattern [X, 2*X, X]
  weight = c(2, 5, 4, 10, 2, 5),                    # Predicted edge weights
  
  # Node labels for each predicted network
  labelnode = rep(c('$v_1$', '$v_2$', '$v_3$'), each = 2),
  
  # Weight labels (converted to fractions for better readability)
  labelweight = as.character(MASS::fractions(c(2, 5, 4, 10, 2, 5))),
  
  # Network labels for faceting (showing prediction values)
  sample = factor(rep(1:2, 3), levels = 1:2, 
                 labels = c('Prediction at X = 2', 'Prediction at X = 5')),
  
  # Node label positions (same as Figure 1)
  xnode = rep(c(-0.1, 2.23, 1), each = 2),         # X-coordinates for node labels
  ynode = rep(c(-0.1, -0.1, sqrt(3) + 0.15), each = 2), # Y-coordinates for node labels
  
  # Weight label positions (adjusted for better visibility)
  xweight = rep(c(0.2, 1, 1.8), each = 2),         # X-coordinates for weight labels
  yweight = rep(c(1, -0.2, 1), each = 2)            # Y-coordinates for weight labels
)
# =============================================================================
# FIGURE 2: PREDICTED NETWORKS AT X = 2 AND X = 5
# =============================================================================

# Create Figure 2 showing the predicted networks using Fréchet regression
# This figure demonstrates the output of Fréchet regression for new covariate values
ggplot(data = df) +
  # Draw edges as line segments with thickness proportional to predicted weight
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, linewidth = weight / 2), 
               color = '#619CFF') +  # Blue color for predicted edges
  
  # Add node labels (v1, v2, v3) using LaTeX formatting
  geom_text(aes(x = xnode, y = ynode, label = TeX(labelnode, output = 'character')), 
            size = 8, parse = TRUE) +
  
  # Add predicted weight labels on each edge
  geom_text(aes(x = xweight, y = yweight, label = labelweight), size = 6) +
  
  # Create separate panels for each prediction (X = 2, X = 5)
  facet_wrap(vars(sample), nrow = 1, strip.position = 'bottom') +
  
  # Use line width as specified (no legend for line width)
  scale_linewidth_identity(guide = 'none') +
  
  # Set axis labels (these are placeholders for the actual plot)
  labs(x = 'Number of nodes m', y = 'Run time (s)') +
  
  # Set plot limits to ensure proper display
  xlim(-0.3, 2.3) +
  ylim(-0.3, sqrt(3) + 0.2) +
  
  # Use fixed aspect ratio to maintain triangle shape
  coord_fixed(ratio = 1) +
  
  # Remove default theme elements for clean appearance
  theme_void() +
  
  # Customize text size and strip positioning (slightly different from Figure 1)
  theme(text = element_text(size=19), 
        strip.text = element_text(margin = margin(0.2,0,0.08,0, "cm")))
