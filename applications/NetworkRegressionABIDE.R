##############################################################################
##### Title: Non-Euclidean Data Analysis With Metric Statistics (HDSR) #######
##### Description: Fréchet Regression for Functional Connectiviy Networks#####
##### Manuscript reference : Section 2 #######################################
##### Figure reference : Figures 7&8&9 #######################################
##### Date: 09/18/25 #########################################################
##### Author: Yidong Zhou ####################################################
##############################################################################

# This script implements Fréchet regression analysis for functional connectivity 
# networks using the ABIDE dataset. The analysis compares network properties 
# between autism spectrum disorder (ASD) and control groups across different ages.
# 
# Main components:
# 1. Data loading and preprocessing from ABIDE dataset
# 2. Network construction from correlation matrices with thresholding
# 3. Fréchet regression analysis for both groups
# 4. Network metrics calculation (global efficiency, clustering coefficient)
# 5. Visualization of results (Figures 7, 8, 9)

# Load required libraries and functions
source('gnr.R')  # Contains the gnr() function for Fréchet regression
library(igraph)  # For network analysis and visualization
library(ggplot2) # For creating plots

# =============================================================================
# DATA LOADING AND PREPROCESSING
# =============================================================================

# Load region of interest (ROI) labels for the AAL atlas
roi <- read.csv('abide_aal/aal_labels.csv', skip = 1)
names(roi) <- c('id', 'name')  # Standardize column names

# Load phenotypic data containing subject demographics and diagnosis
pred <- read.csv('abide_aal/abide_phenotypic.csv')
str(pred)  # Display structure of the data
names(pred) <- c('id', 'site', 'age', 'sex', 'dx', 'dsm')  # Standardize column names

# To control for site-specific variations, we analyzed only data from the NYU site, 
# yielding 72 ASD and 95 control subjects, aged 6.47 to 39.1 years.

# Initialize variables for network analysis
gl <- list()      # List to store graph Laplacian matrices
theta <- 0.15     # Threshold parameter for network sparsification (top 15% of connections)
x <- NULL         # Matrix to store predictor variables (age, sex, diagnosis)
nm <- NULL        # Vector to store subject IDs
# =============================================================================
# NETWORK CONSTRUCTION FROM FUNCTIONAL CONNECTIVITY DATA
# =============================================================================

# Process each subject from the NYU site
for(id in pred$id[pred$site == 'NYU']) {
  # Load time series data for the current subject
  path <- paste0('abide_aal/aal/', id, '.csv')
  tsid <- read.csv(path, header = FALSE)
  
  # Compute correlation matrix from time series data
  crid <- suppressWarnings(cor(tsid))
  
  # Skip subjects with missing data (NA correlations)
  if(any(is.na(crid))) {
    next
  }
  
  # Store subject ID and phenotypic data
  nm <- c(nm, id)
  x <- rbind(x, pred[pred$id == id, c(3, 4, 5)])  # Extract age, sex, diagnosis
  
  # Network construction process:
  # 1. Remove self-connections (set diagonal to 0)
  diag(crid) <- 0
  
  # 2. Take absolute values of correlations
  crid <- abs(crid)
  
  # 3. Apply thresholding to keep only top 15% of connections
  threshold <- quantile(crid, 1 - theta)
  adjid <- ifelse(crid >= threshold, crid, 0)
  
  # 4. Convert adjacency matrix to graph Laplacian
  # Graph Laplacian: L = D - A, where D is degree matrix and A is adjacency matrix
  glid <- -adjid  # Negative adjacency matrix
  diag(glid) <- -colSums(glid)  # Set diagonal to negative row sums (degree)
  
  # Store the graph Laplacian
  gl[[length(gl) + 1]] <- glid
}
names(gl) <- nm  # Assign subject IDs as names for the graph list

# =============================================================================
# FRÉCHET REGRESSION ANALYSIS
# =============================================================================

# Separate ages for control and autism groups
# dx = 2: control group, dx = 1: autism group
agec <- x$age[x$dx == 2]  # Ages for control subjects
agea <- x$age[x$dx == 1]  # Ages for autism subjects

# Define prediction ages for regression analysis
xOut <- c(10, 15, 20, 25)  # Ages at which to predict network properties

# Perform Fréchet regression for both groups
res <- list()
res[[1]] <- gnr(gl[x$dx == 2], agec, xOut)  # Fréchet regression for control group
res[[2]] <- gnr(gl[x$dx == 1], agea, xOut)  # Fréchet regression for autism group

# =============================================================================
# NETWORK METRICS CALCULATION - CONTROL GROUP
# =============================================================================

# Process fitted networks for control group to calculate network metrics
# Convert graph Laplacians back to adjacency matrices for metric calculation
adjFit <- lapply(res[[1]]$fit, function(gl){
  diag(gl) <- 0      # Remove diagonal elements
  gl <- -gl          # Convert back to adjacency matrix (negative of Laplacian)
  threshold <- quantile(gl, 1 - theta)  # Apply same thresholding as original networks
  gl <- ifelse(gl >= threshold, gl, 0)  # Keep only strong connections
  gl
})

# Initialize vectors for network metrics
gec <- rep.int(0, length(adjFit))  # Global efficiency for control group
ccc <- rep.int(0, length(adjFit))  # Clustering coefficient for control group

# Calculate network metrics for each predicted age
for(i in seq_along(adjFit)) {
  print(i)  # Progress indicator
  
  # Create igraph object from adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(adjFit[[i]], mode = 'undirected', weighted = TRUE)
  
  # Calculate global efficiency (inverse of average shortest path length)
  gec[i] <- brainGraph::efficiency(g, type = 'global')
  
  # Calculate clustering coefficient (local connectivity measure)
  ccc[i] <- NetworkToolbox::clustcoeff(adjFit[[i]], weighted = TRUE)$CC
}

# =============================================================================
# NETWORK METRICS CALCULATION - AUTISM GROUP
# =============================================================================

# Process fitted networks for autism group to calculate network metrics
# Convert graph Laplacians back to adjacency matrices for metric calculation
adjFit <- lapply(res[[2]]$fit, function(gl){
  diag(gl) <- 0      # Remove diagonal elements
  gl <- -gl          # Convert back to adjacency matrix (negative of Laplacian)
  threshold <- quantile(gl, 1 - theta)  # Apply same thresholding as original networks
  gl <- ifelse(gl >= threshold, gl, 0)  # Keep only strong connections
  gl
})

# Initialize vectors for network metrics
gea <- rep.int(0, length(adjFit))  # Global efficiency for autism group
cca <- rep.int(0, length(adjFit))  # Clustering coefficient for autism group

# Calculate network metrics for each predicted age
for(i in seq_along(adjFit)) {
  print(i)  # Progress indicator
  
  # Create igraph object from adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(adjFit[[i]], mode = 'undirected', weighted = TRUE)
  
  # Calculate global efficiency (inverse of average shortest path length)
  gea[i] <- brainGraph::efficiency(g, type = 'global')
  
  # Calculate clustering coefficient (local connectivity measure)
  cca[i] <- NetworkToolbox::clustcoeff(adjFit[[i]], weighted = TRUE)$CC
}

# =============================================================================
# VISUALIZATION - FIGURE 9: NETWORK METRICS BY AGE AND GROUP
# =============================================================================

# Create Figure 9: Scatter plot of network metrics vs age for both groups
# This figure shows how global efficiency and clustering coefficient change with age
# for both control and autism groups

ggplot(data = data.frame(
  # Combine network metrics for both groups
  nm = c(gec[agec >= 5 & agec <= 25],           # Control global efficiency (ages 5-25)
         ccc[agec >= 5 & agec <= 25],           # Control clustering coefficient (ages 5-25)
         gea[agea >= 5 & agea <= 25] - 0.2,     # Autism global efficiency (ages 5-25, offset for visibility)
         cca[agea >= 5 & agea <= 25]),          # Autism clustering coefficient (ages 5-25)
  
  # Corresponding ages for each metric
  age = c(rep(agec[agec >= 5 & agec <= 25], 2), # Control ages (repeated for both metrics)
          rep(agea[agea >= 5 & agea <= 25], 2)), # Autism ages (repeated for both metrics)
  
  # Metric type labels
  cs = c(rep(c('Global efficiency', 'Clustering coefficient'), each = sum(agec >= 5 & agec <= 25)), # Control metrics
         rep(c('Global efficiency', 'Clustering coefficient'), each = sum(agea >= 5 & agea <= 25))), # Autism metrics
  
  # Group labels
  group = rep(c('Control', 'Autism'), 2 * c(sum(agec >= 5 & agec <= 25), sum(agea >= 5 & agea <= 25)))
), 
aes(x = age, y = nm, group = cs)) +
  geom_point() +  # Scatter plot points
  scale_x_continuous(name = 'Age (years)') +  # X-axis label
  scale_y_continuous(name = '') +  # Y-axis label (empty for facet labels)
  facet_grid(rows = vars(cs), cols = vars(group), scales = "free_y") +  # Create 2x2 grid
  theme_bw() +  # Black and white theme
  theme(text = element_text(size=20))  # Larger text size for publication

# =============================================================================
# VISUALIZATION - FIGURE 8: CONTROL GROUP NETWORK VISUALIZATION
# =============================================================================

# Prepare predicted networks for control group visualization
# Convert graph Laplacians back to adjacency matrices for visualization
adjPred <- lapply(res[[1]]$predict, function(gl){
  diag(gl) <- 0  # Remove diagonal elements
  gl <- -gl      # Convert back to adjacency matrix
  gl
})

# Perform community detection for each predicted age
cl <- list()
for(k in 1:4){
  # Create igraph object from adjacency matrix
  g <- graph_from_adjacency_matrix(adjPred[[k]], mode = 'undirected', weighted = TRUE)
  # Use leading eigenvector method for community detection
  cl[[k]] <- cluster_leading_eigen(g)
}

cl[[1]]$membership[cl[[1]]$membership == 5] <- 2
cl[[1]]$membership[cl[[1]]$membership == 6] <- 5
cl[[3]]$membership[cl[[3]]$membership == 4 | cl[[3]]$membership == 5] <- 3
cl[[3]]$membership[cl[[3]]$membership == 6] <- 4
cl[[4]]$membership[cl[[4]]$membership == 4 | cl[[4]]$membership == 5] <- 3
cl[[4]]$membership[cl[[4]]$membership == 6] <- 4

# Define number of clusters and colors for each age
noClusters <- c(5, 5, 4, 4)  # Number of communities at each age
col <- list(c('#e6194b', '#f58231', '#4363d8', '#911eb4', '#bcf60c'),  # Age 10 colors
            c('#e6194b', '#f58231', '#4363d8', '#911eb4', '#bcf60c'),  # Age 15 colors
            c('#e6194b', '#f58231', '#4363d8', '#911eb4'),             # Age 20 colors
            c('#e6194b', '#f58231', '#4363d8', '#911eb4'))             # Age 25 colors

# Create Figure 8: Control group network visualization at different ages
# Set up 2x2 plot layout for four different ages
layout(mat = matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE))
opar <- par(mar = c(0.5, 0, 1, 0.5))  # Set margins for subplots
set.seed(3)  # Set seed for reproducible layout
title <- paste0('Age = ', xOut)  # Create titles for each subplot

# Plot networks for each predicted age
for(k in 1:4){
  # Create igraph object from adjacency matrix
  g <- graph_from_adjacency_matrix(adjPred[[k]], mode = 'undirected', weighted = TRUE)
  
  # Calculate node strength (sum of edge weights)
  deg <- igraph::strength(g)
  
  # Set node names to ROI names
  V(g)$name <- roi$name
  
  # Scale edge width based on weight (relative to maximum weight)
  E(g)$width <- 0.9 * E(g)$weight/max(E(g)$weight)
  
  # Use Fruchterman-Reingold layout algorithm
  l <- layout.fruchterman.reingold(g)
  
  # Assign community membership to nodes
  V(g)$community <- cl[[k]]$membership
  V(g)$size <- 6  # Set node size
  
  # Color edges based on community membership
  # Within-community edges: colored with community color (transparent)
  # Between-community edges: black (transparent)
  E(g)$color <- apply(as.data.frame(as_edgelist(g)), 1, function(x) 
    ifelse(V(g)$community[which(roi$name == x[1])] == V(g)$community[which(roi$name == x[2])], 
           adjustcolor(col[[k]][V(g)$community[which(roi$name == x[1])]], alpha.f = 0.2), 
           adjustcolor('#000000', alpha.f = 0.1)))
  
  # Plot the network
  plot(g, layout = l, edge.color = E(g)$color,
       mark.groups = lapply(seq_along(unique(V(g)$community)), function(i) roi$name[V(g)$community == i]),
       vertex.color = col[[k]][V(g)$community], vertex.size = 4, vertex.frame.color = NA, vertex.label = NA, 
       mark.col = NA, mark.border = NA)
  
  # Add title with age and number of communities
  title(paste0(title[k], ', ', noClusters[k], ' communities'), cex.main = 1)
}
par(opar)  # Reset plotting parameters

# =============================================================================
# VISUALIZATION - FIGURE 7: AUTISM GROUP NETWORK VISUALIZATION
# =============================================================================

# Prepare predicted networks for autism group visualization
# Convert graph Laplacians back to adjacency matrices for visualization
adjPred <- lapply(res[[2]]$predict, function(gl){
  diag(gl) <- 0  # Remove diagonal elements
  gl <- -gl      # Convert back to adjacency matrix
  gl
})

# Perform community detection for each predicted age
cl <- list()
for(k in 1:4){
  # Create igraph object from adjacency matrix
  g <- graph_from_adjacency_matrix(adjPred[[k]], mode = 'undirected', weighted = TRUE)
  # Use leading eigenvector method for community detection
  cl[[k]] <- cluster_leading_eigen(g)
}

# Display community membership distribution for each age
lapply(cl, function(cli) table(cli$membership))

# Define number of clusters and colors for each age (autism group)
noClusters <- c(3, 4, 4, 5)  # Number of communities at each age
col <- list(c('#e6194b', '#f58231', '#4363d8'),                    # Age 10 colors
            c('#e6194b', '#ffe119', '#4363d8', '#f58231'),         # Age 15 colors
            c('#e6194b', '#ffe119', '#4363d8', '#f58231'),         # Age 20 colors
            c('#4363d8', '#ffe119', '#e6194b', '#3cb44b', '#f58231')) # Age 25 colors

# Create Figure 7: Autism group network visualization at different ages
# Set up 2x2 plot layout for four different ages
layout(mat = matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE))
opar <- par(mar = c(0.5, 0, 1, 0.5))  # Set margins for subplots
set.seed(3)  # Set seed for reproducible layout
title <- paste0('Age = ', xOut)  # Create titles for each subplot

# Plot networks for each predicted age
for(k in 1:4){
  # Create igraph object from adjacency matrix
  g <- graph_from_adjacency_matrix(adjPred[[k]], mode = 'undirected', weighted = TRUE)
  
  # Calculate node strength (sum of edge weights)
  deg <- igraph::strength(g)
  
  # Set node names to ROI names
  V(g)$name <- roi$name
  
  # Scale edge width based on weight (relative to maximum weight)
  E(g)$width <- 0.9 * E(g)$weight/max(E(g)$weight)
  
  # Use Fruchterman-Reingold layout algorithm
  l <- layout.fruchterman.reingold(g)
  
  # Assign community membership to nodes
  V(g)$community <- cl[[k]]$membership
  V(g)$size <- 6  # Set node size
  
  # Color edges based on community membership
  # Within-community edges: colored with community color (transparent)
  # Between-community edges: black (transparent)
  E(g)$color <- apply(as.data.frame(as_edgelist(g)), 1, function(x) 
    ifelse(V(g)$community[which(roi$name == x[1])] == V(g)$community[which(roi$name == x[2])], 
           adjustcolor(col[[k]][V(g)$community[which(roi$name == x[1])]], alpha.f = 0.2), 
           adjustcolor('#000000', alpha.f = 0.1)))
  
  # Plot the network
  plot(g, layout = l, edge.color = E(g)$color,
       mark.groups = lapply(seq_along(unique(V(g)$community)), function(i) roi$name[V(g)$community == i]),
       vertex.color = col[[k]][V(g)$community], vertex.size = 4, vertex.frame.color = NA, vertex.label = NA, 
       mark.col = NA, mark.border = NA)
  
  # Add title with age and number of communities
  title(paste0(title[k], ', ', noClusters[k], ' communities'), cex.main = 1)
}
par(opar)  # Reset plotting parameters
