# Functions for clustering

#_____________________________________________
# load functions
#_____________________________________________
# get mean intensities for each cluster
get_mean_intensity <- function(expr, cell_clustering){
  
  expr2 <- as.data.frame(expr)
  expr2$cluster <- cell_clustering[match(rownames(expr2),names(cell_clustering))]
  expr_mean <- aggregate(expr[,-dim(expr2)[2]], list(expr2$cluster), mean)
  
  colnames(expr_mean)[1] <- "cell_clustering"
  rownames(expr_mean) <- expr_mean$cell_clustering
  
  return(expr_mean)
}
# get median intensities for each cluster
get_median_intensity <- function(expr, cell_clustering){
  
  expr2 <- as.data.frame(expr)
  expr2$cluster <- cell_clustering[match(rownames(expr2),names(cell_clustering))]
  expr_median <- aggregate(expr[,-dim(expr2)[2]], list(expr2$cluster), median)
  
  colnames(expr_median)[1] <- "cell_clustering"
  rownames(expr_median) <- expr_median$cell_clustering
  
  return(expr_median)
}
# Merge and rename clusters assignment (character vectors)
merge_clusters_alpha = function(cl, maplist, show.warnings = TRUE){
  
  if(!checkmate::checkList(maplist, types = c("integer","character"))){
    stop("maplist must be a named list (character or integer)")
  }
  
  if(sum(duplicated(unlist(maplist)))>0){
    dup.target = unlist(maplist)[duplicated(unlist(maplist))]
    stop(paste(c("maplist is ambiguous:",dup.target),collapse = " "))
  }
  
  targets = unique(unlist(maplist))
  if(length(setdiff(cl,targets))>0 & show.warnings){
    message(paste(c("No mapping for:", setdiff(cl,targets)), collapse=" "))
  }
  
  levels.input = levels(as.factor(cl))
  levels.new = names(maplist)
  levels.merged = unlist(maplist)
  levels.output = sort(union(setdiff(levels.input, levels.merged),levels.new))
  
  cl2 = as.character(cl)
  for(group in names(maplist)){
    cl2[cl %in% maplist[[group]]] = group
  }
  
  factor(cl2,levels = levels.output)
}
# Returns data.frame with clust and dist. M can be cluster centers or cluster ids
get_nearest_center_with_dist = function(X, M, distance = "euclidean", exclude = NULL, include.next.dist=FALSE){
  
  if(is.null(dim(M)) & length(M) == nrow(X)){
    M = get_cluster_centers(X,M)
  }
  
  if(ncol(X) != ncol(M)){
    stop("X and M don't have the same number of rows")
  }
  
  if(is.null(rownames(M))){
    warning("rownames(M) == NULL. Using row number as name.")
    rownames(M) = 1:nrow(M)
  }
  
  if(!is.null(exclude)){
    M = M[!rownames(M) %in% exclude,]
  }
  
  DIST = get_dist_to_centers(X, M, distance = distance)
  UNDEF = (is.nan(DIST)|is.na(DIST))
  if(sum(UNDEF)>0){
    warning("Undefined distances were replaced by Inf")
  }
  DIST = replace(DIST,UNDEF,Inf)
  icol = max.col(-DIST)
  min.dist = matrixStats::rowMins(DIST)
  nc = rownames(M)[icol]
  df = data.frame(row.names = rownames(X),
                  cluster = nc,
                  dist = min.dist,
                  norm = get_row_norm(X),
                  stringsAsFactors = FALSE)
  
  if(include.next.dist){
    df$dist.next = matrixStats::rowOrderStats(DIST,which = 2)
    df$dist.next.ratio =  df$dist/df$dist.next
  }
  df
}
# Returns a matrix M of cluster centers (median or mean)
get_cluster_centers = function(X, group, fun.centers = mean){
  D = ncol(X)
  
  if(class(group)!="factor"){
    group = factor(group)
  }
  
  K = nlevels(group)
  
  M = matrix(nrow = K, ncol = D)
  for(d in 1:D){
    M[,d] = tapply(X = X[,d],INDEX = group, FUN = fun.centers)
  }
  colnames(M) = colnames(X)
  rownames(M) = levels(group)
  M
}
# Return the distances from each point to each center
get_dist_to_centers = function(X, M, distance = "euclidean"){
  
  if(ncol(X) != ncol(M)){
    stop("X and M don't have the same number of rows")
  }
  
  if(distance=="euclidean"){
    euclidean_dist2(X,M)
  }else{
    stop("Unrecognized distance. Valid values are: euclidean, cosine and correlation")
  }
}
# Euclidean dist between rows in X and rows in M
euclidean_dist2 = function(X, M){
  X = as.matrix(X)
  nbx = nrow(X) ; K = nrow(M)
  SSE = matrix(nrow = nbx, ncol = K)
  tX = t(X)
  for(k in 1:K){
    SSE[,k] = matrixStats::colSums2((tX - M[k,])^2) # fastest
  }
  colnames(SSE) = rownames(M)
  rownames(SSE) = rownames(X)
  sqrt(SSE)
}

# get row norm
get_row_norm = function(X) sqrt(rowSums(X^2))

# get sample IDs
# M is a list of cell IDs
get_sample_ids <- function(C){
  S <- sub(':.*$','',C)
  return(S)
}

# plot clustering heatamap
plot_clustering_heatmap <- function(expr, cell_clustering, nbcut=NA, annotation_row = NA, annotation_colors = NA, title = "",
                                    filename = NA, fun.center = mean, fontsize = 18, width = 12, height = 7, legend = FALSE,
                                    display_numbers = TRUE, annotation_legend = FALSE){
  
  require(pheatmap)
  require(RColorBrewer)
  
  expr2 <- as.data.frame(expr)
  expr2$cluster <- cell_clustering[match(rownames(expr2),names(cell_clustering))]
  expr_centers <- aggregate(expr[,-dim(expr2)[2]], list(expr2$cluster), fun.center)
  
  expr_heat <- as.matrix(expr_centers[, 2:ncol(expr_centers)])
  rownames(expr_heat) <- expr_centers$Group.1
  labels_col <- colnames(expr)
  labels_row <- paste0(names(table(cell_clustering)),' (',table(cell_clustering),')')
  
  if(is.na(filename)){
    pheatmap::pheatmap(expr_heat,cluster_cols = FALSE,labels_col = as.character(labels_col),
                       color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                       labels_row = labels_row, annotation_row = annotation_row, annotation_colors = annotation_colors,
                       display_numbers = display_numbers, number_color = 'black', fontsize_number = 12, angle_col = 45,
                       fontsize = fontsize, cutree_rows = nbcut, main = title, legend = legend, annotation_legend = annotation_legend)
  }else{
    pheatmap::pheatmap(expr_heat,cluster_cols = FALSE,labels_col = as.character(labels_col),
                       color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                       labels_row = labels_row, annotation_row = annotation_row, annotation_colors = annotation_colors,
                       display_numbers = display_numbers, number_color = "black", fontsize_number = 12, angle_col = 45,
                       fontsize = fontsize, cutree_rows = nbcut, main = title, legend = legend,
                       filename = filename, width = width, height = height, annotation_legend = annotation_legend)
  }
}
