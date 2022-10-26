# Main functions for spatial analysis of cell data
# Author: Mathieu Lajoie (mathieu.lajoie2@mcgill.ca)
# Date: January 1, 2022

# Public functions ####
NG_create_dt = function(x, y, obj.id, sample.id = NA, group = NA, max.dist = 30, nb.obj.max = 2e4){
  
  # check input
  checkmate::assertNumeric(x, finite = TRUE, any.missing = FALSE, min.len = 1)
  checkmate::assertNumeric(y, finite = TRUE, any.missing = FALSE, len = length(x))
  checkmate::assertCharacter(obj.id, len = length(x), unique = TRUE, any.missing = FALSE)
  checkmate::assertNumeric(max.dist, finite = TRUE, lower = 0, len = 1)
  checkmate::assertNumeric(nb.obj.max, finite = TRUE, lower = 0, len = 1)
  
  D = data.frame(x, y, obj.id, group, sample.id)
  
  spD = split.data.frame(D, D$sample.id)
  LNG = list()
  
  for(sel.sample in names(spD)){
    print(paste("Creating neighborhood graph for:", sel.sample))
    d = spD[[sel.sample]]
    ng = create_neighborhood_graph_dt(
      x = d$x,
      y = d$y,
      obj.id = d$obj.id,
      sample.id = sel.sample,
      max.dist =  max.dist,
      nb.obj.max = nb.obj.max)
    ng$group = d$group
    
    LNG[[sel.sample]] = ng
  }
  
  NG = LNG_merge(LNG)
  NG = NG[match(obj.id,NG$obj.id),]  # Reorder as input
  NG = NG_update_index(NG)
  NG
}

# Split a multi-sample NG into a list of single sample NGs
NG_split = function(NG, f = NG$sample.id){
  
  NG_check(NG)
  checkmate::assertVector(f, any.missing = FALSE, len = nrow(NG))
  message("splitting NG object...")
  LNG = split.data.frame(NG, f)
  message("updating NG objects...")
  for(i in seq_len(length(LNG))){
    LNG[[i]] <- NG_update_index(LNG[[i]])
  }
  LNG
}

# Merge a list of single sample NGs into a multi-sample NG
LNG_merge = function(LNG){
  
  checkmate::assertList(LNG,"DataFrame")
  checkmate::assertTRUE(all(unlist(lapply(LNG,NG_check))))
  
  UNG = do.call(rbind, LNG)
  
  checkmate::assertInteger(UNG$index)
  checkmate::assertCharacter(UNG$sample.id)
  checkmate::assertClass(UNG$nn,"IntegerList")
  checkmate::assertClass(UNG$dist,"NumericList")
  checkmate::assertCharacter(UNG$obj.id,unique = TRUE)
  
  # Update index
  tb = table(UNG$sample.id)
  tb = table(UNG$sample.id)[as.character(unique(UNG$sample.id))]
  cs = c(0,as.integer(cumsum(tb)))[-(length(tb)+1)]
  index.offset = as.integer(cs)[match(UNG$sample.id, names(tb))]
  UNG$index = UNG$index + index.offset
  UNG$nn = UNG$nn + index.offset
  names(UNG$nn) <- UNG$index
  names(UNG$dist) <- UNG$index
  NG_check(UNG)
  UNG
}

# Estimate contact enrichments (X <-> Y) using the permutation approach of Keren at al. (2018) PMID: 30193111
# Specific groups can be held fixed using the fixed.groups parameter.
# Value: A data.frame with z-scores, empirical p-values and other information.
NG_contact_enrichment = function(NG, vgroup = NULL, fixed.groups = NULL, nb.rand = 1e3){
  
  if(!is.null(vgroup)){ # update group
    checkmate::assertVector(vgroup, len = nrow(NG), any.missing = FALSE)
    NG$group <- vgroup
  }
  
  NG$group <- as.factor(NG$group) # use factor to get consistent table orderings
  
  if(!is.null(fixed.groups)){
    unknown.group = setdiff(fixed.groups, levels(NG$group))
    if(length(unknown.group)>0){
      stop(paste("fixed.groups entrie(s) not in group:\n",paste(unknown.group,collapse = " ")))
    }
  }
  
  # Deal with multi-sample using recursive call #
  if(length(unique(NG$sample.id)) > 1){
    LNG <- NG_split(NG)
    res = do.call(rbind, lapply(LNG, NG_contact_enrichment,
                                fixed.groups = fixed.groups,
                                nb.rand = nb.rand))
    rownames(res) <- NULL
    return(res)
  }
  
  # Start single sample
  sample.id = NG$sample.id[[1]]
  message(sample.id)
  NG_check(NG)
  
  # Indexes
  u.nn.int = unlist(NG$nn, use.names = TRUE)
  u.center.int = as.integer(names(u.nn.int))
  
  # Object types
  obj.type = NG$group
  center.type = obj.type[u.center.int]
  nn.type = obj.type[u.nn.int]
  
  # Observed contacts
  tb = table(center.type, nn.type, exclude = NULL)
  tb = tb - (diag(tb)/2 * diag(nrow = nrow(tb), ncol = ncol(tb))) # Remove duplicates on diagonal
  
  # Connectivity with non fixed groups (used to identify undefined scores)
  nb.free.nei = rowSums(tb[,! colnames(tb) %in% fixed.groups])
  
  # Permutations
  tb.rand = NULL
  if(!is.null(fixed.groups)){
    tb.rand = replicate(nb.rand, expr = {
      obj.type.shuffled = partial_shuffle(obj.type, keep = fixed.groups)
      center.type.shuffled = obj.type.shuffled[u.center.int]
      nn.type.shuffled = obj.type.shuffled[u.nn.int]
      tbs = table(center.type.shuffled, nn.type.shuffled, exclude = NULL)
      tbs = tbs - (diag(tbs)/2 * diag(nrow = nrow(tbs), ncol = ncol(tbs))) # Remove duplicates on diagonal
    })}else{
      tb.rand = replicate(nb.rand, expr = {
        obj.type.shuffled = sample(obj.type, replace = FALSE)
        center.type.shuffled = obj.type.shuffled[u.center.int]
        nn.type.shuffled = obj.type.shuffled[u.nn.int]
        tbs = table(center.type.shuffled, nn.type.shuffled, exclude = NULL)
        tbs = tbs - (diag(tbs)/2 * diag(nrow = nrow(tbs), ncol = ncol(tbs))) # Removes duplicates on diagonal
      })
    }
  
  # tb.rand dims : nb.group x nb.group x nb.rand
  tb.rep <- rep(tb, nb.rand) # repeat table of observed values (for vectorized operations below)
  attributes(tb.rep) <- attributes(tb.rand)
  
  # permute dimensions to have replicates in dim 1
  # tb.rand.perm dims : nb.rand x nb.group x nb.group
  tb.rand.perm <- aperm(tb.rand, c(3,1,2))
  tb.rep.perm <- aperm(tb.rep, c(3,1,2))
  
  # Two-pass SD:
  tb.mean.rand <- colMeans(tb.rand.perm, dims = 1)
  tb.mean.rand.rep = rep(tb.mean.rand, nb.rand)
  attributes(tb.mean.rand.rep) <- attributes(tb.rand)
  tb.mean.rand.rep.perm = aperm(tb.mean.rand.rep, c(3,1,2))
  
  tb.sd.rand = sqrt(colSums((tb.rand.perm - tb.mean.rand.rep.perm)^2,dims = 1)/(nb.rand - 1))
  tb.z.score <- (tb - tb.mean.rand) / tb.sd.rand
  
  pc.pval = 1 # Pseudocounts for p-values
  tb.pval.pos <- (colSums(tb.rep.perm <= tb.rand.perm, dims = 1) + pc.pval)/(nb.rand + pc.pval)
  tb.pval.neg <- (colSums(tb.rep.perm >= tb.rand.perm, dims = 1) + pc.pval)/(nb.rand + pc.pval)
  
  # Return as a data.frame
  melt.counts <- reshape2::melt(tb)
  group.i = melt.counts$center.type
  group.j = melt.counts$nn.type
  
  group.counts = table(NG$group)
  Ni = as.numeric(group.counts[match(as.character(group.i),names(group.counts))])
  Nj = as.numeric(group.counts[match(as.character(group.j),names(group.counts))])
  
  dfm <- data.frame(
    sample.id,
    group.i,
    group.j,
    Ni,
    Nj,
    counts = melt.counts$value,
    mean.rand = reshape2::melt(tb.mean.rand)$value,
    sd.rand =   reshape2::melt(tb.sd.rand)$value,
    z.score =   reshape2::melt(tb.z.score)$value,
    pval.pos =  reshape2::melt(tb.pval.pos)$value,
    pval.neg =  reshape2::melt(tb.pval.neg)$value,
    pair = paste0(group.i, ">", group.j),
    nb.rand = nb.rand,
    stringsAsFactors = FALSE
  )
  
  # Update z-scores ##
  
  # The z-scores calculated above will not be finite if no contact are observed in the permutations for a given pair.
  # In some cases, contacts are simply impossible due to the absence of cells or connectivity in the graph.
  # In other cases, contacts are possible and a finite z-score could be obtained by increasing the number of permutations (nb.rand).
  # The following section identifies these cases and update the z-score when possible (latter cases).
  
  dfm$same.group = dfm$group.i == dfm$group.j
  dfm$fixed.i = dfm$group.i %in% fixed.groups
  dfm$fixed.j = dfm$group.j %in% fixed.groups
  
  dfm$nb.free.nei.i = nb.free.nei[as.character(dfm$group.i)]
  dfm$nb.free.nei.j = nb.free.nei[as.character(dfm$group.j)]
  
  dfm$case.unique.cell = dfm$Ni == 1 & dfm$Nj == 1 & (dfm$group.i == dfm$group.j)
  dfm$case.no.cell = pmin(dfm$Ni, dfm$Nj) == 0
  dfm$case.no.nei = (dfm$fixed.i & dfm$nb.free.nei.i == 0) | (dfm$fixed.j & dfm$nb.free.nei.j == 0)
  dfm$case.both.fixed = dfm$fixed.i & dfm$fixed.j
  
  dfm$valid.score = TRUE
  dfm$valid.score[ (dfm$case.unique.cell | dfm$case.no.cell | dfm$case.no.nei | dfm$case.both.fixed) ] <- FALSE
  
  # Now we update undefined but valid scores
  dfm$z.score.original = dfm$z.score
  
  # Case with counts > 0. Positive contact enrichment (counts/0 = Inf)
  # New score obtained with 1 pseudo-count
  r = c(1, rep(0, nb.rand-1))
  pseudo.rand.sd = sqrt(sum((r-mean(r))^2)/(nb.rand-1))
  pseudo.rand.mean = 1 / nb.rand
  dfm$updated.pc = dfm$valid.score & (dfm$mean.rand == 0) & (dfm$counts > 0)
  dfm$z.score[dfm$updated.pc] = (dfm$counts[dfm$updated.pc] - pseudo.rand.mean) / pseudo.rand.sd
  
  # Case with counts == 0. We can't say if enrichment is positive or negative (0/0 = NaN)
  # New score is set to 0
  dfm$updated.0 = dfm$valid.score & (dfm$mean.rand == 0) & (dfm$counts == 0)
  dfm$z.score[dfm$updated.0] = 0
  
  # Set invalid scores and p-values to NA
  dfm$z.score[!dfm$valid.score] <- NA
  dfm$pval.pos[!dfm$valid.score] <- NA
  dfm$pval.neg[!dfm$valid.score] <- NA
  
  # Last check
  checkmate::assertNumeric(dfm$z.score[dfm$valid.score], finite = TRUE, any.missing = FALSE, min.len = 0)
  
  # add info
  dfm$info = "undefined"
  dfm$info[dfm$valid.score] = ""
  dfm$info[dfm$updated.0|dfm$updated.pc] = "updated"
  
  # subset relevant columns
  sel.cols = c("sample.id","group.i","group.j","Ni","Nj","counts","mean.rand","sd.rand","z.score","pval.pos","pval.neg","pair","nb.rand","info")
  dfm = dfm[,sel.cols]
  dfm
}

# Estimate asymmetrical contact enrichment (X -> Y) for each group, where group X is fixed.
NG_contact_enrichment_AS = function(NG, vgroup = NULL, fixed.groups = NULL, nb.rand = 1e3){
  
  if(is.null(vgroup)){
    checkmate::assertVector(NG$group, len = nrow(NG), any.missing = FALSE)
  }else{
    checkmate::assertVector(vgroup, len = nrow(NG), any.missing = FALSE)
    NG$group <- as.factor(vgroup)
  }
  
  NG$group = as.factor(NG$group)
  groups = levels(NG$group) 
  
  if(!is.null(fixed.groups)){
    checkmate::assertVector(fixed.groups)
  }
  
  lres = lapply(groups, function(g){
    print(paste("Processing",g))
    ce = NG_contact_enrichment(NG, fixed.group = unique(c(g,fixed.groups)), nb.rand = nb.rand)
    ce$fixed.group = g
    ce = ce[ce$group.i == g,]
    ce
  })
  
  dfres = do.call(rbind,lres)
  has.fixed.group = dfres$fixed.group == dfres$group.i # IMPORTANT FILTER
  has.free.group = !dfres$both.fixed
  dfres$keep = has.free.group & has.fixed.group
  dfres = dfres[dfres$keep,]
  dfres$free.group = dfres$group.j
  dfres$log2R = log2((dfres$counts + 1e-3) / (dfres$mean.rand + 1e-3))
  
  dfres$key = paste0(dfres$group.i,">",dfres$group.j," ",dfres$sample.id)
  dfres$rel = paste0(dfres$fixed.group,"->",dfres$free.group) # group1 -> group2 labels
  dfres
}

# Returns distance to nearest cell (bird flight) for each group
NG_nn_dist = function(NG, vgroup = NG$group, median.by.sample = TRUE){
  
  checkmate::assertVector(vgroup, len = nrow(NG))
  NG_check(NG)
  
  xNG = data.frame(x = NG$coord.x, y = NG$coord.y, 
                   sample.id = as.character(NG$sample.id), 
                   group = as.character(vgroup), 
                   stringsAsFactors = FALSE)
  
  LNG = split(xNG, xNG$sample.id) 
  
  # loop through samples
  dfres = do.call(rbind, lapply(names(LNG), function(sel.sample){
    # message(sel.sample)
    ng = LNG[[sel.sample]]
    lng = split.data.frame(ng, ng$group)
    
    # loop through target cell groups
    do.call(rbind, lapply(names(lng), function(sel.group){
      ng_group = lng[[sel.group]]
      ng$nn.group = sel.group
      ng$nn.dist = NA_real_
      
      # When a data node is also a query node, we need the second nearest neighbour 
      if(nrow(ng_group)>=2){ #there are at least two cells of type sel.group
        knn = RANN::nn2(data = ng_group[,c("x","y")], query = ng[,c("x","y")], 2, eps = 0)
        ng$nn.dist = knn$nn.dists[,1]
        idx  = ng$group == sel.group # idx of self edges (i.e. group(x) == group(y))
        ng$nn.dist[idx] <- knn$nn.dists[idx,2] # replace self by next neighbour
      }else if(nrow(ng_group)==1){ # special case when there is a single cell of type sel.group
        knn = RANN::nn2(data = ng_group[,c("x","y")], query = ng[,c("x","y")], 1, eps = 0)
        ng$nn.dist = knn$nn.dists[,1]
        ng$nn.dist[ng$group == sel.group] <- NA_real_ # No nearest neighbour for sel.group
      }
      ng
    }))
  }))
  rownames(dfres) <- NULL
  
  dt = data.table::as.data.table(dfres)
  
  # Use different variable names for consistency
  dt = dt[, .(sample.id, group.i = group, group.j = nn.group, nn.dist)]
  
  if(median.by.sample){ # we expand results here to include every possible pairs of groups
    dt = dt[, list(median.nn.dist = stats::median(nn.dist, na.rm = TRUE)), by = list(group.i, group.j, sample.id)]
    dt = data.table::as.data.table(tidyr::complete(dt, tidyr::expand(dt, group.i, group.j, sample.id)))
    dt$metric = "median.nn.dist"
    dt$value = dt$median.nn.dist
    dt$pair = paste(dt$group.i, dt$group.j)
    dt$cell.group = dt$pair
  }
  
  # Add cell counts
  tb = table(sample.id = xNG$sample.id, group = xNG$group)
  tbc = as.data.frame(tb) # cell counts by group
  tbp = as.data.frame(prop.table(tb, margin = 1)) # cell proportions by group
  key = paste(tbc$sample.id, tbc$group) 
  match.i = match(paste(dt$sample.id, dt$group.i), key)
  match.j = match(paste(dt$sample.id, dt$group.j), key)
  dt$Ni = tbc$Freq[match.i]
  dt$Nj = tbc$Freq[match.j]
  dt$Pi = tbp$Freq[match.i]
  dt$Pj = tbp$Freq[match.j]
  dt
}

# Compute assortativity for each group in each sample (considering all other groups as a single one)
NG_assortativity_by_single_group = function(NG, vgroup = NULL, split.samples = TRUE, merge.samples = FALSE){
  
  NG_check(NG)
  
  if(!is.null(vgroup)){
    checkmate::assertVector(vgroup, len = nrow(NG))
    NG$group = vgroup
  }
  
  if(merge.samples){
    NG$sample.id = "MERGED"
  }
  
  # Recursive call for multiple samples
  checkmate::assertLogical(split.samples, len = 1)
  
  if(split.samples == TRUE){
    LNG = NG_split(NG, NG$sample.id)
    lres = lapply(LNG, NG_assortativity_by_single_group, split.samples = FALSE)
    dfres = do.call(rbind, lres)
    return(dfres)
  }else{
    checkmate::assertTRUE(length(unique(NG$sample.id))==1)
  }
  
  if(nrow(NG)<=1 | length(unique(NG$group))<=1){
    warning("Empty NG or single group in NG_assortativity")
    return(data.frame(r=NA, se=NA))
  }
  
  # cell info
  tbc = table(NG$group, exclude = NULL)
  
  # get edges list ###
  u.nei = unlist(NG$nn, use.names = TRUE)
  de = data.frame(id1 = as.integer(names(u.nei)), id2 = u.nei, stringsAsFactors = FALSE)
  demap1 = match(de$id1, NG$index)
  demap2 = match(de$id2, NG$index)
  de$group1 = NG$group[demap1]
  de$group2 = NG$group[demap2]
  groups = unique(NG$group)
  
  # Table of undirected edge counts
  tbe = table(de$group1, de$group2, exclude = NULL)/2; tbe
  tbe
  
  groups = colnames(tbe)
  mdiag = diag(tbe)
  rsum = rowSums(tbe)
  msum = sum(tbe)
  
  dfres = do.call(rbind, lapply(1:length(groups), function(i){
    sel.group = groups[i]
    n.xx = mdiag[[i]]
    n.yy = msum - 2*rsum[[i]] + n.xx
    n.xy = 2*rsum[[i]] - 2*n.xx
    checkmate::assertTRUE(sum(n.xx,n.yy,n.xy) == msum) # sanity check
    r = assortativity_binary(n.xx, n.yy, n.xy)
    r.se = assortativity_binary_se_jacknife(n.xx, n.yy, n.xy)
    data.frame(group = sel.group, assort.r = r, assort.se = r.se, stringsAsFactors = FALSE)
  }))
  
  checkmate::assertTRUE(all(names(tbc) == dfres$group))
  dfres = cbind(sample.id = NG$sample.id[1], dfres)
  dfres$nb.cell.group = as.numeric(tbc)
  dfres$nb.cell.total = sum(tbc)
  dfres$cell.prop = dfres$nb.cell.group/dfres$nb.cell.total
  rownames(dfres) = NULL
  dfres
}

# Plot a single spatial graph object
NG_plot = function(ng, group = ng$group, pal.node, node.cex = 0.5, same.edge.size = 0.15, inter.edge.size = 0.03, xlim = NULL, ylim = NULL){
  
  NG_check(ng)
  checkmate::assertCharacter(unique(ng$sample.id),len = 1)
  checkmate::assertFactor(group, len = length(ng$index),any.missing = FALSE)
  ng$group = group
  
  checkmate::assertCharacter(pal.node, names = "named")
  
  if(is.null(xlim)){
    xlim = c(0,max(ng$coord.x))
  }
  if(is.null(ylim)){
    ylim = c(0,max(ng$coord.y))
  }
  
  checkmate::assertNumeric(xlim, any.missing = FALSE, len = 2)
  checkmate::assertNumeric(ylim, any.missing = FALSE, len = 2)
  
  df = as.data.frame(ng[,c("coord.x","coord.y","group")])
  
  df$`Cell type` = df$group
  de = get_edges_list(ng)
  de$edge.type = ifelse(de$group1 == de$group2, as.character(de$group2), "inter")
  
  pal.node = pal.node[levels(ng$group)]
  pal.edge = pal.node
  pal.edge["inter"] <- "black"
  pal.edge["highlight"] <- "black"
  
  size.edge = setNames(rep(same.edge.size, length(pal.edge)), names(pal.edge))
  size.edge["inter"] <- inter.edge.size
  
  gg.spg = ggplot() + 
    geom_segment(data = de, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = edge.type, size = edge.type)) +
    scale_color_manual(values = pal.edge) + 
    scale_size_manual(values = size.edge) +
    geom_point(df, mapping = aes(x = coord.x, y = coord.y, fill = `Cell type`, colour = `Cell type`), shape = 21, stroke = 0, cex = node.cex) + 
    scale_fill_manual(values = pal.node) + 
    guides(size = "none", shape = "none", colour = "none", linetype = "none") +
    theme_classic(base_size = 8) + 
    xlab(NULL) + ylab(NULL) + 
    theme(legend.position = "right", legend.justification = "top", legend.text = element_text(size = 7), legend.title = element_text(size = 8)) +
    theme(axis.line = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.ticks = element_blank()) + 
    theme(legend.key.size = unit(0, "cm")) +      
    theme(plot.margin = unit(c(0,0,0,0), "lines")) +
    theme(aspect.ratio = 1) + coord_fixed( xlim = xlim, ylim = rev(ylim) ) + 
    guides(fill = guide_legend(override.aes = list(size = 2.5))) 
    #legend.spacing = unit(0.01, "cm")
  gg.spg
}

# Helper functions ####

# Create neighborhood_graph using Delaunay triangulation for a single sample
create_neighborhood_graph_dt = function(x, y, obj.id, sample.id = NA, max.dist = 30, nb.obj.max = 2e4){
  
  # check input
  checkmate::assertNumeric(x)
  checkmate::assertNumeric(y, len = length(x))
  checkmate::assertCharacter(obj.id, len = length(x))
  checkmate::assertCharacter(sample.id, len = 1, any.missing = FALSE)
  
  if(length(x) > nb.obj.max){
    stop(paste("input size exceeded in one sample \"nb.obj.max\":",length(x), ">", nb.obj.max))
  }
  
  NG = S4Vectors::DataFrame(index = 1:length(x))
  
  spNDT = IRanges::IntegerList(tripack::neighbours(tripack::tri.mesh(x,y)))
  names(spNDT) = 1:length(spNDT)
  u.nei = unlist(spNDT,use.names = TRUE)
  u.node = as.integer(names(u.nei))
  coord = cbind(x,y)
  u.dist = sqrt(rowSums((coord[u.node,] - coord[u.nei,])^2))
  spDist = S4Vectors::splitAsList(u.dist, u.node)
  
  NG$nn = spNDT
  NG$dist = spDist
  NG$deg = IRanges::elementNROWS(NG$dist)
  NG$obj.id = obj.id
  NG$coord.x = x
  NG$coord.y = y
  NG$sample.id = sample.id
  
  if(!is.null(max.dist)){
    NG$nn = NG$nn[NG$dist < max.dist]
    NG$dist = NG$dist[NG$dist < max.dist]
    NG$deg = IRanges::elementNROWS(NG$dist)
  }
  
  NG
}

# Partial shuffle (with values in keep fixed) used in permutation analysis
partial_shuffle = function(values, keep){
  rand.values = values
  rand.values[!rand.values %in% keep] <- sample(rand.values[!rand.values %in% keep])
  rand.values
}

# Check if the spatial graph object is valid
NG_check = function(NG){
  
  # object data
  checkmate::assertTRUE(exists("NG"))
  checkmate::assertClass(NG,"DataFrame")
  checkmate::assertInteger(NG$index, unique = TRUE)
  
  # index must correspond to row number
  checkmate::assertTRUE(all(NG$index == seq_len(nrow(NG))))
  
  checkmate::assertNumeric(NG$coord.x)
  checkmate::assertNumeric(NG$coord.y)
  checkmate::assertCharacter(NG$obj.id, unique = TRUE)
  checkmate::assertCharacter(NG$sample.id)
  
  checkmate::assertClass(NG$nn,"IntegerList")
  checkmate::assertClass(NG$dist,"NumericList")
  
  # nn and dist lists are named
  checkmate::assertCharacter(names(NG$nn))
  checkmate::assertCharacter(names(NG$dist))
  
  # nn and dist lists have identical names
  checkmate::assertTRUE(all(names(NG$nn) == NG$index))
  checkmate::assertTRUE(all(names(NG$dist)== NG$index))
  checkmate::assertTRUE(all(IRanges::elementNROWS(NG$dist) == IRanges::elementNROWS(NG$nn)))
  
  # all nn indexes are valid
  checkmate::assertTRUE(all(sort(unique(unlist(NG$nn))) %in% NG$index))
  
  # check degree validity
  checkmate::assertTRUE(all(NG$deg == IRanges::elementNROWS(NG$nn)))
}

# Relabel node indexes after NG is changed (row sub-setting or reordering) to make it valid
NG_update_index = function(NG){
  checkmate::assertInteger(NG$index, unique = TRUE, any.missing = FALSE, len = nrow(NG))
  checkmate::assertClass(NG$nn, classes = "CompressedIntegerList")
  checkmate::assert(length(unique(NG$index))==nrow(NG))

  valid_edges =  S4Vectors::`%in%`(NG$nn,NG$index) 
  NG$nn = NG$nn[valid_edges]
  NG$dist = NG$dist[valid_edges]
  NG$deg = S4Vectors::elementNROWS(NG$nn)
  NG$nn = IRanges::match(NG$nn, NG$index)
  NG$index = 1:nrow(NG)
  names(NG$nn)  <- NG$index
  names(NG$dist) <- NG$index
  NG
}

# Return new spatial graph with subsetted nodes
NG_subset_nodes = function(NG, keep.node){
  
  NG_check(NG)
  checkmate::assertLogical(keep.node, len = length(NG$index))
  
  NG <- NG[keep.node,]
  
  # update nei index
  NG$nn <- IRanges::match(NG$nn, NG$index)
  NG$dist <- NG$dist[!is.na(NG$nn)]
  NG$nn <- NG$nn[!is.na(NG$nn)]
  
  # update node index (1:N)
  NG$index <- seq_len(nrow(NG))
  names(NG$nn) <- NG$index
  names(NG$dist) <- NG$index
  NG$deg <- IRanges::elementNROWS(NG$nn)
  NG_check(NG)
  NG
}

# Returns all NG edges in a DataFrame
# If provided, values in the info param will be returned with each edge.
NG_edges = function(NG, info = NULL){
  NG_check(NG)
  id1 = rep(NG$index, times = IRanges::elementNROWS(NG$nn))
  id2 = as.integer(unlist(NG$nn, use.names = FALSE))
  
  flevels = NG$obj.id
  
  df = S4Vectors::DataFrame(sample.id = NG$sample.id[id1],
                            obj.id.1 = factor(NG$obj.id,levels = flevels)[id1],
                            obj.id.2 = factor(NG$obj.id,levels = flevels)[id2],
                            id1,
                            id2,
                            dist = unlist(NG$dist, use.names = FALSE))
  
  if(!missing(info)){
    checkmate::assertVector(info,len = nrow(NG))
    df$info.1 = info[id1]
    df$info.2 = info[id2]
  }
  
  df
}

# Decompose each group into connected components (cc)
# and label each node with its cc id.
NG_add_cc = function(NG, vgroup = NG$group){
  
  NG_check(NG)
  checkmate::assertVector(vgroup, len = nrow(NG),any.missing = FALSE)
  NG$group = vgroup  
  
  # Create adjacency list for nodes of the same groups (EACH EDGE DUPLICATED)
  ADJ = NG$nn[IRanges::relist(NG$group[unlist(NG$nn)], skeleton = NG$nn) == NG$group]
  
  # Create igraph object from ADJ 
  g = igraph::graph_from_adj_list(as.list(ADJ), mode = "all", duplicate = TRUE) # new better
  cc = igraph::components(g, mode = "weak")
  NG$cc.id = cc$membership
  NG$cc.size = cc$csize[cc$membership]
  NG
}

# Returns node.attribute for neighbors as a list
NG_nei_attribute = function(NG, node.attribute = NG$group){
  NG_check(NG)
  if(!checkmate::checkVector(node.attribute, len = nrow(NG))){
    stop("param \"node.attribute\" must be a vector of length nrow(NG)")
  }
  checkmate::assertVector(node.attribute, len = length(NG$index))
  IRanges::relist(node.attribute[unlist(NG$nn,use.names = FALSE)],NG$nn)
}

# Returns nominal assortativity coefficient (r),
# given edge counts for two types of nodes (X and Y)
assortativity_binary <- function(n.xx, n.yy, n.xy){
  # Inputs are undirected edge counts
  n = n.xx + n.yy + n.xy
  p.xx = n.xx/n
  p.yy = n.yy/n
  p.xy = n.xy/n
  sum_a2 = sum((c(p.xx , p.yy) + p.xy/2)^2)
  sum_d = p.xx + p.yy
  r = (sum_d - sum_a2)/(1 - sum_a2)
  r
}

# Returns standard error (SE) for the nominal assortativity coefficient,
# given edge counts for two types of nodes (X and Y)
assortativity_binary_se_jacknife <- function(n.xx, n.yy, n.xy){
  
  r = assortativity_binary(n.xx, n.yy, n.xy)
  
  if(!is.finite(r)){
    return(NA)
  }
  
  delta_r2_sum = 0;
  
  if(n.xx > 0){
    delta_r2_sum = n.xx * (r - assortativity_binary(n.xx - 1, n.yy, n.xy))^2
  }
  if(n.yy > 0){
    delta_r2_sum = delta_r2_sum + n.yy * (r - assortativity_binary(n.xx, n.yy - 1, n.xy))^2
  }
  if(n.xy > 0){
    delta_r2_sum = delta_r2_sum + n.xy * (r - assortativity_binary(n.xx, n.yy, n.xy - 1))^2
  }
  
  sqrt(delta_r2_sum)
} 

# Helper function for NG_plot().
# Assumes NG is a single sample graph object.
# Duplicate edges are not returned.
get_edges_list = function(NG){
  
  if(length(unique(NG$sample.id))>1){
    stop("NG object contains multiple samples")
  }
  
  NG_check(NG)
  nei = NG$nn
  u.nei = unlist(nei)
  u.dist = unlist(NG$dist)
  de = data.frame(id1 = as.integer(names(u.nei)), id2 = u.nei, len = u.dist, stringsAsFactors = FALSE)
  
  demap1 = match(de$id1,NG$index)
  demap2 = match(de$id2,NG$index)
  
  de$x1 = NG$coord.x[demap1]
  de$y1 = NG$coord.y[demap1]
  de$x2 = NG$coord.x[demap2]
  de$y2 = NG$coord.y[demap2]
  
  de$group1 = NG$group[demap1]
  de$group2 = NG$group[demap2]
  
  de$obj.id1 = NG$obj.id[demap1]
  de$obj.id2 = NG$obj.id[demap2]
  
  de = de[de$id1 <= de$id2,]
  de
  
}
