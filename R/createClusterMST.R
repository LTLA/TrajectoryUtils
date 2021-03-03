#' Minimum spanning trees on cluster centroids
#'
#' Build a MST where each node is a cluster centroid and 
#' each edge is weighted by the Euclidean distance between centroids.
#' This represents the most parsimonious explanation for a particular trajectory
#' and has the advantage of being directly intepretable with respect to any pre-existing clusters.
#'
#' @param x A numeric matrix of coordinates where each row represents a cell/sample and each column represents a dimension
#' (usually a PC or another low-dimensional embedding, but features or genes can also be used).
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object
#' containing such a matrix in its \code{\link{assays}}, as specified by \code{assay.type}.
#' This will be transposed prior to use.
#'
#' Alternatively, for \linkS4class{SingleCellExperiment}s, this matrix may be extracted from its \code{\link{reducedDims}},
#' based on the \code{use.dimred} specification.
#' In this case, no transposition is performed.
#'
#' Alternatively, if \code{clusters=NULL}, a numeric matrix of coordinates for cluster centroids,
#' where each row represents a cluster and each column represents a dimension 
#' Each row should be named with the cluster name.
#' This mode can also be used with assays/matrices extracted from SummarizedExperiments and SingleCellExperiments. 
#' @param ... For the generic, further arguments to pass to the specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method
#' (if \code{use.dimred} is specified) or the ANY method (otherwise).
#' @param clusters A factor-like object of the same length as \code{nrow(x)},
#' specifying the cluster identity for each cell in \code{x}.
#' If \code{NULL}, \code{x} is assumed to already contain coordinates for the cluster centroids.
#'
#' Alternatively, a matrix with number of rows equal to \code{nrow(x)}, 
#' containing soft assignment weights for each cluster (column).
#' All weights should be positive and sum to 1 for each row.
#' @param columns A character, logical or integer vector specifying the columns of \code{x} to use.
#' If \code{NULL}, all provided columns are used by default. 
#' @param outgroup A logical scalar indicating whether an outgroup should be inserted to split unrelated trajectories.
#' Alternatively, a numeric scalar specifying the distance threshold to use for this splitting.
#' @param outscale A numeric scalar specifying the number of MADs above the median, used to define the threshold for outgroup splitting - see Details.
#' Only used if \code{outgroup=TRUE}.
#' @param endpoints A character vector of clusters that must be endpoints, i.e., nodes of degree 1 or lower in the MST.
#' @param assay.type An integer or string specifying the assay to use from a SummarizedExperiment \code{x}.
#' @param use.dimred An integer or string specifying the reduced dimensions to use from a SingleCellExperiment \code{x}.
#' @param use.median A logical scalar indicating whether cluster centroid coordinates should be computed using the median rather than mean.
#' @param dist.method A string specifying the distance measure to be used, see Details.
#' @param with.mnn Logical scalar, deprecated; use \code{dist.method="mnn"} instead.
#' @param mnn.k An integer scalar specifying the number of nearest neighbors to consider for the MNN-based distance calculation when \code{dist.method="mnn"}.
#' See \code{\link[BiocNeighbors]{findMutualNN}} for more details.
#' @param BNPARAM A BiocNeighborParam object specifying how the nearest-neighbor search should be performed when \code{dist.method="mnn"},
#' see the \pkg{BiocNeighbors} package for more details.
#' @param BPPARAM A BiocParallelParam object specifying whether the nearest neighbor search should be parallelized when \code{dist.method="mnn"},
#' see the \pkg{BiocNeighbors} package for more details.
#'
#' @section Computing the centroids:
#' By default, the cluster centroid is defined by taking the mean value across all of its cells for each dimension.
#' If \code{clusters} is a matrix, a weighted mean is used instead.
#' This treats the column of weights as fractional identities of each cell to the corresponding cluster.
#' 
#' If \code{use.median=TRUE}, the median across all cells in each cluster is used to compute the centroid coordinate for each dimension.
#' (With a matrix-like \code{clusters}, a weighted median is calculated.)
#' This protects against outliers but is less stable than the mean.
#' Enabling this option is advisable if one observes that the default centroid is not located near any of its points due to outliers.
#' Note that the centroids computed in this manner is not a true medoid, which was too much of a pain to compute.
#'
#' @section Introducing an outgroup:
#' If \code{outgroup=TRUE}, we add an outgroup to avoid constructing a trajectory between \dQuote{unrelated} clusters (Street et al., 2018).
#' This is done by adding an extra row/column to the distance matrix corresponding to an artificial outgroup cluster,
#' where the distance to all of the other real clusters is set to \eqn{\omega/2}.
#' Large jumps in the MST between real clusters that are more distant than \eqn{\omega} will then be rerouted through the outgroup,
#' allowing us to break up the MST into multiple subcomponents (i.e., a minimum spanning forest) by removing the outgroup.
#'
#' The default \eqn{\omega} value is computed by constructing the MST from the original distance matrix,
#' computing the median edge length in that MST, and then adding it to the MAD of the edge lengths multiplied by \code{outscale}.
#' This effectively identifies branches with large \dQuote{outlier} distances and breaks them to form separate graph components.
#' In this manner, we adapt to the magnitude of the distances and the internal structure of the dataset while still allowing for some variation in the distances. 
#' Alternatively, \code{outgroup} can be set to a numeric scalar in which case it is used directly as \eqn{\omega}.
#'
#' @section Forcing endpoints:
#' If certain clusters are known to be endpoints (e.g., because they represent terminal states), we can specify them in \code{endpoints}.
#' This ensures that the returned graph will have such clusters as nodes of degree 1, i.e., they terminate the path.
#' The function uses an exhaustive search to identify the MST with these constraints.
#' If no configuration can be found, an error is raised - this will occur if all nodes are specified as endpoints, for example.
#'
#' If \code{outgroup=TRUE}, the function is allowed to connect two endpoints together to create a two-node subcomponent.
#' This will result in the formation of a minimum spanning forest if there are more than two clusters in \code{x}.
#' Of course, if there are only two nodes and both are specified as endpoints, a two-node subcomponent will be formed regardless of \code{outgroup}.
#'
#' Note that edges involving endpoint nodes will have infinite confidence values (see below).
#' This reflects the fact that they are forced to exist during graph construction.
#'
#' @section Confidence on the edges:
#' For the MST, we obtain a measure of the confidence in each edge by computing the distance gained if that edge were not present.
#' Ambiguous parts of the tree will be less penalized from deletion of an edge, manifesting as a small distance gain.
#' In contrast, parts of the tree with clear structure will receive a large distance gain upon deletion of an obvious edge.
#'
#' For each edge, we divide the distance gain by the length of the edge to normalize for cluster resolution.
#' This avoids overly penalizing edges in parts of the tree involving broad clusters
#' while still retaining sensitivity to detect distance gain in overclustered regions.
#' As an example, a normalized gain of unity for a particular edge means that its removal
#' requires an alternative path that increases the distance travelled by that edge's length.
#'
#' The normalized gain is reported as the \code{"gain"} attribute in the edges of the MST from \code{\link{createClusterMST}}.
#' Note that the \code{"weight"} attribute represents the edge length.
#'
#' @section Distance measures:
#' Distances between cluster centroids may be calculated in multiple ways:
#' \itemize{
#' \item The default is \code{"simple"}, which computes the Euclidean distance between cluster centroids.
#' \item With \code{"scaled.diag"}, we downscale the distance between the centroids by the sum of the variances of the two corresponding clusters (i.e., the diagonal of the covariance matrix).
#' This accounts for the cluster \dQuote{width} by reducing the effective distances between broad clusters.
#' \item With \code{"scaled.full"}, we repeat this scaling with the full covariance matrix.
#' This accounts for the cluster shape by considering correlations between dimensions, but cannot be computed when there are more cells than dimensions.
#' \item The \code{"slingshot"} option will typically be equivalent to the \code{"scaled.full"} option, 
#' but switches to \code{"scaled.diag"} in the presence of small clusters (fewer cells than dimensions in the reduced dimensional space). 
#' \item For \code{"mnn"}, see the more detailed explanation below.
#' }
#' 
#' If \code{clusters} is a matrix with \code{"scaled.diag"}, \code{"scaled.full"} and \code{"slingshot"}, 
#' a weighted covariance is computed to account for the assignment ambiguity.
#' In addition, a warning will be raised if \code{use.median=TRUE} for these choices of \code{dist.method};
#' the Mahalanobis distances will not be correctly computed when the centers are medians instead of means.
#'
#' @section Alternative distances with MNN pairs:
#' While distances between centroids are usually satisfactory for gauging cluster \dQuote{closeness}, 
#' they do not consider the behavior at the boundaries of the clusters.
#' Two clusters that are immediately adjacent (i.e., intermingling at the boundaries) may have a large distance between their centroids
#' if the clusters themselves span a large region of the coordinate space.
#' This may preclude the obvious edge from forming in the MST.
#'
#' In such cases, we can use an alternative distance calculation based on the distance between mutual nearest neighbors (MNNs).
#' An MNN pair is defined as two cells in separate clusters that are each other's nearest neighbors in the other cluster.
#' For each pair of clusters, we identify all MNN pairs and compute the median distance between them.
#' This distance is then used in place of the distance between centroids to construct the MST.
#' In this manner, we focus on cluster pairs that are close at their boundaries rather than at their centers.
#'
#' This mode can be enabled by setting \code{dist.method="mnn"}, while the stringency of the MNN definition can be set with \code{mnn.k}.
#' Similarly, the performance of the nearest neighbor search can be controlled with \code{BPPARAM} and \code{BSPARAM}.
#' Note that this mode performs a cell-based search and so cannot be used when \code{x} already contains aggregated profiles.
#'
#' @return A \link{graph} object containing an MST computed on \code{centers}.
#' Each node corresponds to a cluster centroid and has a numeric vector of coordinates in the \code{coordinates} attribute.
#' The edge weight is set to the Euclidean distance and the confidence is stored as the \code{gain} attribute.
#'
#' @author Aaron Lun
#'
#' @references
#' Ji Z and Ji H (2016).
#' TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.
#' \emph{Nucleic Acids Res.} 44, e117
#'
#' Street K et al. (2018).
#' Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. 
#' \emph{BMC Genomics}, 477.
#'
#' @examples
#' # Mocking up a Y-shaped trajectory.
#' centers <- rbind(c(0,0), c(0, -1), c(1, 1), c(-1, 1))
#' rownames(centers) <- seq_len(nrow(centers))
#' clusters <- sample(nrow(centers), 1000, replace=TRUE)
#' cells <- centers[clusters,]
#' cells <- cells + rnorm(length(cells), sd=0.5)
#'
#' # Creating the MST:
#' mst <- createClusterMST(cells, clusters)
#' plot(mst)
#'
#' # We could also do it on the centers:
#' mst2 <- createClusterMST(centers, clusters=NULL)
#' plot(mst2)
#'
#' # Works if the expression matrix is in a SE:
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(t(cells), colData=DataFrame(group=clusters))
#' mst3 <- createClusterMST(se, se$group, assay.type=1)
#' plot(mst3)
#' 
#' @name createClusterMST
NULL

#################################################

#' @importFrom igraph graph.adjacency minimum.spanning.tree delete_vertices E V V<-
#' @importFrom stats median dist
.create_cluster_mst <- function(x, clusters, use.median=FALSE, outgroup=FALSE, outscale=3, endpoints=NULL, columns=NULL,
    dist.method = c("simple", "scaled.full", "scaled.diag", "slingshot", "mnn"), 
    with.mnn=FALSE, mnn.k=50, BNPARAM=NULL, BPPARAM=NULL) 
{
    if (!is.null(columns)) {
        x <- x[,columns,drop=FALSE]                
    }

    if (!is.null(clusters)) {
        FUN <- if (use.median) rowmedian else rowmean
        centers <- FUN(x, clusters)
    } else if (is.null(rownames(x))) {
        stop("'x' must have row names corresponding to cluster names")
    } else {
        centers <- as.matrix(x)
    }

    dist.method <- match.arg(dist.method)
    if (with.mnn) {
        .Deprecated(old="with.mnn=TRUE", new="dist.method=\"mnn\"")
        dist.method <- "mnn"
    }

    if (dist.method == "simple") {
        dmat <- dist(centers)
        dmat <- as.matrix(dmat)
    } else {
        if (is.null(clusters)) {
            stop("'clusters' must be specified when 'dist.method!=\"simple\"'")
        }

        if (dist.method == "mnn") {
            dmat <- .create_mnn_distance_matrix(x, clusters, mnn.k=mnn.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
        } else {
            if (use.median) {
                # Distances not really intepretable as Mahalanobis distances anymore.
                warning("'use.median=TRUE' with 'dist.method=\"", dist.method, "\"' may yield unpredictable results")
            }
            use.full <- (dist.method == "scaled.full" || (dist.method == "slingshot" && min(table(clusters)) > ncol(x)))
            dmat <- .dist_clusters_scaled(x, clusters, centers=centers, full=use.full)
        }
    }

    # Ensure all off-diagonal distances are positive, as zero weights = no edge. 
    lower.limit <- min(dmat[dmat > 0])
    dmat[] <- pmax(dmat, lower.limit[1] / 1e6)
    diag(dmat) <- 0

    if (!is.null(endpoints)) {
        # If outgroup=TRUE, then we can have multi-component graphs.
        # If there are only two nodes, we don't really have much choice.
        allow.dyads <- !isFALSE(outgroup) || nrow(dmat) == 2
        dmat <- .enforce_endpoints(dmat, endpoints, allow.dyads=allow.dyads)
    }

    if (!isFALSE(outgroup)) {
        if (!is.numeric(outgroup)) {
            g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
            mst <- minimum.spanning.tree(g)
            mst.w <- E(mst)$weight
            outgroup <- median(mst.w) + mad(mst.w) * outscale
        }

        old.d <- rownames(dmat)
        dmat <- rbind(cbind(dmat, outgroup), outgroup) 
        dmat[length(dmat)] <- 0
        special.name <- strrep("x", max(nchar(old.d))+1L)
        rownames(dmat) <- colnames(dmat) <- c(old.d, special.name)
    }

    g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
    mst <- minimum.spanning.tree(g)
    mst <- .estimate_edge_confidence(mst, g)

    if (!isFALSE(outgroup)) {
        mst <- delete_vertices(mst, special.name)
    }

    # Embed vertex coordinates for downstream use.
    coord.list <- vector("list", nrow(centers))
    names(coord.list) <- rownames(centers)
    for (r in rownames(centers)) {
        coord.list[[r]] <- centers[r,]
    }
    V(mst)$coordinates <- coord.list[names(V(mst))]

    mst 
}

#' @importFrom stats median
#' @importFrom Matrix rowSums
.create_mnn_distance_matrix <- function(x, clusters, mnn.k, BNPARAM=NULL, BPPARAM=NULL) {
    if (is.null(BNPARAM)) {
        BNPARAM <- BiocNeighbors::KmknnParam()
    }
    if (is.null(BPPARAM)) {
        BPPARAM <- BiocParallel::SerialParam()
    }

    if (is.matrix(clusters)) {
        cluster.ids <- .choose_colnames(clusters)
        clusters <- cluster.ids[max.col(clusters, ties.method="first")]
    }
    stopifnot(length(clusters)==nrow(x))
    by.cluster <- split(seq_along(clusters), clusters)
    levels <- names(by.cluster)

    # Looping through all of them.
    collated <- indices <- vector("list", length(levels))
    for (i in seq_along(levels)) {
        chosen <- by.cluster[[i]]
        collated[[i]] <- x[chosen,,drop=FALSE]
        indices[[i]] <- BiocNeighbors::buildIndex(collated[[i]], BNPARAM=BNPARAM)
    }

    distances <- matrix(0, length(levels), length(levels), dimnames=list(levels, levels))

    for (f in seq_along(levels)) {
        left <- collated[[f]]
        lefti <- indices[[f]]

        for (s in seq_along(levels)) {
            if (f==s) break
            right <- collated[[s]]
            righti <- indices[[s]]

            stuff <- BiocNeighbors::findMutualNN(left, right, k1=mnn.k, 
                BNINDEX1=lefti, BNINDEX2=righti, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
            dist2 <- rowSums((left[stuff$first,,drop=FALSE] - right[stuff$second,,drop=FALSE])^2)
            distances[f,s] <- sqrt(median(dist2))
        }
    }

    # Just making it symmetric.
    (distances + t(distances))
}

#' @importFrom S4Vectors head
.enforce_endpoints <- function(dmat, endpoints, allow.dyads=FALSE) {
    available <- dmat[as.character(unique(endpoints)),,drop=FALSE]
    best.stats <- new.env()
    best.stats$distance <- Inf

    SEARCH <- function(path=character(0), distance=0) {
        i <- length(path) + 1L
        if (i > nrow(available)) {
            if (distance < best.stats$distance) {
                best.stats$distance <- distance
                best.stats$path <- path
            }
            return(NULL)
        } else if (distance > best.stats$distance) {
            return(NULL)
        } 
        
        current <- rownames(available)[i]
        self.used <- which(path == current)

        if (length(self.used) == 1) { 
            # Endpoint-to-endpoint dyads should be reciprocated,
            # with no distance added (if they are allowed to exist).
            reciprocal <- rownames(available)[self.used]
            if (!reciprocal %in% path && allow.dyads) {
                SEARCH(c(path, reciprocal), distance)
            }
        } else {
            used.endpoints <- c(current, # currently in use.
                head(rownames(available), length(path)), # endpoints connected from in previous steps.
                intersect(path, rownames(available))) # endpoints connected to in previous steps.
            allowed <- setdiff(colnames(available), used.endpoints)
            for (j in allowed) {
                SEARCH(c(path, j), distance + available[i,j])
            }
        }
    }

    SEARCH()
    if (is.infinite(best.stats$distance)) {
        stop("no solvable tree for specified 'endpoints'")
    }

    for (a in seq_len(nrow(available))) {
        current <- rownames(available)[a]
        others <- setdiff(colnames(dmat), best.stats$path[a])
        dmat[current,others] <- 0
        dmat[others,current] <- 0
    }

    dmat
}

#' @importFrom igraph minimum.spanning.tree E E<- ends get.edge.ids delete.edges V
.estimate_edge_confidence <- function(mst, g) {
    edges <- E(mst)
    ends <- ends(mst, edges)
    reweight <- numeric(length(edges))
    to.skip <- names(V(g))[degree(g) <= 1]

    for (i in seq_along(edges)) {
        cur.ends <- ends[i,]
        if (any(cur.ends %in% to.skip)) {
            reweight[i] <- Inf
        } else {
            id <- get.edge.ids(g, cur.ends)
            g.copy <- delete.edges(g, id)
            mst.copy <- minimum.spanning.tree(g.copy)
            reweight[i] <- sum(E(mst.copy)$weight)
        }
    }

    W <- edges$weight
    total <- sum(W)
    offset <- min(W)
    E(mst)$gain <- (reweight - total)/(W + offset/1e8)
    mst
}

#' @importFrom Matrix t crossprod
.dist_clusters_scaled <- function(x, clusters, centers, full) {
    nclust <- nrow(centers)
    output <- matrix(0, nclust, nclust, dimnames=list(rownames(centers), rownames(centers)))

    # Computing the covariances (possibly with weights).
    all.cor <- vector("list", nclust)
    names(all.cor) <- rownames(centers)

    if (is.matrix(clusters)) {
        # Treating the weights as effective frequencies.
        clusters <- clusters/rowSums(clusters)
        for (i in seq_along(all.cor)) {
            curweight <- clusters[,i]
            out <- t(t(x) - centers[i,]) * sqrt(curweight)
            all.cor[[i]] <- crossprod(out)/(sum(curweight) - 1)
        }
    } else {
        for (i in seq_along(all.cor)) {
            all.cor[[i]] <- cov(x[which(clusters==names(all.cor)[i]),, drop = FALSE])
        }
    }

    for (i in seq_len(nclust)) {
        mu1 <- centers[i,]
        clus1 <- rownames(centers)[i]
        s1 <- all.cor[[i]]
        if (!full) {
            s1 <- diag(diag(s1))
        }

        for (j in seq_len(i - 1L)) {
            mu2 <- centers[j,]
            clus2 <- rownames(centers)[j]
            s2 <- all.cor[[j]]
            if (!full) {
                s2 <- diag(diag(s2))
            }

            diff <- mu1 - mu2
            d <- sqrt(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
            output[i,j] <- output[j,i] <- d
        }
    }

    output
}

#################################################

#' @export
#' @rdname createClusterMST
setGeneric("createClusterMST", function(x, ...) standardGeneric("createClusterMST"))

#' @export
#' @rdname createClusterMST
setMethod("createClusterMST", "ANY", .create_cluster_mst)

#' @export
#' @rdname createClusterMST
#' @importFrom Matrix t
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("createClusterMST", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .create_cluster_mst(t(assay(x, assay.type)), ...)
})

#' @export
#' @rdname createClusterMST
#' @importFrom SingleCellExperiment reducedDim colLabels
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("createClusterMST", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"), ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .create_cluster_mst(reducedDim(x, use.dimred), clusters=clusters, ...)
    } else {
        callNextMethod(x, clusters=clusters, ...)
    }
})
