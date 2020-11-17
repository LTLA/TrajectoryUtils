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
#' @param columns A character, logical or integer vector specifying the columns of \code{x} to use.
#' If \code{NULL}, all provided columns are used by default. 
#' @param outgroup A logical scalar indicating whether an outgroup should be inserted to split unrelated trajectories.
#' Alternatively, a numeric scalar specifying the distance threshold to use for this splitting.
#' @param outscale A numeric scalar specifying the scaling to apply to the median distance between centroids
#' to define the threshold for outgroup splitting.
#' Only used if \code{outgroup=TRUE}.
#' @param assay.type An integer or string specifying the assay to use from a SummarizedExperiment \code{x}.
#' @param use.dimred An integer or string specifying the reduced dimensions to use from a SingleCellExperiment \code{x}.
#' @param use.median A logical scalar indicating whether cluster centroid coordinates should be computed using the median rather than mean.
#' @param dist.method A string specifying the distance measure to be used, see Details.
#' @param with.mnn Logical scalar, deprecated; use \code{dist.method="mnn"} instead.
#' @param mnn.k An integer scalar specifying the number of nearest neighbors to consider for the MNN-based distance calculation when \code{dist.method="mnn"}.
#' See \code{\link[batchelor]{findMutualNN}} for more details.
#' @param BNPARAM A BiocNeighborParam object specifying how the nearest-neighbor search should be performed when \code{dist.method="mnn"},
#' see the \pkg{BiocNeighbors} package for more details.
#' @param BPPARAM A BiocParallelParam object specifying whether the nearest neighbor search should be parallelized when \code{dist.method="mnn"},
#' see the \pkg{BiocNeighbors} package for more details.
#'
#' @section Introducing an outgroup:
#' If \code{outgroup=TRUE}, we add an outgroup to avoid constructing a trajectory between \dQuote{unrelated} clusters (Street et al., 2018).
#' This is done by adding an extra row/column to the distance matrix corresponding to an artificial outgroup cluster,
#' where the distance to all of the other real clusters is set to \eqn{\omega/2}.
#' Large jumps in the MST between real clusters that are more distant than \eqn{\omega} will then be rerouted through the outgroup,
#' allowing us to break up the MST into multiple subcomponents by removing the outgroup.
#'
#' The default \eqn{\omega} value is computed by constructing the MST from the original distance matrix,
#' computing the median edge length in that MST, and then scaling it by \code{outscale}.
#' This adapts to the magnitude of the distances and the internal structure of the dataset
#' while also providing some margin for variation across cluster pairs.
#' Alternatively, \code{outgroup} can be set to a numeric scalar in which case it is used directly as \eqn{\omega}.
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
#' This mode can be enabled by setting \code{method="mnn"}, while the stringency of the MNN definition can be set with \code{mnn.k}.
#' Similarly, the performance of the nearest neighbor search can be controlled with \code{BPPARAM} and \code{BSPARAM}.
#' Note that this mode performs a cell-based search and so cannot be used when \code{x} already contains aggregated profiles.
#'
#' @section Using medians:
#' If \code{use.median=TRUE}, the median across all cells in each cluster is used to compute the centroid coordinate for each dimension.
#' This protects against outliers but is less stable than the mean.
#' Enabling this option is advisable if one observes that the default centroid is not located near any of its points due to outliers.
#' Note that the centroids computed in this manner is not a true medoid, which was too much of a pain to compute.
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
.create_cluster_mst <- function(x, clusters, use.median=FALSE, outgroup=FALSE, outscale=3, columns=NULL, dist.method = c("simple", "scaled.full", "scaled.diag", "slingshot", "mnn"), with.mnn=FALSE, mnn.k=50, BNPARAM=NULL, BPPARAM=NULL) {
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
            dmat <- .create_mnn_distance_matrix(x, clusters, levels=rownames(centers), 
                mnn.k=mnn.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
        } else {
            use.full <- (dist.method == "scaled.full" || (dist.method == "slingshot" && min(table(clusters)) <= ncol(x)))
            dmat <- .dist_clusters_scaled(x, clusters, centers, full=use.full)
        }
    }

    if (!isFALSE(outgroup)) {
        if (!is.numeric(outgroup)) {
            g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
            mst <- minimum.spanning.tree(g)
            med <- median(E(mst)$weight)
            outgroup <- med * outscale
        }

        old.d <- rownames(dmat)

        # Divide by 2 so rerouted distance between cluster pairs is 'outgroup'.
        dmat <- rbind(cbind(dmat, outgroup/2), outgroup/2) 
        diag(dmat) <- 0

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
.create_mnn_distance_matrix <- function(x, clusters, levels, mnn.k, BNPARAM=NULL, BPPARAM=NULL) {
    distances <- matrix(0, length(levels), length(levels), dimnames=list(levels, levels))

    if (is.null(BNPARAM)) {
        BNPARAM <- BiocNeighbors::KmknnParam()
    }
    if (is.null(BPPARAM)) {
        BPPARAM <- BiocParallel::SerialParam()
    }

    # TODO: modify batchelor so that findMutualNN can accept indices.
    for (first in levels) {
        left <- x[clusters==first,,drop=FALSE]
        for (second in levels) {
            if (first==second) break

            right <- x[clusters==second,,drop=FALSE]
            stuff <- batchelor::findMutualNN(left, right, k1=mnn.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
            dist2 <- rowSums((left[stuff$first,,drop=FALSE] - right[stuff$second,,drop=FALSE])^2)
            distances[first,second] <- sqrt(median(dist2))
        }
    }

    # Just making it symmetric.
    (distances + t(distances))
}

#' @importFrom igraph minimum.spanning.tree E E<- ends get.edge.ids delete.edges
.estimate_edge_confidence <- function(mst, g) {
    edges <- E(mst)
    ends <- ends(mst, edges)
    reweight <- numeric(length(edges))

    for (i in seq_along(edges)) {
        id <- get.edge.ids(g, ends[i,])        
        g.copy <- delete.edges(g, id)
        mst.copy <- minimum.spanning.tree(g.copy)
        reweight[i] <- sum(E(mst.copy)$weight)
    }

    W <- edges$weight
    total <- sum(W)
    offset <- min(W)
    E(mst)$gain <- (reweight - total)/(W + offset/1e8)
    mst
}

.dist_clusters_scaled <- function(x, clusters, centers, full) {
    nclust <- nrow(centers)
    output <- matrix(0, nclust, nclust, dimnames=list(rownames(centers), rownames(centers)))

    for (i in seq_len(nclust)) {
        mu1 <- centers[i,]
        clus1 <- rownames(centers)[i]
        s1 <- cov(x[which(clusters==clus1),, drop = FALSE])
        if (!full) {
            s1 <- diag(diag(s1))
        }

        for (j in seq_len(i - 1L)) {
            mu2 <- centers[j,]
            clus2 <- rownames(centers)[j]
            s2 <- cov(x[which(clusters==clus2),, drop = FALSE])
            if (!full) {
                s2 <- diag(diag(s2))
            }

            diff <- mu1 - mu2
            d <- as.numeric(t(diff) %*% solve(s1 + s2) %*% diff)
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
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("createClusterMST", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"), ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .create_cluster_mst(reducedDim(x, use.dimred), clusters=clusters, ...)
    } else {
        callNextMethod(x, clusters=clusters, ...)
    }
})
