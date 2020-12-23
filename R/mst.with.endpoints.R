#' Create a minimum spanning tree with known endpoints
#'
#' Create a minimum spanning tree where certain nodes are forced to be endpoints, i.e., of degree 1.
#' 
#' @param dmat A symmetric adjacency matrix where each row/column is a node.
#' Each entry represents the edge weight between the corresponding nodes;
#' however, no edge is formed when the weight is zero.
#' This should have non-\code{NULL} dimnames.
#' @param endpoints A character vector specifying the nodes to be set as endpoints.
#' @param outgroup.name String specifying the name of the node to use as the outgroup.
#' @param error Logical scalar indicating whether an error should be raised if no tree satisfies the constraints.
#'
#' @return
#' A \link{graph} object containing the minimum spanning tree (or possibly forest, if \code{dyad.penalty} is finite).
#' This may also be \code{NULL} if \code{error=FALSE} and no tree can be found that satisfies the constraints.
#'
#' @details
#' Pretty much as the name suggests, this function will search for the minimum spanning tree with endpoint constraints.
#' For most part, this involves removing the endpoints, identifying the MST from the remaining non-endpoint nodes,
#' and then connecting the endpoints to the closest non-endpoint node to create the full MST.
#'
#' If \code{outgroup.name} is supplied, each endpoint is allowed to connect to the outgroup node in addition to one other node.
#' This makes it possible to form an edge between two endpoints, named here as a \dQuote{dyad}.
#' The idea is that, once the outgroup is removed (see \code{?\link{createClusterMST}}), the connected endpoints form their own subcomponent.
#'
#' The function will perform an exhaustive search for the optimal configuration of edges from endpoints if dyads are allowed.
#' The edge width between the endpoint and the outgroup represents the effective penalty for formation of a dyad.
#' Otherwise, the dyad would eliminate an edge, causing a large drop in the total tree distance and encouraging formation of dyads everywhere.
#' 
#' In some situations, it is impossible to construct a tree, e.g., for an odd number of nodes that are all endpoints.
#' This will result in an error being raised.
#' Users can set \code{error=FALSE} to return a \code{NULL} instead to handle the error state in their own code.
#' 
#' @author Aaron Lun
#'
#' @examples
#' coords <- rbind(A=c(0,0), B=c(1,-1), C=c(1, 1))
#' dmat <- as.matrix(dist(coords))
#'
#' mst.with.endpoints(dmat, endpoints=NULL)
#' mst.with.endpoints(dmat, endpoints="A")
#' mst.with.endpoints(dmat, endpoints="B")
#'
#' mst.with.endpoints(dmat, endpoints=c("A", "B"))
#' try(mst.with.endpoints(dmat, endpoints=c("A", "B", "C")))
#'
#' # Sometimes MSTs are only possible when dyads are allowed to form.
#' # This requires us to pass in a dummy node as an outgroup.
#' dmat2 <- cbind(dummy=10, rbind(dummy=10, dmat))
#' dmat2[1,1] <- 0
#' mst.with.endpoints(dmat2, endpoints=c("A", "B", "C"), outgroup.name="dummy")
#' 
#' @seealso
#' \code{\link{minimum.spanning.tree}}, for the version of this function \emph{sans} the endpoint considerations.
#'
#' \code{\link{createClusterMST}}, where this function gets used.
#'
#' @export
#' @importFrom igraph minimum.spanning.tree add_edges graph.adjacency E
#' @importFrom S4Vectors head
mst.with.endpoints <- function(dmat, endpoints, outgroup.name=NULL, error=TRUE) {
    stopifnot(all(diag(dmat)==0))
    if (length(endpoints)==0) {
        g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
        return(minimum.spanning.tree(g))
    }

    endpoints <- as.character(unique(endpoints))
    if (nrow(dmat)==2L && all(rownames(dmat) %in% endpoints)) {
        return(graph.adjacency(dmat, mode = "undirected", weighted = TRUE))
    }

    # Removing the endpoints before searching for an MST.
    dmat0 <- dmat
    dmat0[endpoints,] <- 0
    dmat0[,endpoints] <- 0
    g0 <- graph.adjacency(dmat0, mode = "undirected", weighted = TRUE)
    mst0 <- minimum.spanning.tree(g0)

    # Connecting the endpoints to the closest non-endpoint. Placeholder dist
    # is 'Inf' so that we get init.weight=Inf for the dyad search.
    chosen.node <- rep(NA_character_, length(endpoints))
    chosen.dist <- rep(Inf, length(endpoints))
    allow.dyads <- !is.null(outgroup.name)

    for (i in seq_along(endpoints)) {
        endpt <- endpoints[i]
        choices <- dmat[endpt,]
        choices <- choices[choices > 0]
        choices <- choices[setdiff(names(choices), endpoints)]

        if (length(choices) > 0) {
            is.min <- which.min(choices)
            chosen.node[i] <- names(choices)[[is.min]]
            chosen.dist[i] <- choices[[is.min]]
        } else if (!allow.dyads) {
            if (error) {
                stop("failed to connect endpoint '", endpt, "' to the MST")
            } else {
                return(NULL)
            }
        }
    }

    if (allow.dyads) {
        # Hot-starting from the greedy solution, to increase the odds of early
        # termination in the recursive search.
        dyad.out <- .explore_dyad_solutions(
            init.path=chosen.node, 
            init.weight=sum(E(mst0)$weight) + sum(chosen.dist),
            dmat=dmat, 
            endpoints=endpoints,
            outgroup.name=outgroup.name
        )

        if (is.null(dyad.out)) {
            if (error) {
                stop("no solvable tree with dyads for specified 'endpoints'")
            } else {
                return(NULL)
            }
        }

        chosen.node <- dyad.out$nodes
        chosen.dist <- dyad.out$distances
    }

    add_edges(mst0, rbind(endpoints, chosen.node), attr=list(weight=chosen.dist))
}

.explore_dyad_solutions <- function(init.path, init.weight, dmat, endpoints, outgroup.name) {
    best.stats <- new.env()
    best.stats$path <- init.path
    best.stats$distance <- init.weight

    # Searching for the optimal configuration of endpoints.
    available <- dmat[endpoints,,drop=FALSE]

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
            # For dyads, we add the distance to the specified outgroup.
            reciprocal <- rownames(available)[self.used]
            if (!reciprocal %in% path) {
                SEARCH(c(path, outgroup.name), distance + available[i,outgroup.name])
            }
        } else {
            choices <- available[i,]

            # Removes self and prevents other disallowed edges (from .estimate_edge_confidence).
            choices <- choices[choices > 0] 

            # Do not connect to other endpoints that have already been used.
            used.endpoints <- c(
                head(rownames(available), length(path)), # endpoints connected from in previous steps.
                intersect(path, rownames(available)) # endpoints connected to in previous steps.
            )

            for (j in setdiff(names(choices), used.endpoints)) {
                SEARCH(c(path, j), distance + choices[j])
            }
        }
    }

    SEARCH()
    if (is.infinite(best.stats$distance)) {
        return(NULL)
    }

    weight <- rep(NA_real_, nrow(available))
    for (a in seq_along(weight)) {
        weight[a] <- available[a,best.stats$path[a]]
    }

    list(nodes=best.stats$path, distances=weight)
}
