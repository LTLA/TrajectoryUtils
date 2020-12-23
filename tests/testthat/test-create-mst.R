# This tests the createClusterMST function.
# library(testthat); library(TrajectoryUtils); source("setup.R"); source("test-create-mst.R")

set.seed(1000)

test_that("MST construction works as expected", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    mst <- createClusterMST(y, cluster=NULL)

    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("A", "D"))
    expect_identical(vertices[igraph::degree(mst)==2], c("B", "C"))

    # For a more complex example.
    y <- rbind(A=c(0, 1), B=c(0, 0), C=c(1, 1), D=c(-1, 1)) 
    mst <- createClusterMST(y, cluster=NULL)

    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C", "D"))
    expect_identical(vertices[igraph::degree(mst)==3], "A")
})

test_that("MST edge and vertex attributes make sense", {
    y <- matrix(rnorm(20), ncol=2)
    rownames(y) <- LETTERS[1:10]
    mst <- createClusterMST(y, cluster=NULL)

    # Check that the distances between edges are what is expected.
    stuff <- Matrix::which(mst[] > 0, arr.ind=TRUE)
    coords <- igraph::V(mst)$coordinates

    for (i in seq_len(nrow(stuff))) {
        M <- stuff[i,1]
        N <- stuff[i,2]
        expect_equal(mst[M, N], sqrt(sum((coords[[M]] - coords[[N]])^2)))
    }

    # Centers must be named!
    y <- matrix(rnorm(20), ncol=2)
    expect_error(createClusterMST(y, cluster=NULL), "must have row names")
})

test_that("MST construction works as expected with clusters", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4))
    mst <- createClusterMST(y, cluster=NULL)

    clusters <- sample(rownames(y), 1000, replace=TRUE)
    y0 <- y[clusters,,drop=FALSE] + runif(1000, -0.1, 0.1)
    mst0 <- createClusterMST(y0, clusters=clusters)
    expect_identical(mst[] > 0, mst0[] > 0)

    # For a more complex example.
    y <- rbind(A=c(0, 1), B=c(0, 0), C=c(1, 1), D=c(-1, 1))
    mst <- createClusterMST(y, cluster=NULL)

    clusters <- sample(rownames(y), 1000, replace=TRUE)
    y0 <- y[clusters,,drop=FALSE] + runif(1000, -0.1, 0.1)
    mst0 <- createClusterMST(y0, clusters=clusters)
    expect_identical(mst[] > 0, mst0[] > 0)
})

test_that("MST construction works as expected with columns", {
    y <- matrix(rnorm(100), ncol=10)
    rownames(y) <- 1:10
    mst <- createClusterMST(y[,1:5], cluster=NULL)
    mst0 <- createClusterMST(y, cluster=NULL, columns=1:5)
    expect_identical(mst[], mst0[])
})

test_that("MST construction works as expected with weights", {
    y <- matrix(rnorm(100), ncol=2)
    clusters <- sample(5, nrow(y), replace=TRUE)
    mat <- factor2matrix(clusters)

    ref <- createClusterMST(y, clusters=clusters)
    mst <- createClusterMST(y, clusters=mat)
    expect_equal(mst[], ref[])

    mst <- createClusterMST(y, clusters=abs(jitter(mat)))
    expect_false(isTRUE(all.equal(mst[], ref[])))

    # And again, with medians. 
    ref <- createClusterMST(y, clusters=clusters, use.median=TRUE)
    mst <- createClusterMST(y, clusters=mat, use.median=TRUE)
    expect_equal(mst[], ref[])

    mst <- createClusterMST(y, clusters=abs(jitter(mat)), use.median=TRUE)
    expect_false(isTRUE(all.equal(mst[], ref[])))
})

test_that("MST construction works as expected with outgroup specification", {
    # No effect with default outgroup settings.
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    ref <- createClusterMST(y, cluster=NULL)
    mst <- createClusterMST(y, outgroup=TRUE, cluster=NULL)
    expect_identical(ref[], mst[])

    # Maximal effect with crazy outgroup settings.
    mst <- createClusterMST(y, outgroup=1e-8, cluster=NULL)
    expect_identical(dim(mst[]), c(4L, 4L))
    expect_identical(sum(mst[]), 0)

    # Default effects.
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(10, 3), D=c(10, 4))
    mst <- createClusterMST(y, outgroup=TRUE, cluster=NULL)
    expect_identical(igraph::components(mst)$no, 2L)
    expect_false(igraph::are_adjacent(mst, "B", "C"))

    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(10, 4), E=c(0,4)) 
    mst <- createClusterMST(y, outgroup=TRUE, cluster=NULL)
    expect_identical(sum(mst[]["D",]), 0)
    expect_identical(sum(mst[][,"D"]), 0)

    # Automated outgroup calculation is done correctly.
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 4), D=c(0, 10.001))
    out <- createClusterMST(y, clusters=NULL, outgroup=TRUE, outscale=3)
    expect_identical(igraph::components(out)$no, 2L)

    y[length(y)] <- 9.999
    out <- createClusterMST(y, clusters=NULL, outgroup=TRUE, outscale=3)
    expect_identical(igraph::components(out)$no, 1L)

    # Gain calculations are *mostly* unaffected. Kind of depends
    # on whether the rerouted path passes through the outgroup.
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    ref <- createClusterMST(y, cluster=NULL)
    mst <- createClusterMST(y, outgroup=TRUE, outscale=10, cluster=NULL)
    expect_identical(igraph::E(ref)$gain, igraph::E(mst)$gain)
})

test_that("MST construction works with endpoint specification", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    out <- createClusterMST(y, endpoint="B", clusters=NULL)
    expect_true(igraph::degree(out, "B")==1L)
    edge.id <- igraph::get.edge.ids(out, c("A", "B"))
    expect_identical(igraph::E(out)$gain[edge.id], 0) # confidence is computed correctly.
    expect_identical(igraph::components(out)$no, 1L)

    out <- createClusterMST(y, endpoint=c("C", "D"), clusters=NULL)
    expect_true(igraph::are_adjacent(out, "B", "C")) 
    expect_true(igraph::are_adjacent(out, "B", "D"))
    expect_identical(igraph::components(out)$no, 1L)

    # Does sensible things when endpoints= doesn't have an effect.
    ref <- createClusterMST(y, clusters=NULL)
    out <- createClusterMST(y, endpoint=c("A", "D"), clusters=NULL)
    expect_identical(ref[], out[])

    # Handles extreme cases correctly.
    out <- createClusterMST(y[1:2,], endpoint=rownames(y)[1:2], clusters=NULL)
    expect_true(igraph::are_adjacent(out, "A", "B")) # forms dyads.
    expect_identical(igraph::E(out)$weight, 1)

    expect_error(createClusterMST(y, endpoint=rownames(y), clusters=NULL), "failed to connect")
})

test_that("MST construction doesn't fail with identical points", {
    y <- rbind(A=c(0, 1), A2=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    out <- createClusterMST(y, clusters=NULL)
    expect_true(igraph::are_adjacent(out, "A", "A2"))
})

library(SingleCellExperiment)
test_that("MST construction works with SE and SCE inputs", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    ref <- createClusterMST(y, cluster=NULL)

    se <- SummarizedExperiment(t(y))
    mst <- createClusterMST(se, cluster=NULL, assay.type=1)
    expect_identical(ref[], mst[])

    sce <- SingleCellExperiment(t(y))
    mst <- createClusterMST(sce, cluster=NULL, assay.type=1)
    expect_identical(ref[], mst[])

    y2 <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(10, 4))
    reducedDim(sce, "PCA") <- y2
    mst <- createClusterMST(sce, cluster=NULL, use.dimred="PCA")
    expect_false(identical(ref[], mst[]))
})

#############################################3

y.simple <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4))
clusters.simple <- sample(rownames(y.simple), 1000, replace=TRUE)
y0.simple <- y.simple[clusters.simple,,drop=FALSE] + runif(2000, -0.1, 0.1)

# Works for a less obvious example. Batch 3 is situated just at the x-midpoint
# of batches 1 and 2, but offset on the y-axis and separated with some space.
# By contrast, batches 1 and 2 are touching each other and should be connected.
y1 <- matrix(runif(500), ncol=2)
y2 <- matrix(runif(500), ncol=2)
y2[,1] <- y2[,1] + 1
y3 <- matrix(runif(500, 0, 0.1), ncol=2)
y3[,2] <- y3[,2] + 1.1
y3[,1] <- y3[,1] + 0.95

y.complex <- rbind(y1, y2, y3)
clusters.complex <- gl(3, 250)

set.seed(100101)
test_that("MST construction works with MNN-based distances", {
    ref <- createClusterMST(y0.simple, clusters=clusters.simple)
    mst <- createClusterMST(y0.simple, clusters=clusters.simple, dist.method="mnn")
    expect_identical(ref[] > 0, mst[] > 0)
    expect_identical(igraph::V(ref)$coordinates, igraph::V(mst)$coordinates)

    ref <- createClusterMST(y.complex, clusters=clusters.complex)
    expect_false(igraph::are_adjacent(ref, "1", "2"))
    mst <- createClusterMST(y.complex, clusters=clusters.complex, dist.method="mnn")
    expect_true(igraph::are_adjacent(mst, "1", "2"))

    # Quietly ignores weight matrices.
    mat <- factor2matrix(clusters.complex)
    wmst <- createClusterMST(y.complex, clusters=mat, dist.method="mnn")
    expect_identical(mst[], wmst[])

    wmst <- createClusterMST(y.complex, clusters=abs(jitter(mat)), dist.method="mnn")
    expect_identical(mst[], wmst[])
})

set.seed(100102)
test_that("MST construction works with scaled distances", {
    ref <- createClusterMST(y0.simple, clusters=clusters.simple)

    mst1 <- createClusterMST(y0.simple, clusters=clusters.simple, dist.method="scaled.diag")
    expect_identical(ref[] > 0, mst1[] > 0)
    expect_identical(igraph::V(ref)$coordinates, igraph::V(mst1)$coordinates)

    mst2 <- createClusterMST(y0.simple, clusters=clusters.simple, dist.method="scaled.full")
    expect_identical(ref[] > 0, mst2[] > 0)
    expect_identical(igraph::V(ref)$coordinates, igraph::V(mst2)$coordinates)

    mst3 <- createClusterMST(y0.simple, clusters=clusters.simple, dist.method="slingshot")
    expect_identical(ref[] > 0, mst3[] > 0)
    expect_identical(igraph::V(ref)$coordinates, igraph::V(mst3)$coordinates)
    expect_identical(mst2[], mst3[]) # chooses scaled.full.

    # Trying something that requires a bit more... finesse.
    ref <- createClusterMST(y.complex, clusters=clusters.complex)
    expect_false(igraph::are_adjacent(ref, "1", "2"))

    mst <- createClusterMST(y.complex, clusters=clusters.complex, dist.method="scaled.diag")
    expect_true(igraph::are_adjacent(mst, "1", "2"))

    mst <- createClusterMST(y.complex, clusters=clusters.complex, dist.method="scaled.full")
    expect_true(igraph::are_adjacent(mst, "1", "2"))

    mst2 <- createClusterMST(y.complex, clusters=clusters.complex, dist.method="slingshot")
    expect_true(igraph::are_adjacent(mst2, "1", "2"))
    expect_identical(mst2[], mst[]) # chooses scaled.full.

    # Works correctly with weight matrices.
    mat <- factor2matrix(clusters.complex)
    wmst <- createClusterMST(y.complex, clusters=mat, dist.method="slingshot")
    expect_equal(mst[], wmst[])

    wmst <- createClusterMST(y.complex, clusters=abs(jitter(mat)), dist.method="slingshot") # Slightly different result with non-identical weights.
    expect_false(isTRUE(all.equal(mst[], wmst[])))

    # warnings out!
    expect_warning(createClusterMST(y.complex, clusters=clusters.complex, dist.method="slingshot", use.median=TRUE), "unpredictable")
})

