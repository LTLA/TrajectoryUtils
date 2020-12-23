# This tests the createClusterMST function on some complicated graphs.
# library(testthat); library(TrajectoryUtils); source("test-mst-endout.R")

set.seed(1000)
test_that("MST construction endpoints interact properly with outgroups", {
    y <- rbind(A=c(0, 1), B=c(0, 1.5), C=c(0, 3.5), D=c(0, 4)) 
    ref <- createClusterMST(y, endpoint=c("B", "C"), clusters=NULL)
    expect_true(igraph::are_adjacent(ref, "A", "D"))
    expect_identical(igraph::components(ref)$no, 1L)

    out <- createClusterMST(y, endpoint=c("B", "C"), outgroup=TRUE, clusters=NULL)
    expect_true(igraph::are_adjacent(out, "A", "B"))
    expect_true(igraph::are_adjacent(out, "C", "D"))
    expect_identical(igraph::components(out)$no, 2L)

    # This example is carefully designed with a spacing of C and D such that it
    # only JUST causes outgroup formation, so it'll fail if endpoints= screws
    # the outgroup= distance calculation.
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 4), D=c(0, 10.001))
    ref <- createClusterMST(y, clusters=NULL, outgroup=TRUE)
    out <- createClusterMST(y, endpoint=c("A", "D"), clusters=NULL, outgroup=TRUE)
    expect_identical(ref[], out[])
})

test_that("endpoints and outgroup arguments work for a simple branch point", {
    y <- rbind(A=c(0, 0), B=c(1, 1), C=c(1, -1)) 
    mst <- createClusterMST(y, cluster=NULL,
                            endpoints = c("B", "C"))
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
    expect_identical(vertices[igraph::degree(mst)==2], "A")
    
    mst2 <- createClusterMST(y, cluster=NULL,
                             endpoints = c("B", "C"), outgroup = 10)
    expect_identical(mst[], mst2[])
})
    
test_that("endpoints and outgroup work for a simple branch point with closer endpoints", {
    y <- rbind(A=c(0, 0), B=c(1, .5), C=c(1, -.5)) 
    mst <- createClusterMST(y, cluster=NULL,
                             endpoints = c("B", "C"), outgroup = 10)
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
    expect_identical(vertices[igraph::degree(mst)==2], "A")
})
    
test_that("endpoints and outgroup work for a branch point with an outlier", {
    y <- rbind(A=c(0, 0), B=c(1, 1), C=c(1, -1), D=c(6, 0)) 
    mst <- createClusterMST(y, cluster=NULL,
                            endpoints = c("B", "C"), outgroup = TRUE)
    
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==0], "D")
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
    expect_identical(vertices[igraph::degree(mst)==2], "A")
})

test_that("endpoints and outgroup work for a branch point with an outlier, closer endpoints", {
    y <- rbind(A=c(0, 0), B=c(1, .5), C=c(1, -.5), D=c(6, 0)) 
    mst <- createClusterMST(y, cluster=NULL,
                            endpoints = c("B","C"), outgroup = TRUE)
    # again makes a dyad, still not convinced
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
    expect_identical(vertices[igraph::degree(mst)==2], "A")
})

test_that("endpoints and outgroup work for three endpoints with a dyad between the closest two", {
    y <- rbind(A=c(0, 0), B=c(1, .5), C=c(1, -.5)) 
    mst <- createClusterMST(y, cluster=NULL,
                             endpoints = c("B", "C"), outgroup = 10)
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==2], "A")
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
})

test_that("Unsolvable errors show up properly", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    expect_error(createClusterMST(y[1:3,], endpoint=rownames(y)[1:3], outgroup=TRUE, clusters=NULL), "could not determine")
    expect_error(createClusterMST(y, endpoint=rownames(y), outgroup=TRUE, clusters=NULL), "could not determine")
})
