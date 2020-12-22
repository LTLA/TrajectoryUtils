# This tests the createClusterMST function on some complicated graphs.
# library(testthat); library(TrajectoryUtils); source("setup.R"); source("test-create-mst.R")

set.seed(1000)

test_that("endpoints and outgroup arguments work together as expected", {
    # simple branch point
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
    
    # simple branch point, closer endpoints
    y <- rbind(A=c(0, 0), B=c(1, .5), C=c(1, -.5)) 
    mst <- createClusterMST(y, cluster=NULL,
                            endpoints = c("B", "C"))
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
    expect_identical(vertices[igraph::degree(mst)==2], "A")
    # this one makes a dyad, but I'm not sure it should
    mst2 <- createClusterMST(y, cluster=NULL,
                             endpoints = c("B", "C"), outgroup = 10)
    expect_identical(mst[], mst2[])
    
    # branch point plus outlier
    y <- rbind(A=c(0, 0), B=c(1, 1), C=c(1, -1), D=c(6, 0)) 
    mst <- createClusterMST(y, cluster=NULL,
                            endpoints = c("B", "C"), outgroup = TRUE)
    
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==0], "D")
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
    expect_identical(vertices[igraph::degree(mst)==2], "A")
    
    # branch point plus outlier, closer endpoints
    y <- rbind(A=c(0, 0), B=c(1, .5), C=c(1, -.5), D=c(6, 0)) 
    mst <- createClusterMST(y, cluster=NULL,
                            endpoints = c("B","C"), outgroup = TRUE)
    # again makes a dyad, still not convinced
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
    expect_identical(vertices[igraph::degree(mst)==2], "A")

    # three endpoints => dyad between the closest two
    y <- rbind(A=c(0, 0), B=c(1, .5), C=c(1, -.5)) 
    mst <- createClusterMST(y, cluster=NULL,
                            endpoints = c("A", "B", "C"))
    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==0], "A")
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C"))
    
    mst2 <- createClusterMST(y, cluster=NULL,
                             endpoints = c("B", "C"), outgroup = 10)
    expect_identical(mst[], mst2[])
})

