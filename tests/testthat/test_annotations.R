context("annotation functions")
library(ribor)

file.path <- system.file("extdata", "sample.ribo", package = "ribor")
ribo.object <- create_ribo(file.path)

lengths <- get_region_lengths(ribo.object)

UTR5 <- c(2, 1, 0)
UTR5J <- c(5, 5, 5)
CDS <- c(5, 7, 2)
UTR3J <- c(5, 5, 5)
UTR3 <- c(3,4,5)

actual <- unname(unlist(lengths[, "UTR5"]))
expected <- UTR5
test_that("get_region_lengths, UTR5",
          expect_equal(actual, expected))

actual <- unname(unlist(lengths[, "UTR5J"]))
expected <- UTR5J
test_that("get_region_lengths, UTR5J",
          expect_equal(actual, expected))

actual <- unname(unlist(lengths[, "CDS"]))
expected <- CDS
test_that("get_region_lengths, CDS",
          expect_equal(actual, expected))

actual <- unname(unlist(lengths[, "UTR3J"]))
expected <- UTR3J
test_that("get_region_lengths, UTR3J",
          expect_equal(actual, expected))

actual <- unname(unlist(lengths[, "UTR3"]))
expected <- UTR3
test_that("get_region_lengths, UTR5J",
          expect_equal(actual, expected))

UTR5_start <- c(1, 1, NA)
UTR5_stop <- c(2, 1, NA)
UTR5J_start <- c(3, 2, 1)
UTR5J_stop <- c(7, 6, 5)
CDS_start <- c(8, 7, 6)
CDS_stop <- c(12, 13, 7)
UTR3J_start <- c(13, 14, 8)
UTR3J_stop <- c(17, 18, 12)
UTR3_start <- c(18,19,13)
UTR3_stop <- c(20, 22, 17)

coordinates <- get_region_coordinates(ribo.object)

actual <- unname(unlist(coordinates[, c("UTR5_start", "UTR5_stop")]))
expected <- c(UTR5_start, UTR5_stop)
test_that("get_region_coordinates, UTR5", expect_equal(actual, expected))

actual <- unname(unlist(coordinates[, c("UTR5J_start", "UTR5J_stop")]))
expected <- c(UTR5J_start, UTR5J_stop)
test_that("get_region_coordinates, UTR5J", expect_equal(actual, expected))

actual <- unname(unlist(coordinates[, c("CDS_start", "CDS_stop")]))
expected <- c(CDS_start, CDS_stop)
test_that("get_region_coordinates, CDS", expect_equal(actual, expected))

actual <- unname(unlist(coordinates[, c("UTR3J_start", "UTR3J_stop")]))
expected <- c(UTR3J_start, UTR3J_stop)
test_that("get_region_coordinates, UTR3J", expect_equal(actual, expected))

actual <- unname(unlist(coordinates[, c("UTR3_start", "UTR3_stop")]))
expected <- c(UTR3_start, UTR3_stop)
test_that("get_region_coordinates, UTR3", expect_equal(actual, expected))