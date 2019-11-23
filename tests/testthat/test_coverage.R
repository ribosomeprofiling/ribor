context("coverage functions")
library(ribor)

file.path <- system.file("extdata", "sample.ribo", package = "ribor")
ribo.object <- create_ribo(file.path)

tidy_cov_1 <- get_coverage(ribo.object,
                           name = "MYC",
                           range.lower = 2,
                           range.upper = 5,
                           tidy = TRUE,
                           experiments = "Hela_1")

tidy_cov_2   <- get_coverage(ribo.object,
                        name = "VEGFA",
                        range.lower = 2,
                        range.upper = 5,
                        tidy = TRUE,
                        experiments = "Hela_1")

tidy_cov_3   <- get_coverage(ribo.object,
                        name = "GAPDH",
                        range.lower = 2,
                        range.upper = 5,
                        tidy = TRUE,
                        experiments = c("Hela_1"))

actual <- sum(tidy_cov_1[, 3]) +
          sum(tidy_cov_2[, 3]) +
          sum(tidy_cov_3[, 3])

actual <- unname(actual)
expected <- 118

test_that("get_coverage, tidy - total reads of an experiment",
          expect_equal(actual, expected))

tidy_cov_1 <- get_coverage(ribo.object,
                           name = "GAPDH",
                           range.lower = 2,
                           range.upper = 5,
                           tidy = TRUE,
                           length = FALSE,
                           experiments = "Hela_1")

expected <- c(1, 1, 0, 1)
actual   <- unname(unlist(tidy_cov_1[as.integer(tidy_cov_1$position) == 2, 4]))

test_that("get_coverage, tidy - randomly choose a value",
          expect_equal(actual, expected))

tidy_cov_2 <- get_coverage(ribo.object,
                           name = "MYC",
                           range.lower = 2,
                           range.upper = 5,
                           tidy = TRUE,
                           length = FALSE,
                           experiments = "Hela_1")


actual <- unname(unlist(tidy_cov_2[as.integer(tidy_cov_2$position) == 11, 4]))
expected <- c(2, 1, 3, 7)
test_that("get_coverage, tidy - randomly choose a value",
          expect_equal(actual, expected))

cov_1   <- get_coverage(ribo.object,
                        name = "MYC",
                        range.lower = 2,
                        range.upper = 5,
                        experiments = "Hela_1")

cov_2   <- get_coverage(ribo.object,
                        name = "VEGFA",
                        range.lower = 2,
                        range.upper = 5,
                        experiments = "Hela_1")

cov_3   <- get_coverage(ribo.object,
                        name = "GAPDH",
                        range.lower = 2,
                        range.upper = 5,
                        experiments = c("Hela_1"))

actual <- sum(unlist(as.list(cov_1[, -1]))) + 
          sum(unlist(as.list(cov_2[, -1]))) + 
          sum(unlist(as.list(cov_3[, -1])))

expected <- 118

test_that("get_coverage- total reads of an experiment",
          expect_equal(actual, expected))

actual <- unname(unlist(cov_1[, 12]))
expected <- 13

test_that("get_coverage- correct sum across lengths",
          expect_equal(actual, expected))

cov_4   <- get_coverage(ribo.object,
                        name = "GAPDH",
                        range.lower = 2,
                        range.upper = 2,
                        experiments = c("Hela_1"))

actual   <- sum(unlist(as.list(cov_4[, -1])))
expected <- 14 

test_that("get_coverage- test individual read length",
          expect_equal(actual, expected))

cov_5   <- get_coverage(ribo.object,
                        name = "MYC",
                        range.lower = 3,
                        range.upper = 3,
                        experiments = c("Hela_2"))

actual <- sum(unlist(as.list(cov_5[, -1])))
expected <- 3 

test_that("get_coverage- test individual read length",
          expect_equal(actual, expected))

file.path <- system.file("extdata", "CBS.ribo", package = "ribor")
CBS <- create_ribo(file.path)

CBS_name <- names(CBS@transcript.info)[1]

CBS_cov <- get_coverage(CBS,
                        name = CBS_name,
                        range.lower = 15,
                        range.upper = 35)
actual <- sum(unlist(as.list(CBS_cov[ , -1])))
expected <- 846
test_that("get_coverage- edge case transcript only", expect_equal(actual, expected))
