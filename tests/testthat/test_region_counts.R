context("region count functions")
library(ribor)

file.path <- system.file("extdata", "sample.ribo", package = "ribor")
ribo.object<- Ribo(file.path)

file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
green <- Ribo(file.path, rename = rename_default)
total.reads <- get_info(green)[["experiment.info"]][,2][1]

green_rc_alias <- get_region_counts(green,
                                    region = "CDS",
                                    range.lower = 28,
                                    range.upper = 28,
                                    transcript = FALSE,
                                    alias = TRUE,
                                    normalize = FALSE,
                                    experiment = experiments(green)[1])

norm.test1 <- 22 * 1000000 / total.reads

green_rc_original <- get_region_counts(green,
                                       region = "CDS",
                                       range.lower = 28,
                                       range.upper = 28,
                                       transcript = FALSE,
                                       alias = FALSE,
                                       experiment = experiments(green)[1])


test_that("get_region_counts- alias values preserved",
          expect_true(all(as.vector(as.character(green_rc_alias[, 4])) ==
                          as.vector(green_rc_original[, 4]))))

original.names <- unlist(as.character(green_rc_original[, 2]))
actual   <- unlist(green_rc_alias[, 2])
expected <- sapply(original.names, rename_default)

test_that("get_region_counts- name ordering preserved",
          expect_true(all(actual == expected)))

green_rc_alias <- get_region_counts(green,
                                    region = "CDS",
                                    range.lower = 28,
                                    range.upper = 28,
                                    transcript = FALSE,
                                    alias = FALSE,
                                    normalize = TRUE,
                                    experiment = experiments(green)[1])

actual <- green_rc_alias[5, 4]
expected <- norm.test1

test_that("get_region_counts- normalize function",
          expect_true(all(actual == expected)))


green_rc <- get_region_counts(green,
                              region = "CDS",
                              range.lower = 28,
                              range.upper = 28,
                              transcript = TRUE,
                              alias = FALSE,
                              normalize = FALSE,
                              experiment = experiments(green)[1])

green_rc_norm <- get_region_counts(green,
                                   region = "CDS",
                                   range.lower = 28,
                                   range.upper = 28,
                                   transcript = TRUE,
                                   alias = FALSE,
                                   normalize = TRUE,
                                   experiment = experiments(green)[1])

expected <- unlist(green_rc[, 3] * 1000000/total.reads)

actual <- unlist(green_rc_norm[, 3])
threshold <- 0.001

test_that("get_region_counts- normalize function",
          expect_true(all(abs(expected - actual) < threshold)))

green_rc_norm <- get_region_counts(green,
                                   region = "CDS",
                                   range.lower = 28,
                                   range.upper = 28,
                                   transcript = FALSE,
                                   alias = TRUE,
                                   normalize = TRUE,
                                   experiment = experiments(green)[1])

actual <- sum(green_rc_norm[, 4])

test_that("get_region_counts- normalize function",
          expect_true(abs(actual - expected) < threshold))

green_rc_norm <- get_region_counts(green,
                                   region = "CDS",
                                   range.lower = 28,
                                   range.upper = 28,
                                   transcript = FALSE,
                                   length = FALSE,
                                   alias = FALSE,
                                   normalize = TRUE,
                                   experiment = experiments(green)[1])

actual <- sum(green_rc_norm[, 5])
test_that("get_region_counts- normalize function",
          expect_true(abs(actual - expected) < threshold))


rc_1 <- get_region_counts(ribo.object,
                          region = c("UtR5", "cDs", "utr3"),
                          2,
                          5,
                          experiment = "Hela_1")

actual <- c(nrow(rc_1), ncol(rc_1))
expected <- c(3, 3)

test_that("get_region_counts- size",
          expect_equal(actual, expected))


rc_all <- get_region_counts(ribo.object,
                            region = c("UtR5", "UTR5j", "cDs", "utr3j", "utr3"),
                             2,
                             5,
                             experiment = "Hela_1")

actual <- c(nrow(rc_all), ncol(rc_all))
expected <- c(5, 3)

test_that("get_region_counts- size",
          expect_equal(actual, expected))

actual   <- sum(rc_all[, 3])
expected <- 118

test_that("get_region_counts- total count",
           expect_equal(actual, expected))

rc_3 <- get_region_counts(ribo.object,
                          region = c("CDS"),
                          2,
                          5,
                          length = TRUE,
                          transcript = FALSE,
                          experiment = c("Hela_1"))

actual <- c(nrow(rc_3), ncol(rc_3))
expected <- c(3, 4)

test_that("get_region_counts- size",
          expect_equal(actual, expected))

actual <- rc_3[1, ][["count"]]
expected <- 13

test_that("get_region-count- preserving transcripts",
          expect_equal(actual, expected))

rc_4 <- get_region_counts(ribo.object,
                          region = c("UtR5j"),
                          2,
                          5,
                          length = FALSE,
                          transcript = TRUE)

actual <- c(nrow(rc_4), ncol(rc_4))
expected <- c(20, 4)

test_that("get_region-count- size",
          expect_equal(actual, expected))

actual <- sum(rc_4$count)
expected <- 102
test_that("get_region_count- preserving lengths",
          expect_equal(actual, expected))



rc_1 <- get_region_counts(ribo.object,
                          region = c("UtR5", "cDs", "utr3"),
                          2,
                          5,
                          compact = FALSE,
                          experiment = "Hela_1")

actual <- c(nrow(rc_1), ncol(rc_1))
expected <- c(3, 3)

test_that("get_region_counts- noncompact, size",
          expect_equal(actual, expected))


rc_all <- get_region_counts(ribo.object,
                            region = c("UtR5", "UTR5j", "cDs", "utr3j", "utr3"),
                            2,
                            5,
                            compact = FALSE,
                            experiment = "Hela_1")

actual <- c(nrow(rc_all), ncol(rc_all))
expected <- c(5, 3)

test_that("get_region_counts- noncompact, size",
          expect_equal(actual, expected))

actual   <- sum(rc_all[, 3])
expected <- 118

test_that("get_region_counts- noncompact, total count",
          expect_equal(actual, expected))

rc_3 <- get_region_counts(ribo.object,
                          region = c("CDS"),
                          2,
                          5,
                          length = TRUE,
                          transcript = FALSE,
                          compact = FALSE,
                          experiment = c("Hela_1"))

actual <- c(nrow(rc_3), ncol(rc_3))
expected <- c(3, 4)

test_that("get_region_counts- size",
          expect_equal(actual, expected))

actual <- rc_3[1, ][["count"]]
expected <- 13

test_that("get_region-count- noncompact, preserving transcripts",
          expect_equal(actual, expected))

rc_4 <- get_region_counts(ribo.object,
                          region = c("UtR5j"),
                          2,
                          5,
                          compact = FALSE,
                          length = FALSE,
                          transcript = TRUE)

actual <- c(nrow(rc_4), ncol(rc_4))
expected <- c(20, 4)

test_that("get_region-count- noncompact, size",
          expect_equal(actual, expected))

actual <- sum(rc_4$count)
expected <- 102
test_that("get_region_count- noncompact, preserving lengths",
          expect_equal(actual, expected))

ld <- get_length_distribution(ribo.object, "CDS", compact = FALSE)
rc <- get_region_counts(ribo.object, region = "CDS", compact = FALSE, length = FALSE)

test_that("get length distribution- compact, returns data.frame",
          expect_true(is.data.frame(ld)))

test_that("get length distribution- compact, correctness for counts",
          expect_equal(ld$count, rc$count))

ld <- get_length_distribution(ribo.object, "CDS")
rc <- get_region_counts(ribo.object, region = "CDS", compact = FALSE, length = FALSE)
test_that("get length distribution- compact, correctness for counts",
          expect_equal(ld$count, rc$count))

