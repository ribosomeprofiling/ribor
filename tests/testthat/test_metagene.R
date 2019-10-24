context("metagene functions")
library(ribor)
file.path <- system.file("extdata", "sample.ribo", package = "ribor")
ribo.object <- create_ribo(file.path)

file.path <- system.file("extdata", "CBS.ribo", package = "ribor")
CBS         <- create_ribo(file.path)

file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")

green <- create_ribo(file.path, rename = rename_default)

green_meta_alias <- get_metagene(green,
                                 "start",
                                 range.lower = 28,
                                 range.upper = 28,
                                 transcript = FALSE,
                                 alias = TRUE,
                                 experiments = get_experiments(green)[1])

green_meta_original <- get_metagene(green,
                                    "start",
                                    range.lower = 28,
                                    range.upper = 28,
                                    transcript = FALSE,
                                    alias = FALSE,
                                    experiments = get_experiments(green)[1])

test_that("get_metagene- alias values preserved",
          expect_true(all(green_meta_alias[, 4] == 
                          green_meta_original[, 4])))

original.names <- unlist(as.character(green_meta_original[, 2]))
actual   <- unlist(as.character(green_meta_alias[, 2]))

expected <- sapply(original.names, rename_default)

test_that("get_metagene- name ordering preserved",
           expect_true(all(actual == 
                           expected)))

meta_1 <- get_metagene(ribo.object,
                       "start",
                       range.lower = 2,
                       range.upper = 5,
                       length = TRUE,
                       transcript = TRUE,
                       experiments = c("Hela_1", "Hela_2"))

actual   <- c(nrow(meta_1), ncol(meta_1)) 
expected <- c(2, 6)

test_that("get_metagene size",
          expect_equal(actual, expected))

expected <- 17
actual <- meta_1[meta_1$experiment == "Hela_1", ][, 3]
test_that("get_metagene, checks the default case",
          expect_equal(actual, 
                       expected))

meta_2 <- get_metagene(ribo.object,
                       "start",
                       range.lower = 2,
                       range.upper = 3,
                       length = TRUE,
                       transcript = FALSE,
                       experiments = c("Hela_1"))

actual   <- c(nrow(meta_2), ncol(meta_2)) 
expected <- c(3, 7)

test_that("get_metagene- size 1",
          expect_equal(actual, expected))

expected <- 4
actual   <- meta_2[meta_2$transcript == "GAPDH", ][["-1"]]
test_that("get_metagene, only condense lengths",
          expect_equal(actual, 
                       expected))

meta_3 <- get_metagene(ribo.object,
                       "stop",
                       range.lower = 3,
                       range.upper = 4,
                       length = FALSE,
                       transcript = TRUE)

expected <- 4
test_that("get_metagene, only condense transcripts",
          expect_equal(meta_3[meta_3$experiment == "Hela_1" &
                              meta_3$length == 3, ][["-2"]], 
                       expected))

meta_4 <- get_metagene(ribo.object,
                       "stop",
                       range.lower = 4,
                       range.upper = 4,
                       length = FALSE,
                       transcript = FALSE)

actual <- meta_4[meta_4$experiment == "Hela_1" &
                   meta_4$transcript == "GAPDH", ][["1"]]
expected <- 1
test_that("get_metagene, preserve everything",
          expect_equal(actual, 
                       expected))


actual   <- c(nrow(meta_4), ncol(meta_4)) 
expected <- c(15, 8)

test_that("get_metagene- size 2",
          expect_equal(actual, expected))

tidy_meta_1 <- get_tidy_metagene(ribo.object,
                                 "start",
                                 2,
                                 3,
                                 length = TRUE)
actual   <- c(nrow(tidy_meta_1), ncol(tidy_meta_1)) 
expected <- c(25, 3)

test_that("get_metagene- size 3",
          expect_equal(actual, expected))


actual <- tidy_meta_1[tidy_meta_1$experiment == "Hela_1" &
                        tidy_meta_1$position == -1, ][["count"]]
expected <- 9

test_that("get_metagene, preserve everything",
          expect_equal(actual, 
                       expected))



tidy_meta_2 <- get_tidy_metagene(ribo.object,
                                 "start",
                                 2,
                                 3,
                                 length = FALSE)
actual <- c(nrow(tidy_meta_2), ncol(tidy_meta_2))
expected <- c(50, 4)

test_that("get_metagene- size 4",
          expect_equal(actual, 
                       expected))

actual <- tidy_meta_2[tidy_meta_2$experiment == "Hela_2" &
                        tidy_meta_2$length == 3 &
                        tidy_meta_2$position == -1, ][["count"]]

expected <- 5
test_that("get_metagene- preserve everything",
          expect_equal(actual, 
                       expected))

CBS_meta <- get_metagene(CBS, "start", 15, 35)
actual <- as.integer(CBS_meta[, 3])
expected <- 1L
test_that("get_metagene- non-tidy random access", 
          expect_equal(actual,
                       expected))

CBS_tidy <- get_tidy_metagene(CBS, "stop", 15, 35)

actual <- CBS_tidy[CBS_tidy$position == -47, 3]
expected <- 1L
test_that("get_tidy_metagene- check random value", expect_equal(actual, expected))
