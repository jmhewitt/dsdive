context("min.tx.R")

test_that("minimum up/down transitions wrt. stage breaks", {
 
  expect_equal(min.tx(d0 = 12, df = 12, ns = 2), 2)
  expect_equal(min.tx(d0 = 12, df = 11, ns = 2), 3)
  expect_equal(min.tx(d0 = 12, df = 10, ns = 2), 2)
  expect_equal(min.tx(d0 = 12, df = 9, ns = 2), 3)
  
  expect_equal(min.tx(d0 = 12, df = 12, ns = 1), 2)
  expect_equal(min.tx(d0 = 12, df = 11, ns = 1), 1)
  expect_equal(min.tx(d0 = 12, df = 10, ns = 1), 2)
  expect_equal(min.tx(d0 = 12, df = 9, ns = 1), 3)
  
  expect_equal(min.tx(d0 = 12, df = 12, ns = 3), 4)
  expect_equal(min.tx(d0 = 12, df = 11, ns = 3), 3)
  expect_equal(min.tx(d0 = 12, df = 10, ns = 3), 4)
  expect_equal(min.tx(d0 = 12, df = 9, ns = 3), 3)
  expect_equal(min.tx(d0 = 12, df = 8, ns = 3), 4)
  expect_equal(min.tx(d0 = 12, df = 7, ns = 3), 5)
  
})
