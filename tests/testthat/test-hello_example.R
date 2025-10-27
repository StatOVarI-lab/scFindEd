test_that("say_hello returns correct greeting", {
  expect_equal(say_hello("Miles"), "Hello, Miles !")
})

test_that("say_hello handles different names", {
  expect_equal(say_hello("Magda"), "Hello, Magda !")
})
