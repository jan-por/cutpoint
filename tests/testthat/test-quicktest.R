
test_that("data1 available", {
   data("data1")

   expect_true( exists("data1") )
   expect_equal(ncol(data1), 5)
   expect_equal(nrow(data1), 100)
   expect_true(is.integer(data1$event))
   expect_true(is.numeric(data1$biomarker))
   expect_true(is.numeric(data1$covariate_1))
   expect_true(is.numeric(data1$covariate_2))
   expect_true(is.numeric(data1$time))
})

test_that("data2_ushape available", {
   data("data2_ushape")

   expect_true( exists("data2_ushape") )
   expect_equal(ncol(data2_ushape), 4)
   expect_equal(nrow(data2_ushape), 200)
   expect_true(is.integer(data2_ushape$event))
   expect_true(is.numeric(data2_ushape$biomarker))
   expect_true(is.numeric(data2_ushape$covariate_1))
   expect_true(is.numeric(data2_ushape$time))
})
