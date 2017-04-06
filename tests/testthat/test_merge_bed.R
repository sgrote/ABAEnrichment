
# context("merge_bed.R")

# input
chroms = c(1,1,1,2,2,2,2,2,2,3,3,4,4,4,5,5,5,5)
starts = c(2,4,10,2,2,2,7,7,7,1,3,8,2,5,1,7,2,3)
stops = c(3,10,12,5,5,5,9,8,10,4,5,10,7,6,2,9,4,5)
bed = data.frame(chroms,starts,stops)

# expected output
chroms = c(1,1,2,2,3,4,4,5,5)
starts = c(2,4,2,7,1,2,8,1,7)
stops = c(3,12,5,10,5,7,10,5,9)
expected = data.frame(chroms,starts,stops)

test_that("candidate and background regions get correctly merged",{
	expect_that(merge_bed(bed), equals(expected))
})


# error-input
chroms = c(1,3,3,4)
starts = c(2,4,10,2)
stops = c(3,10,7,5)
bed = data.frame(chroms,starts,stops)

test_that("error is thrown when start > stop",{
	expect_that(merge_bed(bed), throws_error("Genomic regions must be defined with start < stop."))
})
