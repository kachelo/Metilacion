context("MethylationExperiment")

test_that(
    "MethylationExperiment",
    {
        expect_error(new("MethylationExperiment"))
    }
)
