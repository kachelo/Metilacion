context("MethylationExperiment")

test_that(
    "MethylationExperiment cannot be created",
    {
        expect_error(new("MethylationExperiment"))
    }
)
