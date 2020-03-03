library(testthat)
library(OmnipathR)

op <- tolower(Sys.info()[1])

if (op != "windows"){
    test_check("OmnipathR")
}