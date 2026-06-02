#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(scrapper)
  library(rhdf5)
  library(DelayedArray)
  library(Matrix)
  library(optparse)
  library(yaml)
  library(data.table)
})
cat("OK\n")
