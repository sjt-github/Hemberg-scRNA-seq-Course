#!/bin/bash

Rscript bookdown::render_book("index.Rmd", "bookdown::gitbook")
# Rscript bookdown::render_book("index.Rmd", "bookdown::pdf_book")
# Rscript bookdown::render_book("index.Rmd", "bookdown::epub_book")
