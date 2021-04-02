# Install all of the packages used in this repository (up to ~30 min run time)

# NOTE: Ubuntu (Linux) dependencies: libxml2, libpng, libfortran, libgdal
#
# - tidyverse and xml2 need libxml2
# - Hmisc, ggmap, funModeling need libpng
# - lme4, VCA, stars, sf need libgfortran
# - rgdal needs libgdal

# Create personal library folder if needed (avoids interactive prompt)
mylib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(mylib)) dir.create(mylib, recursive = T, showWarnings = F)

# Install pacman if needed (which will be used to install other packages)
if (!suppressWarnings(require("pacman"))) install.packages("pacman")

# Install LaTeX environment (for rendering PDFs from RMarkdown, etc.)
pdflatex_ver <- try(system("pdflatex -v", intern = T, wait = T), silent = T)
pdflatex_ver <- grep("^pdfTeX", pdflatex_ver, value = T)
if (!(exists("pdflatex_ver") & length(pdflatex_ver) > 0)) {
  pacman::p_load(tinytex)
  if (!dir.exists(tinytex_root(error = F))) tinytex::install_tinytex()
}

# Install limma (for MKmisc) from Bioconductor using BiocManager
pacman::p_load(BiocManager)
if (!suppressWarnings(require("limma"))) BiocManager::install("limma")

# Install other packages with pacman
pacman::p_load(plyr, reshape2, tictoc, stars, sp, sf, hms, slider, feasts)
pacman::p_load(tidyverse, lubridate, broom, rgdal, foreign, downloader)
pacman::p_load(knitr, formatR, ggrepel, Hmisc, EnvStats, codetools, egg)
pacman::p_load(multcomp, modelr, car, lme4, VCA, parallel, NADA, ggmap)  
pacman::p_load(geoR, maps, scatterplot3d, funModeling, scales, akima)
pacman::p_load(MKmisc, tseries, xts, lubridate, tsibble) 
