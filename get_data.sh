#!/bin/sh

# Get data files from instructor's website and store in Datasets folder

wget -r -np -nc -nH -R "index.html*" --cut-dirs=2 -e robots=off \
  https://faculty.washington.edu/sheppard/envh556/
