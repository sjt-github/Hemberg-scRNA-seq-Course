#!/bin/bash

# configure your name and email if you have not done so
git config --global user.email "vladimir.yu.kiselev@gmail.com"
git config --global user.name "wikiselev"

WORKSPACE=$1

# get the docker
docker pull hemberglab/scrna.seq.course:latest
# run the docker
docker run hemberglab/scrna.seq.course:latest

# copy files from the docker
alias dl='docker ps -l -q'
docker cp `dl`:_book $WORKSPACE/scRNA.seq.course/
rm -r docs
mv _book docs
# mkdir docs/blischak
# cp blischak/umi.rds docs/blischak/
# cp blischak/reads.rds docs/blischak/

# push changes to the website
git add docs/*
git commit -m "update the course website"
git push origin master

# cleanup after docker usage
docker rm $(docker ps -a -q)
docker rmi $(docker images -q)



# 
# # clone the repository to the book-output directory
# git clone -b gh-pages \
#   https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git \
#   book-output
# cd book-output
# mkdir blischak
# cp ../blischak/umi.rds blischak/
# cp ../blischak/reads.rds blischak/
# cp -r ../_book/* ./
# rm -rf _bookdown_files
# git add *
# git commit -m "Build the book"
# git push origin HEAD:gh-pages
