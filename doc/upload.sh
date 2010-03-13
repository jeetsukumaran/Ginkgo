#! /bin/bash

rsync -rv --progress build/html/ phylo.bio.ku.edu:/var/www/html/ginkgo/
