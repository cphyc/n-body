#!/bin/sh

# Perform a synchronization with the mesopsl1 cluster, in the n-body
rsync -ahz --progress . mesopsl1.obspm.fr:"~/n-body" --include="*.f90" --exclude=".git/"
