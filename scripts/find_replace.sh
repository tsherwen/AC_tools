#!/bin/bash
#~/.scripts/git-sub
#Author: Khaja Minhajuddin <minhajuddin@cosmicvent.com>
#script which does a global search and replace in the git repository
#it takes two arguments
#e.g. git sub OLD NEW

old=$1
new=$2

for file in $(git grep $old | cut -d':'  -f 1 | uniq)
do
  echo "replacing '$old' with '$new' in '$file'"
  sed -i -e "s/$old/$new/g" $file
done