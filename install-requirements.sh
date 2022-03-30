#!/bin/bash

MATLAB_LIST="matlab-requirements.txt"
MATLAB_DIR="matlab-ext"

git config -f .gitmodules submodule.path/to/submodule.shallow true

cat "$MATLAB_LIST" | \
  while read repo_url; 
  do 
	pkg_name="$(basename $repo_url)"
	echo -e "\E[7m  \E[0m Installing \E[1m$pkg_name\E[0m from <\E[2;3;4m$repo_url\E[0m> ..."
    # git subtree add --prefix "$MATLAB_DIR/$pkg_name" "$repo_url" master --squash
    git submodule add --depth 1 "$repo_url" "$MATLAB_DIR/$pkg_name" 
	echo -e "\E[7m  \E[0;3m ...  $pkg_name installed!\E[0m"
    echo 
  done


