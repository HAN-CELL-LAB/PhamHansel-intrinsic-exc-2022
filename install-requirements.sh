#!/bin/bash

MATLAB_LIST="matlab-requirements.txt"
MATLAB_DIR="matlab-ext"
MAIN_DIR="$(pwd)"

git config -f .gitmodules submodule.path/to/submodule.shallow true

cat "$MATLAB_LIST" | \
  while read repo_opts; 
  do 
	repo_url="$(echo "$repo_opts" | awk '{print $1}')"
	repo_sha="$(echo "$repo_opts" | awk '{print $2}')"
    pkg_name="$(basename "$repo_url")"
    submod_dir="$MATLAB_DIR/$pkg_name"
    
    [ -d "$submod_dir" ] && echo "$submod_dir also exists, skipping ..." && continue 

	echo -e "\E[7m  \E[0m Installing \E[1m$pkg_name\E[0m from <\E[2;3;4m$repo_url\E[0m> ..."
    # git subtree add --prefix "$MATLAB_DIR/$pkg_name" "$repo_url" master --squash
    git submodule add --depth 1 "$repo_url" "$submod_dir"
    
    cd "$submod_dir" && pwd
    git fetch --depth=1 "$repo_url" "$repo_sha"
    git checkout "$repo_sha"
    
    cd "$MAIN_DIR"

	echo -e "\E[7m  \E[0;3m ...  $pkg_name installed!\E[0m"
    echo 
  done


