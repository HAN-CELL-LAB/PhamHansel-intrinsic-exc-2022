#!/bin/bash


WORK_DIR="$(pwd)"
SRC_DIR="../resources"
TMP_DIR="$WORK_DIR/tmp"
VENV_DIR="$WORK_DIR/venv"


UNUSED_FILE="$TMP_DIR/unused-refs.txt"
BIBSRC_FILE="$SRC_DIR/library.bib"
TMPBIB_FILE="$TMP_DIR/library-src.bib"

[[ -d "$TMP_DIR" ]] || mkdir "$TMP_DIR"

# Compile to get main.aux
cd ..
pdflatex main.tex 2>&1 > /dev/null

# Checkcites for unused references
checkcites --unused main.aux | \
	sed '/^[^(=>)]/d;/^$/d;s/=>\s*//g' | \
	sort > "$UNUSED_FILE"

cd "$WORK_DIR"

# Copy src -> tmp for proc 
cp "$BIBSRC_FILE" "$TMPBIB_FILE"

# Optional: Virtual environemnt if needed to use pybtex for cleaning
if [[ -d "$VENV_DIR" ]]; then
	source "$VENV_DIR/bin/activate"
fi

python clean-bib.py
# deactivate

# rm -rf "$TMP_DIR"

