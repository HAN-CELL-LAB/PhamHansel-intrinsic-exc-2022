# VERSION="$(echo $(git rev-parse --abbrev-ref HEAD) | sed 's/-docx//')"
# echo "$VERSION"

cp main.pdf manuscript-PhamHansel.pdf

find . \
  -type d \( -path ./misc -o \
			-path ./archive -o \
			-path ./output -o \
			-path ./clean-bib -o \
			-path ./resources \) -prune -o \
  -name '*.tex' \
  -exec sed -i -E\
	-e 's/\\autoref\{fig([^\{]*)\}/\\textbf{Fig. \\ref{fig\1}}/g' \
	-e 's/\\autoref\{supp([^\{]*)\}/\\textbf{Fig. S\\ref{supp\1}}/g' \
	-e 's/\\autoref\{eq([^\{]*)\}/\\textbf{Eq. \\ref{eq\1}}/g' \
	{} \;

find figures -type f -name 'Fig*.tex' \
	-exec sed -i -E\
	-e 's/\\(HLCh|HLChFig)//g' \
	{} \;

DOCX_FILE=manuscript-PhamHansel.docx
pandoc --defaults front.yml \
	-s main.tex \
	-o "$DOCX_FILE"
	# --verbose &> ./MS-docx.log

cp "$DOCX_FILE" unedited-"$DOCX_FILE"

