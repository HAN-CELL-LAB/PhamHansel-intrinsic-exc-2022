# Notes on cleaning unused references in `library.bib`

- Requirements:
  - `pdflatex` 
  - [`checkcites`](https://www.ctan.org/pkg/checkcites)
  - [`pybtex`](https://pybtex.org/) (either via `virtualenv`, `conda` or whatever). Here I used `virtualenv` for a quick but contained installation and usage.
- Run:
  
  ``` bash 
  # if have not created venv and install pybtex
  ./create-venv.sh
  # then run this
  ./clean-bib.sh
  ```
