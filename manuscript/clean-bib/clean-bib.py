from pybtex.database import parse_file

bib_data = parse_file('tmp/library-src.bib')

with open('tmp/unused-refs.txt','r') as f: 
    unused_refs = f.read().strip().split('\n')

for ref in unused_refs:
    del bib_data.entries[ref]

bib_data.to_file('../resources/library-proc.bib', bib_format='bibtex')

