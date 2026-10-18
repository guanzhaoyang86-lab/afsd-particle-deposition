"""Extract selected reference pages for inspection; source PDF stays intact."""
import argparse
from pathlib import Path
import fitz

parser = argparse.ArgumentParser()
parser.add_argument('pdf')
parser.add_argument('pages', help='comma-separated one-based pages or start:end')
parser.add_argument('--output', required=True)
args = parser.parse_args()
pages = []
for item in args.pages.split(','):
    parts = [int(x) for x in item.split(':')]
    pages.extend(range(parts[0], parts[-1]+1))
doc = fitz.open(args.pdf)
text = '\n\n'.join(f'PDF_PAGE {p}\n'+doc[p-1].get_text() for p in pages)
Path(args.output).write_text(text, encoding='utf-8')
print(f'Extracted {len(pages)} pages to {args.output}')
