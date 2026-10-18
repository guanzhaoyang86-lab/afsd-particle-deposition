"""Render the existing AFSD consultation Markdown to a Chinese, selectable PDF.
Presentation only: never imports the solver or modifies frozen project files.
"""
from pathlib import Path
from html import escape
import hashlib
import json
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, KeepTogether
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'docs/AFSD_CONTEXT_FOR_GPT.md'
TARGET = ROOT/'output/pdf/AFSD_Context_and_GPT_Prompts.pdf'
QA = ROOT/'tmp/pdfs/afsd_context_qa'
TARGET.parent.mkdir(parents=True, exist_ok=True)
QA.mkdir(parents=True, exist_ok=True)

pdfmetrics.registerFont(TTFont('YaHei', r'C:\Windows\Fonts\msyh.ttc', subfontIndex=0))
pdfmetrics.registerFont(TTFont('YaHeiBold', r'C:\Windows\Fonts\msyhbd.ttc', subfontIndex=0))
pdfmetrics.registerFont(TTFont('Consolas', r'C:\Windows\Fonts\consola.ttf'))

NAVY = colors.HexColor('#17324D')
INK = colors.HexColor('#243342')
MUTED = colors.HexColor('#657687')
PALE = colors.HexColor('#F1F5F9')
BORDER = colors.HexColor('#D6E0E8')
WIDTH, HEIGHT = A4
MARGIN = 45
CONTENT = WIDTH - 2*MARGIN

body = ParagraphStyle('Body', fontName='YaHei', fontSize=9.5, leading=15.1,
    textColor=INK, spaceAfter=6, wordWrap='CJK', splitLongWords=True)
styles = {
    'body': body,
    'title': ParagraphStyle('Title', parent=body, fontName='YaHeiBold',
        fontSize=19, leading=28, textColor=NAVY, spaceAfter=14),
    'h2': ParagraphStyle('H2', parent=body, fontName='YaHeiBold',
        fontSize=12.6, leading=19, textColor=NAVY, spaceBefore=12,
        spaceAfter=7, keepWithNext=True),
    'h3': ParagraphStyle('H3', parent=body, fontName='YaHeiBold',
        fontSize=10.6, leading=16, textColor=NAVY, spaceBefore=8,
        spaceAfter=6, keepWithNext=True),
    'code': ParagraphStyle('Code', parent=body, fontName='Consolas',
        fontSize=8.4, leading=12.8, spaceAfter=0),
    'cell': ParagraphStyle('Cell', parent=body, fontSize=8.5, leading=13,
        spaceAfter=0),
    'th': ParagraphStyle('TH', parent=body, fontName='YaHeiBold',
        fontSize=8.5, leading=13, textColor=colors.white, spaceAfter=0),
    'list': ParagraphStyle('List', parent=body, leftIndent=12, firstLineIndent=-10,
        spaceAfter=5),
}

def normalize(text):
    return text.replace('\u2011', '-').replace('\u2013', '-').replace('\u2014', '-')

def rich(text):
    text = normalize(text)
    chunks = re.split(r'(\[[^\]]+\]\(https?://[^)]+\))', text)
    out = []
    for chunk in chunks:
        match = re.fullmatch(r'\[([^\]]+)\]\((https?://[^)]+)\)', chunk)
        if match:
            out.append('<link href="'+escape(match[2], quote=True)+'" color="#166184">'+escape(match[1])+'</link>')
        else:
            out.append(escape(chunk))
    return ''.join(out)

def code_text(text):
    # Consolas for formulas; embedded CJK for explanatory phrases inside fences.
    return ''.join('<font name="YaHei">'+escape(c)+'</font>' if ord(c)>0x2E7F else escape(c) for c in normalize(text))

def p(text, kind='body'):
    return Paragraph(rich(text), styles[kind])

def code_box(lines):
    rows = [[Paragraph(code_text(line) if line else '&#160;', styles['code'])] for line in lines]
    table = Table(rows, colWidths=[CONTENT], hAlign='LEFT', splitByRow=1, splitInRow=0)
    table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), PALE),
        ('LINEBEFORE', (0,0), (0,-1), 2, colors.HexColor('#7A9CB8')),
        ('LEFTPADDING', (0,0), (-1,-1), 10),
        ('RIGHTPADDING', (0,0), (-1,-1), 10),
        ('TOPPADDING', (0,0), (-1,-1), 2),
        ('BOTTOMPADDING', (0,0), (-1,-1), 2),
    ]))
    # Keep compact formulas intact; longer prompts can split only between rows.
    return KeepTogether([table]) if table.wrap(CONTENT, HEIGHT)[1] <= 270 else table

class Document(SimpleDocTemplate):
    def afterFlowable(self, flowable):
        level = getattr(flowable, '_outline_level', None)
        if level is not None:
            key = 'section_'+str(self.page)+'_'+str(getattr(self, '_outline_count', 0))
            self._outline_count = getattr(self, '_outline_count', 0)+1
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(flowable.getPlainText(), key, level, False)

def furniture(canvas, doc):
    canvas.saveState()
    canvas.setFont('YaHei', 7.5)
    canvas.setFillColor(MUTED)
    canvas.drawString(MARGIN, HEIGHT-28, 'AFSD | 独立咨询上下文与GPT提问指南')
    canvas.setStrokeColor(BORDER)
    canvas.setLineWidth(.5)
    canvas.line(MARGIN, HEIGHT-34, WIDTH-MARGIN, HEIGHT-34)
    canvas.line(MARGIN, 34, WIDTH-MARGIN, 34)
    canvas.drawString(MARGIN, 22, 'APPROXIMATE / 科学审查快照 / 非论文复现结果')
    canvas.drawRightString(WIDTH-MARGIN, 22, '第 '+str(doc.page)+' 页')
    canvas.restoreState()

source_text = SOURCE.read_text(encoding='utf-8')
lines = source_text.splitlines()
flow = []
i = 0
section_count = 0
while i < len(lines):
    line = lines[i].strip()
    if not line:
        i += 1
        continue
    if line.startswith('~~~'):
        code = []
        i += 1
        while i < len(lines) and not lines[i].strip().startswith('~~~'):
            code.append(lines[i].rstrip())
            i += 1
        if flow and isinstance(flow[-1], Paragraph):
            flow[-1].keepWithNext = True
        flow.append(code_box(code))
        flow.append(Spacer(1, 8))
        i += 1
        continue
    if line.startswith('|'):
        raw_rows = []
        while i < len(lines) and lines[i].strip().startswith('|'):
            cells = [s.strip() for s in lines[i].strip().strip('|').split('|')]
            if not all(re.fullmatch(r':?-+:?', s) for s in cells):
                raw_rows.append(cells)
            i += 1
        n = len(raw_rows[0])
        ratios = [0.31, 0.69] if n==2 else [0.15,0.11,0.13,0.11,0.15,0.35]
        rows = [[p(cell, 'th' if r==0 else 'cell') for cell in row] for r,row in enumerate(raw_rows)]
        table = Table(rows, colWidths=[CONTENT*x for x in ratios],
            repeatRows=1, hAlign='LEFT', splitByRow=1)
        table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,0),NAVY),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white, PALE]),
            ('VALIGN',(0,0),(-1,-1),'TOP'),
            ('BOX',(0,0),(-1,-1),.45,BORDER),
            ('LINEBELOW',(0,0),(-1,0),.6,NAVY),
            ('LEFTPADDING',(0,0),(-1,-1),7),
            ('RIGHTPADDING',(0,0),(-1,-1),7),
            ('TOPPADDING',(0,0),(-1,-1),6),
            ('BOTTOMPADDING',(0,0),(-1,-1),6),
        ]))
        flow.extend([table, Spacer(1,8)])
        continue
    if line.startswith('# '):
        flow.append(p(line[2:], 'title'))
        i += 1
        continue
    if line.startswith('## '):
        title = p(line[3:], 'h2')
        title._outline_level = 0
        flow.append(title)
        section_count += 1
        i += 1
        continue
    if line.startswith('### '):
        title = p(line[4:], 'h3')
        title._outline_level = 1
        flow.append(title)
        i += 1
        continue
    if line.startswith('- ') or re.match(r'^\d+\. ', line):
        flow.append(p(('• '+line[2:]) if line.startswith('- ') else line, 'list'))
        i += 1
        continue
    block = [line]
    i += 1
    while i < len(lines):
        nxt = lines[i].strip()
        if not nxt or nxt.startswith(('#','|','~~~','- ')) or re.match(r'^\d+\. ', nxt):
            break
        block.append(nxt)
        i += 1
    flow.append(p(' '.join(block)))

document = Document(str(TARGET), pagesize=A4, rightMargin=MARGIN, leftMargin=MARGIN,
    topMargin=47, bottomMargin=45, title='AFSD项目上下文与GPT提问指南',
    author='AFSD Project', subject='Independent scientific review context and consultation prompts',
    allowSplitting=1)
document.build(flow, onFirstPage=furniture, onLaterPages=furniture)

reader = PdfReader(str(TARGET))
extracted = '\n'.join(page.extract_text() or '' for page in reader.pages)
(QA/'extracted_text.txt').write_text(extracted, encoding='utf-8')
compact = re.sub(r'\s+', '', extracted)
anchors = [
 'negative midpoint friction dissipation', '0.063371356',
 '1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5', '35通过、3失败',
 '2.63e12', '第一轮提问词', '第二轮提问词', '对泛泛回答的追问',
 'Theta = 1 - (1 + exp(-c1*Tstar))^(-1/c2)',
]
missing = [x for x in anchors if re.sub(r'\s+', '', x) not in compact]
heads = [x[3:].strip() for x in lines if x.startswith('## ')]
missing_heads = [x for x in heads if re.sub(r'\s+', '', normalize(x)) not in compact]
report = dict(output=str(TARGET), pages=len(reader.pages), sections=section_count,
    bytes=TARGET.stat().st_size, source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    output_sha256=hashlib.sha256(TARGET.read_bytes()).hexdigest(),
    missing_anchors=missing, missing_sections=missing_heads,
    text_extractable=bool(compact), original_physics_files_modified=False)
(QA/'build_manifest.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=True, indent=2))
if missing or missing_heads or section_count != 18:
    raise SystemExit('PDF content check failed')
