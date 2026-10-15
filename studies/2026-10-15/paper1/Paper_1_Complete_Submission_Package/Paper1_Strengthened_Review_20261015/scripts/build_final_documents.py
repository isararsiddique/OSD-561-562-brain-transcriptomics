#!/usr/bin/env python3
"""Editable journal manuscript and supplement with explicit figure/equation sizing."""
import json
import re
from pathlib import Path
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from PIL import Image
HERE=Path(__file__).resolve().parents[1];R=HERE/'results';eq=json.loads((HERE/'figures/equations/manifest.json').read_text())

def runtext(p,text,bold=False,size=None):
    r=p.add_run(text);r.bold=bold
    if size:r.font.size=Pt(size)
    r.font.color.rgb=RGBColor(0,0,0)
    return r

def mathpic(doc,q,k):
    tag=re.search(r'\((S?[0-9]+)\)\.?$',q['text'])[1]
    x=next(e for e in eq if ('\\tag{'+tag+'}') in e['latex']);p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before=Pt(6);p.paragraph_format.space_after=Pt(8)
    shape=p.add_run().add_picture(str(HERE/'figures/equations'/(x['filename']+'.png')),width=Inches(x['width_inches']))
    shape._inline.docPr.set('descr',q['text']+'; LaTeX: '+x['latex'])

def style(doc):
    s=doc.sections[0];s.page_width=Inches(8.5);s.page_height=Inches(11)
    s.left_margin=s.right_margin=Inches(.8);s.top_margin=Inches(.72);s.bottom_margin=Inches(.72)
    s.header_distance=s.footer_distance=Inches(.3)
    base=doc.styles['Normal'];base.font.name='Times New Roman';base.font.size=Pt(12);base.font.color.rgb=RGBColor(0,0,0)
    base.paragraph_format.line_spacing=1.16;base.paragraph_format.space_after=Pt(6)
    for name,n in [('Title',16),('Heading 1',13),('Heading 2',12)]:
        st=doc.styles[name];st.font.name='Times New Roman';st.font.size=Pt(n);st.font.bold=True;st.font.color.rgb=RGBColor(0,0,0)
        st.paragraph_format.space_before=Pt(12 if name!='Title' else 0);st.paragraph_format.space_after=Pt(6);st.paragraph_format.keep_with_next=True
    for st in doc.styles:
        for border in st.element.xpath('.//w:pBdr'):
            border.getparent().remove(border)
        for fonts in st.element.xpath('.//w:rFonts'):
            for attr in ['asciiTheme','hAnsiTheme','eastAsiaTheme','cstheme']:
                fonts.attrib.pop(qn('w:'+attr),None)
    h=s.header.paragraphs[0];h.alignment=WD_ALIGN_PARAGRAPH.RIGHT;runtext(h,'Spaceflight brain transcriptomics',size=9)
    f=s.footer.paragraphs[0];f.alignment=WD_ALIGN_PARAGRAPH.CENTER
    fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');f._p.append(fld)
    doc.core_properties.author='; '.join(x['name']for x in json.loads((HERE/'authors.json').read_text())['authors'])
    doc.core_properties.subject='Computational secondary analysis for journal review'
    doc.core_properties.keywords='spaceflight, transcriptomics, paired analysis, oblique projection'

def table(doc,q):
    compact=q['title'].startswith('Table 6.')
    cp=doc.add_paragraph();runtext(cp,q['title'],bold=True,size=11);cp.paragraph_format.keep_with_next=True
    n=len(q['headers']); t=doc.add_table(rows=1,cols=n);t.autofit=False
    widths=[.75,.7,1.95,.8,1.25,1.4]if n==6 else [2.2,4.65]
    for col,width in zip(t.columns,widths):col.width=Inches(width)
    for cell,width in zip(t.rows[0].cells,widths):cell.width=Inches(width)
    for i,text in enumerate(q['headers']):t.rows[0].cells[i].text=text
    pr=t.rows[0]._tr.get_or_add_trPr();rep=OxmlElement('w:tblHeader');pr.append(rep)
    for row in q['rows']:
        cells=t.add_row().cells
        for i,txt in enumerate(row):
            cells[i].width=Inches(widths[i]);cells[i].text=str(txt).replace('_','_\u200b').replace('/','/\u200b')
    tblpr=t._tbl.tblPr;borders=OxmlElement('w:tblBorders')
    for side in ['top','left','bottom','right','insideH','insideV']:
        el=OxmlElement('w:'+side);el.set(qn('w:val'),'single');el.set(qn('w:sz'),'4');el.set(qn('w:color'),'D9D9D9');borders.append(el)
    tblpr.append(borders)
    for ir,row in enumerate(t.rows):
        rowpr=row._tr.get_or_add_trPr();cant=OxmlElement('w:cantSplit');rowpr.append(cant)
        for ic,cell in enumerate(row.cells):
            cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
            tcpr=cell._tc.get_or_add_tcPr();m=OxmlElement('w:tcMar')
            for side,v in [('top','50' if compact else '80'),('bottom','50' if compact else '80'),('left','90'),('right','90')]:
                z=OxmlElement('w:'+side);z.set(qn('w:w'),v);z.set(qn('w:type'),'dxa');m.append(z)
            tcpr.append(m)
            if ir==0:
                sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'E6E6E6');tcpr.append(sh)
            for p in cell.paragraphs:
                p.paragraph_format.line_spacing=1.04;p.paragraph_format.space_after=Pt(0)
                p.alignment=WD_ALIGN_PARAGRAPH.CENTER if n==6 and ic not in [0,2]else WD_ALIGN_PARAGRAPH.LEFT
                for r in p.runs:r.font.name='Times New Roman';r.font.size=Pt(10 if compact else 10.5);r.bold=ir==0;r.font.color.rgb=RGBColor(0,0,0)
    doc.add_paragraph().paragraph_format.space_after=Pt(1)

def build(source,target,supplement=False):
    blocks=json.loads((R/source).read_text());doc=Document();style(doc);i_eq=4 if supplement else 0;figs=[]
    for q in blocks:
        kind=q['type']
        if kind=='title':
            p=doc.add_paragraph(style='Title');runtext(p,q['text']);doc.core_properties.title=q['text']
        elif kind=='heading':
            if q['text']=='Figures':continue
            # Plain headings; source scientific terminology and numerical identifiers preserved.
            text=re.sub(r'[^\w\s]',' ',q['text']);text=re.sub(r'\s+',' ',text).strip()
            doc.add_heading(text,level=q['level'])
        elif kind=='paragraph':
            paragraphs=q['text'].split('\n\n')
            for text in paragraphs:
                p=doc.add_paragraph();runtext(p,text)
                if text[:3].rstrip('.').isdigit():
                    p.paragraph_format.line_spacing=1.05;p.paragraph_format.space_after=Pt(5)
                    for r in p.runs:r.font.size=Pt(10.5)
        elif kind=='equation':mathpic(doc,q,i_eq);i_eq+=1
        elif kind=='table':table(doc,q)
        elif kind=='figure':
            p=doc.add_paragraph();p.paragraph_format.page_break_before=True;p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.keep_with_next=True
            p.paragraph_format.space_after=Pt(8)
            image=HERE/'publication_figures'/(q['filename']+'.png');im=Image.open(image)
            w=6.9;h=w*im.height/im.width
            if h>6.65:w*=6.65/h
            shape=p.add_run().add_picture(str(image),width=Inches(w));shape._inline.docPr.set('descr',q['caption'])
            cp=doc.add_paragraph();runtext(cp,q['caption'],size=11);cp.paragraph_format.line_spacing=1.1;cp.paragraph_format.keep_together=True
            figs.append(q['filename'])
    doc.save(HERE/target);print('Created',target,'with',len(figs),'figures')

build('manuscript_content.json','MANUSCRIPT_COMPLETE.docx')
build('supplement_content.json','SUPPLEMENTARY_INFORMATION.docx',True)
