"""Minimal .xlsx reader (stdlib only): prints every sheet as rows of cell values."""
import sys,zipfile,re,xml.etree.ElementTree as ET
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def read(path):
    z=zipfile.ZipFile(path)
    ss=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',NS):
            ss.append(''.join(t.text or '' for t in si.iter('{%s}t'%NS['m'])))
    wb=ET.fromstring(z.read('xl/workbook.xml')); rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    rmap={r.get('Id'):r.get('Target') for r in rels}
    out={}
    for s in wb.find('m:sheets',NS):
        tgt=rmap[s.get('{%s}id'%NS['r'])]; tgt=tgt.lstrip('/'); tgt=tgt if tgt.startswith('xl/') else 'xl/'+tgt
        rows=[]
        for row in ET.fromstring(z.read(tgt)).iter('{%s}row'%NS['m']):
            vals={}
            for c in row.findall('m:c',NS):
                col=re.match(r'[A-Z]+',c.get('r')).group(0); n=0
                for ch in col: n=n*26+ord(ch)-64
                v=c.find('m:v',NS); t=c.get('t')
                if t=='s' and v is not None: val=ss[int(v.text)]
                elif t=='inlineStr': val=''.join(x.text or '' for x in c.iter('{%s}t'%NS['m']))
                else: val=v.text if v is not None else ''
                vals[n-1]=val
            rows.append([vals.get(i,'') for i in range(max(vals)+1)] if vals else [])
        out[s.get('name')]=rows
    return out
if __name__=='__main__':
    for name,rows in read(sys.argv[1]).items():
        print(f"## sheet '{name}': {len(rows)} rows")
        for r in rows[:int(sys.argv[2]) if len(sys.argv)>2 else 12]: print('   ',r[:14])
