"""Turn every supplementary file into table rows (list of cell strings), saved as rows/<key>.json"""
import os,json,glob,subprocess,zipfile,re,sys
os.makedirs('rows',exist_ok=True)
def xlsx(f):
    import openpyxl;R=[]
    wb=openpyxl.load_workbook(f,read_only=True,data_only=True)
    for ws in wb.worksheets:
        for r in ws.iter_rows(values_only=True):
            c=['' if v is None else str(v) for v in r]
            if any(c): R.append([ws.title]+c)
    return R
def xls(f):
    import xlrd;R=[];b=xlrd.open_workbook(f)
    for s in b.sheets():
        for i in range(s.nrows): R.append([s.name]+[str(v) for v in s.row_values(i)])
    return R
def docx_(f):
    import docx;R=[];d=docx.Document(f)
    for t in d.tables:
        for r in t.rows: R.append(['table']+[c.text.strip() for c in r.cells])
    for p in d.paragraphs:
        if p.text.strip(): R.append(['text',p.text.strip()])
    return R
def pdf(f):
    t=subprocess.run(['pdftotext','-layout',f,'-'],capture_output=True,text=True).stdout
    return [['pdf']+re.split(r'\s{2,}',l.strip()) for l in t.splitlines() if l.strip()]
def txt(f):
    R=[]
    for l in open(f,errors='ignore'):
        R.append(['txt']+re.split(r'\t|,(?=\S)',l.rstrip('\n')))
    return R
def conv(f,tmp='/tmp/sofconv'):
    os.makedirs(tmp,exist_ok=True)
    ext='xlsx' if f.lower().endswith('.xls') else 'docx'
    subprocess.run(['soffice','--headless','--convert-to',ext,'--outdir',tmp,f],capture_output=True,timeout=120)
    o=os.path.join(tmp,os.path.splitext(os.path.basename(f))[0]+'.'+ext);return o if os.path.exists(o) else None
def read(f):
    e=f.lower().rsplit('.',1)[-1]
    try:
        if e in('xlsx','xlsm'): return xlsx(f)
        if e=='xls':
            try: return xls(f)
            except Exception: o=conv(f);return xlsx(o) if o else []
        if e=='docx': return docx_(f)
        if e=='doc': o=conv(f);return docx_(o) if o else []
        if e=='pdf': return pdf(f)
        if e in('csv','txt','tsv'): return txt(f)
        if e=='zip':
            R=[];d=f+'_x';os.makedirs(d,exist_ok=True);zipfile.ZipFile(f).extractall(d)
            for g in glob.glob(d+'/**/*',recursive=True):
                if os.path.isfile(g): R+=read(g)
            return R
    except Exception as ex: print('fail',f,ex,file=sys.stderr)
    return []
for d in sorted(glob.glob('files/*')):
    k=os.path.basename(d);R=[]
    for f in sorted(glob.glob(d+'/*')):
        if os.path.isfile(f):
            rr=read(f);R+=[[os.path.basename(f)]+r for r in rr]
    json.dump(R,open('rows/'+k+'.json','w'))
    print(k,len(R))
