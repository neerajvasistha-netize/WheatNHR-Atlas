import json,re,sys,io,contextlib,gzip
sys.path[:0]=['scripts/curation','scripts/curation/sb19']
with contextlib.redirect_stdout(io.StringIO()):
    import stb44
MP=json.load(gzip.open('data/markerpos_all.json.gz'))
SSR=r'(gwm|wmc|barc|cfd|cfa|gdm|wms|psp|cfb|dupw|swes|wmm|gpw|gpw|bcd|cdo|psr|fbb|mag|ksm|glk|ksu|abg|wg|ta|tam|xpsp|edm|cnl|hbg|hbe|dp|bwm|wpt)'
def norm(t):
    t=t.strip().strip('.').strip()
    if re.match(r'^\d[ABD]:\d+$',t): return [t.split(':')[1]]
    m=re.match(r'^(?:X)?('+SSR[1:-1]+r')[-_]?(\d+[a-z]?)$',t,re.I)
    if m:
        p,n=m.group(1).lower(),m.group(2)
        if p=='wpt': return ['wPt-'+n.rjust(4,'0'),'wPt-'+n]
        return ['X'+p+n, p+n, 'X'+p+n.lstrip('0')]
    if re.match(r'^wPt\d',t): return ['wPt-'+t[3:]]
    return [t]
def tokens(s):
    s=s.replace('–',' ').replace('—',' ')
    out=[]
    for t in re.split(r'[,;/ ]+|\(|\)',s):
        if not t or not re.search(r'\d',t) or len(t)<4: continue
        if re.match(r'^(peak|flanking|closest)$',t,re.I): continue
        # SSR-SSR ranges joined by '-'
        if '-' in t and not t.startswith(('wPt-','AX-','tPt-','rPt-','GENE-')) and re.match(r'^[A-Za-z]+\d+[a-z]?-[A-Za-z]+\d+',t):
            out+=t.split('-')
        else: out.append(t)
    return out
def encoded(t):
    m=re.match(r'^S(\d[ABD])_(\d{4,})$',t)
    return (m.group(1),int(m.group(2))/1e6) if m else None
def resolve(t,ch):
    e=encoded(t)
    if e: return (e[1],'marker name encodes RefSeq v1.0 position') if e[0]==ch else (None,'')
    for v in norm(t):
        p,h=stb44.P(v,ch)
        if p: return p,h
    return None,''
