"""Propose positions for atlas records that have a chromosome but no position, from (1) markers with a known RefSeq v1.0
position, (2) in-silico PCR of Somers et al. 2004 SSR primers, (3) interpolation on consensus genetic maps (interp.py)."""
import json,glob,sys,re,os
sys.path.insert(0,'.');sys.path.insert(0,'../../p26')
from nm import n2,split_markers; from look2 import pos; import interp as I
CEN={'1A':215,'1B':237,'1D':166,'2A':326,'2B':344,'2D':264,'3A':316,'3B':346,'3D':237,'4A':264,'4B':303,'4D':185,
     '5A':252,'5B':199,'5D':188,'6A':283,'6B':323,'6D':211,'7A':360,'7B':308,'7D':336}   # approx. RefSeq v1.0 centromeres
EP={}
for f in glob.glob('relax2_*.json'):
    ch=f[7:-5]
    for m,h in json.load(open(f)).items():
        h=[x for x in h if x[3]+x[4]<=4]
        if len(h)==1: EP[(m,ch)]=round(h[0][0]/1e6,4)
def arm_of(r):
    t=' '.join([r.get('name',''),r.get('evidence',''),r.get('arm') or ''])
    m=re.findall(r'(?<![0-9A-Za-z])'+r['chr']+r'([SL])(?![a-z])',t)
    return m[0] if len(set(m))==1 else None
out=[]
for f in sorted(glob.glob('/home/claude/wheatnhr-atlas/web/seed*.json')):
    try: d=json.load(open(f))
    except: continue
    if not isinstance(d,dict): continue
    for r in d.get('reports',[]):
        ch=r.get('chr')
        if r.get('start') is not None or r.get('peak') is not None or ch not in CEN or 'reject' in (r.get('note') or '').lower(): continue
        ms=split_markers(r.get('marker'))
        if not ms: continue
        P=[]
        for m in ms:
            p=pos(m,ch)
            if p: P.append((m,p[1],p[1],p[3]));continue
            if (n2(m),ch) in EP: x=EP[(n2(m),ch)];P.append((m,x,x,'in-silico PCR of the Somers et al. 2004 primer pair'));continue
            q=I.interp(m,ch)
            if q and q[0]!='conflict': P.append((m,q[0],q[1],'interpolated on the '+q[2]))
        if not P: continue
        lo=min(p[1] for p in P);hi=max(p[2] for p in P)
        if hi-lo>30: P=[P[0]];lo,hi=P[0][1],P[0][2]
        arm=arm_of(r); c=CEN[ch]
        if arm and ((arm=='S' and lo>c+15) or (arm=='L' and hi<c-15)):
            out.append(dict(file=f.split('/')[-1],id=r['id'],name=r['name'],chr=ch,status='rejected (arm)',P=P,arm=arm));continue
        out.append(dict(file=f.split('/')[-1],id=r['id'],name=r['name'],chr=ch,status='placed',lo=lo,hi=hi,P=P,arm=arm))
json.dump(out,open('proposals.json','w'),indent=1)
import collections;print(collections.Counter(o['status'] for o in out),collections.Counter(o['file'] for o in out if o['status']=='placed'))
print(collections.Counter(p[3].split(' ')[0]+' '+p[3].split(' ')[3] if p[3].startswith('interp') else p[3][:30] for o in out for p in o['P']))
