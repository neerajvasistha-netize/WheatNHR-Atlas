import json,sys,re,collections
sys.path.insert(0,'.');sys.path.insert(0,'../../p26')
from nm import n2,split_markers; from look2 import pos; import interp as I
import place_records as PR   # reuse CEN, EP, arm_of
d=json.load(open('../all_now.json'))
out=[]
for r in d['reports']:
    ch=r.get('chr')
    if r.get('start') is not None or r.get('peak') is not None or ch not in PR.CEN or 'reject' in (r.get('note') or '').lower(): continue
    if r['id']=='S19-081': continue
    ms=split_markers(r.get('marker'))
    if not ms: continue
    P=[]
    for m in ms:
        p=pos(m,ch)
        if p: P.append((m,p[1],p[1],p[3]));continue
        if (n2(m),ch) in PR.EP: x=PR.EP[(n2(m),ch)];P.append((m,x,x,'in-silico PCR of the Somers et al. 2004 primer pair'));continue
        q=I.interp(m,ch)
        if q and q[0]!='conflict': P.append((m,q[0],q[1],'interpolated on the '+q[2]))
    if not P: continue
    lo=min(p[1] for p in P);hi=max(p[2] for p in P)
    if hi-lo>30: P=[P[0]];lo,hi=P[0][1],P[0][2]
    arm=PR.arm_of(r); c=PR.CEN[ch]
    st='placed'
    if arm and ((arm=='S' and lo>c+15) or (arm=='L' and hi<c-15)): st='rejected (arm)'
    out.append(dict(id=r['id'],name=r['name'],chr=ch,status=st,lo=lo,hi=hi,P=P,arm=arm))
json.dump(out,open('proposals2.json','w'),indent=1)
print(collections.Counter(o['status'] for o in out))
