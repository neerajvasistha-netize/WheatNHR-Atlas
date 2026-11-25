"""Write accepted placements into the records (seed files) and prepare the pangenome conversion input."""
import json,collections
P=[o for o in json.load(open('proposals.json')) if o['status']=='placed' and o['id']!='S19-081']
W='/home/claude/wheatnhr-atlas/web/'
byf=collections.defaultdict(dict)
for o in P: byf[o['file']][o['id']]=o
need=[]
for f,ids in byf.items():
    d=json.load(open(W+f))
    for r in d['reports']:
        o=ids.get(r['id'])
        if not o: continue
        lo,hi=o['lo'],o['hi']
        if hi-lo<0.05: r['start']=r['end']=None; r['peak']=round((lo+hi)/2,4)
        else: r['start'],r['end'],r['peak']=round(lo,4),round(hi,4),None
        r['ref']='RefSeq v1.0'; r['grade']='C'
        r['curation']='auto (in-silico from markers, consensus maps Nov 2026)'
        how='; '.join(f"{m} {a:.2f}" + (f"-{b:.2f}" if b-a>=0.05 else '') + f" Mb ({s})" for m,a,b,s in o['P'])
        note=(r.get('note') or '').replace('No marker of this locus could be placed on '+r['chr']+'; left unplaced.','').replace('Peak marker could not be placed physically.','').strip()
        r['note']=(note+f" Placed in silico on RefSeq v1.0 on 25 Nov 2026 from the markers named by the source: {how}. Interpolated positions come from published consensus genetic maps anchored by markers of known physical position; check by hand.").strip()
        a=dict(r.get('asm') or {}); a.pop('_note',None)
        for k in list(a):
            v=a[k]
            if isinstance(v,dict) and v.get('start') is None and v.get('peak') is None: a.pop(k)
        r['asm']=a
        need.append({'rec':r,'file':f})
json.dump(need,open('/home/claude/wheatnhr-atlas/curation/papers_downloads_nov25/need96.json','w'))
print(len(need),{f:len(v) for f,v in byf.items()})
