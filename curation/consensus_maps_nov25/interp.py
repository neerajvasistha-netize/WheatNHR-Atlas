"""Place markers that lack a physical position by interpolation on published consensus genetic maps.
Per map and chromosome, markers that do have a RefSeq v1.0 position are anchors; anchors out of collinear order are
removed (longest non-decreasing run of Mb along cM). A query marker takes the Mb range of anchors within 1 cM of it,
otherwise the two flanking anchors (each within 3 cM). Placements wider than 20 Mb are not used."""
import pandas as pd,bisect,sys
sys.path.insert(0,'.'); from nm import n2,base
M=pd.read_csv('maps_anchored.tsv',sep='\t'); M['k']=M.marker.map(n2); M['b']=M.k.map(base)
def lis(v):
    tails=[];idx=[];prev=[-1]*len(v)
    for i,x in enumerate(v):
        j=bisect.bisect_right([v[t] for t in tails],x)
        if j==len(tails): tails.append(i)
        else: tails[j]=i
        prev[i]=tails[j-1] if j>0 else -1
    out=[];i=tails[-1] if tails else -1
    while i>=0: out.append(i);i=prev[i]
    return set(out[::-1])
CLEAN={}
for (mp,ch),g in M[M.Mb.notna()].groupby(['map','chr']):
    g=g.sort_values(['cM','Mb']); v=list(g.Mb)
    if len(v)<15: continue
    if pd.Series(v).corr(pd.Series(list(g.cM)),method='spearman')<0: v=[-x for x in v]
    keep=lis(v); c=g.iloc[sorted(keep)]
    CLEAN[(mp,ch)]=(list(c.cM),list(c.Mb))
def place_one(mp,ch,cm):
    if (mp,ch) not in CLEAN: return None
    cms,mbs=CLEAN[(mp,ch)]
    w=[mbs[i] for i,c in enumerate(cms) if abs(c-cm)<=1]
    if len(w)>=2: lo,hi,how=min(w),max(w),f'{len(w)} anchor(s) within 1 cM'
    else:
        i=bisect.bisect_left(cms,cm)
        if i==0 or i==len(cms): return None
        if cm-cms[i-1]>2 or cms[i]-cm>2: return None
        lo,hi=sorted([mbs[i-1],mbs[i]]);how=f'flanking anchors at {cms[i-1]:.1f} and {cms[i]:.1f} cM'
    if hi-lo>10: return None
    return (round(lo,4),round(hi,4),f'{mp} {cm:.1f} cM, {how}')
def interp(marker,ch):
    k=n2(marker); g=M[(M.chr==ch)&(M.Mb.isna())&((M.k==k)|(M.b==base(k)))]
    P=[p for p in (place_one(r.map,ch,r.cM) for r in g.itertuples()) if p]
    if not P: return None
    P.sort(key=lambda p:p[1]-p[0]); best=P[0]
    for p in P[1:]:
        if p[0]>best[1]+20 or p[1]<best[0]-20: return ('conflict',P)
    return best
if __name__=='__main__':
    print(len(CLEAN),'map-chromosome sets'); print(interp('Xpsr546','5D'),interp('Xbcd1434','2B'),interp('Xcdo1008','2D'))
