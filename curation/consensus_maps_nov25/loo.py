import sys; sys.path.insert(0,".")
import interp as I,random,bisect,statistics
random.seed(1); A=I.M[I.M.Mb.notna()].sample(3000,random_state=1)
err=[];inside=0;n=0
for r in A.itertuples():
    key=(r.map,r.chr)
    if key not in I.CLEAN: continue
    cms,mbs=I.CLEAN[key]; pairs=[(c,m) for c,m in zip(cms,mbs) if not (abs(c-r.cM)<1e-9 and abs(m-r.Mb)<1e-9)]
    I.CLEAN[key]=([p[0] for p in pairs],[p[1] for p in pairs]); p=I.place_one(r.map,r.chr,r.cM); I.CLEAN[key]=(cms,mbs)
    if not p: continue
    n+=1; mid=(p[0]+p[1])/2; err.append(abs(mid-r.Mb)); inside+= p[0]-2<=r.Mb<=p[1]+2
print(n,'placed of',len(A),'median err',round(statistics.median(err),2),'Mb; within interval±2Mb',round(inside/n*100,1),'%; err>20Mb',sum(e>20 for e in err))
# test only collinear anchors
err=[];n=0
for r in A.itertuples():
    key=(r.map,r.chr)
    if key not in I.CLEAN: continue
    cms,mbs=I.CLEAN[key]
    if not any(abs(c-r.cM)<1e-9 and abs(m-r.Mb)<1e-9 for c,m in zip(cms,mbs)): continue
    pairs=[(c,m) for c,m in zip(cms,mbs) if not (abs(c-r.cM)<1e-9 and abs(m-r.Mb)<1e-9)]
    I.CLEAN[key]=([p[0] for p in pairs],[p[1] for p in pairs]); p=I.place_one(r.map,r.chr,r.cM); I.CLEAN[key]=(cms,mbs)
    if not p: continue
    n+=1; err.append(abs((p[0]+p[1])/2-r.Mb))
print('collinear test',n,'median',round(statistics.median(err),2),'>20Mb',sum(e>20 for e in err),'>10Mb',sum(e>10 for e in err))
