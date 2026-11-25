import sys,bisect;sys.path.insert(0,'.')
import interp as I
from nm import n2,base
def near(mp,ch,cm,win=3):
    if (mp,ch) not in I.CLEAN: return None
    cms,mbs=I.CLEAN[(mp,ch)]
    i=bisect.bisect_left(cms,cm)
    L=[(cms[j],mbs[j]) for j in range(max(0,i-3),min(len(cms),i+3))]
    w=[x for x in L if abs(x[0]-cm)<=win]
    if not w: return None
    return (round(min(x[1] for x in w),2),round(max(x[1] for x in w),2),[(round(a,1),round(b,1)) for a,b in w])
def R(m,ch):
    k=n2(m);g=I.M[(I.M.chr==ch)&((I.M.k==k)|(I.M.b==base(k))|(I.M.k==base(k)))]
    return [(r.map,round(r.cM,1),near(r.map,ch,r.cM)) for r in g.itertuples()]
for m,ch in [('Xglk301','1B'),('Xglk317','1B'),('Xglk483','1B'),('Xpsr131','2A'),('Xglk687','2A'),('Xglk687','2B'),('Xpsr961','2B'),('Xglk407','2B'),('Xpsr924','2B'),('Xpsr644','2B'),('Xpsr956','2B'),('Xpsr932','2D'),('Xpsr331','2D'),('Xpsr304','3A'),('Xpsr598','3A'),('Xglk118','3A'),('Xglk577','3A'),('Xglk538','3B'),('Xpsr902','3B'),('Xglk331','4A'),('Xgwm710','4A'),('Xglk348','4B'),('Xpsr921','4B'),('Xglk302','4D'),('Xpsr1101','4D'),('Xpsr1194','5A'),('Xpsr918','5A'),('Xpsr128','5B'),('Xglk614','5B'),('Xglk163','5B'),('Xpsr426','5B'),('Xglk165','7A'),('Xpsr952','7B'),('Xmwg710','7B'),('Xglk576','7B')]:
    print(m,ch,R(m,ch))
