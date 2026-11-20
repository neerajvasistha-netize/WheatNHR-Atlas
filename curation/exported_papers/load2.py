import json,re
W='/home/claude/wheatnhr-atlas/web/'
ORDER=["seed_cgsw","seed_ad57","seed_ad58","seed_ad59","seed_ad60","seed_ad61","seed_ad62","seed_ad63","seed_ad65","seed_ad66","seed_ad67","seed_ad68","seed_ad69","seed_ad70","seed_ad71","seed_ad72","seed_live","seed_ad64"]
def load():
    txt=open(W+'index.html').read()
    m=re.search(r'<script type="application/json" id="seed">(.*?)</script>',txt,re.S)
    base=json.loads(m.group(1));R={r['id']:r for r in base['reports']};S=dict(base.get('sources',{}))
    for f in ORDER:
        d=json.load(open(W+f+'.json'))
        for r in d.get('reports',[]): R[r['id']]=r
        for x in d.get('del',[]) if isinstance(d.get('del'),list) else []: R.pop(x,None)
        s=d.get('sources',{})
        if isinstance(s,dict): S.update(s)
    return R,S,base
