import json,re,urllib.request,urllib.parse,unicodedata,sys
from concurrent.futures import ThreadPoolExecutor
def norm(s): return re.sub(r'[^a-z]','',unicodedata.normalize('NFKD',s).encode('ascii','ignore').decode().lower())
R=json.load(open('fhb.json'))['refs']
def one(k):
    t=R[k];y=re.search(r'\((\d{4})\)',t);y=y.group(1) if y else None
    fa=norm(re.split(r',',t)[0].split('.')[-1])
    url='https://api.crossref.org/works?rows=5&query.bibliographic='+urllib.parse.quote(t[:300])
    for _ in range(2):
        try: items=json.load(urllib.request.urlopen(url,timeout=40))['message']['items'];break
        except Exception: items=[]
    for i in items:
        a=i.get('author',[]);iy=i.get('issued',{}).get('date-parts',[[None]])[0][0]
        if a and norm(a[0].get('family',''))==fa and y and iy and abs(int(y)-iy)<=1:
            au=', '.join(x.get('family','')+' '+''.join(z[0] for z in x.get('given','').replace('-',' ').replace('.',' ').split()) for x in a[:6])+(' et al.' if len(a)>6 else '')
            return k,[i['DOI'].lower(),au,iy,re.sub(r'<[^>]+>|\s+',' ',i['title'][0]).strip(),i.get('container-title',[''])[0].replace('®',''),i.get('volume'),i.get('page') or i.get('article-number')]
    return k,None
with ThreadPoolExecutor(6) as ex: out=dict(ex.map(one,list(R)))
json.dump(out,open('fhbcr.json','w'),indent=0)
print(sum(1 for v in out.values() if v),len(out))
for k,v in out.items():
    if not v: print('MISS',k,R[k][:120])
