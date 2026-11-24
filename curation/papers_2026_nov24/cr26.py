import json,urllib.request,urllib.parse,re,time,sys
ex=set(json.load(open('../ms/existing_dois.json')))
qs=["wheat resistance QTL","wheat genome-wide association resistance","wheat rust resistance gene mapping","wheat stripe rust","wheat leaf rust","wheat stem rust","wheat powdery mildew resistance","wheat Fusarium head blight","wheat septoria","wheat tan spot","wheat spot blotch","wheat blast resistance","wheat virus resistance","wheat nematode resistance","wheat bunt smut resistance","wheat aphid resistance","Triticum aestivum disease resistance locus","wheat fine mapping resistance gene","wheat crown rot","wheat sharp eyespot","wheat take-all","durum wheat resistance QTL","wheat Hessian fly","wheat bacterial leaf streak","wheat stagonospora nodorum blotch","wheat yellow mosaic virus","wheat eyespot"]
import os
out=json.load(open('cr_raw.json')) if os.path.exists('cr_raw.json') else {}
for q in qs:
  for off in (0,100):
    u='https://api.crossref.org/works?'+urllib.parse.urlencode({'query.bibliographic':q,'filter':'from-pub-date:2026-01-01,type:journal-article','rows':100,'offset':off,'select':'DOI,title,container-title,issued,author','mailto':'neeraj.vasistha@rgu.ac.in'})
    d=None
    for a in range(5):
      try: d=json.load(urllib.request.urlopen(u,timeout=60)); break
      except Exception as e: print('err',q,e,file=sys.stderr); time.sleep(5*(a+1))
    if d is None: continue
    for it in d['message']['items']:
      doi=it['DOI'].lower(); t=(it.get('title') or [''])[0]
      out[doi]=dict(doi=doi,title=t,journal=(it.get('container-title') or [''])[0],year=it['issued']['date-parts'][0][0],first=(it.get('author') or [{}])[0].get('family',''))
    time.sleep(1.5)
json.dump(out,open('cr_raw.json','w'),indent=0)
kw=re.compile(r'wheat|triticum|aegilops',re.I)
dis=re.compile(r'resist|rust|mildew|blight|blotch|spot|blast|virus|mosaic|nematode|bunt|smut|aphid|hessian|rot|eyespot|take-all|streak|pathogen|disease|susceptib|immune|NLR',re.I)
loc=re.compile(r'QTL|GWAS|genome-wide|association|mapping|locus|loci|gene|map|clon|introgress|haplotype|QTN|marker',re.I)
cand=[v for k,v in out.items() if k not in ex and kw.search(v['title']) and dis.search(v['title']) and loc.search(v['title'])]
json.dump(cand,open('cand.json','w'),indent=1)
print(len(out),len(cand))
