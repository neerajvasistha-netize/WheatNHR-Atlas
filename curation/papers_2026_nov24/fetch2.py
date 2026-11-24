import json,re,os,urllib.request,subprocess,time
sel=json.load(open('sel.json'))
UA={'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36'}
def get(u,t=60,raw=False):
    r=urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=t); b=r.read(); return b if raw else b.decode('utf8','ignore')
meta={}
for v in sel:
    doi=v['doi']; fn='txt/'+re.sub(r'[^a-z0-9]','_',doi)+'.txt'
    cur=open(fn).read() if os.path.exists(fn) else ''
    try: w=json.loads(get('https://api.openalex.org/works/doi:'+doi+'?mailto=neeraj.vasistha@rgu.ac.in'))
    except Exception as e: print('oa err',doi,e); continue
    ab=''
    if w.get('abstract_inverted_index'):
        pos={}
        for word,ps in w['abstract_inverted_index'].items():
            for p in ps: pos[p]=word
        ab=' '.join(pos[i] for i in sorted(pos))
    au=[a['author']['display_name'] for a in w.get('authorships',[])]
    meta[doi]=dict(title=w.get('title'),authors=au,year=w.get('publication_year'),date=w.get('publication_date'),journal=((w.get('primary_location') or {}).get('source') or {}).get('display_name'),biblio=w.get('biblio'),abstract=ab,oa=[l.get('pdf_url') for l in w.get('locations',[]) if l.get('pdf_url')])
    if len(cur)<20000:
        got=''
        for u in meta[doi]['oa'][:3]:
            try:
                b=get(u,90,True)
                if b[:4]==b'%PDF':
                    p='pdf/'+re.sub(r'[^a-z0-9]','_',doi)+'.pdf'; open(p,'wb').write(b)
                    got=subprocess.run(['pdftotext',p,'-'],capture_output=True,text=True).stdout; break
            except Exception as e: pass
        if len(got)>len(cur): cur=got
        if ab and ab[:200] not in cur: cur='ABSTRACT: '+ab+'\n\n'+cur
        open(fn,'w').write(cur)
    print(doi,len(cur),len(meta[doi]['oa']),flush=True); time.sleep(0.3)
json.dump(meta,open('meta.json','w'),indent=0)
