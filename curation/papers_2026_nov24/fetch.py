import json,re,os,urllib.request,urllib.parse,time,html as H,sys
sel=json.load(open('sel.json'))
UA={'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36','Accept':'text/html,application/xhtml+xml'}
def get(u,t=60):
    r=urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=t);return r.read().decode('utf8','ignore'),r.geturl()
def strip(s):
    s=re.sub(r'(?s)<(script|style)[^>]*>.*?</\1>',' ',s); s=re.sub(r'<[^>]+>',' ',s); s=H.unescape(s); return re.sub(r'[ \t\r\f\v]+',' ',re.sub(r'\n\s*\n+','\n',s))
log={}
for v in sel:
    doi=v['doi']; fn='txt/'+re.sub(r'[^a-z0-9]','_',doi)+'.txt'
    if os.path.exists(fn) and os.path.getsize(fn)>20000: continue
    how=[];text=''
    try:
        j=json.loads(get('https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=core&query=DOI:%22'+urllib.parse.quote(doi)+'%22')[0])
        res=j['resultList']['result']
        if res:
            r=res[0]; text=strip(r.get('abstractText','')); how.append('epmc abs')
            if r.get('pmcid'):
                try: x,_=get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+r['pmcid']+'/fullTextXML',90); text=strip(x); how.append('epmc full')
                except Exception as e: how.append('nofull')
    except Exception as e: how.append('epmc err')
    if len(text)<20000:
        try:
            x,url=get('https://doi.org/'+doi); t=strip(x)
            if len(t)>len(text): text=t; how.append('pub '+urllib.parse.urlparse(url).netloc)
        except Exception as e: how.append('pub err '+str(e)[:40])
    open(fn,'w').write(text); log[doi]=(len(text),how); print(doi,len(text),how,flush=True)
    time.sleep(1)
json.dump(log,open('fetchlog.json','w'),indent=0)
