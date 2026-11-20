import csv,json,re,os,io,zipfile,urllib.request,urllib.parse,time,sys
rows=list(csv.DictReader(open('../exp/supplementary_needed.tsv'),delimiter='\t'))
UA={'User-Agent':'Mozilla/5.0 (WheatTraitMap curation)'}
def get(u,t=60):
    r=urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=t);return r.read(),r.geturl(),r.headers.get('Content-Type','')
log={}
for r in rows:
    k=r['key'];doi=r['doi_link'].replace('https://doi.org/','');d='files/'+k;os.makedirs(d,exist_ok=True)
    got=[];why=[]
    try:
        j=json.loads(get('https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=lite&query=DOI:%22'+urllib.parse.quote(doi)+'%22')[0])
        pm=[x.get('pmcid') for x in j['resultList']['result'] if x.get('pmcid')]
        if pm:
            try:
                b,_,ct=get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+pm[0]+'/supplementaryFiles',120)
                z=zipfile.ZipFile(io.BytesIO(b))
                for n in z.namelist():
                    if n.endswith('/') : continue
                    fn=os.path.basename(n)
                    if re.search(r'\.(jpg|gif|png|tif)$',fn,re.I): continue
                    open(d+'/'+fn,'wb').write(z.read(n));got.append(fn)
                why.append('europepmc '+pm[0])
            except Exception as e: why.append('epmc supp fail '+str(e)[:60])
    except Exception as e: why.append('epmc search fail '+str(e)[:60])
    if not got:
        try:
            b,url,ct=get('https://doi.org/'+doi)
            html=b.decode('utf8','ignore')
            links=set(re.findall(r'https?://static-content\.springer\.com/esm/[^"\'<> ]+',html))
            links|=set(urllib.parse.urljoin(url,x) for x in re.findall(r'href="([^"]+\.(?:xlsx?|docx?|pdf|zip|csv|txt))"',html,re.I) if re.search(r'supp|esm|SM\d|DataSheet|Table|media',x,re.I))
            links|=set(re.findall(r'https?://www\.frontiersin\.org/api/v\d/articles/\d+/file/[^"\'<> ]+',html))
            for u in sorted(links)[:15]:
                try:
                    b2,u2,ct2=get(u,120);fn=re.sub(r'[^A-Za-z0-9._-]','_',os.path.basename(urllib.parse.urlparse(u2).path))[:80] or 'file'
                    if 'html' in ct2: continue
                    open(d+'/'+fn,'wb').write(b2);got.append(fn)
                except Exception as e: why.append('dl fail '+str(e)[:40])
            why.append('publisher page '+urllib.parse.urlparse(url).netloc)
        except Exception as e: why.append('page fail '+str(e)[:60])
    log[k]={'files':got,'how':why}
    print(k,len(got),why[-1][:60],flush=True)
json.dump(log,open('dl_log.json','w'),indent=1)
