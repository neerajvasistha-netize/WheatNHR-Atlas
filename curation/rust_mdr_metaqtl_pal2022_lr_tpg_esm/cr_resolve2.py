import json,sys,re,urllib.request,urllib.parse,unicodedata
def norm(s): return unicodedata.normalize('NFKD',s).encode('ascii','ignore').decode().lower()
def q(author,year,extra):
    url='https://api.crossref.org/works?rows=12&filter=from-pub-date:%d,until-pub-date:%d&query.author=%s&query.bibliographic=%s'%(int(year)-1,int(year),urllib.parse.quote(author),urllib.parse.quote(extra))
    try: items=json.load(urllib.request.urlopen(url,timeout=40))['message']['items']
    except Exception as e: return None
    best=None
    for i in items:
        a=i.get('author',[])
        if not a or norm(a[0].get('family','')).replace(' ','').replace('-','')!=norm(author).replace(' ','').replace('-',''): continue
        t=re.sub(r'<[^>]+>','',' '.join(i.get('title',['']))).lower()
        if not re.search(r'rust|striiformis|triticina|graminis',t): continue
        if i.get('type') not in('journal-article',): continue
        best=i;break
    if not best: return None
    i=best;a=i['author']
    return [i['DOI'].lower(),', '.join(x.get('family','')+' '+''.join(y[0] for y in x.get('given','').replace('-',' ').replace('.',' ').split()) for x in a[:6])+(' et al.' if len(a)>6 else ''),i['issued']['date-parts'][0][0],re.sub(r'<[^>]+>|\s+',' ',i['title'][0]).strip(),i.get('container-title',[''])[0].replace('®',''),i.get('volume'),i.get('page') or i.get('article-number')]
if __name__=='__main__':
    jobs=json.load(open(sys.argv[1]));out={}
    for key,(a,y) in jobs.items():
        out[key]=q(a,y,'rust resistance wheat');print(key,out[key] and out[key][:4],flush=True)
    json.dump(out,open(sys.argv[2],'w'),indent=0)
