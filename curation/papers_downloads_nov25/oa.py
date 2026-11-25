import json,subprocess,sys,urllib.parse
M=json.load(open(sys.argv[1])); out={}
try: out=json.load(open(sys.argv[2]))
except: pass
for m in M:
  d=m.get('doi')
  if not d or d in out: continue
  r=subprocess.run(['curl','-s','-m','20','https://api.openalex.org/works/doi:'+urllib.parse.quote(d)+'?select=title,authorships,biblio,publication_year,primary_location'],capture_output=True,text=True).stdout
  try:
    j=json.loads(r); out[d]={'title':j['title'],'authors':[a['author']['display_name'] for a in j['authorships']],'biblio':j['biblio'],'year':j['publication_year'],'journal':(j.get('primary_location') or {}).get('source',{}) and j['primary_location']['source'].get('display_name')}
  except Exception as e: out[d]={'err':r[:100]}
json.dump(out,open(sys.argv[2],'w'))
print(len(out),sum('err' in v for v in out.values()))
