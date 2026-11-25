import json,re,sys,os
M=json.load(open('papers_meta.json'))
todo=[m for m in M if not m.get('skip')]
a,b=int(sys.argv[1]),int(sys.argv[2]); L=int(sys.argv[3]) if len(sys.argv)>3 else 650
for m in todo[a:b]:
  p='txt/'+m['file'][:-4]+'.txt'
  t=open(p,errors='ignore').read() if os.path.exists(p) else ''
  t=re.sub(r'\s+',' ',t)
  k=re.search(r'(?i)\b(a ?b ?s ?t ?r ?a ?c ?t|summary)\b',t)
  s=t[k.end():k.end()+3000] if k else t[:3000]
  sents=re.split(r'(?<=[.;])\s+',s)
  keep=[x for x in sents if re.search(r'(?i)chromosom|\b[1-7][ABD][SL]?\b|arm|cM|located|linked|flank|QTL|designated|mapped',x)]
  o=' '.join(keep)[:L]
  print(f"## {m['i']} {m['au'].split(',')[0]} {m['year']} | {m['title'][:90]} | {m['disease']} | {len(t)}\n{o}\n")
