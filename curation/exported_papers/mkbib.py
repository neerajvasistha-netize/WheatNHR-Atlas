import json,re,collections,unicodedata,csv
from supp import supp,SUP
need=json.load(open('need_fig.json'));placed={(r['h'],r['name']) for r in json.load(open('../fig/fig_res.json'))}
P={x['h']:x for x in json.load(open('s1.json'))};meta=json.load(open('rdfmeta.json'))
L=collections.defaultdict(set)
for h,n,ch in need:
    n=n.replace('‐','-').replace('–','-')
    if (h,n) not in placed: L[h].add(n)
def esc(s): return s.replace('&','\\&').replace('%','\\%').replace('_','\\_').replace('#','\\#').replace('{','').replace('}','')
def desc(h,lab):
    t=re.sub(r'\s+',' ',open('txr/'+h+'.txt',errors='ignore').read())
    for m in re.finditer(re.escape(lab)+r'\s*[:|]\s*([A-Z][^.]{10,200})',t):
        d=m.group(1);n=SUP.search(d)
        d=(d[:n.start()] if n else d).strip(' ;,(')
        if len(d)>=10: return d[:110]
    return ''
keys=set();out=[];rows=[];stat=collections.Counter()
hdr=open('bibhdr.txt').read()
for h in sorted(L,key=lambda h:(meta[P[h]['path']]['au'][0][0] if meta[P[h]['path']]['au'] else '',meta[P[h]['path']]['yr'])):
    m=meta[P[h]['path']];au=m['au'];yr=m['yr'] or P[h]['yr']
    base=unicodedata.normalize('NFKD',au[0][0] if au else 'anon').encode('ascii','ignore').decode()
    k=re.sub(r'[^A-Za-z]','',base)+yr+re.sub(r'[^a-z]','',(m['ti'].split()[0] if m['ti'] else 'x').lower())[:8]
    while k in keys: k+='b'
    keys.add(k);loci=sorted(L[h])
    sc,ctx=supp(h,loci)
    tabs=[(lab,s) for lab,s in sc.most_common() if not re.search(r'Fig',lab)]
    top=[lab for lab,s in tabs if s>=max(8,0.35*tabs[0][1])][:3] if tabs else []
    req=top or ([lab for lab,s in sc.most_common() if re.search(r'Fig',lab)][:1])
    if req:
        reqs='Required: '+'; '.join(lab+(' ('+desc(h,lab)+')' if desc(h,lab) else '') for lab in req)+'.';stat['named']+=1
    else:
        reqs='Required: the paper names no supplementary file (it may have none); check the journal page and download any supplementary table with QTL, marker-trait associations, linkage map or marker positions.';stat['none']+=1
    f=[f"  author = {{{' and '.join(esc(s+', '+g) for s,g in au)}}}",f"  title = {{{{{esc(m['ti'])}}}}}",f"  journal = {{{esc(m['j'])}}}",f"  year = {{{yr}}}"]
    if m['vol']: f.append(f"  volume = {{{m['vol']}}}")
    if m['pg']: f.append(f"  pages = {{{m['pg']}}}")
    if m['doi']: f+= [f"  doi = {{{m['doi']}}}",f"  url = {{https://doi.org/{m['doi']}}}"]
    f.append("  keywords = {WheatTraitMap: supplementary needed, "+P[h]['dis'][0][0]+"}")
    f.append("  note = {"+esc(reqs)+" Direction: open the DOI link, go to the Supplementary material section, download the files above and attach them to this item. Loci still without position: "+esc(', '.join(loci))+".}")
    out.append('@article{'+k+',\n'+',\n'.join(f)+'\n}\n')
    rows.append([k,yr,m['ti'],'https://doi.org/'+m['doi'],reqs.replace('Required: ',''),len(loci),', '.join(loci)])
open('supplementary_needed.bib','w').write(hdr+'\n'+'\n'.join(out))
w=csv.writer(open('supplementary_needed.tsv','w'),delimiter='\t');w.writerow(['key','year','title','doi_link','required_supplementary','n_loci','loci']);w.writerows(rows)
print(len(out),stat)
for r in rows[:16]: print(r[0],'|',r[4][:220])
