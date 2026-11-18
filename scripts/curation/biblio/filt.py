import json,re,collections
S=json.load(open('/tmp/bib/sweep.json'))
W=re.compile(r'\b(wheat|triticum|aegilops|durum|spelt|emmer|einkorn|thinopyrum|secale|rye|haynaldia|dasypyrum)\b',re.I)
D={'Stripe rust':r'stripe rust|yellow rust|striiformis|\bYr\d',
   'Leaf rust':r'leaf rust|brown rust|triticina|recondita|\bLr\d',
   'Stem rust':r'stem rust|black rust|graminis f\.? ?sp\.? ?tritici|Ug99|\bSr\d',
   'Powdery mildew':r'powdery mildew|blumeria|erysiphe|\bPm\d|\bMl[A-Z]',
   'Rusts (general)':r'\brusts?\b',
   'Common bunt':r'common bunt|stinking smut|tilletia (caries|tritici|laevis|foetida)|\bBt\d',
   'Dwarf bunt':r'dwarf bunt|tilletia controversa',
   'Karnal bunt':r'karnal bunt|tilletia indica',
   'Loose smut':r'loose smut|ustilago|\bUt\d',
   'Flag smut':r'flag smut|urocystis'}
G=re.compile(r'\b(QTLs?|quantitative trait|genes?|GWAS|genome[- ]wide|association|mapping|mapped|loc(us|i)|markers?|introgress\w*|translocation|clon\w+|haplotype|allel\w+|KASP|SNPs?|linkage|chromosom\w+|R-gene|NLR|resistan\w+ genes?)\b',re.I)
R=re.compile(r'resist|immun|tolera|suscept',re.I)
X=re.compile(r'\b(fungicide|epidemiolog|yield loss|forecast|population structure of|virulence survey|race survey|pathotype|spore|urediniospore|effector repertoire|transcriptom\w+ of the pathogen|barley|oat\b|rice\b|maize|sorghum|arabidopsis|rye grass|ryegrass)\b',re.I)
keep={}
for doi,it in S.items():
    t=' '.join(it.get('title') or ['']); 
    WS=W.search(t) or re.search(r'\b(Yr|Lr|Sr|Pm)\d|triticina|f\. ?sp\. ?tritici|Ug99|Tilletia|Ustilago tritici|stripe rust|leaf rust|stem rust',t,re.I)
    if not t or not WS or not R.search(t) or not G.search(t): continue
    if re.search(r'\b(barley|hordeum|oat|rice|maize|sorghum|arabidopsis|brachypodium|poplar|coffee|soybean|sugarcane|pine|bean)\b',t,re.I) and not re.search(r'wheat|triticum',t,re.I): continue
    ds=[k for k,p in D.items() if re.search(p,t,re.I)]
    if not ds: continue
    if X.search(t) and not re.search(r'wheat',t,re.I): continue
    if re.search(r'\b(barley|oat|rice|maize)\b',t,re.I) and not re.search(r'wheat',t,re.I): continue
    it['_dis']=ds; keep[doi]=it
print(len(keep)); print(collections.Counter(d for v in keep.values() for d in v['_dis']).most_common())
json.dump(keep,open('/tmp/bib/keep.json','w'))
import random; random.seed(1)
for d in random.sample(list(keep),0): print(' -',' '.join(keep[d]['title'])[:130])
