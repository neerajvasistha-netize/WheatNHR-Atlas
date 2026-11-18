"""Public Crossref metadata search for wheat biotrophic-disease resistance-locus papers (no personal data sent)."""
import json, urllib.request, urllib.parse, time, sys
Q = ["wheat stripe rust resistance QTL", "wheat yellow rust resistance gene mapping", "Puccinia striiformis wheat resistance genome-wide association",
     "wheat leaf rust resistance QTL", "wheat leaf rust resistance gene mapping", "Puccinia triticina wheat resistance genome-wide association",
     "wheat stem rust resistance QTL", "wheat stem rust resistance gene mapping Ug99", "Puccinia graminis wheat resistance genome-wide association",
     "wheat powdery mildew resistance gene mapping", "wheat powdery mildew resistance QTL", "Blumeria graminis wheat resistance genome-wide association",
     "adult plant resistance wheat rust QTL mapping", "wheat rust resistance gene cloning", "wheat powdery mildew resistance gene cloning",
     "wheat common bunt resistance gene", "Tilletia caries wheat resistance QTL", "wheat dwarf bunt resistance", "wheat loose smut resistance gene",
     "Karnal bunt resistance wheat QTL", "Aegilops rust resistance gene introgression wheat", "wild emmer powdery mildew resistance gene",
     "wheat multiple disease resistance rust mildew QTL", "durum wheat rust resistance QTL", "wheat stripe rust resistance gene fine mapping",
     "wheat leaf rust resistance gene fine mapping", "wheat stem rust resistance gene fine mapping", "wheat powdery mildew resistance gene fine mapping",
     "Yr gene wheat chromosome molecular markers", "Lr gene wheat chromosome molecular markers", "Sr gene wheat chromosome molecular markers",
     "Pm gene wheat chromosome molecular markers"]
out = {}
for q in Q:
    for page in range(2):
        u = 'https://api.crossref.org/works?' + urllib.parse.urlencode({'query.bibliographic': q, 'filter': 'type:journal-article,from-pub-date:1985', 'rows': 1000,
            'offset': page * 1000, 'select': 'DOI,title,author,issued,container-title,volume,issue,page,type'})
        d = None
        for t in range(3):
            try:
                d = json.load(urllib.request.urlopen(u, timeout=120)); break
            except Exception as e:
                print('retry', q, e, file=sys.stderr); time.sleep(5)
        if not d: continue
        for it in d['message']['items']: out[it['DOI'].lower()] = it
        print(q, page, len(d['message']['items']), len(out), file=sys.stderr)
        time.sleep(1)
json.dump(out, open('/tmp/bib/sweep.json', 'w'))
print(len(out))
