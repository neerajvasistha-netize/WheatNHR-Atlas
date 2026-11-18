"""Targeted public Crossref lookups for original studies cited by Tong et al. (2024) not found by the sweep."""
import json, urllib.request, urllib.parse, time, sys
E = json.load(open('/tmp/bib/excl.json'))
S = json.load(open('/tmp/bib/sweep.json'))
extra = ["wheat rust resistance genome-wide association study", "wheat yellow rust adult plant resistance QTL", "wheat leaf rust adult plant resistance QTL",
         "wheat stem rust adult plant resistance QTL", "wheat seedling rust resistance genome-wide association", "landraces rust resistance wheat association mapping",
         "wheat powdery mildew adult plant resistance QTL", "nested association mapping wheat rust", "wheat rust resistance haplotype", "wheat stripe rust resistance loci validation KASP"]
n = 0
for a, y in E['missing_tong']:
    a = a.strip(',')
    u = 'https://api.crossref.org/works?' + urllib.parse.urlencode({'query.author': a, 'query.bibliographic': 'wheat rust resistance', 'rows': 8,
        'filter': f'from-pub-date:{y},until-pub-date:{y}', 'select': 'DOI,title,author,issued,container-title,volume,issue,page,type'})
    try:
        d = json.load(urllib.request.urlopen(u, timeout=60))
        for it in d['message']['items']: S[it['DOI'].lower()] = it; n += 1
    except Exception as e: print('err', a, y, e, file=sys.stderr)
    time.sleep(0.5)
for q in extra:
    for page in range(2):
        u = 'https://api.crossref.org/works?' + urllib.parse.urlencode({'query.bibliographic': q, 'filter': 'type:journal-article,from-pub-date:1985', 'rows': 1000,
            'offset': page * 1000, 'select': 'DOI,title,author,issued,container-title,volume,issue,page,type'})
        try:
            d = json.load(urllib.request.urlopen(u, timeout=120))
            for it in d['message']['items']: S[it['DOI'].lower()] = it
        except Exception as e: print('err', q, e, file=sys.stderr)
        time.sleep(1)
json.dump(S, open('/tmp/bib/sweep.json', 'w'))
print(len(S), n)
