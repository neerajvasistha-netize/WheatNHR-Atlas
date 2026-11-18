"""Europe PMC lookup (public REST, no personal data) for every DOI in the biotroph BibTeX list: PMCID, open-access flag, full text in Europe PMC."""
import json, re, urllib.request, urllib.parse, time, sys
B = open('docs/literature/WheatRLA_biotroph_resistance_loci_papers.bib').read()
dois = re.findall(r'doi = \{(.+?)\}', B)
res = {}
for i in range(0, len(dois), 20):
    chunk = dois[i:i + 20]
    q = ' OR '.join(f'DOI:"{d}"' for d in chunk)
    u = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({'query': q, 'format': 'json', 'pageSize': 100, 'resultType': 'lite'})
    for t in range(3):
        try:
            d = json.load(urllib.request.urlopen(u, timeout=60)); break
        except Exception as e:
            print('retry', i, e, file=sys.stderr); time.sleep(3); d = None
    if not d: continue
    for r in d.get('resultList', {}).get('result', []):
        if r.get('doi'): res[r['doi'].lower()] = {k: r.get(k) for k in ('pmid', 'pmcid', 'isOpenAccess', 'inEPMC', 'inPMC', 'hasPDF', 'journalTitle', 'pubYear')}
    time.sleep(0.3)
json.dump(res, open('/tmp/bib/epmc.json', 'w'))
n_oa = sum(1 for v in res.values() if v.get('pmcid') and v.get('inEPMC') == 'Y')
print(len(dois), 'dois;', len(res), 'in Europe PMC;', n_oa, 'with full text in Europe PMC;', sum(1 for v in res.values() if v.get('isOpenAccess') == 'Y'), 'flagged open access')
