"""Download Europe PMC full-text JATS XML for the open-access papers found by oa_check.py."""
import json, os, urllib.request, time, sys
E = json.load(open('/tmp/bib/epmc.json'))
os.makedirs('/tmp/bib/xml', exist_ok=True)
from concurrent.futures import ThreadPoolExecutor
todo = [v['pmcid'] for v in E.values() if v.get('pmcid') and v.get('inEPMC') == 'Y']
def get(p):
    f = f'/tmp/bib/xml/{p}.xml'
    if os.path.exists(f) and os.path.getsize(f) > 1000: return 0
    for t in range(2):
        try:
            x = urllib.request.urlopen(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{p}/fullTextXML', timeout=90).read()
            open(f, 'wb').write(x); return 1
        except Exception as e: err = e; time.sleep(2)
    print('err', p, err, file=sys.stderr, flush=True); return 0
with ThreadPoolExecutor(4) as ex: n = sum(ex.map(get, todo))
print(n, 'downloaded', len(os.listdir('/tmp/bib/xml')), 'files of', len(todo))
