"""ad37b: erratum (Euphytica 200:27) corrects transposed author names of Szabo-Hever et al. 2014; repoint V1478-V1483 to a corrected source key."""
import json, re
K = 'SzaboHever2014Euph'
SRC = {'citation': 'Szabo-Hever A, Lehoczki-Krsjak S, Varga M, Purnhauser L, Pauk J, Lantos C, Mesterhazy A (2014) Differential influence of QTL linked to Fusarium head blight, Fusarium-damaged kernel, deoxynivalenol contents and associated morphological traits in a Frontana-derived wheat population. Euphytica 200:9-26 (author names corrected in erratum, Euphytica 200:27)', 'url': 'https://doi.org/10.1007/s10681-014-1124-2', 'short': 'Szabo-Hever et al. (2014)'}
IDS = ['V1478', 'V1479', 'V1480', 'V1481', 'V1482', 'V1483']
def fix(D):
    D.setdefault('sources', {})[K] = SRC
    for r in D['reports']:
        if r['id'] in IDS: r['source'] = K
M = json.load(open('data/merged_v13f.json')); fix(M); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); fix(Sd)
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
ops = [{'op': 'set', 'path': 'sources/' + K, 'data': SRC}] + [{'op': 'merge', 'path': 'reports/' + i, 'data': {'source': K}} for i in IDS]
json.dump(ops, open('/tmp/ops37b.json', 'w'))
