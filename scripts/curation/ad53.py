"""ad53 (Nov 2026): Roy et al. 2021 FPS Table 3 LOD (mean population) added to the stat field of S19-074 to S19-084 (records otherwise unchanged)."""
import json, re
L = {'S19-074': '2.86', 'S19-075': '6.40', 'S19-076': '3.46', 'S19-077': '2.59', 'S19-078': '3.11', 'S19-079': '7.72', 'S19-080': '4.51', 'S19-081': '2.74', 'S19-082': '3.10', 'S19-083': '3.74', 'S19-084': '2.82'}
def build(rs):
    out = []
    for r in rs:
        if r['id'] in L: r = dict(r); r['stat'] = 'LOD ' + L[r['id']] + ' (mean)'
        out.append(r)
    return out
M = json.load(open('data/merged_v13f.json')); M['reports'] = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'] = build(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
print(json.dumps([{'op': 'merge', 'path': 'reports/' + k, 'ifVersion': 1, 'data': {'stat': 'LOD ' + v + ' (mean)'}} for k, v in L.items()]))
