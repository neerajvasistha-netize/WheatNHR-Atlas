"""ad62 (Nov 2026): open-access full texts of the biotroph literature list (Europe PMC JATS XML) read automatically; records built by
scripts/curation/biblio/oa_build.py (curation 'auto', verify). Compilation records (Tong et al. 2024; Catalogue) of the same paper,
locus and chromosome are replaced by the original-paper records."""
import json, re
A = json.load(open('/tmp/bib/ad62.json'))
json.dump(A, open('docs/literature/ad62_open_access_records.json', 'w'), ensure_ascii=False)
def apply(rs, src):
    drop = set(A['deleted']); have = {r['id'] for r in rs}
    out = [r for r in rs if r['id'] not in drop] + [r for r in A['reports'] if r['id'] not in have]
    src.update(A['sources']); return out
ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
M = json.load(open('data/merged_v13f.json')); M['reports'] = apply(M['reports'], M['sources']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'] = apply(Sd['reports'], Sd['sources'])
open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
json.dump({'reports': A['reports'], 'sources': A['sources'], 'deleted': A['deleted']}, open('web/seed_ad62.json', 'w'), ensure_ascii=False, separators=(',', ':'))
print(len(M['reports']), 'records now;', len(A['reports']), 'added;', len(A['deleted']), 'replaced')
