"""ad5: Pch4 (Arora et al. 2026 preprint) replaces the auto-extracted V1600 record.
The paper gives the locus only on the Ae. tauschii TA1675 assembly (1D 6.23-6.93 Mb);
no Chinese Spring position or CS-anchored marker is printed, so the record stays grade D."""
import json, re
FIX = {
 'name': 'Pch4', 'cls': 'DG', 'disease': 'Eyespot', 'chr': '1D', 'arm': '1DS',
 'start': None, 'end': None, 'peak': None, 'ref': 'RefSeq v1.0',
 'marker': 'Aet.TA1675.r1.1D001090', 'grade': 'D', 'curation': 'hand',
 'note': ('Cloned CC-NB-LRR (NLR) gene for resistance to Oculimacula yallundae, identified by k-mer GWAS '
          'in 151 Ae. tauschii lineage 2 accessions and validated in synthetic hexaploids, a biparental F2 and '
          'transgenic Chinese Spring; dosage-dependent dominant; not effective against O. acuformis; belongs to '
          'the Pm3 NLR lineage. Absent from 300 Watkins landraces. The associated peak is 6.23-6.93 Mb on 1DS of '
          'the Ae. tauschii TA1675 assembly (12 genes); the paper gives no Chinese Spring position, so no physical '
          'position is stored. Preprint, not peer reviewed. Corrected from the paper, Nov 2026 (ad5.py).'),
 'evidence': 'Results: k-mer GWAS peak 1DS 6.23-6.93 Mb (TA1675); Pch4 = Aet.TA1675.r1.1D001090; susceptible AL8/78 allele AET1Gv20020600',
 'asm': {'_note': 'no physical position in source'},
}
def build(reports):
    for r in reports:
        if r['id'] == 'V1600':
            r.update(FIX)
            return r
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); R = M['reports'] if isinstance(M, dict) else M
    r = build(R); json.dump(M, open('data/merged_v13f.json', 'w'))
    h = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', h, re.S)
    S = json.loads(m.group(2)); build(S['reports'])
    h = h[:m.start(2)] + json.dumps(S, separators=(',', ':')) + h[m.end(2):]
    open('web/index.html', 'w').write(h)
    d = {k: v for k, v in r.items() if v not in (None, '', [], {}) or k == 'asm'}
    json.dump([{'op': 'set', 'path': 'reports/V1600', 'data': d}], open('/home/claude/work/ad5.json', 'w'))
    print(json.dumps(d))
