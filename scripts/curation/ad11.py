"""ad11 (Nov 2026): Erginbas-Orakci et al. 2018 IJMS 19:2666 (crown rot, F. culmorum, 126 CIMMYT spring lines, DArT GWAS).
Table 3 has three FDR-significant MTAs. wPt-2193 and wPt-2766 are coincident (3B, 39.1 cM): one record (V0455, wPt-2193 has
no physical position), V0454 removed. Both kept records rest on one placed DArT marker, so grade C. Positions unchanged."""
import json, re
E = ' Corrected from the paper, Nov 2026 (ad11.py).'
PANEL = 'in 126 advanced CIMMYT spring bread wheat lines (107 genotyped with DArT), crown rot caused by Fusarium culmorum'
UPD = {
 'V0455': dict(name='wPt-2766', marker='wPt-2766, wPt-2193', grade='C', curation='hand', pve='11.4', stat='p 6.56e-05 (FDR 0.0123)', cls='MTA',
   note=(f'Marker-trait association for field crown score (GLM, Q matrix) {PANEL}. wPt-2193 and wPt-2766 are coincident at 39.1 cM on 3B '
         '(same p and r2); wPt-2193 has no physical position, so the locus is placed at the RefSeq v1.0 position of wPt-2766 (DArT tag BLAST). '
         'The authors place it away from the known 3BL crown rot QTL and consider it likely novel.' + E),
   evidence='Table 3: FCS (GLM) wPt-2193 and wPt-2766, 3B 39.1 cM, p 6.56e-05, FDR 0.01232, r2 11.4%'),
 'V0456': dict(name='wPt-669517', marker='wPt-669517', grade='C', curation='hand', pve='11.6', stat='p 3.24e-05 (FDR 0.0200)', cls='MTA',
   note=(f'Marker-trait association for greenhouse stem score (MLM, Q + kinship) {PANEL}; 2D, 95.7 cM. Single marker; position is the '
         'RefSeq v1.0 position of wPt-669517 (DArT tag BLAST). The authors note it could be the same as a previously reported 2D crown rot QTL.' + E),
   evidence='Table 3: GHSS (MLM) wPt-669517, 2D 95.7 cM, p 3.24e-05, FDR 0.01997, r2 11.6%'),
}
def build(reports):
    out = []
    for r in reports:
        if r['id'] == 'V0454': continue
        if r['id'] in UPD: r = dict(r, **UPD[r['id']])
        out.append(r)
    return out
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'] = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'] = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    ops = [{'op': 'set', 'path': 'reports/V0454', 'data': {'deleted': True}}] + [{'op': 'merge', 'path': 'reports/' + k, 'data': v} for k, v in UPD.items()]
    print(json.dumps(ops, ensure_ascii=False))
