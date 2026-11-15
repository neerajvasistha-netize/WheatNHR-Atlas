"""ad49 (Nov 2026): Tessmann et al. 2019 Crop Sci Supplementary Table S3 gives 90K SNP names for the panel-internal marker IDs (M####).
Applied to Tessmann et al. 2019 (S177) and to the shared IDs in Tessmann and Van Sanford 2018 (S176, same panel)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
N = {'M12955': 'Kukri_c15912_2330', 'M12957': 'Tdurum_contig47317_100', 'M12960': 'BS00085556_51', 'M9385': 'Excalibur_c32979_1152', 'M5248': 'wsnp_Ex_c15269_23491104',
     'M3268': 'wsnp_Ex_c21092_30220342', 'M9432': 'Kukri_c4349_639', 'M6959': 'RAC875_c88582_131', 'M11423': 'wsnp_Ex_c3025_5587183', 'M11046': 'RAC875_c63209_154'}
IDS = ['S177-001', 'S177-002', 'S177-003', 'S177-004', 'S177-005', 'S177-006', 'S177-007', 'S176-007', 'S176-008', 'S176-010']
EXTRA = {'Excalibur_c32979_1152': ' Excalibur_c32979_1152 places only on 5D (520.721 Mb, Sun et al. 2020 PBJ (90K)), not on 5B; left unplaced.',
         'Kukri_c4349_639': ' Kukri_c4349_639 has no physical position in the atlas marker sources; left unplaced.'}
RX = re.compile(r' Markers are panel-internal SNP IDs; positions as printed[^.]*?, so not placed\.')
def fix(r):
    r = dict(r); ids = re.findall(r'M\d+', r['marker']); names = [N[i] for i in ids]
    ok = [(n,) + resolve(n, r['chr']) for n in names]; ok = [x for x in ok if x[1]]
    for k in ('start', 'end', 'peak'): r.pop(k, None)
    t = ' Supplementary Table S3 of Tessmann et al. (2019) names the SNPs: ' + ', '.join(f'{i} = {n}' for i, n in zip(ids, names)) + '; positions as printed are cM.'
    if len(ok) >= 2:
        lo, hi = min(x[1] for x in ok), max(x[1] for x in ok); r.update(start=lo, end=hi, grade='B' if hi - lo <= 100 else 'C')
        t += ' Marker positions (RefSeq v1.0): ' + '; '.join(f'{n} {p:.3f} Mb ({h})' for n, p, h in ok) + '.'
    elif ok:
        n, p, h = ok[0]; r.update(peak=p, grade='C'); t += f' Placed at {n} ({p:.3f} Mb, {h}).'
    else:
        r['grade'] = 'D'; t += ''.join(EXTRA[n] for n in names)
    r['note'] = RX.sub(t, r['note']).replace('(ad38.py).', '(ad38.py); SNP names from Supplementary Table S3 of Tessmann et al. (2019), Nov 2026 (ad49.py).')
    assert 'Table S3' in r['note'], r['id']
    r['marker'] = ', '.join(f'{i} ({n})' for i, n in zip(ids, names)); r['evidence'] += '; Tessmann et al. (2019) Table S3'
    r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
    return r
def build(rs):
    out, up = [], []
    for r in rs:
        if r['id'] in IDS: r = fix(r); up.append(r)
        out.append(r)
    return out, up
M = json.load(open('data/merged_v13f.json')); M['reports'], upd = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
ops = []
for r in upd:
    d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}; d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
    ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'ifVersion': 1, 'data': d}); print(r['id'], r['grade'], r.get('start'), r.get('end'), r.get('peak'), d['asm'])
json.dump(ops, open('/tmp/ops49.json', 'w'), ensure_ascii=False)
