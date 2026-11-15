"""ad51 (Nov 2026): Zanke et al. 2017 TAG Online Resource 3 (SSR MTAs, -log10P > 3). Places the 2AS region (S185-003) with its SSRs and
adds the SSR MTAs that are significant on BLUEs outside 7D (Pch1 region, already S185-001) and 5A BARC303 (S185-002). Unmapped SSRs skipped."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = ' Read from the paper and Online Resource 3, Nov 2026 (ad51.py).'
BASE = 'GWAS for eyespot (field tests in three environments) in 168 European winter and 3 spring wheat varieties with 732 SSRs and 7,761 mapped 90K SNPs; 108 SSR and 235 SNP MTAs (-log10P 3.0 or more). '
NEW = [('S199-001', '1D', 'WMS4780a', 'Xgwm4780', 88.6, '4.02', '14AHL, 14HEY, BLUEs', '9.7', 'the null allele reduces eyespot (allele 97 increases it)'),
       ('S199-002', '2D', 'WMS0539', 'Xgwm539', 135.5, '3.44', 'BLUEs', '9.2', 'allele 129 increases eyespot'),
       ('S199-003', '3D', 'CFD0036a', 'Xcfd36', 8.1, '4.42', '14AHL, 14HEY, BLUEs', '10.6', 'allele 150 reduces eyespot (allele 142 increases it)'),
       ('S199-004', '5A', 'WMS0156b', 'Xgwm156', 98.4, '3.27', '14HEY, BLUEs', '5.6', 'the null allele increases eyespot'),
       ('S199-005', '5D', 'WMS1527', 'Xgwm1527', 63.2, '3.25', '14HEY, BLUEs', '5.0', 'allele 152 reduces eyespot'),
       ('S199-006', '6A', 'WMC0256b', 'Xwmc256', 89.3, '3.61', '14HEY, BLUEs', '6.3-8.8', 'alleles 132 and 134 increase eyespot'),
       ('S199-007', '7A', 'WMC0488b', 'Xwmc488', 143.7, '3.74', '14AHL, 14HEY, BLUEs', '7.8', 'allele 137 increases eyespot')]
EFF = {'res': 'the minor/null allele reduces eyespot', 'sus': 'the scored allele increases eyespot (the alternative allele is the resistant one)', 'mix': 'alleles 142 and 150 have opposite effects (150 reduces eyespot)'}
def mk(i, ch, raw, m, cm, lp, envs, r2, e):
    r = dict(id=i, cls='MTA', ref='RefSeq v1.0', curation='hand', name=f'Eyespot SSR MTA {ch} ({m[1:]}; Zanke 2017)', disease='Eyespot', chr=ch, source='Zanke2017TAG',
             marker=m, pve=r2, stat=f'-log10P {lp} (BLUEs)', evidence=f'Online Resource 3A/3B: {raw} {ch} {cm} cM')
    p, h = resolve(m, ch)
    t = f'SSR {raw} ({cm} cM on the consensus map) significant in {envs}; R2 {r2}%; {e}.'
    if p: r.update(peak=p, grade='C'); t += f' Placed at {m} ({p:.3f} Mb, {h}); single-marker position.'
    else: r['grade'] = 'D'; t += f' {m} could not be placed on {ch}; left unplaced.'
    r['note'] = BASE + t + T; r['asm'] = lift.lift_record(r) if p else {'_note': 'no physical position in source'}; return r
def upd2A(r):
    r = dict(r); a, ha = resolve('Xgwm636', '2A'); b, hb = resolve('Xgwm359', '2A')
    r.update(start=a, end=b, grade='B', marker='Xgwm830, Xgwm636, Xgwm359', evidence='Results; Online Resource 3A/3B',
             note=BASE + 'SSR and SNP data indicated a QTL region on 2AS supported by numerous markers, most alleles reducing eyespot score (smaller effects than Pch1); proposed as a novel resistance source. Online Resource 3 lists the SSRs WMS0830 (15.4 cM, -log10P up to 4.58, R2 9.8%), WMS0636 (37.3 cM, 4.30, 9.5%) and WMS0359 (52.4 cM, 3.23, 7.9%), each with an allele that reduces eyespot (WMS0830 null, WMS0636 106, WMS0359 215). Xgwm830 could not be placed. Marker positions (RefSeq v1.0): ' + f'Xgwm636 {a:.3f} Mb ({ha}); Xgwm359 {b:.3f} Mb ({hb}).' + T)
    r['asm'] = lift.lift_record(r); return r
def build(rs):
    ids = {n[0] for n in NEW}; out, up = [], []
    for r in rs:
        if r['id'] in ids: continue
        if r['id'] == 'S185-003': r = upd2A(r); up.append(r)
        out.append(r)
    new = [mk(*n) for n in NEW]; return out + new, up + new
M = json.load(open('data/merged_v13f.json')); M['reports'], upd = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
ops = []
for r in upd:
    d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}; d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
    o = {'op': 'set', 'path': 'reports/' + r['id'], 'data': d}
    if r['id'] == 'S185-003': o['ifVersion'] = 1
    ops.append(o); print(r['id'], r['grade'], r.get('start'), r.get('end'), r.get('peak'), d['asm'], d.get('cls'))
json.dump(ops, open('/tmp/ops51.json', 'w'), ensure_ascii=False)
