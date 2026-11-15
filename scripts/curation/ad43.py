"""ad43 (Nov 2026): Jin et al. 2020 TAG 133:2063, crown rot GWAS; all 104 significant SNPs from Table S6 (ESM9; positions on RefSeq v1.0)
grouped into regions (SNPs on one chromosome arm within 25 Mb of each other). The 5DL SNPs are summarised in the existing S139-001 record."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); import lift
S = json.load(open('docs/supplements/jin2020_tableS6.json'))
T = ' Read from Table S6 (ESM9), Nov 2026 (ad43.py).'
BASE = 'GWAS of seedling crown rot resistance (F. pseudograminearum) in 358 Chinese wheat germplasms with the 55K SNP array, three trials; significant SNPs (P < 0.001) from Table S6 with RefSeq v1.0 positions as printed. '
def clusters():
    cl = []
    for s in sorted(S, key=lambda x: (x[1], x[2])):
        if cl and cl[-1][-1][1] == s[1] and s[2] - cl[-1][-1][2] <= 25: cl[-1].append(s)
        else: cl.append([s])
    return cl
def build(reports):
    keep = [r for r in reports if not r['id'].startswith('S196-')]; N = []; n = 0
    for c in clusters():
        arm = c[0][1]
        if arm == '5DL': continue
        n += 1; b = min(c, key=lambda x: x[3]); ch = arm[:2]
        r = dict(id=f'S196-{n:03d}', cls='MTA', ref='RefSeq v1.0', curation='hand', disease='Fusarium crown rot', chr=ch, source='Jin2020TAG', grade='A',
                 pve=f'{max(x[4] for x in c):.2f}', stat=f'P {b[3]:.1e} ({b[0]})', evidence='Table S6 (ESM9)')
        if len(c) == 1:
            r.update(name=f'QFcr {arm} {b[2]:.1f} Mb (Chinese germplasm GWAS)', peak=round(b[2], 6), marker=b[0],
                     note=BASE + f'Single SNP {b[0]} on {arm} at {b[2]:.6f} Mb, P {b[3]:.1e}, R2 {b[4]:.2f}%.' + T)
        else:
            lo, hi = c[0][2], c[-1][2]
            r.update(name=f'QFcr {arm} {lo:.0f}-{hi:.0f} Mb (Chinese germplasm GWAS)', start=round(lo, 6), end=round(hi, 6), marker=', '.join(x[0] for x in c[:6]) + (f' and {len(c) - 6} more' if len(c) > 6 else ''),
                     note=BASE + f'{len(c)} SNPs on {arm} from {lo:.3f} to {hi:.3f} Mb; strongest {b[0]} at {b[2]:.3f} Mb (P {b[3]:.1e}), R2 up to {max(x[4] for x in c):.2f}%.' + (' Adjacent to the 5DL SNP cluster (189.7-246.8 Mb, record S139-001) across the centromere.' if arm == '5DS' else '') + (f' The SNPs span {hi - lo:.0f} Mb: a chromosome region, not a locus.' if hi - lo > 100 else '') + T)
        N.append(r)
    for r in N: r['asm'] = lift.lift_record(r)
    out = []
    for r in keep:
        if r['id'] == 'S139-001' and 'Table S6 (ESM9)' not in r['note']:
            r = dict(r); c5 = [x for x in S if x[1] == '5DL']; b = min(c5, key=lambda x: x[3])
            r['note'] = r['note'].replace('(P < 0.001; listed only in Table S6, not uploaded)', '(P < 0.001; Table S6)').replace(' Read from the paper, Nov 2026 (ad29.py).',
                f' Table S6 lists {len(c5)} significant 5DL SNPs from {min(x[2] for x in c5):.3f} to {max(x[2] for x in c5):.3f} Mb, strongest {b[0]} at {b[2]:.3f} Mb (P {b[3]:.1e}, R2 {b[4]:.2f}%); the other regions are records S196. Read from the paper and Table S6, Nov 2026 (ad29.py, ad43.py).')
            r['evidence'] += '; Table S6 (ESM9)'; r['_upd'] = True
        out.append(r)
    return out + N, N
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new = build(M['reports'])
    upd = [r for r in M['reports'] if r.pop('_upd', None)]
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
    for r in Sd['reports']: r.pop('_upd', None)
    open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
    ops = []
    for r in new + upd:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
        o = {'op': 'set', 'path': 'reports/' + r['id'], 'data': d}
        if r['id'] == 'S139-001': o['ifVersion'] = 1
        ops.append(o)
    json.dump(ops, open('/tmp/ops43.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['asm'].get('RefSeq v2.1'))
    print(len(ops), [r['id'] for r in upd])
