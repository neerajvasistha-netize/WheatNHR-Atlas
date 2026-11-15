"""ad54 (Nov 2026): Zanke et al. 2017 TAG Online Resource 4B (SNP MTAs significant on BLUEs, -log10P > 3). Grouped by chromosome and cM cluster;
2AS cluster refines S185-003; 7D SNPs (613.6-630.5 Mb) are the Pch1 region (S185-001) and are not entered separately."""
import json, re, sys, openpyxl
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = ' Read from the paper and Online Resource 4B, Nov 2026 (ad54.py).'
BASE = 'GWAS for eyespot (field tests in three environments) in 168 European winter and 3 spring wheat varieties with 732 SSRs and 7,761 mapped 90K SNPs; 108 SSR and 235 SNP MTAs (-log10P 3.0 or more). '
wb = openpyxl.load_workbook('docs/supplements/zanke2017_ESM4_SNP_MTAs.xlsx', data_only=True); ws = wb['SNP_BLUEs']
ROWS = [dict(m=r[0], ch=str(r[1]).strip(), cm=float(r[2]), lp=r[5], r2=r[14] * 100, eff=r[15]) for r in ws.iter_rows(values_only=True, min_row=3) if r[0]]
G = [('S185-003', '2A', 0, 30, None), ('S199-008', '1B', 50, 65, 'Eyespot SNP MTAs 1BL (54-61 cM; Zanke 2017)'), ('S199-009', '2A', 100, 110, 'Eyespot SNP MTA 2AL (105 cM; Zanke 2017)'),
     ('S199-010', '2B', 125, 135, 'Eyespot SNP MTA 2B (132 cM; Zanke 2017)'), ('S199-011', '2D', 125, 130, 'Eyespot SNP MTAs 2D (127-129 cM; Zanke 2017)'), ('S199-012', '3B', 50, 55, 'Eyespot SNP MTAs 3B (52.8 cM; Zanke 2017)')]
def grp(ch, lo, hi):
    rs = [x for x in ROWS if x['ch'] == ch and lo <= x['cm'] <= hi]
    for x in rs: x['p'], x['h'] = resolve(x['m'], ch)
    return rs
def summ(rs):
    cms = sorted({x['cm'] for x in rs}); best = max(rs, key=lambda x: x['lp']); neg = sum(1 for x in rs if x['eff'] < 0)
    return (f'Online Resource 4B lists {len(rs)} SNP{"s" if len(rs) > 1 else ""} at {", ".join(f"{c:g}" for c in cms)} cM (strongest {best["m"]}, -log10P {best["lp"]:.2f}, R2 {best["r2"]:.1f}%); ' +
            (f'{neg} of {len(rs)} with the minor allele reducing eyespot.' if len(rs) > 1 else ('its minor allele reduces eyespot.' if neg else 'its minor allele increases eyespot.')), best, f'{min(x["r2"] for x in rs):.1f}-{max(x["r2"] for x in rs):.1f}' if len({round(x["r2"], 1) for x in rs}) > 1 else f'{best["r2"]:.1f}')
def place(r, rs, extra=''):
    ok = [x for x in rs if x['p']]
    for k in ('start', 'end', 'peak'): r.pop(k, None)
    if len(ok) >= 2:
        lo, hi = min(x['p'] for x in ok), max(x['p'] for x in ok); r.update(start=lo, end=hi, grade='B' if hi - lo <= 100 else 'C')
        srcs = sorted({re.sub(r' \(v2\.1 [^)]*\)', '', x['h']) for x in ok}); t = f' {len(ok)} SNPs placed in silico at {lo:.3f}-{hi:.3f} Mb (RefSeq v1.0; {"; ".join(srcs)}).'
    elif ok:
        x = ok[0]; r.update(peak=x['p'], grade='C'); t = f' Placed at {x["m"]} ({x["p"]:.3f} Mb, {x["h"]}); single-marker position.'
    else: r['grade'] = 'D'; t = ' None of the SNPs could be placed on ' + r['chr'] + '; left unplaced.'
    r['asm'] = lift.lift_record(r) if ok else {'_note': 'no physical position in source'}; return t + extra
def build(rs):
    ids = {g[0] for g in G}; out, up = [], []
    for r in rs:
        if r['id'] in ids and r['id'] != 'S185-003': continue
        if r['id'] == 'S185-003':
            r = dict(r); s = grp('2A', 0, 30); t, best, r2 = summ(s); core = [x['p'] for x in s if x['p'] and x['cm'] == 10.5]
            r['note'] = r['note'].split(' Read from the paper and Online Resource 3')[0] + ' ' + t + place(r, s, f' The 10.5 cM block ({sum(1 for x in s if x["cm"] == 10.5)} SNPs, the peak of the region) places at {min(core):.3f}-{max(core):.3f} Mb.') + ' SSR data from Online Resource 3 (ad51.py); SNP data from Online Resource 4B, Nov 2026 (ad54.py).'
            r['marker'] = 'Xgwm830, Xgwm636, Xgwm359; ' + best['m'] + ' and 43 other 90K SNPs'; r['evidence'] = 'Results; Online Resource 3A/3B; Online Resource 4B'; up.append(r)
        out.append(r)
    for i, ch, lo, hi, nm in G[1:]:
        s = grp(ch, lo, hi); t, best, r2 = summ(s)
        r = dict(id=i, cls='MTA', ref='RefSeq v1.0', curation='hand', name=nm, disease='Eyespot', chr=ch, source='Zanke2017TAG', marker=', '.join(x['m'] for x in s), pve=r2,
                 stat=f'-log10P up to {best["lp"]:.2f} (BLUEs)', evidence=f'Online Resource 4B: {ch} {lo}-{hi} cM')
        ex = ' The 2D SNPs have the same statistics as the 2A 10.5 cM block (R2 10.9%), so they may be homoeologous calls of that locus.' if ch == '2D' else ''
        r['note'] = BASE + t + place(r, s, ex) + T; out.append(r); up.append(r)
    return out, up
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], upd = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
    open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
    ops = []
    for r in upd:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}; d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
        o = {'op': 'set', 'path': 'reports/' + r['id'], 'data': d}
        if r['id'] == 'S185-003': o['ifVersion'] = 2
        ops.append(o); print(r['id'], r['grade'], r.get('start'), r.get('end'), r.get('peak'), d['asm'])
    json.dump(ops, open('/tmp/ops54.json', 'w'), ensure_ascii=False)
