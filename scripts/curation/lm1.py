"""lm1 (Nov 2026): store the linkage maps of QTL studies (marker, chromosome, cM) with each marker's physical position.
Physical position: the paper's own RefSeq v1.0 bp where given, otherwise placed in silico (insilico_res); then lifted to RefSeq v2.1.
A marker is kept even when unplaced. Output: data/linkage_maps/<source>__<population>.csv and data/linkage_maps/index.json.
Maps transcribed from figures are entered as cumulative cM built from the interval distances printed beside each chromosome."""
import json, csv, re, sys, time, openpyxl, os
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
U = '/mnt/user-data/uploads/'; OUT = 'data/linkage_maps/'
def chrom(x):
    m = re.match(r'([1-7][ABD])', str(x or '')); return m.group(1) if m else ''
def rows_xlsx(f, sh, start, cols):
    ws = openpyxl.load_workbook(U + f, read_only=True)[sh]
    for i, r in enumerate(ws.iter_rows(values_only=True)):
        if i < start: continue
        yield [r[c] if c is not None and c < len(r) else None for c in cols]
def fig(spec):  # spec: {chr: [marker, d, marker, d, ...]} interval distances
    for ch, seq in spec.items():
        cm = 0.0
        for k, x in enumerate(seq):
            if k % 2: cm += x
            else: yield [x, ch, round(cm, 1), None]
MAPS = []
def add(src, pop, origin, rows):
    MAPS.append((src, pop, origin, list(rows)))
# --- supplementary tables
add('Z_1308', 'D8006W x Superior DH', 'Additional file 1, Table S5 (1,254 SNPs used for QTL analysis; bp on RefSeq v1.0)',
    ([m, lg, cm, bp] for m, lg, cm, bp in rows_xlsx('0e9f5df8-12870_2024_5536_MOESM1_ESM.xlsx', 'Table S5 Genetic map data', 2, [1, 4, 3, 5]) if m))
add('V_Thambugala2020', 'Cutler x AC Barrie RIL', 'ESM 1, Table S5 (10,177 SNPs)',
    ([m, c, cm, None] for m, c, cm in rows_xlsx('d1ee61c7-122_2020_3631_MOESM1_ESM.xlsx', 'Table S5 Linkage Map', 5, [0, 1, 2]) if m))
add('V_Thambugala2020', 'AC Barrie x Reeder DH', 'ESM 1, Table S11 (10,877 SNPs)',
    ([m, c, cm, None] for m, c, cm in rows_xlsx('d1ee61c7-122_2020_3631_MOESM1_ESM.xlsx', 'Table S11 Linkage Map', 5, [0, 1, 2]) if m))
T7 = list(rows_xlsx('cafa9c87-122_2020_3708_MOESM2_ESM.xlsx', 'Table S7', 4, [0, 3, 4, 6, 7, 9, 10, 12, 13]))
for pop, a in [('MP-CON (Muellner 2021 multiparent)', 3), ('MP-BLI (Blizzard)', 5), ('MP-BON (Bonneville)', 7)]:
    add('Muellner2021TAG', pop, 'ESM 2, Table S7 (bp on RefSeq v1.0)',
        ([r[0], r[a], r[a + 1], r[2]] for r in T7 if r[0] and r[a + 1] not in (None, '-')))
# --- figures
add('V_Jia2005', "Wangshuibai x Alondra's DH", 'Figure 2 (chromosomes with QTL only; interval distances in cM)', fig({
 '2D': ['Xgwm455', 2.6, 'Xbarc95', 6.2, 'Xgwm261', 16.4, 'Xgwm484'],
 '3B': ['Xgwm389', 2.9, 'Xbarc75', 3.0, 'Xgwm533', 9.8, 'Xgwm493', 11.8, 'Xbarc87', 7.3, 'Xbarc102', 42.7, 'Xgwm284', 10.4, 'Xgwm72', 58.1, 'Xgwm131'],
 '4B': ['Xbarc20', 8.1, 'Xgwm368', 22.7, 'Xgwm149', 15.2, 'Xgwm251', 6.8, 'Xbarc60', 16.8, 'Xgwm6'],
 '5B': ['Xgwm234', 6.7, 'Xgwm443', 10.8, 'Xbarc32', 3.0, 'Xgwm544', 7.4, 'Xgwm213', 4.8, 'Xgwm335', 19.4, 'Xgwm371', 3.7, 'Xgwm67', 8.2, 'Xbarc89'],
 '7A': ['Xgwm666', 7.0, 'Xgwm350', 29.6, 'Xgwm60', 9.7, 'Xgwm130', 78.5, 'Xbarc49', 9.1, 'Xbarc29', 7.7, 'Xgwm276', 19.3, 'Xgwm282', 6.8, 'Xgwm332']}))
if __name__ == '__main__':
    idx = []
    for src, pop, origin, rows in MAPS:
        out = []; n_pl = 0
        for m, lg, cm, bp in rows:
            ch = chrom(lg); m = str(m).strip(); p, how = (None, '')
            if bp not in (None, '', '-') and str(bp).replace('.', '').isdigit() and ch:
                p, how = float(bp) / 1e6, 'paper (RefSeq v1.0 bp)'
            elif ch:
                p, how = resolve(m, ch)
            v2 = None
            if p is not None:
                a = lift.lift_record({'chr': ch, 'ref': 'RefSeq v1.0', 'peak': p}).get('RefSeq v2.1') or {}
                if a.get('chr') == ch: v2 = a.get('peak')
                n_pl += 1
            out.append([m, str(lg), ch, '' if cm in (None, '') else round(float(cm), 2), '' if p is None else round(p, 4), how, '' if v2 is None else v2])
        fn = re.sub(r'[^A-Za-z0-9]+', '_', f'{src}__{pop}').strip('_') + '.csv'
        with open(OUT + fn, 'w', newline='') as fh:
            w = csv.writer(fh); w.writerow(['marker', 'linkage_group', 'chr', 'cM', 'Mb_RefSeq_v1.0', 'v1.0_source', 'Mb_RefSeq_v2.1']); w.writerows(out)
        idx.append({'source': src, 'population': pop, 'origin': origin, 'file': fn, 'markers': len(out), 'placed': n_pl})
        print(fn, len(out), n_pl, flush=True)
    json.dump(idx, open(OUT + 'index.json', 'w'), indent=1)
