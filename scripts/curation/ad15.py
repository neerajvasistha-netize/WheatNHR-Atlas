"""ad15 (Nov 2026): Neupane et al. 2024 BMC Plant Biol 24:852 (FHB, FDK, DON; DH D8006W/Superior, 90K).
Rebuilds the 30 records (Z0191-Z0194, V0285-V0310; most unplaced, two Z records on the wrong chromosome) as the 29 disease QTL of
Table 3 and Table S3. For each QTL the most significant row is used (LOD, PVE, position, one-LOD-drop CI; ICIM where given, else IM); the CI is
converted to RefSeq v1.0 with the map markers of Table S5 (cM and CS v1.0 bp printed in the supplement). Grade B (placed from the paper's map markers) where the markers
in the CI span <= 100 Mb, C where wider. PHT, DA and AE QTL not added."""
import json, re, sys, statistics, openpyxl
sys.path.insert(0, 'scripts/curation/sb19')
import lift
OLD = ['Z0191', 'Z0192', 'Z0193', 'Z0194'] + ['V%04d' % i for i in range(285, 311)]
XL = '/mnt/user-data/uploads/0e9f5df8-12870_2024_5536_MOESM1_ESM.xlsx'
MAP = {}
for r in list(openpyxl.load_workbook(XL, read_only=True)['Table S5 Genetic map data'].iter_rows(values_only=True))[2:]:
    if r[1] and isinstance(r[5], (int, float)): MAP.setdefault(r[4], []).append((float(r[3]), r[5] / 1e6, r[1]))
# name, LG, trait, site-year, method, pos cM, LOD, PVE, L, R, all detections
Q = [
 ('QFhb.ufg-1A', '1A1', 'FHB incidence and index', 'WPG17 (DI)', 'ICIM', 83, 5.53, 11.45, 80.5, 83.5, 'WPG17 DI and index'),
 ('QFhb.ufg-1B', '1B', 'FHB incidence and severity', 'CMN17 (DI)', 'ICIM', 22, 4.73, 9.55, 15.5, 23.5, 'CMN17 DI, CMN16 DS; FHB index Germany 2017 at 33-34 cM'),
 ('QFhb.ufg-2D.1', '2D', 'FHB incidence, severity and index', 'CMN17 (DI)', 'ICIM', 9, 6.02, 14.39, 1.5, 13.5, 'CMN17 DI, WPG17 DS and index'),
 ('QFhb.ufg-2D.2', '2D', 'FHB incidence, severity and index', 'RDG17 (Index)', 'ICIM', 57, 9.86, 14.07, 55.5, 61.5, 'RDG16 DI, CMN17 DS (PVE 23.93%), RDG17 index'),
 ('QFhb.ufg-4B', '4B', 'FHB severity and index', 'CMN16 (DS)', 'ICIM', 22, 6.66, 17.75, 21.5, 22.5, 'CMN16 DS and index; Superior allele increases FHB'),
 ('QFhb.ufg-5A.1', '5A2', 'FHB incidence, severity and index', 'CMN17 (DI)', 'ICIM', 14, 7.47, 18.37, 13.5, 15.0, 'CMN17 DI, WPG17 DS and index'),
 ('QFhb.ufg-5A.2', '5A3', 'FHB incidence, severity and index', 'WPG18 (DI)', 'ICIM', 50, 11.09, 34.34, 48.5, 51.5, 'WPG18 DI, CMN16 DS and index, Germany 2017 index'),
 ('QFhb.ufg-7A', '7A', 'FHB incidence, severity and index', 'CMN17 (DI)', 'ICIM', 97, 5.94, 14.04, 95.5, 98.5, 'CMN17 DI, WPG17 DS, RDG16 index (PVE 20.87%); Superior allele increases FHB'),
 ('QFhb.ufg-2D', '2D', 'FHB index (Germany)', 'G17', 'IM', 54, 3.06, 11.12, 49.5, 60.5, 'Germany 2017 only'),
 ('QFhb.ufg-6B', '6B', 'FHB index (Germany)', 'G17', 'ICIM', 8, 2.99, 10.89, 7.5, 8.5, 'Germany 2017 only'),
 ('QFdk.ufg-1D', '1D', 'Fusarium-damaged kernels', 'WPG16', 'ICIM', 49, 3.94, 10.22, 46.5, 50.5, 'WPG16'),
 ('QFdk.ufg-2B', '2B', 'Fusarium-damaged kernels', 'CMN17', 'ICIM', 1, 3.02, 6.18, 0, 5.5, 'CMN17'),
 ('QFdk.ufg-2D.1', '2D', 'Fusarium-damaged kernels', 'WPG17', 'ICIM', 8, 4.55, 15.89, 1.5, 13.5, 'WPG17'),
 ('QFdk.ufg-2D.2', '2D', 'Fusarium-damaged kernels', 'CMN16', 'ICIM', 69, 4.60, 21.67, 67.5, 81.5, 'CMN16'),
 ('QFdk.ufg-4B', '4B', 'Fusarium-damaged kernels', 'combined', 'ICIM', 18, 5.51, 13.48, 11.5, 18.5, 'combined environments'),
 ('QFdk.ufg-5A', '5A2', 'Fusarium-damaged kernels', 'WPG16', 'ICIM', 7, 10.27, 30.71, 6.5, 7.5, 'WPG16'),
 ('QFdk.ufg-5D', '5D2', 'Fusarium-damaged kernels', 'RDG17', 'ICIM', 113, 2.51, 13.54, 107.5, 118.5, 'RDG17'),
 ('QFdk.ufg-6B', '6B', 'Fusarium-damaged kernels', 'WPG17', 'ICIM', 11, 3.69, 12.57, 8.5, 14.5, 'WPG17'),
 ('QFdk.ufg-7A', '7A', 'Fusarium-damaged kernels', 'CMN17', 'ICIM', 94, 9.50, 23.12, 93.5, 95.5, 'CMN17'),
 ('QDon.ufg-1B', '1B', 'DON content', 'CMN17', 'ICIM', 45, 5.24, 15.43, 44.5, 45.5, 'CMN17'),
 ('QDon.ufg-1D', '1D', 'DON content', 'WPG16', 'ICIM', 53, 2.72, 8.01, 52.5, 53.0, 'WPG16'),
 ('QDon.ufg-2D.1', '2D', 'DON content', 'WPG17', 'ICIM', 8, 11.14, 25.57, 4.5, 11.5, 'WPG17'),
 ('QDon.ufg-2D.2', '2D', 'DON content', 'WPG16', 'IM', 60, 3.76, 12.00, 55.5, 65.5, 'WPG16'),
 ('QDon.ufg-4B', '4B', 'DON content', 'combined', 'ICIM', 27, 2.95, 5.55, 24.5, 27.5, 'combined environments'),
 ('QDon.ufg-4D', '4D1', 'DON content', 'combined', 'ICIM', 4, 7.76, 16.44, 0, 5.5, 'combined environments'),
 ('QDon.ufg-5A', '5A2', 'DON content', 'combined', 'ICIM', 8, 7.51, 16.95, 7.5, 8.5, 'combined environments (Table S3 only)'),
 ('QDon.ufg-6B', '6B', 'DON content', 'RDG16', 'ICIM', 3, 2.95, 9.01, 0, 4.5, 'RDG16 (Table S3 only)'),
 ('QDon.ufg-7A', '7A', 'DON content', 'combined', 'IM', 99, 3.72, 10.24, 96.5, 105.5, 'combined environments (Table S3 only); Superior allele reduces DON'),
 ('QDon.ufg-7D', '7D', 'DON content', 'WPG17', 'IM', 85, 3.57, 8.93, 77.5, 89.5, 'WPG17 (Table S3 only)'),
]
EX = {'QFhb.ufg-5A.1': ' The authors suggest it is likely Fhb5 (46.6-364.4 Mb, Steiner et al.), as it lies at 0-10 cM of 5A.2; the map markers in its CI, however, sit at about 489-508 Mb.',
      'QFhb.ufg-5A.2': ' The authors place it at about 700.4 Mb, away from Fhb5, and consider it likely novel; co-located with the anther-extrusion QTL QAE.ufg-5A.2.',
      'QFhb.ufg-2D.1': ' Most prominent QTL of the study; co-located with QPht.ufg-2D and QDa.ufg-2D.1 (Ppd-D1 region); D8006W allele reduces all FHB traits, plant height and days to anthesis.',
      'QFhb.ufg-4B': ' Unique to this population (Rht-B1 does not segregate).'}
def phys(lg, L, R, pos):
    ms = sorted(MAP[lg])
    inside = [m for m in ms if L <= m[0] <= R]
    if not inside:
        left = [m for m in ms if m[0] < L][-1:]; right = [m for m in ms if m[0] > R][:1]; inside = left + right
    p = [m[1] for m in inside]
    med = statistics.median(p)
    core = [x for x in p if abs(x - med) <= 50] or p       # drop markers that map far from the rest
    near = min(ms, key=lambda m: abs(m[0] - pos))
    return min(core), max(core), near, len(inside), len(p) - len(core)
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S92-')]
    new = []
    for n, (nm, lg, tr, sy, meth, pos, lod, pve, L, R, det) in enumerate(Q, 1):
        a, b, near, k, drop = phys(lg, L, R, pos); ch = lg[:2]
        g = 'B' if b - a <= 100 else 'C'
        r = dict(id=f'S92-{n:03d}', name=nm, cls='QTL', disease='Fusarium head blight', chr=ch, ref='RefSeq v1.0',
                 start=round(a, 4), end=round(b, 4), peak=None, marker=f'{near[2]} (nearest map marker to the peak)',
                 pve=str(pve), stat=f'LOD {lod} ({meth}, {sy})', source='Z_1308', grade=g, curation='hand')
        r['note'] = (f'QTL for {tr} in the D8006W/Superior soft white winter wheat DH population (105 lines; native FHB nurseries in Manitoba and '
                     f'Ontario 2016-2018, and Germany 2017). Detected in: {det}. Most significant row: {sy}, {meth}, {pos} cM on linkage group '
                     f'{lg[:2]}{"." + lg[2:] if len(lg) > 2 else ""}, LOD {lod}, PVE {pve}%, CI {L}-{R} cM. Physical interval is the RefSeq v1.0 span '
                     f'of the {k} map marker{"s" if k > 1 else ""} in the CI (Table S5 positions)' + (f', {drop} outlying marker(s) excluded' if drop else '') + '.'
                     + EX.get(nm, '') + (f' Interval width {b - a:.0f} Mb: a chromosome region, not a locus.' if g == 'C' else '') + ' Corrected from the paper, Nov 2026 (ad15.py).')
        r['evidence'] = ('Table S3: ' if 'S3 only' in det else 'Table 3: ') + f'{nm} {sy} {lg} {pos} cM, LOD {lod}, PVE {pve}%, CI {L}-{R} cM; Table S5 map positions'
        r['asm'] = lift.lift_record(r); new.append(r)
    return keep + new, new, OLD
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new, rem = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], new2, rem2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    ops = [{'op': 'set', 'path': 'reports/' + i, 'data': {'deleted': True}} for i in rem2]
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/home/claude/work/ad15.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['chr'], r['start'], r['end'], r['grade'], r['marker'][:22], r['asm'].get('RefSeq v2.1'))
