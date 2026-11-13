"""Completion of four papers (Nov 2026): results not yet in the atlas.
- Pundir et al. 2022 Sci Rep (CCN, ITMI DH): digenic epistatic loci of Table S3 (one record per locus; partners listed in the note);
  placed between the flanking 90K/SSR markers with the sb26 resolver (grade B; one marker placed -> C; none -> D).
- Manan et al. 2025 TAG (BLS): GWAS MTAs in Table S5 not in Table 3 (positions printed on durum Svevo v1, so the Chinese Spring
  position of the 90K marker is used, grade B; marker not in the CS tables -> D); BP025 QTL QBls.zhl-6A (peak markers at 0 cM:
  wsnp_Ex_rep_c67563_66193104, CS 1.496 Mb; grade B).
- Muellner et al. 2021 TAG (common and dwarf bunt): QTL intervals of the single-population maps MP-BLI and MP-BON (Tables 2, 3;
  RefSeq v1.0, grade A); MP-CON intervals were already in the atlas.
Writes data/merged_v13f.json, web/index.html (seed) and /home/claude/work/ad4_db.json."""
import json, re, sys, io, contextlib
sys.path.insert(0, 'scripts/curation'); sys.path.insert(0, 'scripts/curation/sb19')
with contextlib.redirect_stdout(io.StringIO()):
    import sb26_fix as F, stb44
import lift, openpyxl

R = []
def base(src, name, ch, dis, cls, note, ev, **kw):
    r = dict(name=name, chr=ch, source=src, note=note, evidence=ev, ref='RefSeq v1.0', curation='hand', cls=cls, disease=dis,
             start=None, end=None, peak=None, grade=None, marker='', pve='', stat='')
    r.update(kw); return r
def placed(r, markers, peak=None):
    F.pos = stb44.P
    try: F.place(r, [m for m in markers if m], r['chr'], peak_marker=peak)
    finally: F.pos = stb44.POS0
    if r['peak'] is None and r['start'] is None: r['grade'] = 'D'
    elif peak and r['peak'] is not None: r['grade'] = 'B'
    for k in ('start', 'end', 'peak'):
        if r[k] is not None: r[k] = round(float(r[k]), 4)
    if r['start'] is not None and r['end'] - r['start'] > 100:
        r['grade'] = 'C'; r['note'] += ' The two flanking markers are far apart physically, so the interval is uncertain.'
    R.append(r)

# ---- Pundir 2022 epistatic loci ----
CCN = 'Cereal cyst nematode'
for e in json.load(open('/home/claude/fx/p22epi.json')):
    pos = '/'.join(str(p) for p in e['pos'])
    r = base('Pundir2022SR', 'Epistatic locus %s %s cM' % (e['chr'], pos), e['chr'], CCN, 'QTL',
             'Locus from the digenic (additive x additive) epistatic analysis of H. avenae resistance, ITMI doubled-haploid population (Table S3); not a main-effect QTL. Interacting with: %s.' % '; '.join(e['inter']),
             'Table S3: %s %s cM, %s-%s' % (e['chr'], pos, e['l'], e['r']), marker='%s, %s' % (e['l'], e['r']))
    placed(r, [e['l'], e['r']])

# ---- Manan 2025 extra MTAs (Table S5) ----
BLS = 'Bacterial leaf streak'
IN_T3 = {'RAC875_c13610_2646', 'BS00021688_51', 'Excalibur_c16470_623', 'Excalibur_c16470_494', 'Excalibur_c16470_278'}
wb = openpyxl.load_workbook('/mnt/user-data/uploads/74514a37-122_2025_5111_MOESM6_ESM.xlsx', data_only=True)
mod, seen = None, {}
for row in wb.active.iter_rows(values_only=True):
    if row[0] in ('Blink', 'MLMM'): mod = row[0]
    if not (row[2] and isinstance(row[4], (int, float))) or row[2] in IN_T3: continue
    m, c, p, pv, tr = row[2], row[3], row[4], row[5], row[1]
    seen.setdefault((m, c, p), []).append('%s %s (p %.2g)' % (mod, tr, pv))
for (m, c, p), hits in seen.items():
    r = base('Manan2025TAG', m, c, BLS, 'MTA', 'MTA for seedling resistance to bacterial leaf streak in the Global Durum Panel, listed in Table S5 only (%s). Position in the paper is on the durum Svevo v1 assembly (%s bp); the Chinese Spring position of the 90K marker is used.' % ('; '.join(hits), format(p, ',')),
             'Table S5: %s %s %s bp (Svevo)' % (m, c, format(p, ',')), marker=m)
    placed(r, [], peak=m)
r = base('Manan2025TAG', 'QBls.zhl-6A (BP025)', '6A', BLS, 'QTL',
         'QTL for bacterial leaf streak (infection type, LOD 11.0 CIM; water-soaked area) in the durum BP025 population, 0-2 cM between wsnp_Ex_c16090 and wsnp_CAP_C1339. The peak is defined by three co-segregating SNPs at 0 cM; position is that of wsnp_Ex_rep_c67563_66193104 (Svevo 1,197,947 bp), the one of them placed on 6A in Chinese Spring (wsnp_Ex_c16090_24522660 maps to 6B in the CS marker tables).',
         'Table 4: QBls.zhl-6A(BPIT) 0-14 cM Ex_c16090-Ex_C280 (SIM, LOD 9.4, R2 21.0); CIM 0-2 cM Ex_c16090-CAP_C1339, LOD 11.0; text: peak SNPs at 0 cM',
         marker='wsnp_Ex_rep_c67563_66193104', pve='21.0 (SIM)', stat='LOD 11.0 (CIM)')
placed(r, [], peak='wsnp_Ex_rep_c67563_66193104')

# ---- Muellner 2021 population-specific intervals ----
CB, DB = 'Common bunt', 'Dwarf bunt'
MU = [(CB, 'MP-BLI', 'QBt.ifa-1BS', '1B', 4.35, 26.19, '41.8', '18.6'), (CB, 'MP-BLI', 'QBt.ifa-1AL', '1A', 492.81, 515.55, '36.9', '17.0'),
      (CB, 'MP-BLI', 'QBt.ifa-7AL', '7A', 721.22, 736.69, '22.8', '11.8'), (CB, 'MP-BON', 'QBt.ifa-1BS', '1B', 17.25, 38.91, '51.1', '19.1'),
      (CB, 'MP-BON', 'QBt.ifa-1AL', '1A', 490.09, 513.89, '28.2', '12.8'), (CB, 'MP-BON', 'QBt.ifa-7AL', '7A', 722.30, 736.69, '15.8', '8.2'),
      (DB, 'MP-BLI', 'QBt.ifa-1AL', '1A', 492.81, 516.67, '28.2', '10.9'), (DB, 'MP-BLI', 'QBt.ifa-7AL', '7A', 722.30, 736.69, '16.5', '7.0'),
      (DB, 'MP-BON', 'QBt.ifa-1AL', '1A', 380.97, 516.67, '28.2', '7.2'), (DB, 'MP-BON', 'QBt.ifa-7AL', '7A', 722.30, 736.69, '12.5', '3.5')]
for dis, pop, q, c, a, b, pv, lod in MU:
    pn = {'MP-BLI': 'Blizzard x Rainer', 'MP-BON': 'Bonneville x Rainer'}[pop]
    R.append(base('Muellner2021TAG', '%s (%s)' % (q, pop), c, dis, 'QTL',
                  'QTL for %s incidence in the %s RIL population (%s), multiple QTL mapping across field experiments; interval printed on IWGSC RefSeq v1.0. The consensus-map (MP-CON) interval is a separate record.' % (dis.lower(), pop, pn),
                  'Table %s: %s %s %s-%s Mbp' % ('2' if dis == CB else '3', pop, c, a, b), start=a, end=b, grade='A', pve=pv, stat='LOD ' + lod))

SID = lambda i: 'S83-%03d' % i
def build(reports):
    new = []
    for i, x in enumerate(R, 1):
        r = dict(x); r['id'] = SID(i); r['note'] = (r['note'] + ' Added from the paper, Nov 2026 (ad4.py).').strip()
        r['asm'] = lift.lift_record(r)
        if (r['asm'].get('RefSeq v2.1') or {}).get('flag') == 'wide_gap':
            import ts52
            r['asm']['RefSeq v2.1'].update(ts52.wide_fix(r))
        new.append(r)
    return [r for r in reports if not r['id'].startswith('S83-')] + new, new

if __name__ == '__main__':
    import collections
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], new = build(M['reports'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</'))
    Sd['reports'], new2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs}, open('/home/claude/work/ad4_db.json', 'w'), ensure_ascii=False)
    print(len(new), collections.Counter((r['source'], r['grade']) for r in new),
          collections.Counter((r['asm'].get('RefSeq v2.1') or {}).get('flag', 'none') for r in new))
    for r in new:
        if r['source'] == 'Manan2025TAG': print(r['name'], r['chr'], r['peak'], r['grade'])
