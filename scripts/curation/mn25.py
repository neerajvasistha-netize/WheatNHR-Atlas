"""Position fix (Nov 2026): Manan et al. 2025 TAG (bacterial leaf streak, durum).
The paper gives all positions on the durum Svevo v1 assembly (Table 3 footnote; Table S5), not Chinese Spring.
The atlas had stored them as RefSeq v1.0, which put the 5A GWAS markers about 41 Mb too proximal.
- 90K markers: Chinese Spring RefSeq v1.0 positions from the marker resolver (grade B).
- Excalibur_c16470 markers (not in the CS marker tables): 0.103 Mb proximal to BS00021688_51 on Svevo; placed by that offset (grade C).
- QBls.zhl-6A.1: CIM interval 6A_1729034-6A_6608308 (Svevo) shifted by the local Svevo-to-CS offset (+0.30 Mb, from
  wsnp_Ex_rep_c67563_66193104 and RAC875_c13610_2646 at Svevo 1.198/1.203 Mb, CS 1.496/1.501 Mb) (grade C).
- QBls.zhl-6A.2 and QBls.zhl-3A: Svevo coordinates at chromosome ends, no CS anchor in the paper; Svevo positions used as an
  approximation (grade C).
- QBls.zhl-5A: start was a RefSeq v2.1 value stored as v1.0; now both flanks on RefSeq v1.0.
- Duplicate MTA records (one per trait) merged.
Writes data/merged_v13f.json, web/index.html (seed) and /home/claude/work/mn25_db.json."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift

SV = ' Position in the paper is on the durum Svevo v1 assembly (Table 3, Table S5).'
FIX = [  # (old ids, new fields)
 (['V3250', 'V3255'], dict(peak=1.5008, grade='B', curation='hand', stat='AvIT p 2.11e-08, q 0.00028; Av%WS p 2.53e-07',
     pve='6.72 (AvIT); 5.65 (Av%WS)',
     note='MTA for seedling resistance to bacterial leaf streak in the Global Durum Panel (average infection type and average % water-soaked area).' + SV +
          ' Chinese Spring RefSeq v1.0 position of RAC875_c13610_2646 (Sun et al. 2020, 90K) used; Svevo 1,202,829 bp. Same gene as the BP025 QTL peak marker wsnp_Ex_rep_c67563_66193104.',
     evidence='Table 3: RAC875_c13610_2646 6A 1,202,829 bp (Svevo)')),
 (['V3251', 'V3256'], dict(peak=698.0032, grade='B', curation='hand', stat='AvIT p 4.81e-06; Av%WS p 6.25e-06', pve='4.42 (AvIT); 4.31 (Av%WS)',
     note='MTA for seedling resistance to bacterial leaf streak in the Global Durum Panel (average infection type and average % water-soaked area).' + SV +
          ' Chinese Spring RefSeq v1.0 position of BS00021688_51 (Sun et al. 2020, 90K) used; Svevo 657,084,054 bp.',
     evidence='Table 3: BS00021688_51 5A 657,084,054 bp (Svevo)')),
] + [([i], dict(peak=697.9006, grade='C', curation='hand',
     note='MTA for seedling resistance to bacterial leaf streak in the Global Durum Panel (average infection type, MLMM).' + SV +
          ' %s is not in the Chinese Spring marker tables; it lies 0.103 Mb proximal to BS00021688_51 on Svevo, so it is placed 0.103 Mb proximal to the CS position of BS00021688_51 (698.003 Mb).' % m,
     evidence='Table 3: %s 5A %s bp (Svevo)' % (m, p))) for i, m, p in
     [('V3252', 'Excalibur_c16470_623', '656,981,441'), ('V3253', 'Excalibur_c16470_494', '656,981,570'), ('V3254', 'Excalibur_c16470_278', '656,981,786')]] + [
 (['V3257'], dict(start=2.0291, end=6.9081, grade='C', curation='hand', marker='6A_1729034-6A_6608308; peak 6A_1875580',
     note='QTL for bacterial leaf streak (infection type and % water-soaked area) in the durum RP336 population. Positions in the marker names are on the durum Svevo v1 assembly.'
          ' CIM interval 6A_1729034-6A_6608308 (peak marker 6A_1875580), shifted by the local Svevo-to-Chinese Spring offset (+0.30 Mb; anchors wsnp_Ex_rep_c67563_66193104 and RAC875_c13610_2646).'
          ' The SIM interval was 6A_8610346-6A_15577414; the paper notes the genetic order there does not match the physical order. Same locus as the GWAS MTA and the BP025 QTL at the 6AS end.',
     evidence='Table 4: QBls.zhl-6A.1(RPIT) SIM 6A_8610346-6A_15577414; CIM 6A_1729034-6A_6608308, LOD 13.0')),
 (['V3258'], dict(start=598.1899, end=614.8722, grade='C', curation='hand',
     note='QTL for bacterial leaf streak (infection type and % water-soaked area) in the durum RP336 population; minor, significant in SIM only. Flanking marker positions are on the durum Svevo v1 assembly;'
          ' no Chinese Spring anchor is given near the 6AL end, so the Svevo positions are used as an approximation.',
     evidence='Table 4: QBls.zhl-6A.2(RPIT) 6A_598189921-6A_614872222; (RPWS) 6A_602678106-6A_614872222')),
 (['V3259'], dict(start=1.2146, end=20.2325, grade='C', curation='hand',
     note='QTL for bacterial leaf streak in the durum RP336 population; minor, significant in SIM only (RPIT 3A_12287401-3A_15345505; RPWS 3A_1214628-3A_20232548).'
          ' Flanking marker positions are on the durum Svevo v1 assembly; no Chinese Spring anchor is given near the 3AS end, so the Svevo positions are used as an approximation.',
     evidence='Table 4: QBls.zhl-3A(RPWS) 3A_1214628-3A_20232548')),
 (['V3260'], dict(start=586.605, end=610.217, ref='RefSeq v1.0', grade='B', curation='hand',
     note='QTL for bacterial leaf streak infection type in the durum BP025 population (SIM only). Interval from the Chinese Spring RefSeq v1.0 positions of the flanking 90K markers'
          ' IWA2743 (586.605 Mb, converted from v2.1 588.454 Mb, Hu et al. 2022) and IWA7044 (610.217 Mb, Sun et al. 2020). The paper gives 548-571 Mb on Svevo for these markers.',
     evidence='Table 4: QBls.zhl-5A(BPIT) IWA2743-IWA7044, 124-136 cM')),
]
SID = lambda i: 'S80-%03d' % i

def build(reports):
    byid = {r['id']: r for r in reports}
    old = {i for ids, _ in FIX for i in ids}
    new = []
    for n, (ids, f) in enumerate(FIX, 1):
        r = dict(byid[ids[0]]); r.update(f); r['id'] = SID(n); r['ref'] = f.get('ref', 'RefSeq v1.0')
        for k in ('start', 'end', 'peak'):
            if k not in f: r[k] = None
        r['note'] += ' Corrected from the paper, Nov 2026 (mn25.py).'
        r['asm'] = lift.lift_record(r)
        new.append(r)
    keep = [r for r in reports if r['id'] not in old and not r['id'].startswith('S80-')]
    return keep + new, new, sorted(old)

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], new, rem = build(M['reports'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</'))
    Sd['reports'], new2, rem2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': rem2}, open('/home/claude/work/mn25_db.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['chr'], r['start'], r['end'], r['peak'], r['grade'], r['asm'].get('RefSeq v2.1'))
    print('removed', rem2)
