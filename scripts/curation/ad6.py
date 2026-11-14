"""ad6 (Nov 2026): Martin et al. 2025 Crop Science (Fusarium crown rot, DH PAN 3471/SST 087).
Rebuilds the 9 auto records (V3423-V3431) from Table 4 and Table S3, which print Chinese Spring (IWGSC 2018 = RefSeq v1.0)
flanking-marker positions, and adds the white-head QTL and the remaining single-environment FCR QTL. Yield QTL are not added.
Grade A where the printed interval is <= 100 Mb; C where it is wider (a chromosome region, not a locus)."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift
OLD = ['V3423','V3424','V3425','V3426','V3427','V3428','V3429','V3430','V3431']
S4 = ' Table S4 lists six Axiom probe sets (AX-94431310 to AX-95179073, 7A 638.62-647.22 Mb, RefSeq v1.0) in the region shared by QFcr.cg-7A (2020), QFcr.cg-7A (2022) and QFcr.usq-7A, intended for conversion to KASP markers.'
# name, trait, year, site, LOD, PVE, parent, peak marker, (m1,chr,pos1), (m2,pos2), table, extra
Q = [
 ('QFcr.cg-3D','FCR severity',2020,'SA',3.3,5.38,'PAN 3471','5411532','5411532',495.342360,'cfd152',538.177121,'Table 4',''),
 ('QFcr.cg-7A','FCR severity',2020,'SA',5.0,8.23,'SST 087','1128642','987757',575.259106,'3532952',690.860873,'Table 4',S4),
 ('QFcr.cg-2A','FCR severity',2022,'SA',3.6,6.76,'PAN 3471','3954358','1105077',42.194565,'3034278',698.281086,'Table 4',''),
 ('QFcr.cg-5B','FCR severity',2022,'SA',4.3,7.67,'SST 087','3952912','3028104',492.093143,'1233307',545.244564,'Table 4',''),
 ('QFcr.cg-7A','FCR severity',2022,'SA',7.8,14.58,'SST 087','2261308','1083319',57.742411,'1112611',670.246544,'Table 4',S4),
 ('QFcr.usq-1D','FCR severity',2020,'AU',3.7,7.19,'SST 087','5358480','gdm126b',455.803257,'1360603',486.211404,'Table 4',''),
 ('QFcr.usq-7A','FCR severity',2020,'AU',3.3,6.07,'SST 087','1068343','1083319',57.742411,'1112611',670.246545,'Table 4',S4),
 ('QWhs.cg-1D','white heads',2019,'SA',3.7,6.87,'SST 087','5358480','gdm126b',455.803257,'1360603',486.211404,'Table 4',''),
 ('QWhs.cg-3D','white heads',2019,'SA',5.3,10.33,'SST 087','3957031','1101226',61.758699,'3957031',411.623083,'Table 4',
   ' The distal flanking marker 4734127 has no BLAST hit on 3D, so the interval is closed at the peak marker 3957031 (3D 411,623,083 bp, from the 2022 row).'),
 ('QWhs.cg-7A','white heads',2020,'SA',5.0,9.23,'SST 087','1109891','1140039',685.563919,'1229379',706.764939,'Table 4',''),
 ('QWhs.cg-3D','white heads',2022,'SA',8.5,18.17,'SST 087','3957031','3957031',411.623083,'cfd152',538.177121,'Table 4',''),
 ('QFcr.cg-2B.1','FCR severity',2019,'SA',6.2,12.26,'SST 087','1107117','1148814',50.741284,'1712708',546.870715,'Table S3',''),
 ('QFcr.cg-1A','FCR severity',2020,'SA',5.7,9.42,'SST 087','992499','1103177',585.341665,'1240234',589.475153,'Table S3',''),
 ('QFcr.cg-2B.2','FCR severity',2020,'SA',5.8,9.83,'SST 087','1146975','3934918',739.568493,'1255973',768.025530,'Table S3',''),
 ('QFcr.usq-6D','FCR severity',2020,'AU',3.1,6.26,'SST 087','1088818','4663982',83.216743,'2250869',422.940919,'Table S3',''),
 ('QFcr.usq-7D','FCR severity',2020,'AU',5.0,9.48,'SST 087','1101653','3935878',55.578270,'1240400',551.815757,'Table S3',''),
 ('QWhs.cg-2D','white heads',2020,'SA',3.9,7.08,'SST 087','3944191','1060658',86.032665,'1215137',568.313236,'Table S3',''),
 ('QWhs.cg-3B','white heads',2020,'SA',5.0,9.27,'SST 087','7940357','5970042',802.735409,'984052',829.877125,'Table S3',''),
]
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S84-')]
    new = []
    for n, (nm, tr, yr, site, lod, pve, par, pk, m1, p1, m2, p2, tab, extra) in enumerate(Q, 1):
        chr_ = re.search(r'-(\d[ABD])', nm).group(1)
        w = p2 - p1; g = 'A' if w <= 100 else 'C'
        env = 'South Africa' if site == 'SA' else 'Australia'
        note = (f'QTL for Fusarium crown rot ({tr}) in the PAN 3471 / SST 087 doubled-haploid population, {yr} field trial ({env}); '
                f'favourable allele from {par}. Interval is the Chinese Spring (IWGSC 2018, RefSeq v1.0) span of the flanking markers '
                f'{m1} and {m2} printed in {tab}' + ('; single-environment QTL.' if tab == 'Table S3' else '.') + extra)
        if g == 'C': note += f' Interval width {w:.0f} Mb: a chromosome region as the source reports it, not a locus.'
        note += ' Corrected from the paper, Nov 2026 (ad6.py).'
        r = dict(id=f'S84-{n:03d}', name=f'{nm} ({yr})', cls='QTL', disease='Fusarium crown rot', chr=chr_, start=round(p1,6), end=round(p2,6),
                 peak=None, ref='RefSeq v1.0', marker=f'{pk} (peak); {m1}-{m2} (flanking)', pve=str(pve), stat=f'LOD {lod}',
                 source='A_crown_rot_015', grade=g, curation='hand', note=note,
                 evidence=f'{tab}: {nm} {yr}, {m1} {chr_}:{p1*1e6:,.0f}; {m2} {chr_}:{p2*1e6:,.0f}; LOD {lod}; PVE {pve}%')
        r['asm'] = lift.lift_record(r); new.append(r)
    return keep + new, new, OLD
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
    ops = [{'op': 'set', 'path': 'reports/' + i, 'data': {'deleted': True}} for i in rem2]
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/home/claude/work/ad6.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['chr'], r['start'], r['end'], r['grade'], r['asm'].get('RefSeq v2.1'))
