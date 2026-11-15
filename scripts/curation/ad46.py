"""ad46 (Nov 2026): (1) Pankaj et al. 2022 Australas Plant Pathol ESM1 Table S3 primer sequences re-run by exact-match in-silico PCR
(sb19/ispcr.py on RefSeq v1.0; hits in docs/supplements/pankaj2022_ispcr_hits.json): wmc173 places on 5B, cfd59 on 1B.
(2) Rahman et al. 2021 Field Crops Res Supplementary Table S3 names the QTL_8 SNP (Excalibur_c65830_82)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
PK = ' Placed by exact-match in-silico PCR of the primer sequences in Supplementary Table S3 (ESM1), Nov 2026 (ad46.py).'
p8, h8 = resolve('Excalibur_c65830_82', '2B')
U = {
 'S19-059': dict(peak=119.1647, grade='C', note='PBW343 x IC252874 DH population; 2017-18 and 2018-19; resistance allele from PBW343. Only wmc173 places on 5B (119.165 Mb); the Cfd71 primers amplify only on 4D. Position is that marker, not a reported peak.' + PK, evidence='Table 3: 2017-18 Qsb_rpcau_5B Cfd71-wmc173 LOD 11.16 PVE 18.75; Table S3 (ESM1)'),
 'S19-068': dict(peak=205.4211, grade='C', note='PBW343 x IC252874 DH population; 2018-19; resistance allele from PBW343. Only cfd59 places on 1B (205.421 Mb; 1,464 bp product, so the placement is tentative); the wmc611 primers give no product. Position is that marker, not a reported peak.' + PK, evidence='Table 3: 2018-19 Qsb_rpcau_1B wmc611-cfd59 LOD 6.71 PVE 12.72; Table S3 (ESM1)'),
 'S157-002': dict(peak=p8, grade='C', marker='Excalibur_c65830_82', note='Crown rot (F. pseudograminearum) GWAS with the 90K SNP array in 196 F4/F5 progeny of the three-way cross AUS29529/2/Syn110/3/Sunco/2*Pastor (207 genotypes with parents and checks), paired inoculated and uninoculated field plots at Narrabri in 2014 and 2015; 43 QTL, of which three for resistance (white stem browning) and 40 for tolerance (yield traits under inoculation; not entered). QTL_8 reduced stem browning (R2 6.2%%, PopSeq/consensus position 317.99); Supplementary Table S3 names its SNP, Excalibur_c65830_82, placed at %.3f Mb (%s). Read from the paper and Supplementary Table S3, Nov 2026 (ad46.py).' % (p8, h8), evidence='Table 6; Table S3'),
}
def build(reports):
    out, upd = [], []
    for r in reports:
        if r['id'] in U:
            r = dict(r); r.update(U[r['id']]); r['asm'] = lift.lift_record(r); upd.append(r)
        out.append(r)
    return out, upd
M = json.load(open('data/merged_v13f.json')); M['reports'], upd = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
ops = []
for r in upd:
    d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
    d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
    ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'ifVersion': 1, 'data': d}); print(r['id'], r['peak'], d['asm'], d['source'])
json.dump(ops, open('/tmp/ops46.json', 'w'), ensure_ascii=False)
