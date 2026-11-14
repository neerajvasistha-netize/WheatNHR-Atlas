"""ad14 (Nov 2026): Alahmad et al. 2020 IJMS 21:5260 (durum, 168 genotypes, DArTseq GWAS, crown rot field nursery).
V1627 (auto 'QTL_6B', unplaced) replaced by qCR-6B (with the co-located whitehead QTL qWH1-6B and qWH2-6B). The paper gives only
Svevo consensus-map cM; the DArTseq tags placed by BLAST on Chinese Spring RefSeq v1.0 give the interval. Non-disease QTL not added."""
import json, re, sys
sys.path[:0] = ['scripts/curation', 'scripts/curation/sb19']
import lift, stb44
P = {m: stb44.P(m, '6B') for m in ['4992135', '1023342', '2254131', '2256390', '1039837']}
pl = {m: v[0] for m, v in P.items() if v[0]}
REC = dict(id='S91-001', name='qCR-6B', cls='QTL', disease='Fusarium crown rot', chr='6B', start=min(pl.values()), end=max(pl.values()), peak=None,
           ref='RefSeq v1.0', marker='1023342 (peak), 4992135, 2254131, 2256390, 1039837 (DArTseq)', pve='', stat='-log10(P) 4.20 (CR severity)',
           source='V_Alahmad2020', grade='B', curation='hand',
           note=('Crown rot QTL in 168 durum wheat genotypes (DArTseq GWAS), field crown rot nursery at Warwick 2017; CR severity index from whitehead '
                 'percentage and stem browning: five SNPs at 100.13-101.26 cM on the Svevo consensus map, peak 1023342 (-log10(P) 4.20). The same '
                 'markers gave the whitehead QTL qWH1-6B (-log10(P) up to 3.88) and qWH2-6B (up to 3.58), and the yield QTL qGY-6B in the CR '
                 'nursery. Haplotype hap1 gave 0.79 t/ha more yield under crown rot, a 4.24 d longer grain-filling period and higher root biomass; '
                 'combined with qSRA-6A (root angle) the yield advantage was 0.90 t/ha. Positions are not printed in the paper; the interval is the '
                 'RefSeq v1.0 span of the DArTseq tags placed by BLAST (' + ', '.join(f'{m} {v:.3f} Mb' for m, v in pl.items()) + '); 1023342 and '
                 '1039837 could not be placed. Corrected from the paper, Nov 2026 (ad14.py).'),
           evidence='Table 2: qCR-6B 4992135 100.13 cM, 1023342 100.50, 2254131 100.50, 2256390 101.26, 1039837 101.26 (Svevo map); qWH1-6B, qWH2-6B same markers')
def build(reports):
    keep = [r for r in reports if r['id'] not in ('V1627', 'S91-001')]
    r = dict(REC); r['asm'] = lift.lift_record(r)
    return keep + [r], r
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], r = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], r2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    d = {k: v for k, v in r2.items() if v not in (None, '') and k != 'asm'}
    d['asm'] = {k: r2['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r2['asm']}
    print(json.dumps(d, ensure_ascii=False))
