"""ad13 (Nov 2026): Liu et al. 2011 Czech J Genet Plant Breed 47:S109 (summary of CSIRO crown rot QTL work).
V1585 (auto 'QTL_3B', unplaced) replaced by Qcrs.cpi-3B as reported here. The Ernie 3BL QTL restated in this paper is the
primary result of Li et al. 2010 (S88-001) and is not duplicated; the barley 3HL QTL is outside the wheat atlas."""
import json, re, sys
sys.path[:0] = ['scripts/curation', 'scripts/curation/sb19']
import lift, stb44
pos, how = stb44.P('gwm0181', '3B')
REC = dict(id='S90-001', name='Qcrs.cpi-3B (CSCR6)', cls='QTL', disease='Fusarium crown rot', chr='3B', start=None, end=None, peak=pos,
           ref='RefSeq v1.0', marker='wPt10505, wPt2277 (flanking); gwm0181 (closest SSR, about 2.2 cM)', pve='up to 49 (43 with Rht1 as cofactor)',
           stat='LOD 10.8', source='V_Liu2011', grade='C', curation='hand',
           note=('Major crown rot QTL on 3BL from the T. spelta accession CSCR6, flanked by DArT markers wPt10505 and wPt2277; detected with both '
                 'F. pseudograminearum and F. graminearum isolates; up to 49% of phenotypic variance (LOD 10.8, interval mapping). Not co-located '
                 'with known plant-height genes, but a strong interaction with Rht1 (4BS) was found; with Rht1 as cofactor the maximum effect fell to '
                 '43%. In four validation populations the QTL alone reduced crown rot severity by 29-42% (mean 33%). The flanking DArT markers have '
                 f'no physical position; the locus is placed at the RefSeq v1.0 position of gwm0181 ({how}), about 2.2 cM away. '
                 'The same locus was later fine-mapped (see Qcrs.cpi-3B, Zheng et al.). Corrected from the paper, Nov 2026 (ad13.py).'),
           evidence='Results: Qcrs.cpi-3B flanked by wPt10505 and wPt2277; up to 49% (LOD 10.8); gwm0181 about 2.2 cM; validation 29-42% (mean 33%)')
def build(reports):
    keep = [r for r in reports if r['id'] not in ('V1585', 'S90-001')]
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
