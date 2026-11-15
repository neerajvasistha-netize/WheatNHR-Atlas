"""ad44 (Nov 2026): Kokhmetova et al. 2021 Front Genet 11:581214, tan spot GWAS in a Kazakhstan hexaploid wheat collection.
Supplementary Table S5 names the 11 race-specific MTA SNPs (DArTseq IDs with chromosome) used for the best allele combinations
(race 1: 8 SNPs; race 5: 3 SNPs). The main paper was not uploaded, so P values and effects are not entered; positions by DArT tag BLAST."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = ' Read from Supplementary Table S5, Nov 2026 (ad44.py).'
BASE = 'Tan spot (P. tritici-repentis) GWAS in a hexaploid wheat collection from Kazakhstan (DArTseq); Supplementary Table S5 lists the SNPs whose favourable alleles were combined in the most resistant Kazakh and Russian entries. P values and effect sizes are in the main paper (not uploaded). '
L = [('Race 1', '1B', ['1406319']), ('Race 1', '3B', ['1147153']), ('Race 1', '4B', ['3958510']), ('Race 1', '5A', ['1134493']), ('Race 1', '5B', ['3946488']),
     ('Race 1', '6A', ['1004240']), ('Race 1', '7B', ['3021234', '1081730']), ('Race 5', '5A', ['1008802']), ('Race 5', '5D', ['1093048']), ('Race 5', '6A', ['1862737'])]
def build(reports):
    keep = [r for r in reports if not r['id'].startswith('S197-')]; N = []
    for n, (race, ch, ms) in enumerate(L, 1):
        pl = [(m,) + tuple(resolve(m, ch)) for m in ms]; ok = [x for x in pl if x[1]]
        r = dict(id=f'S197-{n:03d}', cls='MTA', ref='RefSeq v1.0', curation='hand', disease='Tan spot', chr=ch, source='Kokhmetova2021', marker=', '.join(f'{ch}_{m}' for m in ms),
                 name=f'Tan spot {race.lower()} MTA {ch} (Kazakhstan collection)', evidence='Supplementary Table S5')
        if len(ok) == 2:
            r['start'], r['end'], r['grade'] = min(x[1] for x in ok), max(x[1] for x in ok), 'B'
            t = 'Placed in silico: ' + ', '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in ok) + '.'
        elif ok:
            m, p, h = ok[0]; r['peak'], r['grade'] = p, 'C'; t = f'Placed in silico at {p:.3f} Mb ({h}).'
        else:
            r['grade'] = 'D'; t = 'The DArTseq tag could not be placed.'
        r['note'] = BASE + f'{race} MTA on {ch}: ' + ', '.join(f'{ch}_{m}' for m in ms) + '. ' + t + T
        N.append(r)
    for r in N: r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
    return keep + N, N
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
    open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
    ops = []
    for r in new:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/tmp/ops44.json', 'w'), ensure_ascii=False)
    print(json.dumps(ops, ensure_ascii=False))
