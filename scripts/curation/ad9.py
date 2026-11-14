"""ad9 (Nov 2026): R158 Ra_c9427_300 (7A, adult-plant FCR) was filed under source key Rahman2020FCR, whose citation is
Lin et al. 2023. New source Rahman2020MB with the reference as given in Lin et al. 2023 (Mol Breed 40:28).
The 82.09 Mb quoted by Lin et al. sits on their RefSeq v2.0 frame; the marker's RefSeq v1.0 position (Sun et al. 2020, 90K) is used."""
import json, re, sys
sys.path[:0] = ['scripts/curation', 'scripts/curation/sb19']
import lift, stb44
SRC = {'Rahman2020MB': {'citation': 'Rahman M, Davies P, Bansal U, Pasam R, Hayden M, Trethowan R (2020) Marker-assisted recurrent selection improves the crown rot resistance of bread wheat. Molecular Breeding 40:28',
                        'url': 'https://doi.org/10.1007/s11032-020-1105-1', 'short': 'Rahman et al. (2020)'}}
pos, how = stb44.P('Ra_c9427_300', '7A')
def build(reports):
    keep = [r for r in reports if r['id'] not in ('R158', 'S87-001')]
    r = dict(id='S87-001', name='Ra_c9427_300', cls='MTA', disease='Fusarium crown rot', chr='7A', start=None, end=None, peak=pos,
             ref='RefSeq v1.0', marker='Ra_c9427_300', pve='', stat='', source='Rahman2020MB', grade='C', curation='hand',
             note=('SNP significantly associated with adult-plant Fusarium crown rot resistance in a GWAS by Rahman et al. (2020), as cited by '
                   'Lin et al. (2023, Plant Disease 107:1151): "Ra_c9427_300 on chromosome 7A at 82.09 Mb" (the RefSeq v2.0 frame of that paper). '
                   f'Single marker; position is its RefSeq v1.0 position ({how}). Original paper not yet checked. Corrected Nov 2026 (ad9.py).'),
             evidence='Lin et al. 2023 Discussion: Ra_c9427_300, 7A 82.09 Mb, adult FCR resistance (Rahman et al. 2020)')
    r['asm'] = lift.lift_record(r)
    return keep + [r], r
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], r = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], r2 = build(Sd['reports']); Sd['sources'].update(SRC)
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    d = {k: v for k, v in r2.items() if v not in (None, '') and k != 'asm'}
    d['asm'] = {k: r2['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r2['asm']}
    print(json.dumps(d))
