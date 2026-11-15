"""ad30 (Nov 2026): Kwiatek et al. 2015 Cereal Res Commun 43:638 (Pch1 survey), entered from the abstract supplied by the curator
(full text not available)."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {'Kwiatek2015CRC': {'citation': 'Kwiatek M, Wisniewska H, Kaczmarek Z, Korbas M, Gawlowska M, Majka M, Pankiewicz K, Danielewicz J, Belter J (2015) Using markers and field evaluation to identify the source of eyespot resistance gene Pch1 in the collection of wheat breeding lines. Cereal Res Commun 43:638-648',
                           'url': 'https://doi.org/10.1556/0806.43.2015.024', 'short': 'Kwiatek et al. (2015)'}}
def build(reports):
    keep = [r for r in reports if not r['id'].startswith('S141-')]
    r = dict(id='S141-001', cls='DG', ref='RefSeq v1.0', start=631.678, end=635.512, curation='hand', grade='C', name='Pch1 (Polish breeding lines)', disease='Eyespot', chr='7D',
             source='Kwiatek2015CRC', marker='EpD1b (endopeptidase), XustSSR2001-7DL',
             note='372 wheat genotypes from breeding programmes screened with the endopeptidase marker EpD1b and the STS marker XustSSR2001-7DL and inoculated in field tests over three years (2011-2013, K-index): both markers identified the same 18 genotypes as carrying Pch1 (from Aegilops ventricosa), and these were resistant to eyespot, validating the markers. Survey of a known gene; placed at the Pch1 interval in the atlas. Entered from the abstract (full text not available), Nov 2026 (ad30.py).',
             evidence='Abstract: 18 genotypes with Pch1 by both markers')
    r['asm'] = lift.lift_record(r); return keep + [r], [r]
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new = build(M['reports']); M['sources'].update(SRC); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], new2 = build(Sd['reports'])
    if 'sources' in Sd: Sd['sources'].update(SRC)
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    print(json.dumps(new2[0]['asm']))
