"""ad32 (Nov 2026): Liu et al. 2008 Cereal Res Commun 36 (Suppl B):195 (towards positional cloning of Fhb1), entered from the abstract
supplied by the curator (full text not available)."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {'Liu2008CRC': {'citation': 'Liu S, Pumphrey M, Gill B, Trick H, Zhang J, Dolezel J, Chalhoub B, Anderson J (2008) Toward positional cloning of Fhb1, a major QTL for Fusarium head blight resistance in wheat. Cereal Res Commun 36 (Suppl B):195-201',
                       'url': 'https://doi.org/10.1556/CRC.36.2008.Suppl.B.15', 'short': 'Liu et al. (2008)'}}
def build(reports):
    keep = [r for r in reports if not r['id'].startswith('S145-')]
    r = dict(id='S145-001', cls='QTL', ref='RefSeq v1.0', peak=8.528, curation='hand', grade='C', name='Fhb1 (261 kb BAC interval; UMN10)', disease='Fusarium head blight', chr='3B',
             source='Liu2008CRC', marker='UMN10; BAC-derived markers',
             note='A Chinese Spring 3B BAC contig spanning Fhb1 (3BS) was built and sequenced; new markers narrowed Fhb1 to a 261 kb region with seven putative genes, four expressed in spikes but none differentially expressed after Fusarium inoculation. Sumai 3 cosmids covering all seven genes were transformed into Bobwhite; transgenics for four of five cosmids showed no type II resistance (the fifth was still under test). The codominant diagnostic marker UMN10 was developed for MAS. Placed at the cloned Fhb1 (TaHRC) position in the atlas. Entered from the abstract (full text not available), Nov 2026 (ad32.py).',
             evidence='Abstract: Fhb1 in a 261 kb interval, seven candidate genes; UMN10')
    r['asm'] = lift.lift_record(r); return keep + [r], [r]
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new = build(M['reports']); M['sources'].update(SRC); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], new2 = build(Sd['reports'])
    if 'sources' in Sd: Sd['sources'].update(SRC)
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    print(new2[0]['asm']['RefSeq v2.1'])
