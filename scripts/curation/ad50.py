"""ad50 (Nov 2026): Virdi et al. 2016 G3 Table S1 (AL linkage map, 844 markers). Checked the AL tan spot QTL against the map:
S31-004 Tsc2 lies between Xiwa2304 (8.1 cM, 16.90 Mb) and Xwmc154 (35.2 cM, 36.45 Mb), consistent with the atlas Tsc2 position (note only).
S31-005 4B interval 7.1-38.9 cM spans a map gap (9.7-44.2 cM, no markers); its flanks lie 242 Mb apart, so the grade is corrected from B to C."""
import json, re
U = {'S31-004': dict(note='Altar 84 x Langdon durum RIL (AL, 127 lines); race 5 isolate DW5; peak defined by the Tsc2 locus on 2BS; resistance from Langdon. Flanks Xhbg216 and Xwmc25 could not be placed, so the locus is set at the Tsc2 position already in the atlas (Kokhmetova et al. 2021, 24.10 Mb). Read from the paper, Nov 2026 (sb31.py). In the Table S1 map Tsc2 (18.5 cM) lies between Xiwa2304 (8.1 cM, 16.896 Mb, Sun et al. 2020 PBJ (90K)) and Xwmc154 (35.2 cM, 36.447 Mb, Zhao et al. 2019 J Genet), consistent with that position (ad50.py).'),
     'S31-005': dict(grade='C', note='Altar 84 x Langdon durum RIL (AL, 127 lines); race 2 isolate L13-35 and race 5 isolate DW5; resistance from Altar 84. Flanking markers lie 242 Mb apart; interval is wide. Marker positions (RefSeq v1.0): Xiwa2126 14.118 Mb (Hu et al. 2022 BMC Plant Biol (90K on v2.1) (v2.1 15.856 Mb, converted)); Xgwm113 256.314 Mb (in-silico PCR). Read from the paper, Nov 2026 (sb31.py). In the Table S1 map the 4B group has no markers between 9.7 cM (Xiwa4569, 17.258 Mb) and 44.2 cM (Xbarc1045), so the QTL cannot be narrowed; grade corrected to C for the 242 Mb span (ad50.py).')}
def build(rs):
    out = []
    for r in rs:
        if r['id'] in U: r = dict(r); r.update(U[r['id']])
        out.append(r)
    return out
M = json.load(open('data/merged_v13f.json')); M['reports'] = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'] = build(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
