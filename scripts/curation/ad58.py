"""ad58 (Nov 2026): Singh et al. (2023) Front Plant Sci 13:1036064 (Singh2023FPS), spot blotch GWAS, re-curated from the paper.
The 56 earlier records were one per SNP with generic names ("MTA M2266") and no trait, environment, MAF or recommendation.
Now one record per SNP and trait (Table 2; 70 chromosome-assigned MTAs of the 89 detected by four or five models), with
environments, P, MAF, the Table 4 'most important' flag and the five AUDPC MTAs recommended for MAS, plus the eight
epistatic pairs of Table 3 (16 SNPs; the two AUDPC pairs are the other two recommended MTAs). Positions: SNP tag start, IWGSC RefSeq v1.0."""
import json, re, sys, collections
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T2 = '''AUDPC|E1|M2266|1A|445269309|8.37E-04|0.445545
AUDPC|E1|M1654|2B|797177741|2.36E-04|0.237624
AUDPC|E1|M3799|3B|685318543|3.68E-04|0.065677
AUDPC|E1|M4359|6B|22081357|6.66E-04|0.356304
AUDPC|E1|M5103|2A|73179815|8.98E-04|0.227723
AUDPC|E1|M6031|5A|672843173|5.36E-05|0.095215
AUDPC|E1|M9772|3A|659038152|2.32E-04|0.285462
AUDPC|E1|M10198|7A|35630453|6.17E-04|0.057574
AUDPC|E1|M10876|4A|709908601|7.63E-04|0.235066
AUDPC|BLUP|M763|1A|14161093|4.38E-04|0.462046
AUDPC|E2|M5576|3A|23840923|5.63E-04|0.206254
AUDPC|BLUP|M3226|6B|698461241|5.85E-04|0.321634
AUDPC|E2|M6146|5B|598120554|3.61E-04|0.360264
AUDPC|BLUP|M4359|6B|22081357|4.66E-04|0.356304
AUDPC|E2|M8330|2A|700374249|9.51E-04|0.155116
AUDPC|BLUP|M5019|7B|677430865|4.18E-04|0.326337
AUDPC|E2|M10783|7D|490344136|2.01E-04|0.223531
AUDPC|BLUP|M7745|2D|35039085|2.89E-04|0.431749
AUDPC|E3|M1654|2B|797177741|4.09E-04|0.237624
AUDPC|BLUP|M9772|3A|659038152|1.25E-04|0.285462
AUDPC|BLUP|M10783|7D|490344136|1.82E-04|0.223531
AUDPC|E3|M9772|3A|659038152|1.16E-04|0.285462
AUDPC|BLUP|M11095|7A|700115446|7.88E-04|0.079076
IP|E1|M140|2B|91205274|4.95E-04|0.158416
IP|E1|M650|2B|79317077|4.37E-04|0.19802
IP|E3|M10907|7B|730875968|1.24E-04|0.278696
IP|E1|M11792|4B|647704512|2.46E-04|0.348878
IP|E4|M4648|7A|94143005|2.78E-04|0.285132
IP|E4|M5996|3D|479617743|9.84E-04|0.370875
IP|E2|M11418|6A|443334314|8.30E-04|0.203663
IP|E4|M7433|3B|60372653|3.22E-04|0.17495
IP|E4|M10241|3B|11646443|2.24E-04|0.051139
IP|E3|M119|5B|529608885|5.53E-04|0.064274
IP|E4|M11765|2A|107180183|3.43E-04|0.155116
IP|E3|M225|2A|143209405|8.46E-04|0.432343
IP|E3|M804|5A|582958870|3.90E-04|0.054406
IP|E3|M3107|7B|196255937|5.85E-04|0.05264
IP|BLUP|M876|6B|27649755|8.16E-04|0.498548
IP|E3|M3962|5B|228277790|6.15E-04|0.318399
IP|BLUP|M2039|2A|81660284|1.32E-04|0.062706
IP|E3|M4648|7A|94143005|5.89E-05|0.285132
IP|BLUP|M2947|1B|170842706|7.09E-04|0.395776
IP|BLUP|M4089|6A|593657698|1.81E-04|0.159719
IP|E3|M5551|1D|452208385|8.56E-04|0.256997
IP|BLUP|M9665|3B|746854760|9.63E-04|0.154422
IP|E3|M8110|7B|553167382|7.43E-04|0.074868
IP|BLUP|M11765|2A|107180183|3.43E-04|0.155116
LN|E1|M2326|5A|32743559|5.00E-04|0.235924
LN|E3|M6140|5B|522912127|3.24E-04|0.221023
LN|E1|M3426|2B|73990358|1.85E-04|0.480182
LN|E3|M6730|1A|20893697|5.11E-04|0.094604
LN|E1|M5696|5D|382494021|6.70E-04|0.287228
LN|E1|M8024|6A|585062562|8.90E-04|0.168944
LN|E3|M11585|7A|66298367|5.01E-04|0.17066
LN|E1|M8443|4D|134916786|3.98E-06|0.124752
LN|E3|M4854|4A|200817628|8.44E-04|0.186056
LN|E1|M2010|3B|598995436|9.93E-04|0.359736
LN|E2|M164|4B|538706787|5.74E-04|0.067607
LN|E4|M2771|4A|16997957|1.18E-04|0.308531
LN|E4|M3928|2B|602182720|2.83E-04|0.438614
LN|E2|M3416|6A|596590906|6.52E-04|0.409373
LN|E4|M5313|2A|485122901|6.21E-04|0.135314
LN|E2|M9240|6D|467151830|2.54E-04|0.192327
LN|E4|M6140|5B|522912127|5.34E-04|0.221023
LN|E3|M2771|4A|16997957|6.62E-04|0.308531
LN|E4|M6730|1A|20893697|3.99E-04|0.094604
LN|E3|M3928|2B|602182720|2.12E-04|0.438614
LN|E4|M11585|7A|66298367|5.37E-04|0.17066
LN|E3|M5313|2A|485122901|4.94E-04|0.135314
LN|BLUP|M5928|6A|33188097|5.94E-04|0.250957'''
# Table 2 reads down the left column and then continues in the right column; environment assignments (ENVS) follow that order,
# match the per-environment counts in the trait headings, and agree with Table 4.
T4 = {('AUDPC', 'M1654'): 'E3, E4', ('AUDPC', 'M4359'): 'E4, BLUP', ('AUDPC', 'M9772'): 'E3, E4, BLUP', ('AUDPC', 'M10783'): 'E2, BLUP', ('AUDPC', 'M8330'): 'E2',
      ('IP', 'M4648'): 'E3, E4', ('IP', 'M11765'): 'E4, BLUP', ('LN', 'M2771'): 'E3, E4', ('LN', 'M3928'): 'E3, E4', ('LN', 'M5313'): 'E3, E4',
      ('LN', 'M6140'): 'E3, E4', ('LN', 'M6730'): 'E3, E4', ('LN', 'M11585'): 'E3, E4'}
T4P = {'M1654': '4.09E-04', 'M4359': '6.66E-04', 'M9772': '1.16E-04', 'M10783': '2.01E-04', 'M8330': '9.51E-04', 'M4648': '5.89E-05', 'M11765': '3.55E-04',
       'M2771': '1.18E-04', 'M3928': '2.12E-04', 'M5313': '4.94E-04', 'M6140': '5.34E-04', 'M6730': '5.11E-04', 'M11585': '5.01E-04'}
TR = {'AUDPC': 'AUDPC', 'IP': 'incubation period', 'LN': 'lesion number'}
ENV = 'E1 BHU Varanasi 2017-18, E2 BHU 2018-19, E3 BISA Pusa 2017-18, E4 BISA Pusa 2018-19'
BASE = '303 CIMMYT spring wheat genotypes, DArT-seq (12,196 SNPs, 5,400 unassigned), four environments (' + ENV + ') and BLUP; MTAs detected by four of five GWAS models (CMLM, MLMM, FarmCPU, BLINK; not SUPER), P < 0.001. Table 2 heads its P column "-log10(P)" but the values are P. '
ENVS = {'AUDPC': {'E1': 'M2266 M3799 M5103 M6031 M10876', 'E2': 'M5576 M6146 M8330 M10783', 'E3': 'M1654 M9772', 'E4': 'M1654 M4359 M9772 M10198', 'BLUP': 'M763 M3226 M4359 M5019 M7745 M9772 M10783 M11095'}, 'IP': {'E1': 'M140 M650 M11792', 'E2': 'M11418', 'E3': 'M119 M225 M804 M3107 M3962 M4648 M5551 M10907 M8110', 'E4': 'M4648 M5996 M7433 M10241 M11765', 'BLUP': 'M876 M2039 M2947 M4089 M9665 M11765'}, 'LN': {'E1': 'M2326 M3426 M5696 M8024 M8443', 'E2': 'M164 M3416 M9240', 'E3': 'M2771 M3928 M5313 M6140 M6730 M11585 M4854 M2010', 'E4': 'M2771 M3928 M5313 M6140 M6730 M11585', 'BLUP': 'M5928'}}
g = collections.OrderedDict()
for l in T2.splitlines():
    tr, e, m, ch, bp, p, maf = l.split('|')
    k = (tr, m); d = g.setdefault(k, dict(tr=tr, m=m, ch=ch, bp=int(bp), env=[x for x in ('E1', 'E2', 'E3', 'E4', 'BLUP') if m in ENVS[tr][x].split()], p=[], maf=maf)); d['p'].append(float(p))
EPI = [('AUDPC', 'M238', '2B', 67247748, 'M1164', '7A', 94495522, '1.27E-03'), ('AUDPC', 'M4326', '7A', 640669196, 'M11844', '3B', 581263393, '1.07E-03'),
       ('IP', 'M3711', '1B', 499899181, 'M4228', '7B', 572536722, '4.02E-03'), ('LN', 'M921', '3A', 509934375, 'M8709', '6B', 48479624, '1.57E-03'),
       ('LN', 'M1998', '3A', 12868494, 'M8278', '6A', 606427410, '7.26E-03'), ('LN', 'M7467', '3D', 556059162, 'M12122', '7D', 3854256, '9.14E-03'),
       ('LN', 'M9322', '2A', 16016622, 'M3401', '7A', 27478985, '7.55E-03'), ('LN', 'M9496', '3A', 17528445, 'M9267', '3B', 683764839, '2.00E-03')]
TAIL = ' Read from the paper, Nov 2026 (ad58.py).'
def records(old):
    pool = collections.defaultdict(list)
    for r in old:
        if r['source'] == 'Singh2023FPS': pool[r['marker']].append(r['id'])
    used = set(); out = []; n = 0
    def nid(m):
        nonlocal n
        for i in pool.get(m, []):
            if i not in used: used.add(i); return i
        n += 1; return f'S205-{n:03d}'
    for (tr, m), d in g.items():
        imp = (tr, m) in T4; rec = tr == 'AUDPC' and imp and m != 'M8330' or m == 'M8330'
        P = float(T4P[m]) if imp else min(d['p'])
        nt = BASE + f'Trait {TR[tr]}; significant in {", ".join(d["env"])} ; P {P:.2e}; MAF {float(d["maf"]):.3f}.'
        if imp: nt += f' Listed among the most important MTAs (Table 4).'
        if tr == 'AUDPC' and imp: nt += ' One of the five main-effect AUDPC MTAs recommended for MAS after validation.'
        if m == 'M8330': nt += ' Lies within the QSb.bhu-2A interval (Kumar et al. 2010).'
        r = dict(id=nid(m), cls='MTA', curation='hand', name=f'{m} ({tr})', disease='Spot blotch', chr=d['ch'], peak=round(d['bp'] / 1e6, 6), start=None, end=None,
                 ref='RefSeq v1.0', marker=m, stat=f'P {P:.2e}', pve='', source='Singh2023FPS', grade='A',
                 evidence=f'Table 2: {m} {d["ch"]} {d["bp"]} bp' + ('; Table 4' if imp else ''), note=nt + TAIL)
        out.append(r)
    for tr, m1, c1, p1, m2, c2, p2, p in EPI:
        for a, ca, pa, b, cb, pb in ((m1, c1, p1, m2, c2, p2), (m2, c2, p2, m1, c1, p1)):
            nt = BASE.split('; MTAs detected')[0] + f'. First-order epistatic interaction (PLINK) for {TR[tr]}: {a} ({ca}, {pa} bp) x {b} ({cb}, {pb} bp), P {p} (Table 3).'
            if tr == 'AUDPC': nt += ' One of the two AUDPC epistatic pairs recommended for MAS.'
            out.append(dict(id=nid(a), cls='MTA', curation='hand', name=f'{a} ({tr}, epistatic with {b})', disease='Spot blotch', chr=ca, peak=round(pa / 1e6, 6), start=None, end=None,
                            ref='RefSeq v1.0', marker=f'{a}, {b}', stat=f'P {p} (epistasis)', pve='', source='Singh2023FPS', grade='A', evidence=f'Table 3: {a} {ca}:{pa}', note=nt + TAIL))
    drop = [i for v in pool.values() for i in v if i not in used]
    for r in out: r['asm'] = lift.lift_record(r)
    return out, drop
def apply(rs):
    new, drop = records(rs); ids = {r['id'] for r in new}
    keep = [r for r in rs if r['source'] != 'Singh2023FPS' and r['id'] not in ids]
    return keep + new, new, drop
if __name__ == '__main__':
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    M = json.load(open('data/merged_v13f.json')); M['reports'], new, drop = apply(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _, _ = apply(Sd['reports'])
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    json.dump({'reports': new, 'sources': {}, 'deleted': drop}, open('web/seed_ad58.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    print(len(new), 'records', len(drop), 'dropped', drop, sum(r['id'].startswith('S205') for r in new), 'new ids')
