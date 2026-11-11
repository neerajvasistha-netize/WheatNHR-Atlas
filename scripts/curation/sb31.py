"""Batch 3 (Nov 2026): five more SNB / tan spot sources.
  Cockram et al. 2015 (G3 5:2257; file Cockram2010.pdf), Reddy et al. 2008 (Plant Genome 1:55),
  Virdi et al. 2016 (G3 6:4139), Running 2022 PhD dissertation (NDSU; front matter and chapter 2 only),
  Saini et al. 2022 (Mol Breed 42:11, meta-QTL: not added, the atlas does not keep meta-QTL).
Uses the marker resolver of sb26_fix.py (genes, in-silico PCR incl. Roder 1998 GWM primers, marker tables).
Run from the repository root after sb26_apply.py and r98_apply.py:  python scripts/curation/sb31.py"""
import json, re, sys, io, contextlib
sys.path.insert(0, 'scripts/curation'); sys.path.insert(0, 'scripts/curation/sb19')
with contextlib.redirect_stdout(io.StringIO()):
    import sb26_fix as F
import lift

R = []
def add(**k):
    for f, v in (('ref', 'RefSeq v1.0'), ('curation', 'hand'), ('cls', 'QTL'), ('start', None), ('end', None), ('peak', None),
                 ('pve', ''), ('stat', ''), ('note', ''), ('grade', None), ('marker', ''), ('disease', 'Septoria nodorum blotch')):
        k.setdefault(f, v)
    R.append(k); return k

def gene(rec, g, how):
    c, s, e = F.GENE[g]
    rec['start'], rec['end'], rec['grade'] = s, e, 'B'
    rec['note'] = (rec['note'] + ' ' + how).strip()
    return rec

# ---------- Cockram et al. 2015 G3, eight-founder NIAB MAGIC, SnTox1 infiltration, Table 1 ----------
r = add(name='Snn1 (QSnn.niab-1B; Cockram 2015)', chr='1B', cls='SG', marker='Excalibur_c21898_1423 (peak); BS00026180_51, IAAV5782',
        pve='25.9', stat='-log10P 55.29', source='Cockram2015G3',
        note='NIAB 8-founder MAGIC, 642 lines, SnTox1 sensitivity; peak SNPs Excalibur_c21898_1423 and BS00093078_51 cosegregate at 8.36 cM. '
             'Interval spans the most significant mapped 1BS SNPs; the cloned Snn1 gene (2.36 Mb) lies 0.26 Mb proximal to its end.',
        evidence='Table 1: Excalibur_c21898_1423 -log10P 55.29, 1B 8.361 cM; text: peak explains 25.9%')
F.place(r, ['BS00026180_51', 'IAAV5782'], '1B', peak_marker='Excalibur_c21898_1423')
r = add(name='QSnn.niab-5A.1', chr='5A', marker='BS00068108_51 (peak); Excalibur_c42255_425, Excalibur_c77910_272, RAC875_c32639_395, BS00022299_51',
        pve='4.1', stat='-log10P 7.25', source='Cockram2015G3',
        note='NIAB 8-founder MAGIC, SnTox1 sensitivity; new SnTox1 sensitivity locus on 5AL (11 significant SNPs at 124.6-144.3 cM).',
        evidence='Table 1: BS00068108_51 -log10P 7.25, 5A 144.263 cM; text: 4.1% of variation')
F.place(r, ['RAC875_c32639_395', 'BS00022299_51', 'Excalibur_c42255_425', 'Excalibur_c77910_272'], '5A', peak_marker='BS00068108_51')

# ---------- Reddy et al. 2008 Plant Genome, CS x CS-DIC 1B, 8510 gametes ----------
gene(add(name='Snn1 (Reddy 2008 high-resolution map)', chr='1B', cls='SG', marker='Xksud14.2–XBF293222 (0.46 cM)', source='Reddy2008PG',
         note='Snn1 delimited to a 0.46 cM interval (8510 gametes); SSR flanks Xfcp618 (0.74 cM distal) and Xpsp3000 (0.16 cM proximal, '
              'in-silico PCR 5.14 Mb). Candidate NBS-LRR ESTs BE498831 and BF474204 cosegregate.',
         evidence='Abstract and Fig. 4: Snn1 between Xksud14.2 and XBF293222, 0.46 cM'),
     'Snn1', 'Placed at the cloned Snn1 gene (RefSeq v1.0).')

# ---------- Virdi et al. 2016 G3, Altar 84 x Langdon RIL (AL) and Langdon x LDN-DIC 5B (LD5B) ----------
AL = 'Altar 84 x Langdon durum RIL (AL, 127 lines)'
r = add(name='QTs.fcu-2B (Tsc2; Virdi 2016)', chr='2B', disease='Tan spot', marker='Tsc2 (Xhbg216–Xwmc25)', pve='26', stat='LOD 12.0',
        source='Virdi2016G3', grade='C', peak=24.1012,
        note=AL + '; race 5 isolate DW5; peak defined by the Tsc2 locus on 2BS; resistance from Langdon. Flanks Xhbg216 and Xwmc25 could not be '
             'placed, so the locus is set at the Tsc2 position already in the atlas (Kokhmetova et al. 2021, 24.10 Mb).',
        evidence='Table 3: QTs.fcu-2B Tsc2 18.7 cM LOD 12.00 R2 0.26 (DW5)')
r = add(name='QTs.fcu-4B (Virdi 2016)', chr='4B', disease='Tan spot', marker='Xiwa2126–Xgwm113', pve='11-12', stat='LOD 4.0 (L13-35), 8.6 (DW5)',
        source='Virdi2016G3', note=AL + '; race 2 isolate L13-35 and race 5 isolate DW5; resistance from Altar 84.',
        evidence='Table 3: QTs.fcu-4B Xiwa2126-Xgwm113 7.1-38.9 cM LOD 4.00/8.60 R2 0.11/0.12')
F.place(r, ['Xiwa2126', 'Xgwm113'], '4B')
r = add(name='QTs.fcu-5B (AL; Virdi 2016)', chr='5B', disease='Tan spot', marker='Xgwm234', pve='12', stat='LOD 4.2',
        source='Virdi2016G3', note=AL + '; race 2 isolate L13-35; 5BS, peak at Xgwm234 (0 cM); resistance from Altar 84.',
        evidence='Table 3: QTs.fcu-5B Xgwm234 0.0 cM LOD 4.20 R2 0.12')
F.place(r, [], '5B', peak_marker='Xgwm234')
r = add(name='QTs.fcu-6B (Virdi 2016)', chr='6B', disease='Tan spot', marker='Xiwa5148–Xiwa967', pve='22', stat='LOD 6.9',
        source='Virdi2016G3', note=AL + '; race 2 isolate 86-124; 6BL; resistance from Altar 84.',
        evidence='Table 3: QTs.fcu-6B Xiwa5148-Xiwa967 84.4-103.4 cM LOD 6.90 R2 0.22')
F.place(r, ['Xiwa5148', 'Xiwa967'], '6B')
gene(add(name='QSnb.fcu-5B (AL, Tsn1; Virdi 2016)', chr='5B', marker='Tsn1 (Xiwa7024–Xiwa6915)', pve='38', stat='LOD 13.0', source='Virdi2016G3',
         note=AL + '; SNB isolate Sn2000; QTL peak is the Tsn1 locus (flanked by Xiwa7024 and Xiwa6915 at 0.4 and 0.8 cM).',
         evidence='Table 3: QSnb.fcu-5B Tsn1 110.2 cM LOD 13.00 R2 0.38 (Sn2000)'),
     'Tsn1', 'Placed at the Tsn1 gene (RefSeq v1.0).')
add(name='QTs.fcu-5B (LD5B; Virdi 2016)', chr='5B', disease='Tan spot', marker='Xbcd873–Xabg705', pve='19', stat='LOD 3.8', source='Virdi2016G3', grade='D',
    note='Langdon x LDN-DIC 5B RICL (LD5B, 85 lines); race 2 isolate 86-124; distal 5BS, possibly the same as QTs.fcu-5B in AL. '
         'RFLP flanks, not placed on RefSeq.',
    evidence='Table 5: QTs.fcu-5B Xbcd873-Xabg705 0-13.6 cM LOD 3.80 R2 0.19', ref='RefSeq v1.0')
gene(add(name='QSnb.fcu-5B (LD5B, Tsn1; Virdi 2016)', chr='5B', marker='Tsn1', pve='95', stat='LOD 54.0', source='Virdi2016G3',
         note='Langdon x LDN-DIC 5B RICL; SNB isolate Sn2000; Tsn1 locus (first reported by Faris and Friesen 2009).',
         evidence='Table 5: QSnb.fcu-5B Tsn1 48.2 cM LOD 54.0 R2 0.95'),
     'Tsn1', 'Placed at the Tsn1 gene (RefSeq v1.0).')

# ---------- Running 2022 dissertation: Snn5-B2 is on 2BS (the atlas had it on 4B by analogy) ----------
add(name='Snn5-B2', chr='2B', cls='SG', marker='', source='Running2022Diss', grade='D', ref='arm only',
    note='Second SnTox5 sensitivity locus, mapped to the short arm of chromosome 2B in durum wheat (Running 2022 dissertation, abstract and '
         'chapter 5). Marker positions are in chapter 5, which was not among the pages supplied, so the gene is not placed. '
         'Replaces atlas record V1644, which put Snn5-B2 on 4B by analogy with Snn5-B1.',
    evidence='Abstract: a second SnTox5 sensitivity locus, designated Snn5-B2, was mapped to the short arm of chromosome 2B')

NEW_SOURCES = {
    'Cockram2015G3': {'citation': 'Cockram J, Scuderi A, Barber T, Furuki E, Gardner KA, Gosman N, Kowalczyk R, Phan HP, Rose GA, Tan KC, Oliver RP, Mackay IJ (2015) Fine-mapping the wheat Snn1 locus conferring sensitivity to the Parastagonospora nodorum necrotrophic effector SnTox1 using an eight founder multiparent advanced generation inter-cross population. G3 Genes|Genomes|Genetics 5:2257-2266', 'url': 'https://doi.org/10.1534/g3.115.021584', 'short': 'Cockram et al. (2015)'},
    'Reddy2008PG': {'citation': 'Reddy L, Friesen TL, Meinhardt SW, Chao S, Faris JD (2008) Genomic analysis of the Snn1 locus on wheat chromosome arm 1BS and the identification of candidate genes. The Plant Genome 1:55-66', 'url': 'https://doi.org/10.3835/plantgenome2008.03.0181', 'short': 'Reddy et al. (2008)'},
    'Virdi2016G3': {'citation': 'Virdi SK, Liu Z, Overlander ME, Zhang Z, Xu SS, Friesen TL, Faris JD (2016) New insights into the roles of host gene-necrotrophic effector interactions in governing susceptibility of durum wheat to tan spot and Septoria nodorum blotch. G3 Genes|Genomes|Genetics 6:4139-4150', 'url': 'https://doi.org/10.1534/g3.116.036525', 'short': 'Virdi et al. (2016)'},
    'Running2022Diss': {'citation': 'Running KLD (2022) Genomic characterization of necrotrophic effector sensitivity genes in wheat. PhD dissertation, North Dakota State University, Fargo', 'url': '', 'short': 'Running (2022)'},
}
RM_IDS = ['V1644']
SID = lambda i: 'S31-%03d' % i

def build(reports, sources):
    keep = [r for r in reports if r['id'] not in RM_IDS and not r['id'].startswith('S31-')]
    new = []
    for i, x in enumerate(R, 1):
        r = dict(x); r['id'] = SID(i)
        r['note'] = (r.get('note', '') + ' Read from the paper, Nov 2026 (sb31.py).').strip()
        if r.get('peak') is not None and r.get('start') is not None:
            r['start'], r['end'] = min(r['start'], r['peak']), max(r['end'], r['peak'])
        r['asm'] = lift.lift_record(r)
        new.append(r)
    src = dict(sources); src.update(NEW_SOURCES)
    return keep + new, src, new

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], M['sources'], new = build(M['reports'], M['sources'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    S = json.loads(m.group(2).replace('<\\/', '</'))
    S['reports'], S['sources'], new2 = build(S['reports'], S['sources'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(S) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': RM_IDS, 'sources': NEW_SOURCES}, open('/home/claude/work/sb31_db.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
