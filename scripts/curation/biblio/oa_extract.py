"""Extract resistance loci from Europe PMC open-access full text (JATS XML) of the biotroph literature list.
Tables are read column by column: locus/marker names, chromosome, physical position (only when the header states Mb, bp or a physical
position and the paper states a RefSeq version), PVE/R2, LOD, P. Loci without a printed physical position are placed in silico from their
markers (resolve()). Every record is marked curation 'auto' with a note to verify against the paper."""
import json, re, sys, glob, collections, xml.etree.ElementTree as ET
sys.path.insert(0, 'scripts/curation/sb19')
exec(open('scripts/curation/lib/insilico_res.py').read())
E = json.load(open('/tmp/bib/epmc.json'))
P2D = {v['pmcid']: d for d, v in E.items() if v.get('pmcid')}
BIB = open('docs/literature/WheatRLA_biotroph_resistance_loci_papers.bib').read()
BE = {}
for e in BIB.split('\n@article{')[1:]:
    g = lambda k: (re.search(k + r' = \{(.*?)\},?\n', e) or [None, ''])[1]
    BE[g('doi').lower()] = dict(key=e.split(',')[0], title=g('title'), author=g('author'), year=g('year'), journal=g('journal'), volume=g('volume'), pages=g('pages'), kw=g('keywords'))
CH = r'([1-7][ABD])'
DIS = [('Stripe rust', r'stripe rust|yellow rust|\bYR\b|\bYr|striiformis|\bPst\b'), ('Leaf rust', r'leaf rust|brown rust|\bLR\b|\bLr|triticina|\bPt\b'),
       ('Stem rust', r'stem rust|\bSR\b|\bSr\d|graminis f\. ?sp\. ?tritici|\bPgt\b|Ug99'), ('Powdery mildew', r'powdery mildew|\bPM\b|\bPm\d|Blumeria|\bBgt\b'),
       ('Common bunt', r'common bunt|Tilletia caries|T\. caries|stinking smut'), ('Dwarf bunt', r'dwarf bunt|controversa'), ('Karnal bunt', r'karnal bunt|Tilletia indica'),
       ('Loose smut', r'loose smut|Ustilago')]
def T(e): return re.sub(r'\s+', ' ', ''.join(e.itertext())).strip() if e is not None else ''
def dis_of(s, default):
    hit = [d for d, p in DIS if re.search(p, s)]
    return hit[0] if len(hit) == 1 else (default if default in hit or not hit else hit[0])
def num(s):
    s = s.replace('−', '-').replace('–', '-').replace(',', '')
    m = re.findall(r'-?\d+(?:\.\d+)?(?:[eE]-?\d+)?', s); return [float(x) for x in m]
H = dict(name=r'^(QTLs?|QTL name|Locus|Loci|Gene|Genes?|MTAs?|QTNs?|Region|Chromosome region|Designation|Name|Putative gene)\b',
         marker=r'marker|SNP|flanking|closest|interval|peak|nearest|linked', chr=r'^(chr|chromosome)', pve=r'PVE|R2|R²|r2|PEV|variance explained|phenotypic variation|% ?var',
         lod=r'^LOD|LOD score|LOD value', p=r'^p[- ]?value|^P$|^-?log|^−log|probability', phys=r'\b(Mb|Mbp|bp)\b|physical|position on (the )?(CS|Chinese Spring|RefSeq|IWGSC)')
def assembly(text):
    if re.search(r'RefSeq\s*v?\s*2\.\d|v2\.1|IWGSC.{0,30}v2', text, re.I): return 'RefSeq v2.1'
    if re.search(r'RefSeq\s*v?\s*1\.\d|IWGSC.{0,40}(RefSeq|v1\.[01]|2018)|Ref[Ss]eq v1', text, re.I): return 'RefSeq v1.0'
    return None
def rows_of(tw):
    rows = []
    for tr in tw.iter('tr'):
        cells = []
        for c in tr:
            if c.tag not in ('td', 'th'): continue
            cells += [T(c)] * int(c.get('colspan', '1') or 1)
        rows.append(cells)
    return rows
recs, skipped, papers = [], collections.Counter(), {}
for f in sorted(glob.glob('/tmp/bib/xml/*.xml')):
    pmc = f.split('/')[-1][:-4]; doi = P2D.get(pmc); b = BE.get(doi or '')
    if not b: skipped['not in list'] += 1; continue
    try: root = ET.parse(f).getroot()
    except Exception: skipped['bad xml'] += 1; continue
    if (root.get('article-type') or '') in ('review-article', 'systematic-review') or re.search(r'\breview\b|meta-?analys|meta-?QTL|consensus map', b['title'], re.I): skipped['review/meta'] += 1; continue
    body = T(root.find('.//body')); asm = assembly(body)
    defdis = next((d for d in ['Stripe rust', 'Leaf rust', 'Stem rust', 'Powdery mildew', 'Common bunt', 'Dwarf bunt', 'Karnal bunt', 'Loose smut'] if d in b['kw']), None) or dis_of(b['title'], 'Stripe rust')
    got = []
    for tw in root.iter('table-wrap'):
        cap = T(tw.find('label')) + ' ' + T(tw.find('caption'))
        if re.search(r'primer|ANOVA|analysis of variance|heritab|correlat|population structure|polymorphism|LD decay|linkage disequilibrium|statistics|pedigree|virulence|pathotype|race|isolate|gene postulation|presence|absence|haplotype frequenc|plant height|heading|flowering|yield|grain|kernel weight|combination|previously|reported|literature|candidate gene|gene model|annotation|expression|epista|gene family|family genes|details of selected|genes in (the )?wheat genome|members', cap, re.I): continue
        R = [r for r in rows_of(tw) if any(x for x in r)]
        if len(R) < 2: continue
        # header = first row(s) until one with a chromosome-like value
        hdr = R[0]; i0 = 1
        if len(R) > 2 and not any(re.fullmatch(CH + r'[SL]?', x) for x in R[1]): hdr = [a + ' ' + c for a, c in zip(R[0], R[1] + [''] * len(R[0]))]; i0 = 2
        col = {}
        for k, p in H.items():
            for j, h in enumerate(hdr):
                if k not in col and re.search(p, h, re.I) and not (k == 'marker' and re.search(r'^(No|Number|N)\b|primer|sequence|allele', h, re.I)):
                    if k == 'name' and 'marker' in col and col['marker'] == j: continue
                    col.setdefault(k, j)
        mcols = [j for j, h in enumerate(hdr) if re.search(H['marker'], h, re.I) and not re.search(r'^(No|Number|N)\b|primer|sequence|allele|position|distance|cM|Mb|bp|effect|frequenc', h, re.I)]
        if any(re.search(r'^(Reference|References|Citation|Literature|Author)', h, re.I) for h in hdr): continue
        if 'chr' not in col and 'name' not in col: continue
        if 'marker' not in col and 'phys' not in col and 'name' not in col: continue
        for r in R[i0:]:
            g = lambda k: r[col[k]] if k in col and col[k] < len(r) else ''
            nm, cs = g('name'), g('chr'); mk = ' '.join(r[j] for j in mcols if j < len(r))
            m = re.search(CH, cs) or re.search(r'[.\-]' + CH + r'[SL]?\b', nm) or re.search(r'(?:QTL|Q\w+?)' + CH, nm)
            if not m: continue
            ch = m.group(1)
            if not (nm or mk): continue
            if not (nm or mk.split()): continue
            name = nm or mk.split()[0]
            if len(name) > 60 or re.fullmatch(r'[\d.\-– ]+', name): continue
            r2 = dict(name=name, chr=ch, marker=mk, pve=g('pve'), lod=g('lod'), p=g('p'), ph=(re.sub(r'\s+[a-z]$', '', hdr[col['p']]).strip() if 'p' in col else ''), cap=cap[:120], dis=dis_of(cap + ' ' + ' '.join(r), defdis))
            ph = g('phys')
            if ph and asm and 'phys' in col:
                v = num(ph); hdrp = hdr[col['phys']]
                if v:
                    unit = 1e-6 if (re.search(r'\bbp\b', hdrp) and not re.search(r'Mb', hdrp)) or max(v) > 1e5 else 1
                    v = [x * unit for x in v if x >= 0]
                    if v and max(v) < 900: r2.update(phys=v, asm=asm)
            got.append(r2)
    if got: papers[doi] = dict(b, pmc=pmc, asm=asm, n=len(got)); recs += [dict(x, doi=doi) for x in got]
    else: skipped['no locus table'] += 1
print(len(papers), 'papers with loci;', len(recs), 'rows;', dict(skipped))
json.dump({'papers': papers, 'rows': recs}, open('/tmp/bib/oa_rows.json', 'w'))
