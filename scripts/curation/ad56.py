"""ad56 (Nov 2026): Catalogue of Gene Symbols for Wheat, Supplements 2004-2017 (2012 file uploaded twice, identical).
Disease-resistance entries (designated, temporary and QTL designations; insects and mites excluded) read by hand into
docs/supplements/cgsw/e<year>.txt (name|cls|disease|chr|markers|ref|pve|stat|note). Each entry is attributed to the original
paper cited in the supplement ({ref} resolved from the supplement reference lists); where the reference is not listed, the
supplement itself is the source. Markers placed in silico; where none places and the atlas already holds the same designated
gene, that position is reused (grade C). Entries already in the atlas from the same paper are skipped."""
import json, re, sys, glob, difflib, collections
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
DIS = {'FHB': 'Fusarium head blight', 'CR': 'Fusarium crown rot', 'STB': 'Septoria tritici blotch', 'SNB': 'Septoria nodorum blotch', 'TS': 'Tan spot',
       'SB': 'Spot blotch', 'LR': 'Leaf rust', 'YR': 'Stripe rust', 'SR': 'Stem rust', 'PM': 'Powdery mildew', 'CCN': 'Cereal cyst nematode',
       'EYE': 'Eyespot', 'BYD': 'Barley yellow dwarf virus', 'WSMV': 'Wheat streak mosaic virus', 'WYMV': 'Wheat yellow mosaic virus',
       'SBCMV': 'Soil-borne cereal mosaic virus', 'WSSMV': 'Wheat spindle streak mosaic virus', 'CB': 'Common bunt', 'DB': 'Dwarf bunt',
       'KB': 'Karnal bunt', 'LS': 'Loose smut', 'WB': 'Wheat blast', 'CS': 'Cephalosporium stripe', 'RKN': 'Root-knot nematode',
       'ANT': 'Anthracnose', 'CRR': 'Crown rust'}
REF = json.load(open('docs/supplements/cgsw/refs.json'))
D = 'docs/supplements/cgsw/'
M = json.load(open('data/merged_v13f.json'))
SRC = M['sources']
def ttl(c):
    c = re.sub(r'<[^>]+>', '', c); m = re.search(r'\b(19|20|21)\d\d[a-z]?\)?\.?\s+(.+?)\.\s', c + ' ')
    return re.sub(r'[^a-z]', '', (m.group(2) if m else c).lower())[:90]
STT = {k: ttl(v.get('citation', '')) for k, v in SRC.items()}
def short(c):
    m = re.match(r'(.+?)\s\(?(\d{4})', c)
    if not m: return None
    au, yr = m.groups()
    names = re.findall(r"([A-Z][A-Za-z'\- ]+?) [A-Z]{1,4}\b(?:,|\s*&|\s*and|$)", au.replace('andand', 'and'))
    if re.search(r'et al', au): return f"{re.split(r'[ ,]', au.strip())[0]} et al. ({yr})"
    if not names: return f"{au.split()[0]} ({yr})"
    if len(names) == 1: return f'{names[0]} ({yr})'
    if len(names) == 2: return f'{names[0]} and {names[1]} ({yr})'
    return f'{names[0]} et al. ({yr})'
def _old_short(c):
    m = re.match(r'(.+?)\s(\d{4})[a-z]?\s', c); au, yr = (m.group(1), m.group(2)) if m else (c[:30], '')
    first = re.split(r'[ ,]', au.strip())[0]; n = len(re.findall(r'[A-Z][a-zA-Z\'-]+ [A-Z]{1,3}\b', au))
    return f'{first} et al. ({yr})' if n > 2 or ' & ' not in au and ',' in au else (f"{first} and {re.split(r' & | and ', au)[-1].split()[0]} ({yr})" if (' & ' in au or ' and ' in au) and n == 2 else f'{first} ({yr})')
NEWSRC = {}
FIXC = {'10677': 'Li Q, Chen XM, Wang MN & Jing JX (2011) Yr45, a new wheat gene for stripe rust resistance on the long arm of chromosome 3D. Theoretical and Applied Genetics 122: 189-197.',
        '0301': 'Xie et al. (2003) Theoretical and Applied Genetics 106: 341-345.'}
FIXS = {'10175': 'Singrun et al. (2004)', '10272': 'Uauy et al. (2005)', '10340': 'Xu et al. (2006)', '10498': 'Sukhwinder-Singh et al. (2003)',
        '10499': 'Sukhwinder-Singh et al. (2007)', '10591': 'Marais et al. (2009)', '10717': 'Chu et al. (2009)', '10968': 'Herrera-Foessel et al. (2015)',
        '10039': 'Zahravi et al. (2003)', '10677': 'Li et al. (2011)', '0301': 'Xie et al. (2003)', '10013': 'De Majnik et al. (2003)'}
def clean(ref, c):
    c = FIXC.get(ref, c); c = re.sub(r'^(Update:|Replace with:|Revised reference:|\d{4,5}\.)\s*', '', c.strip())
    c = re.split(r'(?<=communication\.)\s\d{4,5}\.?\s', c)[0] if re.search(r'Personal communication\.\s\d{4,5}', c) else c
    return c.replace('andand', 'and').strip()
def source_for(ref, year):
    c = REF.get(ref); c = clean(ref, c) if c else c
    if not c:
        k = f'CGSW_{year}'
        NEWSRC[k] = {'citation': f'McIntosh RA et al. ({year}) Catalogue of gene symbols for wheat: {year} supplement. Annual Wheat Newsletter.', 'short': f'McIntosh et al. ({year})', 'url': ''}
        return k, False
    t = ttl(c); best = max(((difflib.SequenceMatcher(None, t, v).ratio(), k) for k, v in STT.items() if v), default=(0, None))
    if best[0] >= 0.9: return best[1], True
    k = f'CGS_{ref}'
    NEWSRC[k] = {'citation': re.sub(r'\s(\d{4})\s', r' (\1) ', c, count=1), 'short': FIXS.get(ref) or short(c), 'url': ''}
    if not NEWSRC[k]['short']: del NEWSRC[k]; return source_for(None, year)
    return k, False
def nname(n): return re.sub(r'[^a-z0-9]', '', re.split(r'\s\(', n)[0].lower())
EXIST = collections.defaultdict(list)
for r in M['reports']: EXIST[nname(r.get('name', ''))].append(r)
def base(n): return re.split(r'[\s(/]', n.strip())[0].lower()
GP = collections.defaultdict(list)
for r in M['reports']:
    if r.get('peak') is not None or r.get('start') is not None: GP[(base(r['name']), r['chr'])].append(r)
def gene_pos(name, ch, need=('A', 'B', 'C')):
    xs = [r for r in GP.get((base(name), ch), []) if r.get('grade') in need and r.get('cls') in ('DG', 'SG', 'CG')
          and (r.get('end') or r.get('peak') or 0) - (r.get('start') or r.get('peak') or 0) <= 100]
    return sorted(xs, key=lambda r: ('ABC'.index(r['grade']), r['id']))[0] if xs else None
DUPX = {('stb1', '5B'): 'Adhikari2004', ('qfhsndsu3as', '3A'): 'Z_1325', ('qfhsnau2dl', '2D'): 'Z_1183', ('fhb7', '7E'): 'Guo2015TAG'}
rows = []
for f in sorted(glob.glob(D + 'e20*.txt')):
    yr = re.search(r'e(\d{4})', f).group(1)
    for l in open(f):
        p = l.rstrip('\n').split('|')
        if len(p) == 9: rows.append([yr] + p)
# within-catalogue dedupe: same ref + base name + chr -> keep latest supplement
seen = {}
for i, (yr, name, cls, dis, ch, mk, ref, pve, stat, note) in enumerate(rows): seen[(ref, nname(name), ch)] = i
rows = [r for i, r in enumerate(rows) if seen[(r[6], nname(r[1]), r[4])] == i]
out = []; skipped = []; n = 0
for yr, name, cls, dis, ch, mk, ref, pve, stat, note in rows:
    sk, existing = source_for(ref, yr)
    if existing and any(r['source'] == sk and r['chr'] == ch for r in EXIST.get(nname(name), [])):
        skipped.append((name, sk, 'same paper already in atlas')); continue
    if (nname(name), ch) in DUPX: skipped.append((name, sk, 'same paper already in atlas under ' + DUPX[(nname(name), ch)])); continue
    n += 1
    r = dict(id=f'K{n:04d}', cls=cls, ref='RefSeq v1.0', curation='hand', name=name, disease=DIS[dis], chr=ch, source=sk,
             evidence=f'Catalogue of Gene Symbols for Wheat {yr} Supplement' + (f' {{{ref}}}' if ref else ''))
    if pve: r['pve'] = pve
    if stat: r['stat'] = stat
    ms = [m for m in mk.split(';') if m and m != '-']
    if ms: r['marker'] = ', '.join(ms)
    ok = [(m,) + resolve(m, ch) for m in ms]; ok = [x for x in ok if x[1]]
    alien = 'alien' in note
    srcnote = '' if sk.startswith('CGS_') or existing else f' The original paper {{{ref}}} is not in the supplement reference list, so the supplement is cited.'
    g0 = gene_pos(name, ch, ('A', 'B')) if cls in ('DG', 'SG') else None
    weak = len(ok) < 2 or (max(x[1] for x in ok) - min(x[1] for x in ok) > 100)
    if g0 and weak and ok: ok = []
    if len(ok) >= 2:
        lo, hi = min(x[1] for x in ok), max(x[1] for x in ok); r.update(start=lo, end=hi, grade='B' if hi - lo <= 100 else 'C')
        t = ' Placed in silico (RefSeq v1.0): ' + '; '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in ok) + '.' + (f' The markers span {hi-lo:.0f} Mb: a chromosome region, not a locus.' if hi - lo > 100 else '')
    elif ok:
        m, p, h = ok[0]; r.update(peak=p, grade='C'); t = f' Placed at {m} ({p:.3f} Mb, {h}); single-marker position.'
    else:
        g = gene_pos(name, ch) if cls in ('DG', 'SG') else None
        if g is None and base(name) == 'lr34': g = {'name': 'Lr34 (TraesCS7D02G080300)', 'id': 'RefSeq v1.0 gene model', 'start': 47.410, 'end': 47.422, 'grade': 'A'}
        if g:
            for k in ('start', 'end', 'peak'):
                if g.get(k) is not None: r[k] = g[k]
            r['grade'] = 'C'; t = (' The listed markers place loosely, so' if ms else ' No marker listed;') + f' the position of {g["name"]} ({g["id"]}) is used.'
        else:
            r['grade'] = 'D'; t = (' Alien segment; no wheat marker to place.' if alien else (' No listed marker could be placed on ' + ch + '; left unplaced.' if ms else ' No placeable marker given; left unplaced.'))
    r['note'] = note.rstrip('.') + '.' + t + srcnote + f' Read from the Catalogue of Gene Symbols for Wheat {yr} Supplement, Nov 2026 (ad56.py).'
    r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
    out.append(r)
json.dump({'reports': out, 'sources': NEWSRC, 'skipped': skipped}, open('/tmp/cat/ad56_out.json', 'w'), ensure_ascii=False)
print(len(rows), 'entries;', len(out), 'records;', len(skipped), 'skipped;', len(NEWSRC), 'new sources')
print(collections.Counter(r['grade'] for r in out)); print(collections.Counter(r['disease'] for r in out).most_common())
if __name__ == '__main__':
    NEWDIS = ',["Crown rust","#a0522d","Biotroph"],["Cephalosporium stripe","#556b2f","Necrotroph"],["Anthracnose","#cd5c5c","Hemibiotroph"],["Soil-borne cereal mosaic virus","#4682b4","Virus"],["Wheat spindle streak mosaic virus","#20b2aa","Virus"],["Root-knot nematode","#9acd32","Nematode"]'
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    M['reports'] += out; M['sources'].update(NEWSRC); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'] += out; Sd['sources'].update(NEWSRC)
    H = H[:m.start(2)] + ser(Sd) + H[m.end(2):]
    if '"#9acd32","Nematode"' not in H: H = H.replace('["Cereal cyst nematode","#6b8e23","Nematode"]]', '["Cereal cyst nematode","#6b8e23","Nematode"]' + NEWDIS + ']')
    open('web/index.html', 'w').write(H)
    ops = [{'op': 'set', 'path': 'sources/' + k, 'data': v} for k, v in NEWSRC.items()]
    for r in out:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/tmp/ops56.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    print(len(ops), 'ops')
