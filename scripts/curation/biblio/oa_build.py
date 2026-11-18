"""ad62 (Nov 2026): build atlas records from rows extracted out of open-access full texts (oa_extract.py), place them, and write
data/merged_v13f.json, the page seed and web/seed_ad62.json. Records from a paper that is the original study behind a Tong et al. (2024)
compilation record or a Catalogue (CGS) record replace that record when name and chromosome agree."""
import json, re, sys, collections, unicodedata, difflib
sys.path.insert(0, 'scripts/curation/sb19'); import lift
exec(open('scripts/curation/lib/insilico_res.py').read())
D = json.load(open('/tmp/bib/oa_rows.json'))
M = json.load(open('data/merged_v13f.json'))
def nt(s): s = unicodedata.normalize('NFKD', re.sub(r'<[^>]+>', '', s or '')).lower(); return re.sub(r'[^a-z0-9]', '', s)
def nn(n): return re.sub(r'[^a-z0-9]', '', re.split(r'\s\(', n or '')[0].lower())
TAG = ' Auto-extracted from the open-access full text (Europe PMC), Nov 2026 (ad62); verify against the paper.'
GENE = re.compile(r'^(Q[\w.\-]+|Yr\w+|Lr\w+|Sr\w+|Pm\w+|Bt\w+|Ut\w+|Ml\w+|Ug\w*|QTL\w*)$')
MK = re.compile(r'^(X?(gwm|wmc|barc|cfd|cfa|gdm|wms|psp|cdo|bcd|psr|mag|ksum|cfb|stm|hbg|gpw|swes|wpt|cnl|bf|be|bq|cd|ta)\S*|wPt-?\d+|IW[AB]\d+|BS\d+_\d+|BobWhite\S+|Excalibur\S+|Kukri\S+|RAC875\S+|Tdurum\S+|wsnp\S+|Ra_c\S+|Ku_c\S+|GENE-\S+|tplb\S+|CAP\d+\S*|AX-\d+|S\d[ABD]_\d+|IAAV\d+|D_\S+|Jagger\S+|TA\d+-\d+)$', re.I)
def citation(b):
    au = b['author'].split(' and ')
    def fmt(a):
        f, _, g = a.partition(', '); ini = ''.join(x[0] for x in re.split(r'[\s.\-]+', g) if x)
        return f'{f} {ini}'.strip()
    s = ', '.join(fmt(a) for a in au[:6]) + (' et al.' if len(au) > 6 else '')
    first = au[0].split(', ')[0] if au and au[0] else 'Anon'
    short = f'{first} et al. ({b["year"]})' if len(au) > 2 else (f'{first} and {au[1].split(", ")[0]} ({b["year"]})' if len(au) == 2 else f'{first} ({b["year"]})')
    return f'{s} ({b["year"]}) {b["title"]}. {b["journal"]} {b["volume"]}' + (f':{b["pages"].replace("--", "-")}' if b['pages'] else ''), short, first
# compilation records that an original paper can replace
tong = collections.defaultdict(list)
for r in M['reports']:
    if r['source'] == 'Tong2024TAG':
        m = re.search(r'Original study: ([^;.]+?) (?:et al\.? )?(\d{4})', r.get('note', ''))
        if m: tong[(nt(m.group(1).split()[0]), m.group(2))].append(r)
cgs = collections.defaultdict(list)
for r in M['reports']:
    if r['source'].startswith('CGS'): cgs[r['source']].append(r)
cgs_title = {}
for k, v in M['sources'].items():
    if k.startswith('CGS'):
        m = re.search(r'\(\d{4}[a-z]?\)\.?\s*(.+?)\.\s', v.get('citation', '') + ' ')
        if m: cgs_title[nt(m.group(1))[:50]] = k
direct = {nt((re.search(r'\(\d{4}[a-z]?\)\.?\s*(.+?)\.\s', v.get('citation', '') + ' ') or [None, ''])[1])[:50] for k, v in M['sources'].items() if not k.startswith('CGS')}
SRC, NEW, DROP, log = {}, [], set(), collections.Counter()
n = 0; per = collections.Counter()
by = collections.defaultdict(list)
for r in D['rows']: by[r['doi']].append(r)
for doi, rows in by.items():
    b = D['papers'][doi]
    if re.search(r'meta-?analys|meta-?QTL|consensus map|\breview\b', b['title'], re.I): log['skipped meta/review'] += 1; continue
    if nt(b['title'])[:50] in direct: log['skipped already direct'] += 1; continue
    cit, short, first = citation(b); key = 'OA_' + b['key']
    seen = set(); recs = []
    for x in rows:
        name = re.sub(r'\s*[a-z]$|\s*[*†‡§]+$|[\-.]+$', '', x['name']).strip()
        if re.match(r'^(Os|Bradi|Sb|GRMZM|AT\dG|Traes)', name) or len(name) < 2: continue
        if re.fullmatch(r'(MLM|GLM|CMLM|MLMM|FarmCPU|BLINK|SUPER|EMMAX|GAPIT|TASSEL|BLUE|BLUP|E\d|Env\d*|Mean|Average|Total)', name, re.I): continue
        for pat, dd in ((r'^(Yr|QYr|Qyr|QYR|qYr|QPst)', 'Stripe rust'), (r'^(Lr|QLr|Qlr|QLR|qLr)', 'Leaf rust'), (r'^(Sr|QSr|Qsr|QSR|qSr)', 'Stem rust'), (r'^(Pm|QPm|Qpm|QPM|qPm|Ml)', 'Powdery mildew'),
                         (r'^(Bt\d|QCbt|Qcbt|QBt)', 'Common bunt'), (r'^(QKb|Qkb|QKB)', 'Karnal bunt'), (r'^(Ut\d|QLs|QUt)', 'Loose smut')):
            if re.match(pat, name): x = dict(x, dis=dd); break
        ch = x['chr']; k0 = (nn(name), ch, x['dis'])
        if k0 in seen: continue
        seen.add(k0)
        mks = [m for m in re.split(r'[\s,;/]+|(?<=\w)[\-–](?=[Xxw])', re.sub(r'\([^)]*\)', ' ', x['marker'] or '')) if m and MK.match(m)]
        if not mks and MK.match(name): mks = [name]
        isgene = bool(GENE.match(name.split()[0])) and not MK.match(name)
        disp = name   # the locus name as reported (gene, QTL, MTA marker, QTN)
        r = dict(cls='QTL' if re.match(r'^Q', name) else ('DG' if re.match(r'^(Yr|Lr|Sr|Pm|Bt|Ut)\d', name) else 'MTA'), curation='auto', name=disp, disease=x['dis'], chr=ch,
                 source=key, marker=', '.join(mks) if mks else (x['marker'] or ''), pve=(x['pve'] or '')[:20], stat=('LOD ' + x['lod'][:12]) if x['lod'] else (((x.get('ph') or 'P') + ' ' + x['p'][:14]) if x['p'] else ''),
                 evidence=x['cap'], start=None, end=None, peak=None)
        if x.get('phys'):
            v = x['phys']; r['ref'] = x['asm']; r['grade'] = 'A'
            if len(v) >= 2: r.update(start=min(v[:2]), end=max(v[:2]))
            else: r['peak'] = v[0]
            how = f'Physical position as printed ({x["asm"]}).'
            if r.get('end') is not None and r['end'] - r['start'] > 100: r['grade'] = 'C'; how += f' Interval {r["end"] - r["start"]:.0f} Mb wide: a chromosome region, not a locus.'
        else:
            ok = [(m,) + tuple(resolve(m, ch)) for m in mks]; ok = [o for o in ok if o[1]]
            r['ref'] = 'RefSeq v1.0'
            if len(ok) >= 2:
                lo, hi = min(o[1] for o in ok), max(o[1] for o in ok); r.update(start=lo, end=hi, grade='B' if hi - lo <= 100 else 'C')
                how = 'Placed in silico (RefSeq v1.0): ' + '; '.join(f'{m} {p:.3f} Mb' for m, p, h in ok) + '.'
            elif ok:
                r.update(peak=ok[0][1], grade='C'); how = f'Placed in silico at {ok[0][0]} ({ok[0][1]:.3f} Mb, RefSeq v1.0); single marker.'
            elif isgene:
                r['grade'] = 'D'; how = 'No listed marker could be placed.'
            else:
                log['row skipped: no position, not a named locus'] += 1; continue
        r['note'] = f'From {x["cap"][:100]}. ' + how + TAG
        recs.append(r)
    if not recs: log['papers with no placeable row'] += 1; continue
    SRC[key] = {'citation': cit, 'short': short, 'url': f'https://doi.org/{doi}'}
    for r in recs:
        n += 1; r['id'] = f'OA{n:04d}'; r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}; NEW.append(r)
    per[key] = len(recs)
    # replace compilation records of the same paper
    names = {(nn(r['name'].split(' (')[0]), r['chr']) for r in recs}
    for t in tong.get((nt(first), b['year']), []):
        if (nn(t['name']), t['chr']) in names: DROP.add(t['id'])
    ck = cgs_title.get(nt(b['title'])[:50])
    for t in cgs.get(ck, []):
        if (nn(t['name']), t['chr']) in names: DROP.add(t['id'])
print(len(SRC), 'papers;', len(NEW), 'records;', collections.Counter(r['grade'] for r in NEW), ';', len(DROP), 'compilation records replaced;', dict(log))
print(collections.Counter(r['disease'] for r in NEW))
json.dump({'reports': NEW, 'sources': SRC, 'deleted': sorted(DROP), 'per_paper': per}, open('/tmp/bib/ad62.json', 'w'), ensure_ascii=False)
