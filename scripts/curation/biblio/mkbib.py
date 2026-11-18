import json, re, unicodedata, collections
K = json.load(open('/tmp/bib/keep.json'))
M = json.load(open('/home/claude/wheatnhr-atlas/data/merged_v13f.json'))
BIO = {'Leaf rust', 'Stripe rust', 'Stem rust', 'Powdery mildew', 'Multiple rusts', 'Crown rust', 'Common bunt', 'Dwarf bunt', 'Loose smut', 'Karnal bunt'}
def nt(s): s = unicodedata.normalize('NFKD', re.sub(r'<[^>]+>', '', s or '')).lower(); return re.sub(r'[^a-z0-9]', '', s)
# direct atlas sources (curated from the paper itself): not Catalogue-derived (CGS_), not compilations
direct_doi, direct_t, cat_t = set(), set(), {}
for k, v in M['sources'].items():
    c = v.get('citation', ''); m = re.search(r'\(\d{4}[a-z]?\)\.?\s*(.+?)\.\s', c + ' ')
    t = nt(m.group(1))[:60] if m else None
    if k.startswith('CGS'):
        if t: cat_t[t] = k
        continue
    u = (v.get('url') or '').lower(); md = re.search(r'10\.\d{4,9}/\S+', u)
    if md: direct_doi.add(md.group(0).rstrip('.'))
    if t: direct_t.add(t)
# Tong 2024 original studies (first author, year)
tong = set()
for r in M['reports']:
    if r['source'] == 'Tong2024TAG':
        m = re.search(r'Original study: ([^;.]+?) (?:et al\.? )?(\d{4})', r.get('note', ''))
        if m: tong.add((nt(m.group(1).split()[0]), m.group(2)))
def kind(t):
    k = []
    if re.search(r'genome[- ]wide association|GWAS|association (mapping|analysis|study)|marker[- ]trait association', t, re.I): k.append('GWAS')
    if re.search(r'\bQTLs?\b|quantitative trait', t, re.I): k.append('QTL mapping')
    if re.search(r'\bclon(ed|ing)\b|encodes|map-based', t, re.I): k.append('gene cloning')
    if re.search(r'fine[- ]map', t, re.I): k.append('fine mapping')
    if re.search(r'introgress|translocation|alien|addition line|substitution line|aegilops|thinopyrum|dasypyrum|haynaldia|secale|\brye\b|elymus|agropyron', t, re.I): k.append('alien introgression')
    if re.search(r'\b(mapping|mapped|map|location|locat\w+|linkage|tagging|identification of .*gene)\b', t, re.I) and not k: k.append('gene mapping')
    if re.search(r'detection|presence|postulat|screening|distribution|diversity|molecular characteri[sz]ation of .* (cultivars|germplasm|lines)', t, re.I): k.append('marker screening / postulation')
    if re.search(r'overexpress|silenc|knock|positively|negatively|regulat|transcriptom|express', t, re.I): k.append('functional / expression')
    return k or ['other genetics']
def esc(s): return re.sub(r'([&%$#_{}])', r'\\\1', re.sub(r'<[^>]+>', '', s or '')).replace('\n', ' ').strip()
used = set(); out = []; stats = collections.Counter()
for doi, it in sorted(K.items(), key=lambda x: ((x[1].get('issued', {}).get('date-parts') or [[0]])[0][0] or 0, x[0])):
    t = ' '.join(it.get('title') or ['']).strip(); tn = nt(t)[:60]
    if doi in direct_doi or tn in direct_t: stats['excluded: already curated directly'] += 1; continue
    au = it.get('author') or []
    y = (it.get('issued', {}).get('date-parts') or [[None]])[0][0]
    fam = (au[0].get('family') or au[0].get('name') or 'Anon') if au else 'Anon'
    key = re.sub(r'[^A-Za-z]', '', unicodedata.normalize('NFKD', fam).encode('ascii', 'ignore').decode()) + str(y or '') + (re.sub(r'[^a-z]', '', nt(t))[:6] or 'x')
    while key in used: key += 'b'
    used.add(key)
    tags = list(it['_dis']) + kind(t)
    st = []
    if tn in cat_t: st.append('atlas: loci entered via Catalogue of Gene Symbols')
    if (nt(fam), str(y)) in tong: st.append('atlas: loci entered via Tong et al. 2024 compilation')
    if not st: st.append('atlas: not yet in atlas')
    tags += st; stats[st[0]] += 1
    tags = [x.replace('Rusts (general)', 'rust') for x in tags]
    authors = ' and '.join(esc((a.get('family') or a.get('name') or '') + (', ' + a['given'] if a.get('given') else '')) for a in au)
    f = [f'  title = {{{esc(t)}}}', f'  author = {{{authors}}}' if authors else None, f'  year = {{{y}}}' if y else None,
         f'  journal = {{{esc(" ".join(it.get("container-title") or []))}}}' if it.get('container-title') else None,
         f'  volume = {{{esc(it["volume"])}}}' if it.get('volume') else None, f'  number = {{{esc(it["issue"])}}}' if it.get('issue') else None,
         f'  pages = {{{esc(it["page"]).replace("-", "--")}}}' if it.get('page') else None, f'  doi = {{{doi}}}', f'  url = {{https://doi.org/{doi}}}',
         f'  keywords = {{{", ".join(dict.fromkeys(tags))}}}']
    out.append('@article{' + key + ',\n' + ',\n'.join(x for x in f if x) + '\n}')
open('/tmp/bib/WheatRLA_biotroph_resistance_loci_papers.bib', 'w').write('% WheatRLA: wheat resistance-locus papers for biotrophic diseases (rusts, powdery mildew, bunts, smuts)\n% Built 18 Nov 2026 from Crossref metadata; papers already curated directly into the atlas are left out.\n% keywords become Zotero tags: disease, study type, and atlas status.\n\n' + '\n\n'.join(out) + '\n')
print(len(out), stats)
print(collections.Counter(k for e in out for k in re.findall(r'keywords = \{(.*)\}', e)[0].split(', ')).most_common(25))
