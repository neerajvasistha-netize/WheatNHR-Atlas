"""Batch 4 (Nov 2026): Septoria tritici blotch papers from Archive.zip. Records are defined in stb44_recs.py;
markers are placed with the sb26 resolver (genes, exact-match in-silico PCR incl. Roder 1998 GWM primers, marker tables).
Writes data/merged_v13f.json, web/index.html (seed) and work/stb44_db.json. Run after sb31.py."""
import json, re, sys, io, contextlib
sys.path.insert(0, 'scripts/curation'); sys.path.insert(0, 'scripts/curation/sb19')
with contextlib.redirect_stdout(io.StringIO()):
    import sb26_fix as F
import lift, stb44_recs as RR

RM_SOURCES = ('A_STB_042 V_Adhikari2015 V_Arraiano2006 A_STB_045 A_STB_003 V_Cobertera2022 V_Pb V_Ferjaoui2021 A_STB_005 Yang2020Phyto '
              'A_STB_013 IHAR_STB Z_1223 A_STB_006 Agron2019STB Agron2021Mandub A_STB_024 V_103390agronomy1106 A_STB_011 X5_1871 '
              'A_STB_021 Teagasc_STB').split()
KEYMAP = {'V_Arraiano2006': 'Arraiano2007PP', 'A_STB_005': 'Ghaffary2012TAG', 'A_STB_006': 'Liu2013Eu', 'V_Ferjaoui2021': 'Ferjaoui2022BMC', 'Z_1223': 'LanglandsPerry2022'}
C = lambda cit, doi, short: {'citation': cit, 'url': ('https://doi.org/' + doi) if doi else '', 'short': short}
NEW_SOURCES = {
 'Arraiano2007PP': C('Arraiano LS, Chartrain L, Bossolini E, Slatter HN, Keller B, Brown JKM (2007) A gene in European wheat cultivars for resistance to an African isolate of Mycosphaerella graminicola. Plant Pathology 56:73-78', '10.1111/j.1365-3059.2006.01499.x', 'Arraiano et al. (2007)'),
 'Ghaffary2012TAG': C('Tabib Ghaffary SM, Faris JD, Friesen TL, Visser RGF, van der Lee TAJ, Robert O, Kema GHJ (2012) New broad-spectrum resistance to septoria tritici blotch derived from synthetic hexaploid wheat. Theoretical and Applied Genetics 124:125-142', '10.1007/s00122-011-1692-7', 'Tabib Ghaffary et al. (2012)'),
 'Liu2013Eu': C('Liu Y, Zhang L, Thompson IA, Goodwin SB, Ohm HW (2013) Molecular mapping re-locates the Stb2 gene for resistance to Septoria tritici blotch derived from cultivar Veranopolis on wheat chromosome 1BS. Euphytica 190:145-156', '10.1007/s10681-012-0796-8', 'Liu et al. (2013)'),
 'Ferjaoui2022BMC': C('Ferjaoui S, Aouini L, Slimane RB, Ammar K, Dreisigacker S, Schouten HJ, Sapkota S, Bahri BA, Ben MBarek S, Visser RGF, Kema GHJ, Hamza S (2022) Deciphering resistance to Zymoseptoria tritici in the Tunisian durum wheat landrace accession Agili39. BMC Genomics 23:372', '10.1186/s12864-022-08560-2', 'Ferjaoui et al. (2022)'),
 'LanglandsPerry2022': C('Langlands-Perry C, Cuenin M, Bergez C, Krima SB, Gelisse S, Sourdille P et al. (2022) Resistance of the wheat cultivar Renan to Septoria leaf blotch explained by a combination of strain specific and strain non-specific QTL mapped on an ultra-dense genetic map. Genes 13:100', '10.3390/genes13010100', 'Langlands-Perry et al. (2022)'),
 'Chartrain2004Phyto': C('Chartrain L, Brading PA, Widdowson JP, Brown JKM (2004) Partial resistance to Septoria tritici blotch (Mycosphaerella graminicola) in wheat cultivars Arina and Riband. Phytopathology 94:497-504', '10.1094/PHYTO.2004.94.5.497', 'Chartrain et al. (2004)'),
 'Dreisigacker2015TAG': C('Dreisigacker S, Wang X, Martinez Cisneros BA, Jing R, Singh PK (2015) Adult-plant resistance to Septoria tritici blotch in hexaploid spring wheat. Theoretical and Applied Genetics 128:2317-2329', '10.1007/s00122-015-2587-9', 'Dreisigacker et al. (2015)'),
 'Ghaffary2011TAG': C('Tabib Ghaffary SM, Robert O, Laurent V, Lonnet P, Margale E, van der Lee TAJ, Visser RGF, Kema GHJ (2011) Genetic analysis of resistance to septoria tritici blotch in the French winter wheat cultivars Balance and Apache. Theoretical and Applied Genetics 123:741-754', '10.1007/s00122-011-1623-7', 'Tabib Ghaffary et al. (2011)'),
 'Goodwin2015CS': C('Goodwin SB, Cavaletto JR, Hale IL, Thompson I, Xu SS, Adhikari TB, Dubcovsky J (2015) A new map location of gene Stb3 for resistance to Septoria tritici blotch in wheat. Crop Science 55:35-43', '10.2135/cropsci2013.11.0766', 'Goodwin et al. (2015)'),
 'Ilincic2019MB': C('Tamburic-Ilincic L, Barcellos Rosa S (2019) QTL mapping of Fusarium head blight and Septoria tritici blotch in an elite hard red winter wheat population. Molecular Breeding 39:94', '10.1007/s11032-019-0999-y', 'Tamburic-Ilincic and Barcellos Rosa (2019)'),
 'Karlstedt2019Eu': C('Karlstedt F, Kopahnke D, Perovic D, Jacobi A, Pillen K, Ordon F (2019) Mapping of quantitative trait loci (QTL) for resistance against Zymoseptoria tritici in the winter spelt wheat accession HTRI1410 (Triticum aestivum subsp. spelta). Euphytica 215:108', '10.1007/s10681-019-2432-3', 'Karlstedt et al. (2019)'),
 'Kelm2012MB': C('Kelm C, Tabib Ghaffary SM, Bruelheide H, Roder MS, Miersch S, Weber WE, Kema GHJ, Saal B (2012) The genetic architecture of seedling resistance to Septoria tritici blotch in the winter wheat doubled-haploid population Solitaer x Mazurka. Molecular Breeding 29:813-830', '10.1007/s11032-011-9592-8', 'Kelm et al. (2012)'),
 'Mergoum2013CRC': C('Mergoum M, Harilal VE, Singh PK, Adhikari TB, Kumar A, Ghavami F, Elias E, Alamri MS, Kianian SF (2013) Genetic analysis and mapping of seedling resistance to Septoria tritici blotch in Steele-ND/ND 735 bread wheat population. Cereal Research Communications 41:199-210', '10.1556/CRC.2013.0005', 'Mergoum et al. (2013)'),
 'Miedaner2012MB': C('Miedaner T, Risser P, Paillard S, Schnurbusch T, Keller B, Hartl L, Holzapfel J, Korzun V, Ebmeyer E, Utz HF (2012) Broad-spectrum resistance loci for three quantitatively inherited diseases in two winter wheat populations. Molecular Breeding 29:731-742', '10.1007/s11032-011-9586-6', 'Miedaner et al. (2012)'),
 'Naz2015PB': C('Naz AA, Klaus M, Pillen K, Leon J (2015) Genetic analysis and detection of new QTL alleles for Septoria tritici blotch resistance using two advanced backcross wheat populations. Plant Breeding 134:514-519', '10.1111/pbr.12301', 'Naz et al. (2015)'),
 'Raman2009MB': C('Raman R, Milgate AW, Imtiaz M, Tan MK, Raman H, Lisle C, Coombes N, Martin P (2009) Molecular mapping and physical location of major gene conferring seedling resistance to Septoria tritici blotch in wheat. Molecular Breeding 24:153-164', '10.1007/s11032-009-9280-0', 'Raman et al. (2009)'),
 'SanchezVidana2022': C('Sanchez-Vidana MR, Vargas-Hernandez M, He X, Sansaloni CP, Kumar-Singh P, Hernandez-Anguiano AM (2022) Mapping of QTL in the Mutus#1 wheat line resistant to leaf spot. Revista Mexicana de Ciencias Agricolas 13(7)', '', 'Sanchez-Vidana et al. (2022)'),
 'Simon2010EJPP': C('Simon MR, Khlestkina EK, Castillo NS, Borner A (2010) Mapping quantitative resistance to septoria tritici blotch in spelt wheat. European Journal of Plant Pathology 128:317-324', '10.1007/s10658-010-9640-y', 'Simon et al. (2010)'),
 'Stadlmeier2019G3': C('Stadlmeier M, Jorgensen LN, Corsi B, Cockram J, Hartl L, Mohler V (2019) Genetic dissection of resistance to the three fungal plant pathogens Blumeria graminis, Zymoseptoria tritici, and Pyrenophora tritici-repentis using a multiparental winter wheat population. G3 Genes|Genomes|Genetics 9:1745-1757', '10.1534/g3.119.400068', 'Stadlmeier et al. (2019)'),
}

POS0 = F.pos
def P(m, c):
    if m.lower().lstrip('x') in ('barc008', 'barc8') and c == '1B':   # cosegregates with Xwmc230 on the Dreisigacker et al. 2015 map
        p, h = POS0('Xwmc230', '1B')
        return p, 'position of Xwmc230, with which Xbarc008 cosegregates (Dreisigacker et al. 2015 map)'
    for v in (m, ('wsnp_' + m) if m.startswith(('Ku_', 'Ex_', 'Ra_')) else None, m.replace('Tag_', ''), m.replace('wPt_', 'wPt-')):
        if not v: continue
        p, h = POS0(v, c)
        if p is not None: return p, h
    return None, ''

def riaz():
    T = json.load(open('scripts/curation/stb44_riaz_table5.json'))
    fixed = {}
    for q, rows in T.items():
        for x in rows:
            name = q
            if q.endswith('5D.1') and x['chr'] == '6A': name = 'QStb.teagasc-6A.1'
            fixed.setdefault(name, []).append(x)
    fixed['QStb.teagasc-6B.3'] = [dict(ds='2018_T1_flag', chr='6B', s=None, e=None, peak='BS00048295_51', pp='513.69', lod='14.3', r2='17.94')]
    out = []
    for q, rows in sorted(fixed.items()):
        b = max(rows, key=lambda x: float(x['lod']))
        pp = float(b['pp']) if b['pp'] not in ('-', '') else None
        fx = dict(grade='A')
        if b['s'] is not None: fx.update(start=b['s'], end=b['e'])
        if pp is not None and (b['s'] is None or b['s'] <= pp <= b['e']): fx['peak'] = pp
        note = ('NIAB Elite MAGIC (Teagasc), field STB; %d dataset(s); interval and peak are the physical positions printed in Table 5 for the '
                'dataset with the highest LOD (%s, LOD %s).' % (len(rows), b['ds'], b['lod']))
        if pp is not None and 'peak' not in fx: note += ' The printed peak position (%.2f Mb) lies outside that interval and is not used.' % pp
        out.append(dict(source='A_STB_021', name=q, chr=b['chr'], markers=[], peak_marker=None, fixed=fx, pve=', '.join(x['r2'] for x in rows),
                        stat='LOD ' + ', '.join(x['lod'] for x in rows), note=note, evidence='Table 5: %s peak %s' % (q, b['peak']), extra={}))
    return out

def make(x):
    rec = dict(name=x['name'], chr=x['chr'], source=KEYMAP.get(x['source'], x['source']), pve=x['pve'], stat=x['stat'], note=x['note'],
               evidence=x['evidence'], ref='RefSeq v1.0', curation='hand', cls=x['extra'].get('cls', 'QTL'),
               disease=x['extra'].get('disease', RR.STB), marker=', '.join(x['markers'] + ([x['peak_marker']] if x['peak_marker'] else [])),
               start=None, end=None, peak=None, grade=None)
    fx = x['fixed']
    if fx:
        for k in ('start', 'end', 'peak', 'grade'):
            if fx.get(k) is not None: rec[k] = fx[k]
        if x['markers']:
            got = [(m,) + P(m, x['chr']) for m in x['markers']]
            got = [g for g in got if g[1] is not None]
            if got: rec['note'] += ' Marker positions (RefSeq v1.0): ' + '; '.join('%s %.3f Mb (%s)' % g for g in got) + '.'
    else:
        F.pos = P
        try: F.place(rec, x['markers'], x['chr'], peak_marker=x['peak_marker'])
        finally: F.pos = POS0
    for k in ('start', 'end', 'peak'):
        if rec[k] is not None: rec[k] = round(float(rec[k]), 4)
    if rec['peak'] is not None and rec['start'] is not None:
        rec['start'], rec['end'] = min(rec['start'], rec['peak']), max(rec['end'], rec['peak'])
    rec['grade'] = rec['grade'] or 'D'
    return rec

R = [make(x) for x in RR.R + riaz()]
SID = lambda i: 'S44-%03d' % i

def build(reports, sources):
    keep = [r for r in reports if r['source'] not in RM_SOURCES and not r['id'].startswith('S44-')]
    removed = [r['id'] for r in reports if r['source'] in RM_SOURCES]
    new = []
    for i, x in enumerate(R, 1):
        r = dict(x); r['id'] = SID(i)
        r['note'] = (r['note'] + ' Read from the paper, Nov 2026 (stb44.py).').strip()
        r['asm'] = lift.lift_record(r)
        new.append(r)
    src = dict(sources); src.update(NEW_SOURCES)
    return keep + new, src, new, removed

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], M['sources'], new, rem = build(M['reports'], M['sources'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    S = json.loads(m.group(2).replace('<\\/', '</'))
    S['reports'], S['sources'], new2, rem2 = build(S['reports'], S['sources'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(S) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': rem2, 'sources': NEW_SOURCES}, open('/home/claude/work/stb44_db.json', 'w'), ensure_ascii=False)
    import collections
    print('new', len(new), 'removed', len(rem), len(rem2), collections.Counter(r['grade'] for r in new),
          collections.Counter((r['asm'].get('RefSeq v2.1') or {}).get('flag', 'none') for r in new))
    for r in new: print(r['id'], r['grade'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['name'][:40])
