"""ad60 (Nov 2026): audit fixes after ad56-ad59.
1. 83 Catalogue records (K) reused the position of the same gene already in the atlas but kept ref 'RefSeq v1.0' although the
   donor record's coordinates were on RefSeq v2.1 (or v1.1/v2.0). Ref now follows the donor record; re-lifted.
2. Source attribution re-checked against each supplement's own reference list (the merged reference table had collisions):
   K0122 -> Srnic et al. (2005); K0270-K0272, K0275, K0429 -> Friesen et al. (2007) Plant J; K0307-K0308 -> Schmolke et al. (2008);
   K0610-K0613 -> Powell et al. (2013).
3. Duplicates of records already in the atlas from the same paper, now that the source is right: K0430 (Snn4 = H031), K0492 (Fhb4 = H005;
   supplement cites {10884}), K0143/K0144 (2-49 x Janz crown rot QTL = H011/H012; the 2006 supplement cites {10132}, a Kanyuka et al. reference,
   for Collard et al. 2005), K0291 (Pch1 = S144-001), K0176 (QSbv.ksu-5D, kept as K0290). Removed."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {'CGS_2006_10274': {'citation': 'Srnic G, Murphy JP, Lyerly JH, Leath S, Marshall DS (2005) Inheritance and chromosomal assignment of powdery mildew resistance genes in two winter wheat germplasm lines. Crop Science 45: 1578-1586', 'short': 'Srnic et al. (2005)', 'url': ''},
       'CGS_10507': {'citation': 'Friesen TL, Meinhardt SW, Faris JD (2007) The Stagonospora nodorum-wheat pathosystem involves multiple proteinaceous host-selective toxins and corresponding host sensitivity genes that interact in an inverse gene-for-gene manner. The Plant Journal 51: 681-692', 'short': 'Friesen et al. (2007)', 'url': ''},
       'CGS_10588': {'citation': 'Schmolke M, Zimmermann G, Schweizer G, Miedaner T, Korzun V, Ebmeyer E, Hartl L (2008) Molecular mapping of quantitative trait loci for field resistance to Fusarium head blight in a European winter wheat population. Plant Breeding 127: 459-464', 'short': 'Schmolke et al. (2008)', 'url': ''},
       'CGS_10962': {'citation': 'Powell NM, Lewis CM, Berry ST, MacCormack R, Boyd LA (2013) Stripe rust resistance genes in the UK winter wheat cultivar Claire. Theoretical and Applied Genetics 126: 1599-1612', 'short': 'Powell et al. (2013)', 'url': ''}}
RESRC = {'K0122': 'CGS_2006_10274', 'K0270': 'CGS_10507', 'K0271': 'CGS_10507', 'K0272': 'CGS_10507', 'K0275': 'CGS_10507', 'K0429': 'CGS_10507',
         'K0307': 'CGS_10588', 'K0308': 'CGS_10588', 'K0610': 'CGS_10962', 'K0611': 'CGS_10962', 'K0612': 'CGS_10962', 'K0613': 'CGS_10962'}
DEL = ['K0430', 'K0492', 'K0143', 'K0144', 'K0291', 'K0176']
DROPSRC = ['CGS_10275']
T = ' Audit, Nov 2026 (ad60.py):'
# positions beyond the chromosome end
REREF = {'V0448': 'RefSeq v2.1'}   # printed 5158017xx bp is a v2.1 coordinate
REPLACE = {'V1768': ('Xgwm102', 72.3878), 'V1769': ('Xgwm102', 72.3878), 'V2216': ('Xgwm538', 641.4823), 'V2495': ('psp3103', 104.7441)}
UNPLACE = ['V0259', 'V0671', 'V1693', 'V1770', 'V1821', 'V2392', 'V2566', 'R164']
def apply(rs):
    R = {r['id']: r for r in rs}; out = []; changed = []
    for r in rs:
        if r['id'] in DEL: continue
        if r['id'].startswith('K'):
            n = None
            m = re.search(r'\(([A-Za-z0-9_-]+)\) is used\.', r.get('note', ''))
            if m and m.group(1) in R and R[m.group(1)]['ref'] != r['ref']:
                d = R[m.group(1)]; n = dict(r); n['ref'] = d['ref']
                for k in ('start', 'end', 'peak'): n[k] = d.get(k)
                n['note'] = n['note'] + f'{T} coordinates are those of {d["id"]} on {d["ref"]}; reference assembly corrected from RefSeq v1.0.'
            if r['id'] in RESRC:
                n = n or dict(r); old = n['source']; n['source'] = RESRC[r['id']]
                rn = re.search(r'\{(\d+)\}', n['evidence']); rn = rn.group(1) if rn else '?'
                n['note'] = n['note'] + f'{T} source corrected (supplement reference {{{rn}}} is {SRC[RESRC[r["id"]]]["short"]}, not the paper first assigned).'
            if n:
                n['asm'] = lift.lift_record(n) if n.get('grade') != 'D' else n.get('asm'); out.append(n); changed.append(n); continue
        if r['id'] in REREF or r['id'] in REPLACE or r['id'] in UNPLACE:
            n = dict(r); was = '-'.join(f'{n[k]}' for k in ('start', 'end', 'peak') if n.get(k) is not None)
            if r['id'] in REREF:
                n['ref'] = REREF[r['id']]; n['note'] = (n.get('note') or '') + f'{T} the printed coordinate is on RefSeq v2.1; reference assembly corrected.'
            elif r['id'] in REPLACE:
                mk, p = REPLACE[r['id']]; n.update(ref='RefSeq v1.0', start=None, end=None, peak=p, grade='C')
                n['note'] = (n.get('note') or '') + f'{T} the source position ({was} Mb, {r["ref"]}) lies beyond the end of {r["chr"]}; replaced by the in-silico position of {mk} ({p:.3f} Mb, RefSeq v1.0).'
            else:
                n.update(start=None, end=None, peak=None, grade='D')
                n['note'] = (n.get('note') or '') + f'{T} the source position ({was} Mb, {r["ref"]}) lies beyond the end of {r["chr"]} and no listed marker could be placed; left unplaced.'
            n['asm'] = lift.lift_record(n) if n['grade'] != 'D' else {'_note': 'no physical position in source'}
            out.append(n); changed.append(n); continue
        out.append(r)
    return out, changed
if __name__ == '__main__':
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    M = json.load(open('data/merged_v13f.json')); M['reports'], ch = apply(M['reports']); M['sources'].update(SRC)
    for k in DROPSRC: M['sources'].pop(k, None)
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = apply(Sd['reports']); Sd['sources'].update(SRC)
    for k in DROPSRC: Sd['sources'].pop(k, None)
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    json.dump({'reports': ch, 'sources': SRC, 'deleted': DEL}, open('web/seed_ad60.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    print(len(ch), 'changed', sum('reference assembly corrected' in r['note'] for r in ch), 'reref', len(DEL), 'deleted', len(M['reports']))
