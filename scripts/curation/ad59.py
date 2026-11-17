"""ad59 (Nov 2026): Singh et al. (2023) DArT-seq key (SNP__TAG_GBS.xlsx, supplied by the user; M-code -> DArT CloneID and 69-bp tag).
Tags of the 90 M-codes in Tables 2-3 matched to IWGSC RefSeq v1.0 (exact 24-mer seed, then whole-tag mismatch count).
- 69 of 73 records: tag lies at the printed position. Clone ID and tag added to every record.
- M1654: printed 797,177,741 bp; tag at 797,277,741 (1 mismatch): a one-digit misprint, corrected.
- M119: printed 529,608,885; tag at 529,987,697 (1 mismatch); corrected, printed value kept in the note.
- M238 (epistasis): printed 2B 67,247,748; tag at 2B 682,939,231 (no mismatch); corrected.
- M5103: printed 2A 73,179,815; tag matches only 2D 73,592,164 (no mismatch); moved to 2D with the printed 2A noted.
- 17 MTAs on unassigned SNPs (Table 2) added: six placed (tag with at most 6 mismatches at one locus), eleven left unplaced."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); import lift
K = json.load(open('docs/supplements/singh2023_dartseq_key.json'))
T = ' DArT-seq key from the authors (M-code to clone ID and tag), Nov 2026 (ad59.py).'
FIX = {'M1654': ('2B', 797.277741, 'Table 2 prints 797,177,741 bp; the tag lies at 797,277,741 bp (1 mismatch in 69), a one-digit misprint; corrected.'),
       'M119': ('5B', 529.987697, 'Table 2 prints 529,608,885 bp; the tag lies at 529,987,697 bp (1 mismatch in 69); the tag position is used.'),
       'M238': ('2B', 682.939231, 'Table 3 prints 2B 67,247,748 bp; the tag lies at 2B 682,939,231 bp (exact match); the tag position is used.'),
       'M5103': ('2D', 73.592164, 'Table 2 assigns it to 2A (73,179,815 bp), but the tag matches only 2D at 73,592,164 bp (exact, no 2A hit); placed on 2D.')}
BASE = '303 CIMMYT spring wheat genotypes, DArT-seq (12,196 SNPs, 5,400 unassigned), four environments (E1 BHU Varanasi 2017-18, E2 BHU 2018-19, E3 BISA Pusa 2017-18, E4 BISA Pusa 2018-19) and BLUP; MTAs detected by four of five GWAS models (CMLM, MLMM, FarmCPU, BLINK; not SUPER), P < 0.001. Table 2 heads its P column "-log10(P)" but the values are P. '
TR = {'AUDPC': 'AUDPC', 'IP': 'incubation period', 'LN': 'lesion number'}
# unassigned MTAs (Table 2, chromosome NA): trait, environments with P, MAF, placement (chr, bp, mismatches) or candidate loci
UN = [('M7254', 'AUDPC', 'E1 3.09E-04', 0.11495, ('6A', 101362751, 2)), ('M10592', 'AUDPC', 'E1 7.72E-04', 0.062706, 'seed hits 1A 21,638,515 and 1B 32,775,840, whole tag 16 and 14 mismatches'),
      ('M8026', 'AUDPC', 'E3 4.99E-04; E4 4.53E-04', 0.236584, ('6D', 469075133, 5)), ('M11870', 'AUDPC', 'E4 2.70E-04', 0.156617, 'with no match in RefSeq v1.0'),
      ('M976', 'AUDPC', 'E3 3.36E-04; E4 5.54E-05', 0.257228, 'seed hit 1A 323,192,054, whole tag 20 mismatches'),
      ('M7025', 'AUDPC', 'BLUP 7.40E-04', 0.055941, 'seed hit 7D 5,376,610, whole tag 19 mismatches'), ('M9397', 'AUDPC', 'BLUP 5.26E-04', 0.116386, ('4B', 536105824, 1)),
      ('M11338', 'AUDPC', 'BLUP 9.06E-04', 0.459703, 'with no match in RefSeq v1.0'), ('M8531', 'IP', 'E2 9.08E-04', 0.108812, 'with no match in RefSeq v1.0'),
      ('M4837', 'IP', 'E3 1.86E-04', 0.224422, 'seed hit 1A 559,780,823, whole tag 20 mismatches'), ('M4958', 'IP', 'E3 9.23E-04', 0.113449, ('2A', 753061760, 4)),
      ('M9264', 'IP', 'E3 7.62E-04', 0.174719, 'seed hits 6A 4,614,210 and 4,773,825, whole tag 26 and 22 mismatches'), ('M12439', 'IP', 'E4 8.90E-04', 0.183828, 'with no match in RefSeq v1.0'),
      ('M7868', 'IP', 'BLUP 3.70E-04', 0.156155, ('7D', 20237435, 4)), ('M9994', 'LN', 'E1 4.57E-04', 0.306931, 'seed hit 2A 769,631,158, whole tag 18 mismatches'),
      ('M992', 'LN', 'E2 7.62E-05', 0.145017, ('1A', 58188105, 6)), ('M9263', 'LN', 'E3 7.11E-04', 0.440545, 'seed hits 1B 587,982,742, 1B 588,284,359 and 1D 434,594,788, whole tag 26-31 mismatches')]
def seqs(marker):
    return [{'m': f'{m} (DArT {K[m][0]})', 's': K[m][1]} for m in [x.strip() for x in marker.split(',')] if m in K]
def apply(rs):
    out = []; upd = []
    for r in rs:
        if r['source'] != 'Singh2023FPS': out.append(r); continue
        r = dict(r); m = r['marker'].split(',')[0].strip()
        r['snpseq'] = seqs(r['marker'])
        nt = r['note'].replace(' Read from the paper, Nov 2026 (ad58.py).', '').replace(T, '')
        if m in FIX:
            ch, pk, why = FIX[m]; r['chr'] = ch; r['peak'] = pk; r['grade'] = 'B'; nt += ' ' + why
        else:
            nt += ' The DArT tag lies at the printed position (IWGSC RefSeq v1.0).'
        r['note'] = nt + ' Read from the paper, Nov 2026 (ad58.py).' + T
        r['asm'] = lift.lift_record(r); out.append(r); upd.append(r)
    have = {r['id'] for r in out}; n = max([int(i.split('-')[1]) for i in have if i.startswith('S205-')] + [0])
    for m, tr, ev, maf, pl in UN:
        n += 1; i = f'S205-{n:03d}'
        if i in have: continue
        P = min(float(x.split()[1]) for x in ev.split('; '))
        r = dict(id=i, cls='MTA', curation='hand', name=f'{m} ({tr}, unassigned SNP)', disease='Spot blotch', ref='RefSeq v1.0', marker=m, stat=f'P {P:.2e}', pve='',
                 source='Singh2023FPS', evidence=f'Table 2: {m} chromosome NA ({ev})', start=None, end=None, peak=None, snpseq=seqs(m))
        nt = BASE + f'Trait {TR[tr]}; significant in {ev.replace(" ", " (P ").replace(";", ");")}); MAF {maf:.3f}. Chromosome not assigned by the vendor (NA in Table 2).'
        if isinstance(pl, tuple):
            ch, bp, mm = pl; r.update(chr=ch, peak=bp / 1e6, grade='B' if mm <= 2 else 'C')
            nt += f' Placed in silico from the DArT tag: {ch} {bp:,} bp ({mm} mismatches in 69 bp, single locus).'
        else:
            r.update(chr='unm', grade='D'); nt += f' Not placed: DArT tag {pl}.'
        r['note'] = nt + T
        r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
        out.append(r); upd.append(r)
    return out, upd
if __name__ == '__main__':
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    M = json.load(open('data/merged_v13f.json')); M['reports'], upd = apply(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = apply(Sd['reports'])
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    json.dump({'reports': upd, 'sources': {}, 'deleted': []}, open('web/seed_ad59.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    print(len(upd), sum(1 for r in upd if r['id'] >= 'S205-018'), [ (r['id'],r['name'],r['chr'],r.get('peak'),r['grade']) for r in upd if 'unassigned' in r['name'] or r['grade']=='B'])
