"""ad61 (Nov 2026): Tomar et al. (2021) Front Genet 11:613217 (Tomar2021FG) re-curated from the paper (user query on QSb.bisa-1A).
Problems: names lacked the dot used in the paper (Q.Sb.bisa-); R091-R094 had no marker, P or effect; R093 (2B) and R094 (6D) used the
abstract positions (421.92 and 6.84 Mb) instead of Table 2; R0238-R0241 carried 'Effect 10X' (allele effect) in the PVE field (the paper
reports no PVE). Now each of the eight QTL spans its Table 2 SNPs (RefSeq v1.0, GBS), with the peak SNP, environments, P, MAF and allele
effects; nearest genes from Table 3."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); import lift
B = ('139 advanced CIMMYT spring wheat lines, GBS SNPs called on RefSeq v1.0, FarmCPU GWAS; Env1 Pusa 2019, Env2 Pusa 2020 (India), '
     'Env3 BARI Jamalpur 2017 (Bangladesh), natural infection; P < 0.003 in Table 2. The paper gives allele effects on disease severity '
     '(marker effect and "Effect 10X"), not PVE. ')
T = ' Re-read from the paper, Nov 2026 (ad61.py).'
Q = [  # id, chr, snps(list of (snp, Mb)), peak snp, envs, best P, maf, effect text, gene note, extra
 ('R091', '1A', [('S1A_497201550', 497.201550), ('S1A_497201682', 497.201682)], 'S1A_497201550', 'Env1, Env2', 0.001232, '0.28',
  'C allele +1.88 (Effect 10X 18.83) in Env1, +1.83 in Env2', 'Both SNPs lie 9.4 kb from TraesCS1A01G303600 (LURP-one-like protein); TraesCS1A01G303700 (GTP cyclohydrolase 1) 23 kb away.', 'Reported as a new QTL.'),
 ('R0238', '1B', [('S1B_636840957', 636.840957)], 'S1B_636840957', 'Env1, Env2', 0.001754, '0.37',
  'A allele +1.61 (Effect 10X 16.14) in Env1, +1.63 in Env2', 'Nearest genes TraesCS1B01G409800 (60S ribosomal protein L35a-like) and TraesCS1B01G409900 (transmembrane protein).', ''),
 ('R092', '1D', [('S1D_89835681', 89.835681)], 'S1D_89835681', 'Env1, Env2', 0.001896, '0.13',
  'A allele +2.40 (Effect 10X 23.95) in Env1, +2.40 in Env2', 'Nearest genes TraesCS1D01G101800 (pre-rRNA-processing protein TSR2, 238 kb) and TraesCS1D01G101900 (glucan 1,3-beta-glucosidase, 313 kb).', 'Reported as a new QTL; minor allele effective.'),
 ('R0239', '2A', [('S2A_703111105', 703.111105), ('S2A_703358397', 703.358397), ('S2A_703391915', 703.391915), ('S2A_703391992', 703.391992), ('S2A_703427639', 703.427639), ('S2A_704446408', 704.446408)],
  'S2A_703111105', 'Env1, Env2', 0.002881, '0.16', 'C allele +2.28 (Effect 10X 22.83) in Env2, +2.22 in Env1', 'Six SNPs among senescence-associated family protein genes (TraesCS2A01G453900 to 454700); S2A_704446408 near a lectin receptor kinase (TraesCS2A01G455600).', ''),
 ('R093', '2B', [('S2B_419320960', 419.320960), ('S2B_419456700', 419.456700), ('S2B_420723687', 420.723687), ('S2B_421708152', 421.708152), ('S2B_422031949', 422.031949), ('S2B_422190984', 422.190984),
                 ('S2B_422191033', 422.191033), ('S2B_422507900', 422.507900), ('S2B_422983662', 422.983662), ('S2B_423094651', 423.094651), ('S2B_423836280', 423.836280)],
  'S2B_422983662', 'Env1, Env2, Env3', 0.001657, '0.11 (S2B_422983662); 0.45-0.46 (Env3 SNPs)', 'A allele of S2B_422983662 +3.01 (Effect 10X 30.07) in Env1, +3.07 (30.68) in Env2, the largest effect in the study; ten Env3 SNPs about 2.3',
  'Genes in the interval include a MADS-box transcription factor region, isocitrate lyase (TraesCS2B01G300200) and an auxin response factor (TraesCS2B01G300300).', 'Detected in all three environments; reported as a new QTL. The abstract gives 421.92 Mb.'),
 ('R0240', '4A', [('S4A_725538462', 725.538462), ('S4A_725660945', 725.660945), ('S4A_725880148', 725.880148)], 'S4A_725660945', 'Env1, Env2', 0.000855, '0.21-0.24',
  'A allele of S4A_725660945 +2.02 (Effect 10X 20.18) in Env1, +1.89 in Env2', 'S4A_725660945 near TraesCS4A01G460100 (cytochrome P450).', ''),
 ('R0241', '5B', [('S5B_682958475', 682.958475), ('S5B_683240735', 683.240735)], 'S5B_683240735', 'Env1, Env3', 0.002373, '0.07 (S5B_683240735); 0.46 (S5B_682958475)',
  'S5B_683240735 G allele -3.14 (Effect 10X -31.37) in Env1; S5B_682958475 A allele +2.14 in Env3', '', 'Minor allele effective.'),
 ('R094', '6D', [('S6D_6395796', 6.395796), ('S6D_6944636', 6.944636)], 'S6D_6395796', 'Env1, Env2', 0.002335, '0.27',
  'T allele -1.78 (Effect 10X -17.77) in Env2, -1.73 in Env1', '', 'Reported as a new QTL. The abstract gives 6.84 Mb; Table 2 gives 6.40 and 6.94 Mb.')]
def recs(old):
    O = {r['id']: r for r in old}; out = []
    for i, ch, snps, pk, env, p, maf, eff, gene, extra in Q:
        r = dict(O[i]); pos = dict(snps)
        r.update(name=f'Q.Sb.bisa-{ch}', cls='MTA', chr=ch, ref='RefSeq v1.0', grade='A', curation='hand', pve='', stat=f'P {p}', effect=eff,
                 marker=', '.join(s for s, _ in snps), peak=round(pos[pk], 6), start=round(min(pos.values()), 6) if len(snps) > 1 else None,
                 end=round(max(pos.values()), 6) if len(snps) > 1 else None, evidence=f'Table 2: {len(snps)} SNP(s), peak {pk} at {pos[pk]:.2f} Mb; Table 3 genes')
        r['note'] = B + f'Significant in {env}; peak SNP {pk}; MAF {maf}. Allele effect: {eff}. ' + (gene + ' ' if gene else '') + (extra + ' ' if extra else '') + 'Earlier record lacked markers and effects' + (' and stored the allele effect as PVE' if i.startswith('R02') else '') + '.' + T
        r['note'] = r['note'].replace('  ', ' ')
        r['asm'] = lift.lift_record(r); out.append(r)
    return out
def apply(rs):
    new = recs(rs); N = {r['id']: r for r in new}
    return [N.get(r['id'], r) for r in rs], new
if __name__ == '__main__':
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    M = json.load(open('data/merged_v13f.json')); M['reports'], new = apply(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = apply(Sd['reports'])
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    json.dump({'reports': new, 'sources': {}, 'deleted': []}, open('web/seed_ad61.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new: print(r['id'], r['name'], r['chr'], r.get('start'), r.get('end'), r['peak'], r['stat'])
