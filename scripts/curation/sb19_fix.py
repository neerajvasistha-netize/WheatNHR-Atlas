"""Rebuild the loci of 19 spot blotch / leaf spot source papers from the papers themselves (Nov 2026).

Every record below was read from the paper's own results table. Physical positions:
  grade A  position printed in the paper (RefSeq v1.0 or v2.0, as stated by the paper)
  grade B  position derived from the paper's own markers/genes on the stated chromosome
           (exact-match in-silico PCR on RefSeq v1.0, 90K/DArT positions, gene coordinates)
  grade C  only one flanking/linked marker could be placed, or the placement relies on
           another paper's projection of the same marker
  grade D  no marker of the locus could be placed on the stated chromosome
SSR markers are placed by exact-match in-silico PCR of their primers on RefSeq v1.0
(work/ispcr_res.json); a marker is used only if it amplifies on the chromosome the paper states.
Run from repository root after work/ispcr_res.json exists:  python scripts/curation/sb19_fix.py
"""
import json, re, gzip, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift
W = 'scripts/curation/sb19/'
PCR = json.load(open(W + 'ispcr_res.json'))
MK = json.load(gzip.open('data/markerpos_all.json.gz'))
LOG = []
# positions from an earlier atlas session (primer search with mismatches); exact-match PCR here gives no product
LEGACY = {('barc175', '6D'): 411.884, ('barc121', '7D'): 531.586}

def nk(m):
    k = re.sub(r'[^a-z0-9]', '', m.lower())
    if k.startswith('x') and not k.startswith('xfcp'): k = k[1:]
    if k.startswith('xfcp'): k = k[1:]
    k = re.sub(r'^wms', 'gwm', k)
    return re.sub(r'^(gwm|wms|barc|wmc|cfd|cfa|gdm|swm)0+(\d)', r'\1\2', k)

def pos(m, ch, near=None):
    """RefSeq v1.0 Mb of marker m on chromosome ch; returns (Mb, how) or (None, why)."""
    k = nk(m)
    hits = [h for h in PCR.get(k, []) if h[0] == ch]
    if hits:
        if near is not None: hits.sort(key=lambda h: abs(h[1] - near))
        if len(hits) > 1 and near is None and max(h[1] for h in hits) - min(h[1] for h in hits) > 5:
            return None, '%s amplifies at %d sites on %s' % (m, len(hits), ch)
        return hits[0][1], 'in-silico PCR'
    if PCR.get(k):   # primers amplify, but not on this chromosome: do not trust a dictionary copy
        return None, '%s amplifies only on %s' % (m, ','.join(sorted({h[0] for h in PCR[k]})))
    ent = MK.get(k, [])
    prim = {c for c, p, ref, src in ent if src == 'Zhao et al. 2019 J Genet'}
    for c, p, ref, src in ent:
        if c == ch and ref == 'RefSeq v1.0' and 'atlas marker dictionary' not in src:
            if '(Table S1)' in src and prim and ch not in prim:
                return None, '%s: primary locus on %s; %s copy is secondary' % (m, ','.join(sorted(prim)), ch)
            return round(p, 4), src
    if (k, ch) in LEGACY: return LEGACY[(k, ch)], 'earlier atlas in-silico PCR, not reproduced by exact match'
    for c, p, ref, src in MK.get(k, []):
        if c == ch and ref == 'RefSeq v2.1':
            y = lift.lift_record({'chr': ch, 'peak': p, 'ref': 'RefSeq v2.1'})['RefSeq v1.0'].get('peak')
            return y, '%s (RefSeq v2.1 %.3f Mb, converted to v1.0)' % (src, p)
    return None, 'not placed on ' + ch

def place(rec, markers, ch, peak_marker=None):
    """Fill start/end (and peak if a peak marker is named) from marker positions."""
    got = [];
    for m in markers:
        p, how = pos(m, ch)
        if p is not None: got.append((m, p, how))
        else: LOG.append((rec['name'], m, how))
    hows = sorted({h for _, _, h in got})
    if len(got) >= 2:
        ps = sorted(p for _, p, _ in got); rec['start'], rec['end'] = ps[0], ps[-1]; rec['grade'] = rec.get('grade') or 'B'
    elif len(got) == 1:
        rec['peak'] = got[0][1]; rec['grade'] = 'C'
        rec['note'] = (rec.get('note', '') + ' Only %s could be placed on %s; position is that marker, not a reported peak.' % (got[0][0], ch)).strip()
    else:
        rec['grade'] = 'D'
        rec['note'] = (rec.get('note', '') + ' No flanking marker could be placed on %s; left unplaced.' % ch).strip()
    if peak_marker:
        p, how = pos(peak_marker, ch)
        if p is not None and len(got) >= 2: rec['peak'] = p
    if got:
        rec['note'] = (rec.get('note', '') + ' Marker positions (RefSeq v1.0): ' + '; '.join('%s %.3f Mb (%s)' % (m, p, h) for m, p, h in got) + '.').strip()
    return rec

R = []   # corrected records
def add(**k):
    k.setdefault('ref', 'RefSeq v1.0'); k.setdefault('curation', 'hand'); k.setdefault('cls', 'QTL')
    k.setdefault('start', None); k.setdefault('end', None); k.setdefault('peak', None)
    k.setdefault('pve', ''); k.setdefault('marker', ''); k.setdefault('stat', ''); k.setdefault('note', ''); k.setdefault('grade', None)
    R.append(k); return k

# ---------------- Gahtyari et al. 2021 Plants 10:973, Table 4 (Mb on CS RefSeq v1.0 printed) ----------------
G = [('WC','1A','76.5-88.5','2.8-8.0','1218247','1096735',500.8,512.5,'4.9-9.9','WUYA','co-located with PH QTL'),
     ('WC','1B','175.5-178','2.8-4.8','1137809','5411162',670.6,673.7,'5.1-5.4','WUYA','matches Lr46 region'),
     ('WC','2A','37.6-40.0','2.8-6.2','3020873','1132194',28.8,30.9,'2.4-2.6','WUYA','co-located with PH QTL'),
     ('WC','2A','90.4-98.0','2.7-4.3','4991898','3064660',73.5,101.0,'2.8-8.2','WUYA','differed between resistant and susceptible bulks'),
     ('WC','2D','71.5-73.5','3.5-5.2','2261906','1085831',61.8,65.7,'5.0-5.3','WUYA','differed between resistant and susceptible bulks'),
     ('WC','4B','39.5-41.5','4.8-5.7','985312','1241652',37.5,65.1,'5.5-13.7','CIANO T79','PH QTL / Rht-B1 region'),
     ('WC','4B','57.5-60.5','2.5-3.4','1159447','13375761',558.6,572.1,'3.7-4.8','CIANO T79','DRSB; near Sb4'),
     ('WC','4B','92.5-96.5','4.1-4.4','1195526','100167524',658.0,659.5,'4.5-4.9','WUYA',''),
     ('WC','5A','139.3-156.1','3.1-7.1','1102419','5411712',570.2,584.6,'5.8-8.5','WUYA','DH QTL / Vrn-A1 region'),
     ('WC','5D','38.5-51.5','7.1','1058378','1048778',458.8,542.6,'15.0','WUYA',''),
     ('WC','7B','8.5-12.5','2.7-4.2','7341261','1007925',18.1,18.1,'3.0-5.1','WUYA',''),
     ('KC','2A','48.1-48.6','2.5-3.3','1127593','3034005',7.9,8.0,'3.2-5.3','CIANO T79','DRSB'),
     ('KC','2A','224.8-334.0','4.5','1098973','1013641',16.8,65.1,'7.6','KATH',''),
     ('KC','2B','13.7-22.0','2.6-2.8','1218621','4992694',24.1,27.2,'4.3-5.4','CIANO T79','DRSB'),
     ('KC','2B','37.1-43.0','2.6-3.0','1066619','1278607',56.9,454.8,'4.2-6.5','CIANO T79','DH QTL; 398 Mb span as printed'),
     ('KC','4B','31.8-48.0','3.1-5.7','3222467','1233562',21.6,483.8,'7.8-11.5','CIANO T79','PH QTL / Rht-B1; 462 Mb span as printed'),
     ('KC','4B','56.5-57.0','3.8','1863994','1132777',535.1,548.2,'6.6','CIANO T79','near Sb4'),
     ('KC','4D','0-30.5','2.9-4.7','Rht-D1','BS00036421_51',18.8,32.3,'2.3-9.3','KATH','PH QTL / Rht-D1'),
     ('KC','4D','85.5-100.5','4.1','2257171','2256312',457.5,502.6,'7.8','CIANO T79','')]
popn = {'WC': 'WUYA x CIANO T79', 'KC': 'KATH x CIANO T79'}
seen = {}
for p, ch, cm, lod, l, r, s, e, pve, src, rem in G:
    seen[(p, ch)] = seen.get((p, ch), 0) + 1
    tag = '' if sum(1 for x in G if x[0] == p and x[1] == ch) == 1 else '.%d' % seen[(p, ch)]
    add(name='QSb-%s-%s%s (Gahtyari 2021)' % (p, ch, tag), disease='Spot blotch', chr=ch, start=s, end=e,
        marker='%s–%s' % (l, r), pve=pve, stat='LOD ' + lod, source='Gahtyari2021', grade='A',
        note='%s RIL, %s cM; resistance allele from %s. %s Physical range printed in Table 4 (RefSeq v1.0).' % (popn[p], cm, src, (rem + '.') if rem else ''),
        evidence='Table 4: %s %s %s cM LOD %s %s–%s %s–%s Mb PVE %s' % (p, ch, cm, lod, l, r, s, e, pve))

# ---------------- He et al. 2020 Phytopathology 110:1980-1987, Table 3 (RefSeq v1.0 printed) ----------------
H = [('CC','1B','','210.8-215.6',619.6,622.7,'1042822','1117717','2.2-3.5','CIANO T79',''),
     ('CC','2B','','121.5-143.1',52.7,73.3,'1228793','1167654','2.8-3.0','CIANO T79',''),
     ('CC','4A','','25.5-35.3',620.8,624.1,'3023781','992456','2.8-4.6','CASCABEL',''),
     ('CC','4B','','72.6-98.9',30.9,622.4,'Rht-B1','2279882','2.4-6.3','CASCABEL','associated with DH and PH; 592 Mb span as printed'),
     ('CC','5A','','44.2-48.6',582.6,587.4,'2341646','Vrn-A1','14.1-22.0 (mean 19.4)','CASCABEL','Vrn-A1 region, associated with DH'),
     ('CC','5B','','73.9-77.7',412.7,431.9,'996745','10592866','8.9-21.4 (mean 17.2)','CASCABEL','major QTL; text projects 406.9-431.9 Mb'),
     ('KC','1B','','149.0-151.8',673.6,674.5,'2263671','100193437','3.1-4.6','KATH',''),
     ('KC','2B','','66.2-71.2',65.4,68.6,'1017820','3947227','4.2-15.0 (mean 8.4)','CIANO T79',''),
     ('KC','4B','.1','88.7-91.3',609.5,619.1,'4989130','1114716','3.4-8.0','CIANO T79',''),
     ('KC','4B','.2','94.4-115.4',613.1,652.4,'998971','1095421','4.5-5.6','KATH',''),
     ('KC','4D','','0.0-17.5',32.3,336.8,'BS00036421_51','1059032','2.3-4.1','KATH','associated with PH'),
     ('KC','5A','','174.6-188.2',570.5,591.2,'987242','IWA4449','13.5-25.9 (mean 25.6)','KATH','Vrn-A1 region, associated with DH')]
hp = {'CC': 'CASCABEL x CIANO T79', 'KC': 'KATH x CIANO T79'}
for p, ch, t, cm, s, e, l, r, pve, src, rem in H:
    add(name='QSb-%s-%s%s (He 2020)' % (p, ch, t), disease='Spot blotch', chr=ch, start=s, end=e, marker='%s–%s' % (l, r),
        pve=pve, source='He2020Phyto', grade='A',
        note='%s RIL, %s cM; resistance allele from %s. %sPhysical range printed in Table 3 (RefSeq v1.0).' % (hp[p], cm, src, (rem + '. ') if rem else ''),
        evidence='Table 3: %s %s %s cM %s–%s Mb %s %s PVE %s' % (p, ch, cm, s, e, l, r, pve))
add(name='5BL spot blotch QTL cluster (He 2020)', disease='Spot blotch', chr='5B', start=404.2, end=418.8,
    marker='wsnp_Ku_c50354_55979952–gwm213', source='He2020Phyto', grade='A', cls='QTL',
    note='Cluster of 5BL spot blotch QTL (CASCABEL x CIANO T79, Yangmai 6, Ning 8201) delimited by He et al. (2020) in the text; RefSeq v1.0 as printed.',
    evidence='Text: hotspot delimited by wsnp_Ku_c50354_55979952 and gwm213, 404.2 to 418.8 Mb')

# ---------------- Kaur et al. 2021 Front Plant Sci 12:650400, Table 2 (Mb printed, RefSeq v1.0) ----------------
for nm, mk, ch, p, gs, lod, pve, donor in [
        ('QSb.pau-2A', 'S2A_755774702', '2A', 755.7747, 'HDS', '3.12-3.18', '18.44-18.77', 'Ae. speltoides'),
        ('QSb.pau-2B', 'S2B_673595704', '2B', 673.5957, 'AUDPC, HDS', '3.04-5.02', '16.03-25.56', 'PDW274'),
        ('QSb.pau-3B', 'S3B_104700872', '3B', 104.7009, 'AUDPC, HDS', '3.22-4.58', '17.86-25.27', 'Ae. speltoides'),
        ('QSb.pau-5B', 'S5B_703858864', '5B', 703.8589, 'DS', '3.23-3.32', '18.97-19.48', 'Ae. speltoides'),
        ('QSb.pau-6A', 'S6A_131743987', '6A', 131.744, 'FS', '3.07-3.08', '16.38-16.42', 'PDW274')]:
    n = 'T. durum-Ae. speltoides introgression lines (PDW274 background); trait %s; resistance allele from %s. Peak SNP position printed in Table 2.' % (gs, donor)
    if ch == '3B': n += ' Table 2 gives S3B_104700872, Table 3 gives S3B_104700839 for the same QTL (33 bp apart).'
    add(name=nm, disease='Spot blotch', chr=ch, peak=p, marker=mk, pve=pve, stat='LOD ' + lod, source='Kaur2021FPS', grade='A', note=n,
        evidence='Table 2: %s %s %s %.2f Mb LOD %s PVE %s' % (nm, mk, ch, p, lod, pve))

# ---------------- Kumar et al. 2009 TAG 118:783-792, Table 4 (Yangmai 6 x Sonalika; SSR) ----------------
for nm, l, r, arm, pve, lod in [('QSb.bhu-2A', 'Xbarc353', 'Xgwm445', '2A', '13.79-15.58', '2.8-3.8'),
                                ('QSb.bhu-2B', 'Xgwm148', 'Xgwm374', '2B', '8.72-20.50', '2.4-5.8'),
                                ('QSb.bhu-5B', 'Xgwm067', 'Xgwm371', '5B', '11.76-41.10', '3.4-12.8'),
                                ('QSb.bhu-6D', 'Xbarc175', 'Xgwm732', '6D', '14.61-22.50', '2.5-3.3')]:
    place(add(name=nm, disease='Spot blotch', chr=arm, marker='%s–%s' % (l, r), pve=pve, stat='LOD ' + lod, source='Kumar2009TAG',
              note='Yangmai 6 x Sonalika RIL; CIM.', evidence='Table 4: %s %s-%s' % (nm, l, r)), [l, r], arm)

# ---------------- Kumar et al. 2010 Mol Breed 26:477-491, Tables 4 and 6 ----------------
for pop, nm, l, r, ch, pve, lod in [
        ('Ning 8201 x Sonalika', 'QSb.bhu-2A', 'Xgwm425', 'Xbarc159', '2A', '15.2-22.7', '2.9-7.2'),
        ('Ning 8201 x Sonalika', 'QSb.bhu-2B', 'Xgwm148', 'Xbarc91', '2B', '15.1-23.7', '2.8-4.3'),
        ('Ning 8201 x Sonalika', 'QSb.bhu-5B', 'Xgwm067', 'Xgwm213', '5B', '10.7-16.9', '2.5-2.7'),
        ('Ning 8201 x Sonalika', 'QSb.bhu-7D', 'Xgwm111', 'Xgwm1168', '7D', '25.3-51.3', '3.5-6.9'),
        ('Chirya 3 x Sonalika', 'QSb.bhu-2B', 'Xgwm148', 'Xgwm129', '2B', '13.1-18.8', '3.1-3.7'),
        ('Chirya 3 x Sonalika', 'QSb.bhu-2D', 'Xgwm455', 'Xgwm815', '2D', '6.5-10.7', '2.2-3.5'),
        ('Chirya 3 x Sonalika', 'QSb.bhu-3B', 'Xgwm533', 'Xgwm1037', '3B', '9.5-9.7', '3.0-3.1'),
        ('Chirya 3 x Sonalika', 'QSb.bhu-7B', 'Xgwm263', 'Xgwm255', '7B', '7.5-13.1', '2.3-4.1'),
        ('Chirya 3 x Sonalika', 'QSb.bhu-7D', 'Xgwm111', 'Xswm008', '7D', '10.4-17.1', '3.5-4.8')]:
    tag = ' (Ning 8201)' if pop.startswith('Ning') else ' (Chirya 3)'
    place(add(name=nm + tag, disease='Spot blotch', chr=ch, marker='%s–%s' % (l, r), pve=pve, stat='LOD ' + lod, source='Kumar2010MB',
              note='%s RIL; CIM.' % pop, evidence='Table %s: %s %s-%s' % ('4' if tag == ' (Ning 8201)' else '6', nm, l, r)), [l, r], ch)

# ---------------- Kumar et al. 2015 Mol Breed 35:218: Sb2 ----------------
place(add(name='Sb2', cls='DG', disease='Spot blotch', chr='5B', marker='Xgwm639–Xgwm1043 (0.62 cM)', source='Kumar2015MB',
          note='QSb.bhu-5B Mendelized in YS116 x Sonalika and designated Sb2; 5BL1-0.55-0.75 deletion bin. Xgwm1043 primers not available to the atlas.',
          evidence='Abstract: Sb2 flanked by Xgwm639 (0.28 cM) and Xgwm1043'), ['Xgwm639', 'Xgwm1043'], '5B')
# ---------------- Kumar et al. 2016 Mol Breed 36:95 (NDVI) ----------------
place(add(name='Sb2 region QTL, NDVI and %DS (Kumar 2016)', disease='Spot blotch', chr='5B', marker='Xgwm639–Xgwm1043',
          pve='51.3-56.1', stat='LOD up to 40.5', source='Kumar2016MB',
          note='YS116 x Sonalika; the same Sb2 locus mapped with NDVI and visual %DS; Tsn1-linked Xfcp1 lies 1.1 cM away (Xfcp1 at 546.15 Mb).',
          evidence='Results: QTL between Xgwm639 and Xgwm1043, PVE up to 54.9% (NDVI) and 56.1% (%DS)'), ['Xgwm639', 'Xgwm1043'], '5B')

# ---------------- Lillemo et al. 2013 TAG 126:711-719 ----------------
add(name='Sb1 (Lr34 locus)', cls='DG', disease='Spot blotch', chr='7D', start=47.412, end=47.4245, marker='gwm1220–swm10; csLV34–gwm295',
    pve='17.1-55.2', stat='LOD 4.4-18.3', source='Lillemo2013', grade='B',
    note='Avocet x Saar RIL; major QTL at Lr34 given the designation Sb1. Placed at the Lr34 gene TraesCS7D02G080300 (RefSeq v1.0 7D:47,412,062-47,424,490, Ensembl).',
    evidence='Table 3: 7DS 81 cM gwm1220-swm10 LOD 18.3 R2 55.2 (mean)')
add(name='QSb.7DS minor, Lr34-susceptible subpopulation (Lillemo 2013)', disease='Spot blotch', chr='7D', marker='wPt-7654–gdm88',
    pve='22.6-40.4', stat='LOD 3.12-5.81', source='Lillemo2013', note='Detected in the 52 RILs carrying the susceptible Lr34 allele.',
    evidence='Table 4: 7DS 42-45 cM wPt-7654-gdm88')
place(R[-1], ['gdm88'], '7D')
place(add(name='QSb.7DS minor 2 (Lillemo 2013)', disease='Spot blotch', chr='7D', marker='gdm145–barc352', pve='21.9', stat='LOD 4.59',
          source='Lillemo2013', note='Lr34-susceptible subpopulation, 2007 date 1 only.', evidence='Table 4: 2007 date 1 7DS 59 cM gdm145-barc352 LOD 4.59 R2 21.9'),
      ['gdm145', 'barc352'], '7D')
add(name='QSb.1BL at Lr46 (Lillemo 2013)', disease='Spot blotch', chr='1B', start=670.6, end=673.7, marker='hbe248–ncw1-V (Lr46)', pve='15.1', stat='LOD 3.78',
    source='Lillemo2013', grade='C',
    note='Minor QTL at the Lr46 locus in the Lr34-susceptible subpopulation. hbe248 and ncw1-V could not be placed; interval is the Lr46 region as projected by Gahtyari et al. (2021, 670.6-673.7 Mb RefSeq v1.0).',
    evidence='Table 4: 2007 date 1 1BL 120 cM hbe248-ncw1-V LOD 3.78 R2 15.1')

# ---------------- Lu et al. 2016 TAG 129:577-589: Sb3 ----------------
sb3 = add(name='Sb3', cls='DG', disease='Spot blotch', chr='3B', marker='XWGGC5911–XWGGC4320 (co-segregating XWGGC3957); fine map 0.15 cM, 602 kb',
          source='Lu2016', note='621-7-1; dominant gene on 3BS.', evidence='Abstract: Sb3 mapped in a 0.15 cM interval spanning a 602 kb region of CS 3BS')
place(sb3, ['wggc12530', 'wggc9893'], '3B', peak_marker='wggc3957')

# ---------------- Pankaj et al. 2022 Australas Plant Pathol 51:441-452, Table 3 (SSR) ----------------
PK = [('2017-18','Qsb_rpcau_2A','2A','Barc212','wmc382','2.61','15.32','PBW343'),
      ('2017-18','Qsb_rpcau_5B','5B','Cfd71','wmc173','11.16','18.75','PBW343'),
      ('2017-18','Qsb_rpcau_1B','1B','gwm131','gwm498','4.37','12.18','PBW343'),
      ('2017-18','Qsb_rpcau_5A','5A','Gwm304','barc141','7.30','9.81','PBW343'),
      ('2017-18','Qsb_rpcau_1D','1D','cfa2129','cfd48','8.04','8.04','PBW343'),
      ('2017-18','Qsb_rpcau_2D','2D','wmc503','gwm261','7.42','12.01','IC252874'),
      ('2017-18 and 2018-19','Qsb_rpcau_2B','2B','Wmc109','gwm312','6.04-6.87','10.32-15.34','IC252874'),
      ('2018-19','Qsb_rpcau_4A','4A','Cfd257','wmc718','8.21','10.26','PBW343'),
      ('2018-19','Qsb_rpcau_4D','4D','Wmc48','cfd106','4.53','7.43','PBW343'),
      ('2018-19','Qsb_rpcau_2A','2A','Gwm359','wmc728','5.09','8.04','IC252874'),
      ('2017-18 and 2018-19','Qsb_rpcau_5B','5B','Cfd71','wmc173','15.39','21.12','PBW343'),
      ('2018-19','Qsb_rpcau_1B','1B','wmc611','cfd59','6.71','12.72','PBW343'),
      ('2018-19','Qsb_rpcau_3B','3B','Barc68','gwm72','4.2','11.53','PBW343'),
      ('2018-19','Qsb_rpcau_5A','5A','Wmc492','barc40','7.25','16.25','PBW343'),
      ('2018-19','Qsb_rpcau_6A','6A','Gwm570','wmc553','8.04','12.52','PBW343'),
      ('2018-19','Qsb_rpcau_2A','2A','Gwm614','gwm497','9.67','22.13','PBW343'),
      ('2018-19','Qsb_rpcau_4D','4D','Gwm213','gwm608','10.36','14.63','IC252874')]
done = set()
for yr, nm, ch, l, r, lod, pve, al in PK:
    if (nm, l) in done:  # same QTL, same interval in both years: merge into the first record
        for x in R:
            if x['source'] == 'Z_1329' and x['name'].startswith(nm) and x['marker'].startswith(l):
                x['pve'] = x['pve'] + ', ' + pve; x['stat'] = x['stat'] + ', ' + lod; x['note'] = x['note'].replace('2017-18;', '2017-18 and 2018-19;')
        continue
    done.add((nm, l))
    k = sum(1 for x in R if x['source'] == 'Z_1329' and x['name'].startswith(nm))
    place(add(name=nm + ('' if k == 0 else ' (%s–%s)' % (l, r)), disease='Spot blotch', chr=ch, marker='%s–%s' % (l, r), pve=pve, stat='LOD ' + lod,
              source='Z_1329', note='PBW343 x IC252874 DH population; %s; resistance allele from %s.' % (yr, al),
              evidence='Table 3: %s %s %s-%s LOD %s PVE %s' % (yr, nm, l, r, lod, pve)), [l, r], ch)

# ---------------- Roy et al. 2021 Front Plant Sci 12:641324, Table 3 (DArTseq) ----------------
RY = [('BC','1A','153.79-155.61','4989967','1026215','2.31-4.01','BARTAI'),
      ('BC','4D','70.49-90.31','BS00036421_51','1119387','9.53-12.23','BARTAI'),
      ('BC','4D','112.48-131.93','12002205','1072422','4.49','BARTAI'),
      ('BC','5B','135.84-138.59','9724385','2267710','4.75','BARTAI'),
      ('CC','1B','261.53-263.82','1168776','1037914','2.79-3.11','CIANO T79'),
      ('CC','5A','331.49-332.06','1067537','2257572','8.93-12.62 (mean 10.32)','CASCABEL'),
      ('CC','5A','472.56-481.59','1218172','1683258','4.13-5.89','CIANO T79'),
      ('CC','5B','522.08-525.15','3958735','1137742','4.02-5.77','CASCABEL'),
      ('CC','6A','5.08-61.14','1125980','100193832','5.24-6.48','CASCABEL'),
      ('CC','7A','328.34-367.52','1126352','1208614','5.19-5.45','CIANO T79'),
      ('CC','7B','134.52-173.55','1125523','1007745','4.18-4.56','CIANO T79')]
OLD_ROY = {}
DS = 'Flanking DArTseq tag sequences from Roy et al. (2021) Data Sheet 1 Table S4 located on RefSeq v1.0 (earlier atlas curation).'
ROYDS = {('BC','1A','153.79-155.61'): ({'start': 512.492, 'end': 516.584, 'peak': None, 'grade': 'B'}, DS),
         ('BC','4D','112.48-131.93'): ({'start': None, 'end': None, 'peak': 425.973, 'grade': 'C'}, DS + ' Only one flank placed.'),
         ('BC','5B','135.84-138.59'): ({'start': 455.741, 'end': 458.387, 'peak': None, 'grade': 'B'}, DS),
         ('CC','6A','5.08-61.14'): ({'start': None, 'end': None, 'peak': 17.659, 'grade': 'C'}, DS + ' Only one flank placed.'),
         ('CC','5B','522.08-525.15'): ({'start': None, 'end': None, 'peak': None, 'grade': 'D'}, 'The two flanking tags place 193 Mb apart on 5B (3958735 at 406.93 Mb, 1137742 at 600.13 Mb) for a 3 cM interval, so at least one is a paralogous hit; left unplaced.'),
         ('CC','7B','134.52-173.55'): ({'start': None, 'end': None, 'peak': 0.642, 'grade': 'C'}, DS + ' Only one flank placed.')}   # positions the atlas previously derived from Roy et al. Data Sheet 1 (tag BLAST)
rp = {'BC': 'BARTAI x CIANO T79', 'CC': 'CASCABEL x CIANO T79'}
for p, ch, cm, l, r, pve, src in RY:
    rec = add(name='QSb-%s-%s %s cM (Roy 2021)' % (p, ch, cm), disease='Spot blotch', chr=ch, marker='DArTseq %s–%s' % (l, r), pve=pve,
              source='Roy2021FPS', note='%s RIL; resistance allele from %s.' % (rp[p], src), evidence='Table 3: %s %s %s cM %s %s PVE %s' % (p, ch, cm, l, r, pve))
    OLD_ROY[(ch, cm)] = rec
    place(rec, [l, r], ch)
    ov = ROYDS.get((p, ch, cm))
    if ov:
        rec.update(ov[0]); rec['note'] = rec['note'].split(' Only ')[0].split(' No flanking')[0].split(' Marker positions')[0] + ' ' + ov[1]

# ---------------- Singh et al. 2016 Mol Biol Rep 43:1293-1303, Tables 7-8 ----------------
place(add(name='QSb.iiwbr-7B', disease='Spot blotch', chr='7B', marker='Xwmc758–Xwmc335 (peak Xwmc758)', pve='7.2-12.3', stat='LOD 3.4-6.4',
          source='Singh2016MBR', note='Sonalika x BH 1146 RIL; resistance allele from BH 1146. The paper calls this 7BL, but both flanking markers lie on 7BS (centromere near 300 Mb).',
          evidence='Table 7: wmc758-wmc335 8.6 cM, mean LOD 6.3 R2 11.4'), ['wmc758', 'wmc335'], '7B', peak_marker='wmc758')
place(add(name='QSb.iiwbr-7D', disease='Spot blotch', chr='7D', marker='Xwmc653–Xbarc121 / Xwmc653–Xcfa2040 (peak Xbarc121)', pve='5.6-18.3', stat='LOD 2.9-5.4',
          source='Singh2016MBR', note='Sonalika x BH 1146 RIL; resistance allele from BH 1146; 7DL.',
          evidence='Table 8: wmc653-barc121 26.6 cM / wmc653-cfa2040 47.7 cM'), ['wmc653', 'barc121', 'cfa2040'], '7D', peak_marker='barc121')

# ---------------- Singh et al. 2018 IJMS 19:4054, Table 3 (DArTseq) ----------------
S18 = [('BC','1B','0.4-2.8','995296','5410777','5.4-7.1','CIANO T79',''),
       ('BC','1B','209.4-218.6','5324988','1110815','5.5-8.9','BARTAI',''),
       ('BC','1D','48.4-53.3','100142243','1037975','2.9-3.5','CIANO T79',''),
       ('BC','3A','103.2-146.0','1109808','990692','3.2-5.1','CIANO T79',''),
       ('BC','4B','111.3-117.1','1247021','1702494','2.2-5.2','BARTAI',''),
       ('BC','4D','0-9.8','BS00036421_51','1119387','3.6','BARTAI','associated with PH'),
       ('BC','5A','175.9-179.4','Vrn-A1','3064415','6.9-21.9 (mean 12.5)','BARTAI','associated with DH'),
       ('BC','6D','21.9-33.2','1239681','1095962','2.4-2.7','BARTAI',''),
       ('WC','1B','232.3-252.9','1283166','2263671','2.9-8.3','WUYA',''),
       ('WC','2D','2.0-3.2','1085831','1098973','4.7-9.0','WUYA','associated with DH'),
       ('WC','4B','47.9-54.1','Rht-B1','1241652','2.3-2.7','CIANO T79','associated with PH'),
       ('WC','4B','99.2-107.4','1217659','997944','4.5-6.0','WUYA',''),
       ('WC','5A','147.5-148.4','1135154','2260918','21.9-27.1 (mean 25.1)','WUYA','associated with DH')]
GENE = {'Vrn-A1': ('5A', 587.4), 'Rht-B1': ('4B', 30.86), 'Rht-D1': ('4D', 18.78)}
sp = {'BC': 'BARTAI x CIANO T79', 'WC': 'WUYA x CIANO T79'}
for p, ch, cm, l, r, pve, src, rem in S18:
    rec = add(name='QSb-%s-%s %s cM (Singh 2018)' % (p, ch, cm), disease='Spot blotch', chr=ch, marker='%s–%s' % (l, r), pve=pve,
              source='V_Singh2018', note='%s RIL; resistance allele from %s.%s' % (sp[p], src, (' ' + rem + '.') if rem else ''),
              evidence='Table 3: %s %s %s %s %s PVE %s' % (p, ch, cm, l, r, pve))
    ms = [m for m in (l, r) if m not in GENE]
    place(rec, ms, ch)
    for m in (l, r):
        if m in GENE:
            g = GENE[m][1]
            if rec['peak'] is not None and rec['start'] is None:
                rec['start'], rec['end'] = sorted((rec['peak'], g)); rec['peak'] = None; rec['grade'] = 'B'
            elif rec['start'] is None:
                rec['peak'] = g; rec['grade'] = 'C'
            rec['note'] += ' %s placed at the gene (%.2f Mb RefSeq v1.0).' % (m, g)

# ---------------- Vasistha et al. 2024 Sci Rep 14:13083, Table 1 (16 MQTLs) ----------------
VG = {'MQTL3': ('1B', 'TraesCS1B02G035600', 'TraesCS1B02G038100'), 'MQTL4': ('2A', 'TraesCS2A02G193300', 'TraesCS2A02G194400'),
      'MQTL5': ('2A', 'TraesCS2A02G012800', 'TraesCS2A02G014900'), 'MQTL8': ('5A', 'TraesCS5A02G487900', 'TraesCS5A02G489800')}
GV = {}
for l in gzip.open(W + 'genes_v10.tsv.gz', 'rt'):
    g, c, s, e = l.split('\t'); GV[g] = (c, int(s) / 1e6, int(e) / 1e6)
MQ = [('MQTL1','1A','73.8/57.8','14.0 (5.4)','tplb0025b13_2792','1115313','SBR, GLMS'),
      ('MQTL2','1B','1.2/1.5','16.0 (3.4)','Xgwm374','1047726','SBR'),
      ('MQTL3','1B','97.2/10.0','14.2 (6.9)','RFL_Contig2826_548','IACX886','SBR, GLAD, GLMS'),
      ('MQTL4','2A','26.7/4.5','10.7 (5.5)','Xpsr934','Xfba385','GLMS'),
      ('MQTL5','2A','52.3/2.4','21.8 (4.7)','wms210','wmc474','SBR, BPR, GLMS, GLAD'),
      ('MQTL6','2D','85.5/8.3','16.0 (5.4)','AX-109422256','AX-110072786','SBR, BPR'),
      ('MQTL7','4B','71.8/4.6','5.4 (5.5)','BS00062892_51','wsnp_Ex_c296_574790','SBR, SG, GLMS'),
      ('MQTL8','5A','86.8/3.6','13.6 (4.8)','12772364','Tdurum_contig56267_180','SBR, GLAD, GLMS'),
      ('MQTL9','5A','108.3/3.5','18.0 (4.5)','Bobwhite_c23972_1007','5411879','GLAD, FLS'),
      ('MQTL10','5B','130.3/42.1','14.4 (6.1)','Xcfd5','Xdupw115','SG, GLAD, GLMS'),
      ('MQTL11','5B','303.3/20.7','19.0 (3.3)','QSb.bhu-5B','Q.sb.bisa-5B','SBR'),
      ('MQTL12','6A','22.7/10.3','13.0 (5.8)','RAC875_c21158_99','1124657','SBR, GLMS'),
      ('MQTL13','6A','61.5/36.1','10.5 (5.2)','D_contig80449_91','BobWhite_c20926_313','SBR, FLS'),
      ('MQTL14','7A','63.6/19.1','10.8 (11.5)','RAC875_c17331_79','BobWhite_c27734_355','FLS, GLAD, GLMS'),
      ('MQTL15','7B','43.2/1.8','8.5 (3.8)','Excalibur_c29698_76','fbb226','SBR'),
      ('MQTL16','7D','37.9/38.7','9.0 (3.2)','XgbxR799a','Xgdm145','SBR; Sb1 lies within the interval')]
for nm, ch, cm, pv, l, r, tr in MQ:
    rec = add(name='%s-%s (Vasistha 2024)' % (nm, ch), cls='MQTL', disease='Spot blotch', chr=ch, marker='%s–%s' % (l, r),
              pve=pv.split(' ')[0], stat='mean LOD ' + pv.split('(')[1].rstrip(')'), source='Vasistha2024SR',
              note='Meta-QTL, consensus map peak/CI %s cM; traits %s.' % (cm, tr), evidence='Table 1: %s-%s %s cM PVE %s %s-%s' % (nm, ch, cm, pv, l, r))
    if nm in VG:
        c, g1, g2 = VG[nm]
        a, b = GV[g1], GV[g2]
        rec['start'], rec['end'] = round(min(a[1], b[1]), 3), round(max(a[2], b[2]), 3); rec['grade'] = 'B'
        rec['note'] += ' Placed on the span of the candidate genes the paper lists for this MQTL (Table 2; genes were taken from 1 Mb either side of the MQTL peak), %s-%s, RefSeq v1.0 coordinates from Ensembl.' % (g1, g2)
        place_chk = dict(rec, note='', start=None, end=None, peak=None); place(place_chk, [l, r], ch)
        if place_chk['start'] is not None or place_chk['peak'] is not None:
            fl = [x for x in (place_chk['start'], place_chk['end'], place_chk['peak']) if x is not None]
            if min(abs(f - rec['start']) for f in fl) > 20:
                rec['note'] += ' CONFLICT: the flanking markers in Table 1 place on %s at %s Mb, far from the candidate-gene window; check the MQTL table or the gene list.' % (ch, ', '.join('%.1f' % f for f in fl))
                rec['conflict'] = True
    elif nm == 'MQTL11':
        rec['start'], rec['end'] = 315.534, 504.302; rec['grade'] = 'C'
        rec['note'] += ' Flanks are QTL names, not markers; placed from QSb.bhu-5B (Xgwm067, 315.5 Mb) to Sb2/Xgwm639 (504.3 Mb).'
    elif nm == 'MQTL16':
        place(rec, [r], ch)
        rec['start'], rec['end'], rec['peak'] = rec['peak'], 47.4245, None; rec['grade'] = 'C'
        rec['note'] += ' XgbxR799a not placed; interval runs from Xgdm145 to Sb1/Lr34 (47.42 Mb), which the paper states lies inside the MQTL.'
    else:
        place(rec, [l, r], ch)

# ---------------- Zhang et al. 2020 TAG 133:2451-2459: Sb4 ----------------
a, b = GV['TraesCS4B02G295400'], GV.get('TraesCS4B02G296700') or (None, None, 582.11)
add(name='Sb4', cls='DG', disease='Spot blotch', chr='4B', start=round(a[1], 3), end=582.11, marker='YK12928 co-segregating; 1.19 cM, 1.34 Mb; candidate genes TraesCS4B01G295400–TraesCS4B01G296700',
    source='Zhang2020TAG', grade='A', note='GY17; fine-mapped to 1.34 Mb on 4BL. Start from TraesCS4B02G295400 (Ensembl RefSeq v1.0), end 581.07 + 1.34 Mb as stated.',
    evidence='Abstract: Sb4 mapped in a 1.19 cM interval corresponding to 1.34 Mb of CS 4BL; 21 predicted genes')

# ---------------- Zhu et al. 2014 Mol Breed 34:1215-1228, Table 4 ----------------
for nm, l, r, ch, pve, lod in [('QSb.cim-1B', 'Xwmc128', 'Xgwm374', '1B', '5.8-8.5', '2.1-3.1'),
                               ('QSb.cim-3B', '990937', '1123330', '3B', '7.1-26.0', '2.6-10.2'),
                               ('QSb.cim-5A', '1086218', '982608', '5A', '8.4-12.3', '2.6-4.0'),
                               ('QSb.cim-1A', '1026215', '1088769', '1A', '5.3', '2.8'),
                               ('QSb.cim-14', '1016317', '1104599', '6B', '4.9', '2.8')]:
    rec = add(name=nm, disease='Spot blotch', chr=ch, marker='%s–%s' % (l, r), pve=pve, stat='LOD ' + lod, source='Zhu2014MB',
              note='SYN1 x Ocoroni DH; ICIM.', evidence='Table 4: %s %s-%s PVE %s' % (nm, l, r, pve))
    if nm == 'QSb.cim-14':
        rec['note'] += ' Linkage group 14 was not anchored to a chromosome in the paper; chromosome 6B is inferred here from the position of DArT 1104599 only.'
    place(rec, [l, r], ch)

# ---------------- Bemister et al. 2019 Crop Sci 59:650-658, Table 3 ----------------
LIU = {'d3384830': ('4B', 27.7987), 's1047233': ('7D', 71.5695), 's1022494': ('7D', 84.8045)}
for tr, dis, nm, ch, cmpos, l, r, lod, r2, src in [
        ('Ls', 'Leaf spot complex', 'QLs.dms-2A', '2A', '75 (74.5-75.5)', 'd1117218', 'd1089684', '5.7', '10', 'Carberry'),
        ('Ls', 'Leaf spot complex', 'QLs.dms-4B', '4B', '43 (41.5-44.5)', 'd3384830', 's4991673', '4.3', '8', 'Peace'),
        ('Ls', 'Leaf spot complex', 'QLs.dms-7D', '7D', '13 (11.5-14.5)', 's1047233', 's1022494', '5.9', '10', 'Carberry'),
        ('Lr', 'Leaf rust', 'QLr.dms-4A', '4A', '64 (63.5-64.5)', 'd1164232', 's1118950', '3.4', '6', 'Carberry'),
        ('Lr', 'Leaf rust', 'QLr.dms-3D', '3D', '45 (44.5-45.5)', 'd977229', 'd1269275', '2.9', '5', 'Carberry'),
        ('Yr', 'Stripe rust', 'QYr.dms-3A', '3A', '70 (66.5-72.5)', 's3948706', 'd1103494', '3.4', '6', 'Carberry'),
        ('Yr', 'Stripe rust', 'QYr.dms-4B', '4B', '21 (20.5-21.5)', 'd1111841', 'd1087899', '7.9', '16', 'Carberry')]:
    rec = add(name=nm, disease=dis, chr=ch, marker='%s–%s' % (l, r), pve=r2, stat='LOD ' + lod, source='Bemister2019CS',
              note='Peace x Carberry RIL, combined environments; %s cM; resistance allele from %s.' % (cmpos, src) +
                   (' Leaf spot complex scored in the field (tan spot, Septoria blotches and spot blotch not separated).' if tr == 'Ls' else ''),
              evidence='Table 3: %s %s %s cM %s-%s LOD %s R2 %s%%' % (nm, ch, cmpos, l, r, lod, r2))
    got = []
    for m in (l, r):
        p = None
        for c, x, ref, s in MK.get(m[1:], []):
            if c == ch and ref == 'RefSeq v1.0': p = (round(x, 4), s); break
        if p is None and m in LIU and LIU[m][0] == ch: p = (LIU[m][1], 'Liu et al. 2020 Table S5')
        if p: got.append((m, p[0], p[1]))
        else: LOG.append((nm, m, 'not placed'))
    if len(got) == 2:
        rec['start'], rec['end'] = sorted([got[0][1], got[1][1]]); rec['grade'] = 'B'
    elif got:
        rec['peak'] = got[0][1]; rec['grade'] = 'C'; rec['note'] += ' Only %s placed; position is that marker.' % got[0][0]
    else:
        rec['grade'] = 'D'
    if nm == 'QYr.dms-4B':
        rec['ref'] = 'RefSeq v2.0'; rec['peak'] = 4.7; rec['grade'] = 'C'
        rec['note'] += ' d1111841 at 4.7 Mb on 4B (RefSeq v2.0) per Iqbal et al. (2023); d1087899 maps to 4D (1.8 Mb), so the right flank is on the wrong chromosome.'
        got = [('d1111841', 4.7, 'Iqbal et al. 2023')]
    if got: rec['note'] += ' Marker positions: ' + '; '.join('%s %.3f Mb (%s)' % g for g in got) + '.'

# ---------------- Iqbal et al. 2023 Crop Sci 63:2310-2328, Table 1 (RefSeq v2.0 printed) ----------------
IQ = """PAC Yr QYr.dms-2A.4 2A 657.9 654.7 661.2 3.0 3.5
PAC Yr QYr.dms-3A.4 3A 667.6 666.4 667.9 5.1 4.9
PAC Yr QYr.dms-3D.1 3D 27.9 26.8 28.1 3.9 5.1
PAC Yr QYr.dms-3D.2 3D 533.7 517.5 550.4 4.4 4.4
PAC Yr QYr.dms-4B.1 4B 2.1 1.8 2.2 11.7 12.1
PAC Yr QYr.dms-5A 5A 557.0 556.2 557.4 3.3 3.1
PAC Yr QYr.dms-5B.2 5B 579.7 579.4 580.8 6.1 5.9
PAC Yr QYr.dms-6B.2 6B 725.2 725.1 725.3 2.5 2.4
PAC Lr QLr.dms-1B 1B 5.1 5.0 5.4 3.0 2.1
PAC Lr QLr.dms-2D.2 2D 71.6 64.2 79.0 7.8 5.7
PAC Lr QLr.dms-3D.1 3D 26.3 24.3 27.7 3.6 3.2
PAC Lr QLr.dms-3D.2 3D 489.0 481.8 496.6 5.9 9.8
PAC Lr QLr.dms-7B 7B 751.7 744.1 759.4 5.2 4.2
PAC Ls QLs.dms-1A 1A 368.2 367.7 369.2 6.0 1.4
PAC Ls QLs.dms-2A.2 2A 641.6 628.1 654.7 10.6 2.6
PAC Ls QLs.dms-2D.2 2D 616.1 616.0 616.1 4.0 0.7
PAC Ls QLs.dms-3A 3A 53.7 49.9 57.9 3.4 0.5
PAC Ls QLs.dms-3B 3B 54.2 52.0 56.7 2.8 0.5
PAC Ls QLs.dms-4B.1 4B 3.5 2.5 5.0 52.1 13.9
PAC Ls QLs.dms-5B 5B 561.9 561.4 563.1 4.2 0.6
PAC Ls QLs.dms-7D 7D 69.8 67.5 72.3 14.5 3.2
PAC Cbt QCbt.dms-1B 1B 21.4 21.0 21.5 13.7 20.2
PAC Cbt QCbt.dms-6D 6D 7.6 7.4 7.6 14.8 20.9
ACG Yr QYr.dms-3A.2 3A 53.8 53.5 53.9 3.8 4.6
ACG Yr QYr.dms-4A.1 4A 323.5 322.8 325.8 3.2 7.7
ACG Yr QYr.dms-5B.2 5B 547.4 544.7 548.7 2.8 3.3
ACG Lr QLr.dms-2B.1 2B 25.6 21.9 26.5 2.5 1.7
ACG Lr QLr.dms-2D.2 2D 99.8 98.9 100.9 2.6 10.8
ACG Ls QLs.dms-2D.1 2D 34.2 34.0 35.3 4.7 12.8
ACG Ls QLs.dms-5B 5B 553.7 553.5 555.7 2.8 6.2
ACG Ls QLs.dms-6B 6B 66.5 66.4 66.7 3.1 6.5
ACG Cbt QCbt.dms-2B 2B 244.0 243.5 244.2 15.4 14.7"""
DIS = {'Yr': 'Stripe rust', 'Lr': 'Leaf rust', 'Ls': 'Leaf spot complex', 'Cbt': 'Common bunt'}
IP = {'PAC': 'Peace x Carberry', 'ACG': 'Attila x CDC Go'}
for line in IQ.splitlines():
    p, tr, nm, ch, pk, s, e, lod, r2 = line.split()
    add(name=nm + ' (%s)' % p, disease=DIS[tr], chr=ch, start=float(s), end=float(e), peak=float(pk), ref='RefSeq v2.0', pve=r2, stat='LOD ' + lod,
        source='Iqbal2023CS', grade='A', cls='QTL',
        note='%s RIL, means of all environments. Positions printed in Table 1 (IWGSC RefSeq v2.0).' % IP[p] +
             (' Leaf spot complex scored in the field.' if tr == 'Ls' else ''),
        evidence='Table 1: %s %s %s %s %s-%s Mb LOD %s R2 %s' % (p, nm, ch, pk, s, e, lod, r2))

if __name__ == '__main__':
    print(len(R), 'records')
    for x in R: print(x['source'], '|', x['name'], '|', x['chr'], x['start'], x['end'], x['peak'], x['grade'], '|', x['marker'][:40])
    print('\nunplaced markers:'); [print(' ', *l) for l in LOG]
