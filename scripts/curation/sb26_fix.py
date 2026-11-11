"""Batch 2 (11 Nov 2026): loci of 26 Septoria nodorum blotch / tan spot papers re-read from the papers.
Same rules as sb19_fix.py: grade A = position printed in the paper; B = from the paper's markers or genes;
C = one marker only, or another paper's projection; D = unplaced. SSRs by exact-match in-silico PCR on
RefSeq v1.0 (sb19/ispcr_res*.json), used only on the chromosome the paper states."""
import json, re, gzip, sys, os
sys.path.insert(0, 'scripts/curation/sb19')
import lift
W = 'scripts/curation/sb19/'
PCR = {}
for f in ('ispcr_res.json', 'ispcr_res3.json', 'ispcr_roder.json'):
    if os.path.exists(W + f):
        for k, v in json.load(open(W + f)).items(): PCR.setdefault(k, []).extend(v)
MK = json.load(gzip.open('data/markerpos_all.json.gz'))
RODER = json.load(open(W + 'roder1998_primers.json'))  # Roder et al. 1998 Genetics 149:2007, mapped chromosomes per GWM primer set
LOG = []
LEGACY = {('barc175', '6D'): 411.884}
# genes placed from Ensembl RefSeq v1.0 coordinates or cloned-gene records already in the atlas
GENE = {'Tsn1': ('5B', 546.698, 546.843), 'Snn1': ('1B', 2.356, 2.357), 'Snn3-B1': ('5B', 6.648, 6.975),
        'Vrn-A1': ('5A', 587.40, 587.42), 'Rht1': ('4B', 30.86, 30.86), 'Ppd-D1': ('2D', 33.95, 33.96), 'Ppd-A1': ('2A', 36.93, 36.94),
        'Lr46': ('1B', 670.6, 673.7)}

def nk(m):
    k = re.sub(r'[^a-z0-9]', '', m.lower())
    if k.startswith('x'): k = k[1:]
    k = re.sub(r'^wms', 'gwm', k)
    k = re.sub(r'^(gwm|barc|wmc|cfd|cfa|gdm|swm|edm)0+(\d)', r'\1\2', k)
    k = re.sub(r'^(gwm|barc|wmc|cfd|cfa|gdm|glk|fba|edm|psr|mwg)(\d+)[a-z]$', r'\1\2', k)   # locus suffix a/b/c
    k = re.sub(r'^wpt(\d)', r'wpt\1', k)
    return k

def pos(m, ch):
    if m in GENE and GENE[m][0] == ch: return GENE[m][1], 'gene'
    k = nk(m)
    hits = [h for h in PCR.get(k, []) if h[0] == ch]
    if hits:
        if len(hits) > 1 and max(h[1] for h in hits) - min(h[1] for h in hits) > 5:
            return None, '%s amplifies at %d sites on %s' % (m, len(hits), ch)
        return hits[0][1], 'in-silico PCR'
    if PCR.get(k) and ch not in RODER.get(k, {}).get('chrs', []):
        return None, '%s amplifies only on %s' % (m, ','.join(sorted({h[0] for h in PCR[k]})))
    if (k, ch) in LEGACY: return LEGACY[(k, ch)], 'earlier atlas in-silico PCR'
    keys = [k, 'wpt' + k[3:] if k.startswith('wpt') else k]
    for kk in keys:
        ent = MK.get(kk, [])
        prim = {c for c, p, ref, src in ent if src == 'Zhao et al. 2019 J Genet'}
        for c, p, ref, src in ent:
            if c == ch and ref == 'RefSeq v1.0' and 'atlas marker dictionary' not in src:
                if '(Table S1)' in src and prim and ch not in prim: return None, 'secondary copy'
                return round(p, 4), src
        for c, p, ref, src in ent:
            if c == ch and ref == 'RefSeq v2.1':
                y = lift.lift_record({'chr': ch, 'peak': p, 'ref': 'RefSeq v2.1'})['RefSeq v1.0'].get('peak')
                return y, '%s (v2.1 %.3f Mb, converted)' % (src, p)
    return None, 'not placed on ' + ch

def place(rec, markers, ch, peak_marker=None):
    got = []
    for m in markers:
        if not m: continue
        p, how = pos(m, ch)
        if p is not None: got.append((m, p, how))
        else: LOG.append((rec['name'], m, how))
    if peak_marker:
        p, how = pos(peak_marker, ch)
        if p is not None:
            rec['peak'] = p; rec['grade'] = rec.get('grade') or 'B'
            rec['note'] += ' Peak marker %s at %.3f Mb (%s).' % (peak_marker, p, how)
    if len(got) >= 2:
        ps = sorted(p for _, p, _ in got)
        rec['start'], rec['end'] = ps[0], ps[-1]; rec['grade'] = rec.get('grade') or 'B'
        if ps[-1] - ps[0] > 150: rec['note'] += ' Flanking markers lie %.0f Mb apart; interval is wide.' % (ps[-1] - ps[0])
    elif len(got) == 1 and rec.get('peak') is None:
        rec['peak'] = got[0][1]; rec['grade'] = 'C'
        rec['note'] += ' Only %s could be placed on %s; position is that marker, not a reported peak.' % (got[0][0], ch)
    elif not got and rec.get('peak') is None:
        rec['grade'] = 'D'; rec['note'] += ' No marker of this locus could be placed on %s; left unplaced.' % ch
    if got: rec['note'] += ' Marker positions (RefSeq v1.0): ' + '; '.join('%s %.3f Mb (%s)' % g for g in got) + '.'
    rec['note'] = rec['note'].strip()
    return rec

R = []
def add(**k):
    for f, v in (('ref', 'RefSeq v1.0'), ('curation', 'hand'), ('cls', 'QTL'), ('start', None), ('end', None), ('peak', None),
                 ('pve', ''), ('stat', ''), ('note', ''), ('grade', None), ('marker', ''), ('disease', 'Septoria nodorum blotch')):
        k.setdefault(f, v)
    R.append(k); return k

# ---------- Lin et al. 2020 TAG 133:785-808 (NIAB Elite MAGIC), Table 3: v1.0 bp printed ----------
for name, rows in json.load(open('/home/claude/batch2/ex/Lin2020.json')).items():
    s, e = min(r['s'] for r in rows), max(r['e'] for r in rows)
    best = max(rows, key=lambda r: r['lp'])
    ch = re.search(r'-(\d[ABD])', name).group(1)
    nt = 'NIAB Elite MAGIC (8 founders). Interval spans the flanking-marker positions of all %d table rows: %s.' % (
        len(rows), '; '.join('%s %.2f-%.2f Mb R2 %s%%' % (r['trait'], r['s'], r['e'], r['r2']) for r in rows))
    if name.endswith('7D.3'):
        # Table 3 prints 629325776-724128709 on 7D, but 7D is only 638.7 Mb long in RefSeq v1.0:
        # 724.13 Mb is the 7A position of RAC875_c10022_23. Keep the 7D flank JD_c2708_1512 only.
        s, e = 629.3258, 629.3258
        nt += ' The printed end (724.13 Mb) is beyond the end of 7D and is the 7A position of RAC875_c10022_23; the locus is placed at the 7D flank JD_c2708_1512.'
    if name.endswith('5B.1'): nt += ' Co-locates with Snn3-B1.'
    if name.endswith('5B.2'): nt += ' Co-locates with Tsn1.'
    add(name=name, chr=ch, start=s, end=e, marker=best['fl'], pve=', '.join(r['r2'] for r in rows), stat='-log10p ' + ', '.join(str(r['lp']) for r in rows),
        source='A_SNB_007', grade='B' if name.endswith('7D.3') else 'A', note=nt, evidence='Table 3: %s %s %s%.4f-%.4f Mb' % (name, ch, 'printed ' if name.endswith('7D.3') else '', min(r['s'] for r in rows), max(r['e'] for r in rows)))

# ---------- Lin et al. 2021 TAG 134:125-142 (BMWpop MAGIC), Table 2: v1.0 bp printed ----------
for name, ch, s, e, r2, tr, fl in [
        ('QSnb.nmbu-2A.1', '2A', 507.691472, 718.885511, '5.4 (2018), 5.2 (2016)', 'leaf blotch 2018 and 2016 (2016 interval 688.62-693.29 Mb)', 'Kukri_c11327_977–Tdurum_contig33398_106'),
        ('QSnb.nmbu-2B.1', '2B', 572.591268, 648.083739, '6.9', 'leaf blotch 2017', 'Excalibur_rep_c66577_159–AX-94522698'),
        ('QSnb.nmbu-2D.1', '2D', 14.636197, 15.115231, '7.0', 'leaf blotch 2016', 'D_F1BEJMU01A0OMY_356–BS00067698_51'),
        ('QInf.nmbu-5A.1', '5A', 391.548987, 413.418597, '6.2', 'culture filtrate infiltration, isolate 203649', 'BS00040623_51–BobWhite_rep_c64913_315'),
        ('QSnb.nmbu-5A.1', '5A', 558.69278, 571.683217, '6.7 (2016), 5.0 (2018)', 'leaf blotch 2016 and 2018', 'Excalibur_c472_914–Tdurum_contig54785_216'),
        ('QTox3.nmbu-5B.1', '5B', 12.324613, 16.421984, '36.0', 'SnTox3 infiltration; Snn3-B1', 'wsnp_Ku_c5308_9450093–wsnp_Ex_c2459_4591587'),
        ('QInf.nmbu-5B.1', '5B', 508.795354, 516.041245, '6.0', 'culture filtrate infiltration, isolate 203649', 'BobWhite_c7070_196–BobWhite_rep_c50822_462'),
        ('QInf.nmbu-6A.1', '6A', 0.796322, 2.40251, '6.9', 'culture filtrate infiltration, isolate 203649', 'wsnp_Ex_c21633_30782312–RAC875_c6135_95'),
        ('QTox3.nmbu-6A.1', '6A', 545.83235, 574.479996, '5.5', 'SnTox3 infiltration', 'wsnp_RFL_Contig1871_1020122–IAAV151'),
        ('QInf.nmbu-7B.1', '7B', 679.800093, 700.551772, '17.1', 'culture filtrate infiltration, isolate 203649', 'GENE-4442_121–BS00057323_51')]:
    add(name=name, chr=ch, start=round(s, 4), end=round(e, 4), marker=fl, pve=r2, source='Lin2021TAG', grade='A',
        note='BMWpop German MAGIC (8 founders); %s. RefSeq v1.0 range printed in Table 2.' % tr, evidence='Table 2: %s %s %.3f-%.3f Mb R2 %s' % (name, ch, s, e, r2))

# ---------- Singh et al. 2019 IJMS 20:5432 (PBW343 x Kenya Nyangumi), Tables 3-4 ----------
for nm, ch, cm, l, r, lod, pve, dis in [
        ('QTs.cim-5B.1', '5B', '112', 'wPt-8942', 'wPt-4936', '4.33', '2.98', 'Tan spot'), ('QTs.cim-5B.2', '5B', '180', 'XFCP393', 'wPt-1733', '17.73', '23.32', 'Tan spot'),
        ('QTs.cim-7B', '7B', '14', 'wPt-0266', 'wPt-5892', '35.27', '11.75', 'Tan spot'), ('QTs.cim-1B', '1B', '52.8', 'CSLV46-1B', 'wPt-1070', '4.83', '3.00', 'Tan spot'),
        ('QTs.cim-7D', '7D', '43', 'wPt-663918', 'wPt-733104', '8.86', '5.78', 'Tan spot'), ('QTs.cim-1A.1', '1A', '0', 'wPt-0595', 'tPt-7214', '6.02', '3.69', 'Tan spot'),
        ('QTs.cim-1A.2', '1A', '25', 'wPt-741704', 'tPt-7559', '23.18', '19.00', 'Tan spot'),
        ('QSnb.cim-5B.1', '5B', '168', 'wPt-3661', 'wPt-3457', '4.88', '8.00', 'Septoria nodorum blotch'), ('QSnb.cim-5B.2', '5B', '180', 'XFCP393', 'wPt-1733', '11.38', '19.51', 'Septoria nodorum blotch'),
        ('QSnb.cim-1A.1', '1A', '47', 'wPt-8172', 'wPt-6139', '5.70', '5.24', 'Septoria nodorum blotch'), ('QSnb.cim-1A.2', '1A', '104', 'wPt-729972', 'wPt-2527', '12.05', '11.45', 'Septoria nodorum blotch'),
        ('QSnb.cim-5A', '5A', '10', 'wPt-1370', 'wPt-9513', '11.24', '10.69', 'Septoria nodorum blotch')]:
    rec = add(name=nm, disease=dis, chr=ch, marker='%s–%s' % (l, r), pve=pve, stat='LOD ' + lod, source='A_SNB_019',
              note='PBW343 x Kenya Nyangumi RIL, seedling; %s cM.' % cm, evidence='Table %s: %s %s %s cM %s %s LOD %s PVE %s' % ('3' if dis == 'Tan spot' else '4', nm, ch, cm, l, r, lod, pve))
    if l == 'XFCP393':   # Tsn1-linked marker
        rec['note'] += ' XFCP393 is a Tsn1 marker; Tsn1 gene position used for that flank.'
        place(rec, ['Tsn1', r], ch)
    elif l == 'CSLV46-1B':
        rec['note'] += ' csLV46 marks Lr46 (Lr46 region as projected by Gahtyari et al. 2021).'
        place(rec, [r], ch)
        if rec['peak'] is not None: rec['start'], rec['end'], rec['peak'] = sorted([rec['peak'], 670.6])[0], max(rec['peak'], 673.7), None; rec['grade'] = 'C'
    else:
        place(rec, [l, r], ch)

# ---------- Schnurbusch et al. 2003 TAG 107:1226-1234 (Arina x Forno), Tables 2-3, glume blotch ----------
for nm, ch, l, r, cm, pve, al in [
        ('QSng.sfr-2AL', '2A', 'cfd276a', 'cfa2086', '188', '9.1', 'Arina'), ('QSng.sfr-2BS', '2B', 'OA102', 'psr933b', '0', '3.6', 'Forno'),
        ('QSng.sfr-2BL', '2B', 'glk600', 'psr644b', '158', '7.8', 'Forno'), ('QSng.sfr-3AS', '3A', 'cfd79a', 'gwm369', '0', '7.1-8.3', 'Arina'),
        ('QSng.sfr-3BS', '3B', 'gwm389', 'cfd79c', '0', '12.8-31.2', 'Arina'), ('QSng.sfr-3BL', '3B', 'cfa2134b', 'gwm131b', '86', '9.4', 'Arina'),
        ('QSng.sfr-4BL', '4B', 'gwm165', 'glk335', '14-16', '7.0-22.6', 'Forno'), ('QSng.sfr-5AL', '5A', 'glk317a', 'gwm639b', '68', '15.7', 'Arina'),
        ('QSng.sfr-5BL', '5B', 'gwm371', 'gwm639a', '58-64', '7.5-14.4', 'Forno'), ('QSng.sfr-5DL', '5D', 'cfd81', 'cfd266', '40', '7.3', 'Arina'),
        ('QSng.sfr-6BL', '6B', 'fba81', 'psr924', '46', '6.9', 'Arina'), ('QSng.sfr-6Dc', '6D', 'cfd19', 'gwm55b', '108', '9.2', 'Arina'),
        ('QSng.sfr-7BS', '7B', 'cfa2174b', 'gwm46', '26', '10.5', 'Forno')]:
    place(add(name=nm, chr=ch, marker='%s–%s' % (l, r), pve=pve, source='Z_1366',
              note='Arina x Forno RIL; Stagonospora glume blotch; peak %s cM; resistance allele from %s.' % (cm, al),
              evidence='Table 3: %s %s–%s peak %s cM R2 %s' % (nm, l, r, cm, pve)), [l, r], ch)

# ---------- Shankar et al. 2008 Phytopathology 98:886-893 (6HRWSN125 x WAWHT2074), Table 3 ----------
SH = [('QSns.daw-2D', '2D', 'cfd50', 'gwm382', 'wPt-8524', '3.12', '5', 'seedling 2003'),
      ('QSns.daw-4D', '4D', 'cfd54', 'wPt-7009', 'cfd54', '3.45', '5', 'seedling 2003'),
      ('QSns.daw-5B.1', '5B', 'wPt-6136', 'wPt-1420', 'gwm234', '5.40', '8', 'seedling 2003'),
      ('QSns.daw-5B.2', '5B', 'wPt-9613', 'wPt-5896', 'wPt-0498', '7.43', '13', 'seedling 2003'),
      ('QSns.daw-6D', '6D', 'gdm46', 'edm148', 'cfd76a', '4.22', '6', 'seedling 2003'),
      ('QSng.daw-2A', '2A', 'wPt-7026', 'gwm71c', 'wmc258', '4.71', '5', 'glume 2003'),
      ('QSng.daw-2B (2003)', '2B', 'gwm148', 'wPt-8072', 'gwm319', '3.01', '3', 'glume 2003'),
      ('QSng.daw-2B (2004)', '2B', 'wPt-4199', 'wPt-9336', 'gwm191c', '7.15', '13', 'glume 2004'),
      ('QSng.daw-4B', '4B', 'Rht1', 'gwm149', 'barc340', '6.35-15.83', '8-19', 'glume 2003 (Rht1-gwm495) and 2004 (Rht1-gwm149)'),
      ('QSng.daw-5A', '5A', 'gwm617a', 'wPt-0373', 'wPt-0373', '10.43', '8', 'glume 2003'),
      ('QSng.daw-5B', '5B', 'wPt-5896', 'gwm604', 'wPt-5896', '3.01', '3', 'glume 2003'),
      ('QSng.daw-7B', '7B', 'wPt-4025', 'wPt-6498', 'wPt-2572', '4.17', '5', 'glume 2003'),
      ('QSng.daw-1B', '1B', 'cfd59b', 'wPt-3566', 'wPt-0170', '7.51', '10', 'glume 2004'),
      ('QSng.daw-2D', '2D', 'cfd11', 'gwm157', 'wPt-2544', '3.50', '4', 'glume 2004'),
      ('QSnl.daw-1B', '1B', 'wPt-3566', 'cfa2129b', 'wPt-8279', '4.00', '5', 'flag leaf 2003'),
      ('QSnl.daw-2D (2003)', '2D', 'wPt-0638', 'gwm30', 'gwm157', '5.46', '8', 'flag leaf 2003'),
      ('QSnl.daw-2D (2004)', '2D', 'cfd11', 'gwm157', 'wPt-2544', '8.69', '17', 'flag leaf 2004'),
      ('QSnl.daw-5A', '5A', 'gwm617a', 'wPt-0373', 'wPt-0373', '10.39', '21', 'flag leaf 2003'),
      ('QSnl.daw-4B', '4B', 'wPt-3991', 'barc340', 'wPt-7569', '2.94-8.80', '4-11', 'flag leaf 2004'),
      ('QSnl.daw-5B', '5B', 'wPt-6136', 'wPt-1420', 'gwm234', '2.72', '4', 'flag leaf 2004'),
      ('QSnr.daw-2A', '2A', 'wPt-7026', 'wmc296', 'wmc453', '2.70-6.56', '3-6', 'relative grain weight 2003, 2004'),
      ('QSnr.daw-2B', '2B', 'wPt-4199', 'wPt-9336', 'gwm191c', '3.40', '8', 'relative grain weight 2004'),
      ('QSnr.daw-2D.1', '2D', 'cfd11', 'wPt-8330', 'gwm157', '2.18-3.06', '3-5', 'relative grain weight 2003, 2004'),
      ('QSnr.daw-2D.2 (2003)', '2D', 'wPt-8524', 'wPt-8319', 'wPt-3757', '2.60', '3', 'relative grain weight 2003'),
      ('QSnr.daw-2D.2 (2004)', '2D', 'gwm311', 'barc197', 'cfd44', '4.95', '8', 'relative grain weight 2004'),
      ('QSnr.daw-4B', '4B', 'Rht1', 'gwm149', 'barc340', '8.74', '8', 'relative grain weight 2003'),
      ('QSnr.daw-5A', '5A', 'gwm617a', 'wPt-0373', 'wPt-0373', '5.06-17.39', '9-21', 'relative grain weight 2003, 2004'),
      ('QSnr.daw-5B', '5B', 'wPt-0498', 'gwm604', 'wPt-5896', '8.44-13.06', '10-13', 'relative grain weight 2003, 2004')]
for nm, ch, l, r, pk, lod, r2, tr in SH:
    rec = add(name=nm, chr=ch, marker='%s–%s (peak %s)' % (l, r, pk), pve=r2, stat='LOD ' + lod, source='Y_166',
              note='6HRWSN125 x WAWHT2074 DH; %s.' % tr, evidence='Table 3: %s %s-%s max LOD %s (%s) R2 %s' % (nm, l, r, lod, pk, r2))
    place(rec, [l, r], ch, peak_marker=pk)

# ---------- Rivera-Burgos et al. 2022 PLoS ONE 17:e0268546 (AGS2033 x GA03185-12LE29 DH), Tables 2-3 ----------
for nm, ch, cm, ci, lod, pve, tr, extra in [
        ('Qncb.snl-1A', '1A', '25', '22-28', '5.0', '7.3', 'leaf severity, CIM', None), ('Qncb.sng-1A', '1A', '0', '0-3', '4.7', '12.7', 'glume severity, CIM', None),
        ('Qncb.sng-2A', '2A', '34', '30-36', '7.5', '32.1', 'glume severity, CIM', None), ('Qncb.sng-5B', '5B', '153', '143-182', '4.5', '14.0', 'glume severity, CIM', None),
        ('Ppd-A1 (leaf severity)', '2A', '13', '10-16', '8.0', '27.4', 'leaf severity, CIM; the Ppd-A1 heading-date locus', 'Ppd-A1'),
        ('Ppd-D1 (leaf and glume severity)', '2D', '86', '83-89', '7.8-13.3', '31.9-37.2', 'leaf and glume severity, CIM; the Ppd-D1 heading-date locus', 'Ppd-D1'),
        ('Qncb.snl-1Ac', '1A', '4', '0-36', '5.7', '10.2', 'leaf severity, SIM with phenology covariates', None),
        ('Qncb.snl-1Bc', '1B', '11', '0-30', '3.3', '6.7', 'leaf severity, SIM with covariates', None),
        ('Qncb.snl-5Ac', '5A', '201', '10-201', '3.5', '4.0', 'leaf severity, SIM with covariates; text places it at 579.7-579.8 Mb (RefSeq v2.0) with Vrn-A1', 'v2'),
        ('Qncb.snl-6Bc', '6B', '43', '35-279', '3.7', '11.8', 'leaf severity, SIM with covariates', None),
        ('Qncb.sng-1Ac', '1A', '36', '0-48', '7.2', '11.8', 'glume severity, SIM with covariates', None),
        ('Qncb.sng-1Bc', '1B', '20', '0-51', '4.6', '5.6', 'glume severity, SIM with covariates', None)]:
    rec = add(name=nm, chr=ch, pve=pve, stat='LOD ' + lod, source='V_RiveraBurgos2022',
              note='AGS2033 x GA03185-12LE29 DH; %s; %s cM (CI %s cM). Flanking markers are given only in Fig. 2.' % (tr, cm, ci),
              evidence='Table %s: %s %s %s cM LOD %s PVE %s' % ('3' if nm.endswith('c') else '2', nm, ch, cm, lod, pve))
    if extra == 'v2':
        rec.update(ref='RefSeq v2.0', start=579.7, end=579.8, grade='A')
    elif extra:
        g = GENE[extra]; rec.update(start=g[1], end=g[2], grade='B', marker=extra); rec['note'] += ' Placed at the %s gene (RefSeq v1.0).' % extra
    else:
        rec['grade'] = 'D'; rec['note'] += ' Left unplaced.'

# ---------- Czembor et al. 2019 Eur J Plant Pathol 155:1001-1016 (Liwilla x Begra DH), Table 4 ----------
for nm, ch, ms, lod, pve, ds in [('QSnl.ihar-1B', '1B', ['wPt-7160'], '3.34', '12.4', 'f02 leaf AUDPC'),
                                 ('QSnl.ihar-7D', '7D', ['gdm46'], '3.12-4.42', '11.4-20.0', 'f03, f06 leaf AUDPC'),
                                 ('QSng.ihar-3A.1', '3A', ['wPt-6012', 'wPt-3041'], '4.69', '9.5', 't07 head necrosis'),
                                 ('QSng.ihar-3A.2', '3A', ['wPt-9049', 'wPt-2698'], '3.98', '10.3', 'f02 head AUDPC'),
                                 ('QSng.ihar-4A', '4A', ['gwm610', 'barc5'], '4.61', '9.1', 'f06 head AUDPC')]:
    place(add(name=nm, chr=ch, marker=' – '.join(ms), pve=pve, stat='LOD ' + lod, source='Z_1021', note='Liwilla x Begra DH; %s.' % ds,
              evidence='Table 4: %s %s %s LOD %s R2 %s' % (nm, ch, '-'.join(ms), lod, pve)), ms, ch)

# ---------- Abeysekara et al. 2012 Phytopathology 102:94-104 (Salamouni x Katepwa RIL), Table 5 ----------
add(name='QSnb.fcu-1AS (Snn4)', cls='QTL', chr='1A', marker='Snn4', pve='23', stat='LOD 7.03', source='V_Abeysekara2012', grade='C',
    peak=7.227, note='Salamouni x Katepwa RIL, isolate Sn99CH 1A7a; QTL defined by the Snn4 locus. Placed at the Snn4 position already in the atlas (Abeysekara et al. 2009, 7.23 Mb).',
    evidence='Table 5: 1AS QSnb.fcu-1AS Snn4 R2 0.23 LOD 7.03')
place(add(name='QSnb.fcu-7AS', chr='7A', marker='XBF293121', pve='16', stat='LOD 4.64', source='V_Abeysekara2012',
          note='Salamouni x Katepwa RIL; peak at EST-STS XBF293121 (primers from the paper\'s Table 2).', evidence='Table 5: 7AS QSnb.fcu-7AS XBF293121 R2 0.16 LOD 4.64'), ['bf293121'], '7A')

# ---------- Chu et al. 2010 TAG 120:1451-1459 (TA4152-60 x ND495 DH), Table 3 ----------
add(name='Snn1 locus (Chu 2010)', cls='QTL', chr='1B', start=2.356, end=2.357, marker='Xfcp566–Xfcp596', pve='22-50', stat='LOD 9.98-20.26', source='Z_1090', grade='B',
    note='TA4152-60 x ND495 DH; QTL at Snn1 (Snn1-SnTox1). Placed at the cloned Snn1 gene (atlas record from Shi et al. 2016).', evidence='Table 3: Snn1 Xfcp566-Xfcp596 R2 0.22/0.50')
add(name='Tsn1 locus (Chu 2010)', cls='QTL', chr='5B', start=546.698, end=546.843, marker='Xgwm260–Xbarc73', pve='28', stat='LOD 12.61', source='Z_1090', grade='B',
    note='TA4152-60 x ND495 DH; QTL at Tsn1 (Tsn1-SnToxA), isolate Sn2000. Placed at the Tsn1 gene.', evidence='Table 3: Tsn1 Xgwm260-Xbarc73 R2 0.28 LOD 12.61')
for nm, ch, l, r, pve, lod in [('QSnb.fcu-3D.1', '3D', 'Xgwm161', 'Xbarc6', '8', '4.57'), ('QSnb.fcu-3D.2', '3D', 'Xcfa2134', 'Xfcp597', '4-9', '2.52-3.39'),
                               ('QSnb.fcu-4A', '4A', 'Xcfa2121', 'Xcfa2173', '6', '3.04'), ('QSnb.fcu-5A.1', '5A', 'Xcfa2250', 'Xfcp497', '5-8', '2.83-3.09'),
                               ('QSnb.fcu-5A.2', '5A', 'Xbarc1061', 'Xcfa2185', '6', '3.04')]:
    place(add(name=nm, chr=ch, marker='%s–%s' % (l, r), pve=pve, stat='LOD ' + lod, source='Z_1090', note='TA4152-60 x ND495 DH; resistance from TA4152-60.',
              evidence='Table 3: %s %s-%s R2 %s LOD %s' % (nm, l, r, pve, lod)), [l, r], ch)

# ---------- Francki et al. 2011 Phytopathology 101:1278-1284, Table 4 ----------
for nm, ch, l, r, pk, lod, r2, pop in [('QSnl.daw-2A (2004)', '2A', 'wPt2448', 'wPt7056', 'wPt9432', '6.4', '11.2', 'P9819RB1, flag leaf 2004'),
                                       ('QSnl.daw-2A (2005)', '2A', 'gwm614a', 'wPt9432', 'barc124', '13.3', '21.8', 'P9819RB1, flag leaf 2005'),
                                       ('QSnl.daw-1B', '1B', 'wPt8949', 'gwm264', 'wPt8949', '10.8', '15.8', '05Y001, flag leaf 2007'),
                                       ('QSnl.daw-4B', '4B', 'barc163', 'wPt0391', 'wmc349', '5.1', '8.2', '05Y001, flag leaf 2007'),
                                       ('QSnl.daw-5B', '5B', 'wPt4628', 'wPt1733', 'fcp001', '8.8', '13.9', '05Y001, flag leaf 2007')]:
    rec = add(name=nm, chr=ch, marker='%s–%s (peak %s)' % (l, r, pk), pve=r2, stat='LOD ' + lod, source='Z_1039', note='%s.' % pop,
              evidence='Table 4: %s %s-%s max LOD %s (%s) R2 %s' % (nm, l, r, lod, pk, r2))
    place(rec, [l, r], ch, peak_marker=pk)

# ---------- Lu and Lillemo 2014 TAG 127:2635-2644 (SHA3/CBRD x Naxos), Table 3 ----------
for ch, mk, r2, src in [('1B', 'wmc619', '2.3-7.1', 'Naxos'), ('3A', 'gwm2', '4.4-9.1', 'Naxos'), ('3B', 'wPt-4127', '2.6-11.7', 'SHA3/CBRD'),
                        ('3B', 'wPt-4933', '2.5-11.7', 'Naxos'), ('5B', 'wPt-5346', '3.7', 'SHA3/CBRD'), ('5B', 'fcp1', '2.6-6.5', 'SHA3/CBRD'),
                        ('7A', 'wmc603', '3.3-8.3', 'Naxos'), ('7B', 'wPt-0963', '2.9-8.8', 'Naxos')]:
    arm = {'wPt-4933': '3BL', 'wPt-5346': '5BS', 'fcp1': '5BL', 'gwm2': '3AS'}.get(mk, ch)
    rec = add(name='QSnb-%s (%s; Lu and Lillemo 2014)' % (arm, mk), chr=ch, marker=mk, pve=r2, source='A_SNB_022',
              note='SHA3/CBRD x Naxos RIL, adult-plant leaf blotch; resistance from %s; closest marker only.' % src,
              evidence='Table 3: %s %s R2 %s source %s' % (arm, mk, r2, src))
    if mk == 'fcp1':
        rec.update(peak=546.15, grade='C'); rec['note'] += ' Xfcp1 (Tsn1-linked) at 546.15 Mb (Kumar et al. 2016 / atlas marker table).'
    else:
        place(rec, [mk], ch)

# ---------- Li et al. 2021 IJMS 22:5580 (EGA Blanco x Millewa), Table 2: v1.0 bp printed ----------
for nm, ch, mk, p, lod, pve in [('QSnl07.daw-1B', '1B', 'IWB64368', 4.320856, '18.6', '25.8'), ('QSnl07.daw-1B', '1B', '1BR_NRC_5A', 1.782029, '13.6', '20.3'),
                                ('QSnl08.daw-1B', '1B', 'X3476283 (Snn1)', 2.357713, '10.1', '17.1'), ('QSnl08.daw-1B', '1B', '1BR_RGA_4B', 2.004576, '10.0', '16.9'),
                                ('QSnl08.daw-5B', '5B', '5B_RGA_SNP_3A1', 562.793027, '4.3', '6.9')]:
    add(name='%s, %s' % (nm, mk), cls='MTA', chr=ch, peak=round(p, 4), marker=mk, pve=pve, stat='LOD ' + lod, source='Z_1239', grade='A',
        note='EGA Blanco x Millewa DH; R-gene marker associated with the flag-leaf QTL; physical position printed in Table 2 (RefSeq v1.0).',
        evidence='Table 2: %s %s %.6f Mb LOD %s PVE %s' % (nm, mk, p, lod, pve))
add(name='QSnl07.daw-5B, Xfcp623 (Tsn1)', cls='MTA', chr='5B', start=546.698, end=546.843, marker='Xfcp623 (Tsn1)', pve='10.1', stat='LOD 7.5', source='Z_1239', grade='B',
    note='EGA Blanco x Millewa DH; Tsn1 marker (no position printed); placed at the Tsn1 gene.', evidence='Table 2: QSnl07.daw-5B Xfcp623 (Tsn1) LOD 7.5 PVE 10.1')

# ---------- Thapa et al. 2016 BMC Res Notes 9:58: QSng.sfr-3BS order relative to Fhb1 and Sr2 ----------
add(name='QSng.sfr-3BS (fine order, Thapa 2016)', chr='3B', marker='Xsnp3BS-3 / g01130 (about 1 cM distal); Xsnp3BS-2, Xumn10, Xsnp3BS-9', pve='0.9-1.4', source='Z_1440', grade='C',
    peak=0.806, note='Arina x (Alsen x Ocoroni86) F2, 1600 plants; QSng.sfr-3BS lies between Sr2 (distal) and Fhb1 (proximal) on 3BS. Position taken from Xgwm389 (0.81 Mb), the QTL marker of Schnurbusch et al. 2003; the KASP markers could not be placed here.',
    evidence='Abstract and Table 3: QSng.sfr-3BS located between Sr2 and Fhb1')

# ---------- Hu et al. 2019 Crop J 7:809-818, Tables 2-3 (DArTseq) ----------
for pop, dis, ch, cm, l, r, pve, src in [
        ('BxC', 'Tan spot', '3B', '83.9-93.8', '1056249', '100141919', '2.5-3.7', 'CIANO T79'), ('BxC', 'Tan spot', '4D', '0-9.8', 'BS00036421_51', '1119387', '2.8-4.1', 'CIANO T79'),
        ('BxC', 'Tan spot', '5A', '175.9-178.7', 'Vrn-A1', '3064415', '11.2-14.5', 'BARTAI'), ('BxC', 'Tan spot', '5B', '112.3-115.1', '1062126', '1209358', '3.7', 'CIANO T79'),
        ('BxC', 'Tan spot', '5B', '137.6-146.1', '4911220', '1108969', '16.7-25.9', 'BARTAI'),
        ('BxC', 'Septoria nodorum blotch', '2B', '12.9-16.7', '1068352', '1024275', '3.4-3.6', 'CIANO T79'), ('BxC', 'Septoria nodorum blotch', '4D', '33.5-48.7', '1072422', '1069862', '6.9-7.1', 'CIANO T79'),
        ('BxC', 'Septoria nodorum blotch', '5A', '173.9-175.9', '1204040', 'Vrn-A1', '14.6-16.8', 'BARTAI'), ('BxC', 'Septoria nodorum blotch', '5B', '137.6-146.1', '4911220', '1108969', '2.2-6.4', 'BARTAI'),
        ('CxC', 'Tan spot', '3A', '196.4-202.4', '1076524', '990692', '3.7-6.8', 'CASCABEL'), ('CxC', 'Tan spot', '5A', '108.0-111.4', '2278040', '7340557', '4.4-5.1', 'CASCABEL'),
        ('CxC', 'Tan spot', '5A', '44.2-46.3 (LG 5A.2)', '2341646', '1141498', '2.5-6.2', 'CASCABEL'), ('CxC', 'Tan spot', '5B', '46.7-50.6', '100168878', '1090946', '3.7-6.5', 'CIANO T79'),
        ('CxC', 'Tan spot', '5D', '0-1.1', '2333933', '1273057', '3.3-6.6', 'CASCABEL'),
        ('CxC', 'Septoria nodorum blotch', '1B', '186.0-191.4', '100184803', '100182947', '4.1-5.3', 'CIANO T79'), ('CxC', 'Septoria nodorum blotch', '2A', '0-28.0', '100142576', '4989077', '2.7-8.3', 'CASCABEL'),
        ('CxC', 'Septoria nodorum blotch', '2B', '247.7-250.8', '1094999', '1124944', '3.9-6.9', 'CASCABEL'), ('CxC', 'Septoria nodorum blotch', '5A', '43.1-46.3 (LG 5A.2)', '4911066', '1141498', '5.9-9.7', 'CASCABEL'),
        ('CxC', 'Septoria nodorum blotch', '6B', '87.0-97.9', '1080914', '1228394', '4.2-8.3', 'CASCABEL'), ('CxC', 'Septoria nodorum blotch', '6D', '19.6-22.6', '994751', '979610', '4.4-5.6', 'CIANO T79')]:
    pre = 'QTs' if dis == 'Tan spot' else 'QSnb'
    rec = add(name='%s-%s-%s %s cM (Hu 2019)' % (pre, pop, ch, cm.split(' ')[0]), disease=dis, chr=ch, marker='%s–%s' % (l, r), pve=pve, source='V_Hu2019',
              note='%s RIL (%s), seedling; resistance from %s.' % ({'BxC': 'BARTAI x CIANO T79', 'CxC': 'CASCABEL x CIANO T79'}[pop], cm, src),
              evidence='Table %s: %s %s %s cM %s %s PVE %s' % ('2' if pop == 'BxC' else '3', dis, ch, cm, l, r, pve))
    ms = [m for m in (l, r) if m != 'Vrn-A1']
    place(rec, ms, ch)
    if 'Vrn-A1' in (l, r):
        g = GENE['Vrn-A1'][1]
        if rec['peak'] is not None and rec['start'] is None:
            rec['start'], rec['end'], rec['peak'] = min(rec['peak'], g), max(rec['peak'], g), None; rec['grade'] = 'B'
            rec['note'] += ' Vrn-A1 placed at the gene (587.40 Mb).'
        elif rec['start'] is None:
            rec['peak'] = g; rec['grade'] = 'C'
            other = [m for m in ms][0]
            rec['note'] = re.sub(r' No marker of this locus could be placed on 5A; left unplaced\.', ' Placed at the Vrn-A1 gene (587.40 Mb); %s could not be placed on 5A.' % other, rec['note'])
    if ch == '5A' and float(cm.split('-')[0]) > 170:
        rec['note'] += ' The paper places the common 5AL TS/SNB region at 582.6-589.3 Mb (RefSeq v1.0).'

# ---------- Francki et al. 2018 Genome 61:145-149 (90K SNP fine mapping), glume QTL on 2D named in text ----------
add(name='QSng.daw-2D (P92201D5/P91193D1, USA)', chr='2D', marker='IWB31450–IWB21124 (0.3 cM)', source='V_Francki2018',
    note='P92201D5/P91193D1; glume resistance co-located across three US environments within this SNP interval. Other QTL of the paper are shown only in Figs 1-2.',
    evidence='Text: glume resistance on 2D within the SNP interval IWB31450-IWB21124 (0.3 cM)')
place(R[-1], ['IWB31450', 'IWB21124'], '2D')

# ---------- Ruud et al. 2017 TAG 130:1361-1374 (SHA3/CBRD x Naxos), Tables 5 and 7 ----------
for ch, mk, pve, src, kind in [('1A', 'RAC875_c10083_800', '11.7', 'Naxos', 'seedling'), ('1B', 'psp3000', '10.4', 'SHA3/CBRD', 'seedling'),
                               ('2D', 'wsnp_RFL_Contig3960_4401914', '11.1', 'Naxos', 'seedling'), ('5B', 'BS00091518_51', '27.5-51.8', 'SHA3/CBRD', 'seedling, Snn3'),
                               ('7B', 'wsnp_BE498662B_Ta_2_5', '15.5', 'Naxos', 'seedling'),
                               ('1A', 'wsnp_Ex_c25734_34995416', '2.4-10.3', 'SHA3/CBRD', 'field'), ('1B', 'SCM9', '5.2-8.1', 'Naxos', 'field, 1RS'),
                               ('3A', 'gwm2', '3.7-11.5', 'Naxos', 'field, 3AS.1'), ('3A', 'Ku_c41007_116', '2.2-9.4', 'SHA3/CBRD', 'field, 3AS.2'),
                               ('3B', 'BS00030534_51', '5.7', 'SHA3/CBRD', 'field, 3BS.1'), ('3B', 'wBE445348B_Ta_2_1', '6.9', 'SHA3/CBRD', 'field, 3BS.2'),
                               ('3B', 'wPt-4933', '3.5-11.2', 'Naxos', 'field, 3BL'), ('5B', 'BS00091518_51', '4.7-24.0', 'SHA3/CBRD', 'field, 5BS (Snn3)'),
                               ('5B', 'wPt-5914', '2.4-5.6', 'SHA3/CBRD', 'field, 5B.2'), ('7A', 'RAC875_c14195_1155', '2.9-6.5', 'Naxos', 'field'),
                               ('7B', 'BobWhite_rep_c50229_413', '2.7-8.4', 'Naxos', 'field')]:
    rec = add(name='QSnb-%s %s (%s; Ruud 2017)' % (ch, kind.split(',')[0], mk), chr=ch, marker=mk, pve=pve, source='Z_1360',
              note='SHA3/CBRD x Naxos RIL; %s; resistance from %s; closest marker only.' % (kind, src), evidence='Table %s: %s %s %s' % ('5' if kind.startswith('seedling') else '7', ch, mk, pve))
    place(rec, [mk], ch)

# ---------- Arseniuk et al. 2004 Euphytica 137:225-231 (Alba x Begra DH) ----------
add(name='QSnl.ihar-6A', chr='6A', peak=579.126, marker='Xgwm570 (Xmwg934)', pve='14 (INC), 36 (DIS)', stat='LOD 4.08 (INC), 7.26 (DIS)', source='Arseniuk2004', grade='C',
    note='Leaf blotch, partial resistance components; peak at Xgwm570 (Fig. 2A); placed at Xgwm570.', evidence='Results: QSnl.ihar-6A LOD 4.08 INC, 7.26 DIS; 14% and 36%')
place(add(name='Putative 6DL QTL (Arseniuk 2004)', chr='6D', marker='gdm98–gdm108–Xbarc175', pve='10 (INC), 8 (DIS)', stat='LOD 2.85, 2.23 (below threshold)', source='Arseniuk2004',
          note='Putative QTL on distal 6DL, below the LOD 3 threshold.', evidence='Results: putative QTL on 6DL, LOD 2.85 and 2.23, 10% and 8%'), ['gdm98', 'gdm108', 'barc175'], '6D')

# ---------- Liu et al. 2004 Phytopathology 94:1061-1067 (ITMI W-7984 x Opata 85), Table 3 ----------
add(name='snn1 (Liu 2004)', cls='QTL', chr='1B', start=2.356, end=2.357, marker='snn1', pve='27.3-58.3', source='Liu2004Phyto', grade='B',
    note='ITMI RI population, seedling, isolate Sn2000; snn1 explained most of the variation. Placed at the cloned Snn1 gene.', evidence='Table 3: snn1 1BS R2 0.583/0.480/0.273')
for mk, ch, arm, r2, src in [('Xcdo1312', '4B', '4BL', '6-9', 'W-7984'), ('XksuG12', '4A', '4AL', '4-7', 'W-7984'), ('Xksu912 (Prp)', '3A', '3AS', '8', 'W-7984'),
                             ('Xcdo1508', '5D', '5DL', '8-9', 'W-7984'), ('Xbcd361', '3D', '3DL', '4-7', 'W-7984'), ('Xbcd758', '6A', '6AL', '4-5', 'W-7984'),
                             ('Xksu927 (Grp94)', '7B', '7BL', '4-6', 'Opata 85')]:
    add(name='QSnb %s (%s; Liu 2004)' % (arm, mk.split(' ')[0]), chr=ch, marker=mk, pve=r2, source='Liu2004Phyto', grade='D',
        note='Minor locus from single-marker regression; RFLP marker, not placed on RefSeq. Resistance from %s.' % src, evidence='Table 3: %s %s R2 %s' % (mk, arm, r2))

# ---------- Gonzalez-Hernandez et al. 2009 Euphytica 166:199-206 (LDN x LDN(Dic-5B) RICL) ----------
add(name='QSnb.ndsu-5B', chr='5B', start=504.302, end=546.698, marker='Xbcd9–Xfbb237 (proximal to tsn1)', pve='37.6', source='GonzalezHernandez2009', grade='C',
    note='T. turgidum var. dicoccoides 5B; resistance from dicoccoides. RFLP flanks not placed; interval bracketed by Xgwm639 (504.30 Mb, proximal on the paper\'s map) and Tsn1 (546.70 Mb, distal).',
    evidence='Results: QSnb.ndsu-5B in Xbcd9-Xfbb237 on 5BL, proximal to tsn1, 37.6%')

# ---------- Reszka et al. 2007 Plant Pathol Bull 16:161-167 (triticale Bogo DH) ----------
for nm, ch, ms, pve in [('QTL 3B (Reszka 2007)', '3B', ['barc132', 'barc344'], '10 (DIS), 14 (LAT)'),
                        ('QTL 4B (Reszka 2007)', '4B', ['barc193', 'gwm149', 'barc106', 'gwm495', 'gwm513'], '19 (INC), 30 (DIS), 14 (LAT)'),
                        ('QTL 6A (Reszka 2007)', '6A', ['barc107', 'barc171', 'barc195', 'barc206', 'gwm570'], '6 (INC)')]:
    place(add(name=nm, chr=ch, marker=', '.join(ms), pve=pve, source='Reszka2007PPB',
              note='Winter triticale cv. Bogo DH; partial resistance components. Wheat SSRs on the triticale A/B chromosome placed on wheat RefSeq v1.0; interval spans the associated markers.',
              evidence='Results: QTL on %s with %s; %s%%' % (ch, ', '.join(ms), pve)), ms, ch)

# ---------- Zhang et al. 2009 Mol Breed 23:? (BR34 x Grandin): Snn2 and Tsn1 markers ----------
place(add(name='Snn2 (Zhang 2009)', cls='SG', chr='2D', marker='XTC253803–Xcfd51 (Xcfd56 distal)', source='Zhang2009MB',
          note='BR34 x Grandin RIL; SnTox2 sensitivity gene delineated to a 4.0 cM interval on 2DS.', evidence='Abstract and Fig.: Snn2 between XTC253803 and Xcfd51, 4.0 cM'),
      ['cfd56', 'cfd51'], '2D')
add(name='Tsn1 (Zhang 2009 markers)', cls='SG', chr='5B', start=546.698, end=546.843, marker='Xfcp1–Xfcp620–Tsn1–Xfcp394 (1.0 cM)', source='Zhang2009MB', grade='B',
    note='BR34 x Grandin; Tsn1 delineated to 1.0 cM with fcp markers; placed at the Tsn1 gene.', evidence='Abstract: SSRs delineated Tsn1 to a 1.0 cM interval')
# ---------- Shi et al. 2016 Mol Genet Genomics 291:107-119: high-resolution mapping of Snn3-B1 ----------
add(name='Snn3-B1 (Shi 2016 high-resolution map)', cls='SG', chr='5B', start=6.648, end=6.975, marker='Sumai 3 F2 saturation and high-resolution maps', source='Shi2016MGG', grade='B',
    note='Snn3-B1 high-resolution mapping in two Sumai 3 F2 populations; placed on the Snn3-B1 interval already in the atlas (Downie et al. 2018, RefSeq v1.0).',
    evidence='Abstract: saturation and high-resolution mapping of Snn3-B1')

if __name__ == '__main__':
    import collections
    print(len(R), collections.Counter(r['grade'] for r in R))
    for x in R: print(x['source'][:14].ljust(14), x['name'][:42].ljust(42), x['chr'], x['start'], x['end'], x['peak'], x['grade'])
    print('\nunplaced:'); [print(' ', *l) for l in LOG]
