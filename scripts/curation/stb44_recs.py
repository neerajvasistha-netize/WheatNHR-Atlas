"""Records re-read from the STB Archive.zip papers (Nov 2026). One tuple per locus:
  (source_key, name, chr, markers, peak_marker, fixed, pve, stat, note, evidence, extra)
markers: placed with the sb26 resolver; fixed: dict of printed positions (start/end/peak, ref, grade) that override markers;
extra: dict of other fields (cls, disease, gene)."""
STB = 'Septoria tritici blotch'
R = []
def r(src, name, ch, markers=(), peak=None, fixed=None, pve='', stat='', note='', ev='', **extra):
    R.append(dict(source=src, name=name, chr=ch, markers=list(markers), peak_marker=peak, fixed=fixed or {}, pve=pve, stat=stat,
                  note=note, evidence=ev, extra=extra))

# ---- Adhikari et al. 2004 Phytopathology 94:1198 (Stb4, file Adhikari2012.pdf)
r('A_STB_042', 'Stb4', '7D', ['Xgwm111'], cls='DG', note='Tadinia x (Yecora Rojo x UC554) RIL; Stb4 on 7DS, Xgwm111 0.7 cM and AFLP EcoRI-ACTG/MseI-CAAA5 2.1 cM.',
  ev='Abstract: Xgwm111 and EcoRI-ACTG/MseI-CAAA5 linked to Stb4 at 0.7 and 2.1 cM')
# ---- Adhikari et al. 2015 Euphytica 205:699
r('V_Adhikari2015', 'QStb.2A (Adhikari 2015)', '2A', ['wPt_3896', 'wPt_4197'], pve='10.8', note='Largo x ND495 RIL (118 lines); resistance from Largo; CIM peak at 12 cM, flanked by wPt-3896 and wPt-4197.',
  ev='Table 1 and text: QStb.2A flanked by wPt_3896 and wPt_4197, R2 10.8%')
# ---- Arraiano et al. 2007 Plant Pathol 56:73 (Stb15)
r('V_Arraiano2006', 'Stb15 (Arraiano 2007)', '6A', [], cls='DG', fixed=dict(start=48.527, end=48.527, grade='C'),
  note='Arina x Forno RIL; Stb15 on 6AS (MQM, RFLP Xpsr904 14 cM). Placed at the cloned Stb15 position already in the atlas (Hafeez et al. 2025, 48.53 Mb).',
  ev='Abstract: gene for resistance to IPO88004 located on 6AS, named Stb15')
# ---- Chartrain et al. 2004 Phytopathology 94:497
r('Chartrain2004Phyto', 'QStb.psr-6B-1', '6B', ['Xgwm133'], note='Arina x Riband DH; partial adult-plant resistance from Riband; closest SSR Xgwm133 (t test P 5e-4).',
  ev='Table 3: Xwgm133 Diff 0.83 P 5e-4; abstract: QStb.psr-6B-1 from Riband on 6B')
r('Chartrain2004Phyto', 'Stb6 (Arina; Chartrain 2004)', '3A', [], cls='DG', fixed=dict(peak=26.2035, grade='C'),
  note='Arina x Riband DH; resistance to IPO323 linked to Xgwm369 (10% recombination), probably Stb6; placed at the cloned Stb6 gene.',
  ev='Results: Stb6 and Xgwm369 recombination fraction 10%')
# ---- Chartrain et al. 2005 TAG 110:1138 (TE 9111)
r('A_STB_045', 'Stb11 (TE 9111)', '1B', ['Xbarc008'], cls='DG', pve='50', note='TE 9111 x Exeter F4 SSD; Stb11 on 1BS, tightly linked to Xbarc008.',
  ev='Abstract and results: Stb11 1BS, closest marker Xbarc008')
r('A_STB_045', 'Stb6 (TE 9111)', '3A', [], cls='DG', fixed=dict(peak=26.2035, grade='C'), note='TE 9111 carries Stb6 (IPO323), linked to Xgwm369; placed at the cloned Stb6 gene.',
  ev='Abstract: the Stb6 gene for resistance to isolate IPO323')
r('A_STB_045', 'Stb7 (probable; TE 9111)', '4A', ['Xwmc313'], cls='DG', pve='11', note='Resistance to IPO87019 on distal 4AL near Xwmc313, probably Stb7 (the atlas had called it Stb12).',
  ev='Results: TE 9111 allele of Xwmc313 closest to the IPO87019 resistance gene; LOD 1.96, 11%')
# ---- Chartrain et al. 2009 Plant Pathol 58:547 (Stb9)
r('A_STB_003', 'Stb9', '2B', ['Xwmc317', 'Xbarc0129', 'Xgpw1214', 'Xgpw4044', 'Xgpw5039'], cls='DG',
  note='Courtot x Chinese Spring; Stb9 on 2BL between RFLPs Xfbb226 (3.6 cM) and XksuF1b (9 cM); SSRs within 10 cM: Xgpw1214, Xgpw4044, Xgpw5039, Xwmc317, Xbarc0129.',
  ev='Abstract and discussion: Stb9 between Xfbb226 and XksuF1b; SSRs within 10 cM listed')
# ---- Cobertera et al. 2022 (Research Square preprint) Foote x Madsen; Table 3 bp positions (IWGSC RefSeq v1)
for q, c, s, e in [('QStb-2A', '2A', 8735846, 8735909), ('QStb-2B', '2B', 87673429, 87673492), ('QStb-4B', '4B', 597979606, 597979669),
                   ('QStb-5A', '5A', 509569472, 509569535), ('QStb-6B', '6B', 2609365, 2609428), ('QStb-6D', '6D', 61003540, 61003598),
                   ('QStb-7D', '7D', 614103520, 614103583)]:
    r('V_Cobertera2022', q + ' (Cobertera 2022)', c, [], fixed=dict(peak=round(s / 1e6, 4), grade='A'),
      note='Foote x Madsen RIL, natural STB in Oregon 2016-2019; physical position of the peak marker printed in Table 3 (IWGSC RefSeq v1).',
      ev='Table 3: %s %s %d-%d bp' % (q, c, s, e))
# ---- Dreisigacker et al. 2015 TAG 128:2317 (spelt-derived, Table 4)
for q, c, m, pv, ld, pop in [('QStb.cim-1BS', '1B', ['Xwmc230', 'Ku_c16938_25916260'], '11.7', '6.2', 'NASMA x IAS20*5/H567.71'),
                             ('QStb.cim-5AL-1', '5A', ['Ex_rep_c68441_67261679', 'Xgwm186'], '5.4', '3.1', 'NASMA x IAS20*5/H567.71'),
                             ('QStb.cim-3BL', '3B', ['Ex_c19778_28779907', 'Xwmc291'], '22.7', '12.2', 'NASMA x RPB709.71/COC'),
                             ('QStb.cim-5AL-2', '5A', ['Ex_c621_1231444', 'Ex_c5998_10513766'], '18.5', '10.3', 'NASMA x RPB709.71/COC'),
                             ('QStb.cim-7AS', '7A', ['Xbarc174', 'Ra_c2094_4066675'], '5.1', '3.0', 'NASMA x RPB709.71/COC')]:
    r('Dreisigacker2015TAG', q, c, m, pve=pv, stat='LOD ' + ld, note=pop + ' RIL (spelt-derived); adult-plant STB in at least three of four environments.',
      ev='Table 4: %s %s-%s LOD %s PVE %s%%' % (q, m[0], m[1], ld, pv))
# ---- Ferjaoui et al. 2022 BMC Genomics 23:372 (durum Agili39 x Khiar)
for q, c, m, pv, ld in [('Qstb1A', '1A', ['Tag_1694925', 'Tag_1127081'], '9.1-16.0', '4.1-4.6'), ('Qstb7A', '7A', ['Tag_2277193', 'Tag_100009953'], '9.0', '4.0'),
                        ('Qstb2B_1', '2B', ['Tag_1056626', 'Tag_1117572'], '10.0-61.6', '3.2-33.9'), ('Qstb2B_2', '2B', ['Tag_100031118', 'Tag_3027184'], '14.9-54.3', '4.8-32.8'),
                        ('Qstb5B', '5B', ['Tag_1125523', 'Tag_490992'], '', 'epistatic')]:
    r('V_Ferjaoui2021', q + ' (Ferjaoui 2022)', c, m, pve=pv, stat='LOD ' + ld,
      note="Durum 'Agili39' x Khiar RIL; DArTseq flanking tags; physical positions only in the paper's additional tables.",
      ev='Table 6 / text: %s %s-%s' % (q, m[0], m[1]))
# ---- Tabib Ghaffary et al. 2011 TAG 123:741 (Apache x Balance DH)
r('Ghaffary2011TAG', 'Stb18', '6D', ['Xgpw5176', 'Xgpw3087'], cls='DG', pve='8.9-67.5', stat='LOD 3.1-23.2',
  note='Apache x Balance DH, seedling; new gene from Balance, flanked by Xgpw5176 and Xgpw3087 (6DS).', ev='Table 4 and abstract: Stb18 6DS Xgpw5176-Xgpw3087')
r('Ghaffary2011TAG', 'Stb6 (Balance)', '3A', [], cls='DG', fixed=dict(peak=26.2035, grade='C'), pve='27.7-73.7', stat='LOD 7.3-25.6',
  note='3AS QTL at wPt-0836 for IPO323, most likely Stb6; placed at the cloned Stb6 gene.', ev='Table 4: IPO323 wPt-0836 3AS')
r('Ghaffary2011TAG', 'Stb11 (Apache; 1BS QTL)', '1B', ['wPt-2019'], cls='QTL', pve='59.3-68.3', stat='LOD 17.8-21.1',
  note='1BS QTL at wPt-2019 for IPO87016, most likely Stb11.', ev='Table 4: IPO87016 wPt-2019 1BS')
r('Ghaffary2011TAG', 'QStb 7DS (Apache; Stb4/Stb5 region)', '7D', ['Xgwm111'], pve='5.9-27.5', stat='LOD 2.2-9.5',
  note='7DS QTL at Xgwm111 (IPO98022, IPO98046), most probably Stb4 or Stb5.', ev='Table 4: Xgwm111 7DS')
r('Ghaffary2011TAG', 'QStb 7DL (Apache)', '7D', ['wPt-1859', 'Xgpw313'], pve='8-20.5', stat='LOD 4.5-6.8',
  note='7DL QTL (wPt-1859, Xgpw313) for IPO89011 and IPO98046; possibly a new Stb gene.', ev='Table 4: wPt-1859 7DL; Xgpw313 7DL')
# ---- Tabib Ghaffary et al. 2012 TAG 124:125 (Kulm x M3; Stb16q, Stb17)
r('A_STB_005', 'Stb16q', '3D', [], cls='DG', fixed=dict(peak=588.0, grade='C'), pve='28-71', stat='LOD 7.2-27.0',
  note='Kulm x SHW M3 RIL; broad-spectrum seedling and adult resistance on 3DL (Xbarc125-Xbarc128); placed at the cloned Stb16q position already in the atlas (588.0 Mb).',
  ev='Tables 5-6: Stb16q 3DL Xbarc125-Xbarc128')
r('A_STB_005', 'Stb17', '5A', ['Xgwm617', 'Xhbg247'], cls='DG', pve='12-32', stat='LOD 3.0-8.9',
  note='Kulm x SHW M3 RIL; adult-plant resistance on 5AL (Xgwm617-Xhbg247).', ev='Table 6: Stb17/5AL Xgwm617-Xhbg247')
# ---- Goodwin et al. 2015 Crop Sci (Stb3)
r('Goodwin2015CS', 'Stb3', '7A', ['Xwmc83', 'Xcfa2028', 'Xbarc222'], cls='DG',
  note='Israel 493 x Chinese Spring; Stb3 on 7AS, completely linked to Xwmc83, flanked by Xcfa2028 (12.4 cM distal) and Xbarc222 (2.1 cM proximal); not on 6DS as reported before.',
  ev='Abstract: Stb3 7AS, Xwmc83, Xcfa2028-Xbarc222')
# ---- He et al. 2021 Phytopathology 111:1001 (Murga x Huirivis#1)
for q, c, m, pv in [('2B', '2B', ['988088', '7335682'], '1.6-3.9'), ('2D.1', '2D', ['3961111', '2245461'], '2.1-2.3'), ('2D.2', '2D', ['1091926', '4540601'], '1.7-3.6'),
                    ('2D.3', '2D', ['3534347', '1076603'], '6.6-7.2'), ('3A', '3A', ['1094665', '1073870'], '4.2'), ('3B', '3B', ['1211784', '3955846'], '5.3'),
                    ('3D', '3D', ['1073418', '3956057'], '27.5-62.5'), ('5B.1', '5B', ['2263706', '100045912'], '2.4-12.4'), ('5B.2', '5B', ['1209730', '4989083'], '1.7')]:
    fx = dict(start=588.3, end=590.1, grade='A') if q == '3D' else {}
    r('Yang2020Phyto', 'QStb-%s (He 2021)' % q, c, m, fixed=fx, pve=pv, note='Murga x Huirivis#1 RIL, adult plants in Mexico and Uruguay.' + (' Interval printed in the text (588.3-590.1 Mb), within the Stb16 region.' if q == '3D' else ''),
      ev='Table 4: %s %s-%s PVE %s%%' % (q, m[0], m[1], pv))
# ---- Tamburic-Ilincic et al. 2019 Mol Breed 39:94 (Maxine x FTHP Redeemer), STB only
for c, m, pv, ld, srcp in [('7A', ['1270330', '1372059'], '12.6', '3.37', 'Maxine'), ('4B', ['1242543', '3020840'], '12', '3.23', 'FTHP Redeemer'), ('2D', ['1103191', '988130'], '24.3', '7.58', 'FTHP Redeemer')]:
    r('Ilincic2019MB', 'QStb-%s (Tamburic-Ilincic 2019)' % c, c, m, pve=pv, stat='LOD ' + ld, note='Maxine x FTHP Redeemer DH; resistance from %s.' % srcp,
      ev='Table 4: %s %s-%s LOD %s %s%%' % (c, m[0], m[1], ld, pv))
# ---- Radecka-Janusik and Czembor 2014 Euphytica 200:109 (Liwilla x Begra), Table 3 union of flanks
for q, c, m, pv, ld in [('QStb.ihar-3A.2', '3A', ['wPt1688', 'gwm666'], '23.9-50.3', '7.15-9.69'), ('QStb.ihar-6A', '6A', ['wPt9679', 'wPt7027'], '9.5', '3.26'),
                        ('QStb.ihar-7A.2', '7A', ['wPt0114', 'wPt6184', 'wPt0040', 'wPt5153', 'wPt0429', 'wPt5590'], '13.5-30.3', '2.76-7.24'),
                        ('QStb.ihar-1B', '1B', ['wPt8779', 'wPt7188'], '15.8-40.0', '3.43-7.17'), ('QStb.ihar-2B.2', '2B', ['wPt8460', 'wPt7859'], '14.4', '3.80'),
                        ('QStb.ihar-3B', '3B', ['wPt1171', 'wPt1804', 'wPt7186', 'wPt9912'], '13.1-17.8', '3.13-3.62'), ('QStb.ihar-5D', '5D', ['gwm205', 'barc110', 'gwm174'], '26.3-35.1', '3.63-4.14')]:
    r('A_STB_013', q, c, m, pve=pv, stat='LOD ' + ld, note='Liwilla x Begra DH; interval spans the flanking markers of all Table 3 rows.', ev='Table 3: %s %s' % (q, '-'.join(m)))
# ---- Karlstedt et al. 2019 Euphytica 215:108 (HTRI1410 spelt), Table 5 Mbp (IWGSC RefSeq)
for q, c, s, e, p, pv, ld in [('QTL_5A', '5A', 388.1, 503.1, 494.8, '17.5', '6.83'), ('QTL_4B', '4B', 558.5, 641.5, 641.5, '8.8', '3.25'),
                              ('QTL_7B_1', '7B', 583.4, 605.9, 593.1, '8.5', '3.12'), ('QTL_7B_2', '7B', 608.7, 646.0, 614.2, '15.3', '6.11')]:
    r('Karlstedt2019Eu', q + ' (Karlstedt 2019)', c, [], fixed=dict(start=s, end=e, peak=p, grade='A'), pve=pv, stat='LOD ' + ld,
      note='DH lines from HTRI1410 (spelt) crossed with Solitaer, Alcedo and Jenga; field STB; physical positions of flanking and peak markers printed in Table 5.',
      ev='Tables 4-5: %s %.1f-%.1f Mb, peak %.1f Mb' % (q, s, e, p))
# ---- Kelm et al. 2012 Mol Breed 29:813 (Solitaer x Mazurka DH, seedling), Table 3
for q, c, m, gene, pv in [('QStb.4AL (Stb7/Stb12)', '4A', ['Xwmc313', 'Xgwm160'], '', '2.9-24.8'), ('QStb.3AS (Stb6)', '3A', ['Xgwm369'], 'Stb6', '19.1-42.7'),
                          ('QStb.1B.a (Stb11)', '1B', ['Xgwm752'], '', '3.7-17.8'), ('QStb.1B.b', '1B', ['Xgwm1078'], '', '2.4'), ('QStb.3DS', '3D', ['Xgwm1243'], '', '3.7-8.5'),
                          ('QStb.7DS (Stb4/Stb5)', '7D', ['Xgwm885'], '', '1.9-3.9'), ('QStb.6B', '6B', ['Xgwm1076', 'Xgwm219'], '', '4.0-6.5'),
                          ('QStb.3B (Stb2/Stb14)', '3B', [], '', '9.0'), ('QStb.1A', '1A', [], '', '7.1')]:
    fx = dict(peak=26.2035, grade='C') if gene == 'Stb6' else {}
    r('Kelm2012MB', q + ' (Kelm 2012)', c, m, fixed=fx, pve=pv, note='Solitaer x Mazurka DH, seedling tests with seven isolates.' + (' Placed at the cloned Stb6 gene.' if gene else '') + ('' if m or gene else ' Peak marker is an AFLP; not placed.'),
      ev='Table 3: %s %s' % (q, ', '.join(m) or 'AFLP'))
# ---- Langlands-Perry et al. 2022 Genes 13:100 (Renan x Chinese Spring), extents of the annotated genes in Table 6 (RefSeq v1.0 bp)
for q, c, s, e, n in [('Qstb.renan-1D', '1D', 7.277, 11.494, 'QTL 4.28 Mb long'), ('Qstb.renan-5D', '5D', 42.926, 341.193, 'QTL about 312 Mb long'),
                      ('Qstb.renan-7B', '7B', 720.131, 723.900, 'QTL 50.86 Mb long; Table 6 lists only two of its genes')]:
    r('Z_1223', q, c, [], fixed=dict(start=s, end=e, grade='B' if c != '7B' else 'C'),
      note='Renan x Chinese Spring RIL, ultra-dense map; %s. Interval spans the candidate genes listed for the QTL (RefSeq v1.0 coordinates in Table 6). Replaces 21 per-gene records.' % n,
      ev='Text and Table 6: %s' % q)
# ---- Liu et al. 2013 Euphytica 190:145 (Stb2)
r('A_STB_006', 'Stb2', '1B', ['Xwmc406', 'Xbarc008', 'Xwmc230'], cls='DG', note='Veranopolis x RAC875-2 DH; Stb2 re-located to 1BS, flanked by Xwmc406 (6 cM) and Xbarc008 (5 cM, Xwmc230 close).',
  ev='Abstract: Stb2 on 1BS flanked by Xwmc406 and Xbarc008')
# ---- Mergoum et al. 2013 Cereal Res Commun 41:199 (Steele-ND x ND 735), Table 3
for c, m, pv, ld in [('5B', ['XwPt7101', 'XwPt-250', 'X377410'], '7.1-10.2', '2.2-3.8'), ('1D', ['XwPt1799'], '13.0', '3.5'), ('7A', ['XwPt8043'], '5.5', '2.1')]:
    r('Mergoum2013CRC', 'QStb-%s (Mergoum 2013)' % c, c, m, pve=pv, stat='LOD ' + ld, note="'Steele-ND' x 'ND 735' RIL; seedling STB." + (' Stable 5BL QTL in all experiments.' if c == '5B' else ''),
      ev='Table 3: %s %s' % (c, '/'.join(m)))
# ---- Miedaner et al. 2012 Mol Breed 29:731, Table 1 (STB QTL)
for pop, c, m, ld, pv in [('Arina/Forno', '2B', ['Xpsr540', 'Xcfd276'], '5.5', '9.6'), ('Arina/Forno', '3B', ['Xcfa2134', 'Xgwm131'], '7.3', '14.2'), ('Arina/Forno', '5B', ['Xwmc473', 'Xpsr574'], '4.7', '6.7'),
                          ('Arina/Forno', '6D', ['Xcfd19', 'Xgdm14'], '6.6', '9.5'), ('Arina/Forno', '7B', ['XksuD2', 'Xgwm146'], '3.7', '3.0'), ('History/Rubens', '4D', ['Rht-D1', 'Xbarc105'], '16.6', '24.4'),
                          ('History/Rubens', '5B.1', ['XwPt-4996', 'Xgwm274'], '11.5', '12.0'), ('History/Rubens', '5B.2', ['Xbarc142'], '4.5', '11.3'), ('History/Rubens', '6B', [], '7.9', '8.2'),
                          ('History/Rubens', '7B', ['Xgwm263'], '10.4', '7.1')]:
    r('Miedaner2012MB', 'QStb-%s (%s; Miedaner 2012)' % (c, pop), c[:2], m, pve=pv, stat='LOD ' + ld, note='%s; adult-plant STB; broad-spectrum resistance study.' % pop + ('' if m else ' Flanks are AFLPs; not placed.'),
      ev='Table 1: %s %s %s LOD %s' % (pop, c, '-'.join(m) or 'AFLP flanks', ld))
# ---- Naz et al. 2015 Plant Breed 134:514 (advanced backcross B22, Z86), Table 3
for q, c, m in [('QStb.B22-1B', '1B', 'wsnp_Ex_c5101_9053178'), ('QStb.B22-1D.a', '1D', 'GENE_0014_822'), ('QStb.B22-1D.b', '1D', 'wsnp_Ku_c53270_57959459'), ('QStb.B22-1D.c', '1D', 'Xgdm111'),
                ('QStb.B22-2D', '2D', 'Excalibur_c2311_2260'), ('QStb.B22-3A', '3A', 'BS00049032_51'), ('QStb.B22-5A.a', '5A', 'RAC875_c30566_230'), ('QStb.B22-5A.b', '5A', 'Kukri_c49033_52'),
                ('QStb.B22-7B.a', '7B', 'RAC875_rep_c105494_293'), ('QStb.B22-7B.b', '7B', 'Kukri_c51101_351'), ('QStb.Z86-3A', '3A', 'BobWhite_c28950_147'), ('QStb.Z86-4B.a', '4B', 'Xwmc349'),
                ('QStb.Z86-4B.b', '4B', 'BS00022988_51'), ('QStb.Z86-5A', '5A', 'BobWhite_c51109_415'), ('QStb.Z86-7A', '7A', 'RAC875_c3450_836')]:
    r('Naz2015PB', q, c, [], peak=m, note='Advanced backcross population %s (synthetic wheat donor), adult-plant field STB; linked marker only.' % q.split('.')[1][:3],
      ev='Table 3: %s %s' % (q, m))
# ---- Odilbekov et al. 2019 Agronomy 9:510 (Nimbus x Stigg), positions printed in the text (RefSeq v1.0)
for q, c, s, e, m, pv in [('3DL', '3D', 593.4, 614.4, None, '6.7-10.5'), ('2BL', '2B', 789.4, 799.3, None, '13.9-22.7'), ('3AS', '3A', 7.4, 8.7, None, '5.9-10.2'),
                          ('1BS', '1B', None, None, ['AX-94509078', 'AX-94436250', 'AX-95162604'], '4.7-13.2')]:
    r('Agron2019STB', 'QStb-%s (Odilbekov 2019)' % q, c, m or [], fixed=dict(start=s, end=e, grade='A') if s else {}, pve=pv,
      note="Nimbus x Stigg DH; chlorosis, necrosis and pycnidia." + (' Physical interval printed in the text.' if s else ''), ev='Table 1 and discussion: %s' % q)
# ---- Piaskowska et al. 2021 Agronomy 11:1108 (Mandub x Begra), Table 4
for q, c, m, pv, fx in [('QStb.ihar-1B.2', '1B', [], '55.9-70.0', dict(start=48.2, end=49.6, grade='A')), ('QStb.ihar-2B.3', '2B', ['1221257', '1067532', '3947887'], '11.3-18.2', {}),
                        ('QStb.ihar-2B.4', '2B', ['5323643', '5410831', '3533738'], '29.3-31.4', {}), ('QStb.ihar-2B.5', '2B', ['3958573', '4009747'], '13.9-16.6', {}),
                        ('QStb.ihar-5B', '5B', ['4411210', '1117248'], '11.9', {})]:
    r('A_STB_024', q, c, m, fixed=fx, pve=pv, note='Mandub x Begra DH; seedling and adult tests.' + (' Interval printed in the text (48.2-49.6 Mbp).' if fx else ''),
      ev='Table 4: %s %s' % (q, ', '.join(m) or '1041519/1034868'))
# ---- Raman et al. 2009 Mol Breed 24:153 (StbWW)
r('Raman2009MB', 'StbWW (WW2449, WW1842, WW2451)', '1B', ['Xwmc230', 'Xbarc119', 'Xbarc008'], cls='DG', pve='up to 38 (adult)',
  note='Three DH populations; single seedling resistance gene on 1BS, within 7 cM of Xwmc230, Xbarc119b and Xksum045; distal bin 1BS.sat-4.',
  ev='Abstract: StbWW on 1BS, Xwmc230, Xbarc119b, Xksum045, Xbarc008')
# ---- Risser et al. 2011 Phytopathology 101:1209, Table 3
for pop, c, m, ld, pv in [('Florett/Biscay', '1A', ['Xwmc312', 'XwPt-7030'], '4.4', '3.0'), ('Florett/Biscay', '1B', ['XwPt-0260', 'Xwmc419'], '4.3', '2.9'), ('Florett/Biscay', '2B', ['Xwmc344'], '7.1', '4.8'),
                          ('Florett/Biscay', '3B', [], '20.7', '12.6'), ('Florett/Biscay', '4B', ['XwPt-8092'], '14.7', '9.3'), ('Florett/Biscay', '5B', ['Xgwm408', 'XwPt-4577'], '3.4', '2.3'),
                          ('Florett/Biscay', '6D', ['Xgwm469', 'Xcfd13'], '19.4', '11.9'), ('Florett/Biscay', '7A', ['Xbarc108', 'Xwmc9'], '7.6', '5.1'), ('Florett/Biscay', '7D', ['XwPt-7842', 'XwPt-7368'], '4.3', '3.0'),
                          ('Tuareg/Biscay', '1A', ['Xwmc24', 'XwPt-3904'], '3.1', '3.1'), ('Tuareg/Biscay', '4A', ['XwPt-5434', 'Xwmc219'], '5.8', '5.8'), ('Tuareg/Biscay', '4B', ['Xwmc471', 'Xwmc238'], '19.3', '17.3'),
                          ('Tuareg/Biscay', '4D', ['Xcfd71', 'Xgwm129'], '5.2', '5.2'), ('Tuareg/Biscay', '6B', ['XwPt-6286'], '12.5', '11.8'), ('Tuareg/Biscay', '7B', ['Xwmc517'], '8.2', '8.0')]:
    r('A_STB_011', 'QStb-%s (%s; Risser 2011)' % (c, pop), c, m, pve=pv, stat='LOD ' + ld, note='%s; adult-plant STB, means across environments.' % pop + (' Flanks are Xstb10 and an AFLP; the 3B QTL is the Stb10-linked region, not placed.' if not m else ''),
      ev='Table 3: %s %s %s LOD %s' % (pop, c, '-'.join(m) or 'XP2553-237-Xstb10', ld))
# ---- Sanchez-Vidana et al. 2022 Rev Mex Cienc Agric 13 (Mutus#1 x Huirivis#1), Mb printed in the text
for c, q, s, e, pv, src in [('1B', '1B', 40.8, 563.0, '6.4', 'Mutus#1'), ('4A', '4A', 39.5, 378.5, '8.4', 'Mutus#1'), ('4B', '4B', 96.2, 652.9, '1.7', 'Huirivis#1'),
                            ('5B', '5B.1', 47.4, 604.0, '4.6', 'Mutus#1'), ('5B', '5B.2', 489.8, 583.6, '5.1', 'Huirivis#1')]:
    r('SanchezVidana2022', 'QStb-%s (Sanchez-Vidana 2022)' % q, c, [], fixed=dict(start=s, end=e, grade='A'), pve=pv,
      note='Mutus#1 x Huirivis#1 RIL, adult plants at Toluca; resistance from %s; physical interval printed in the text.' % src, ev='Table 4 and text: %s %.1f-%.1f Mb' % (q, s, e))
# ---- Simon et al. 2010 Eur J Plant Pathol 128:317 (CS x CS(spelt 7D))
r('Simon2010EJPP', 'QStb.ipk-7D1', '7D', ['Xgwm111'], note='CS x CS(T. spelta 7D) single-chromosome recombinant DH; seedling resistance near the 7D centromere at Xgwm111 (Stb4/Stb5 region).',
  ev='Results: QStb.ipk-7D1 linked to Xgwm0111')
r('Simon2010EJPP', 'QStb.ipk-7D2', '7D', ['Xgwm1220'], note='Adult-plant resistance on 7DS near Xgwm1220, in the Lr34/Yr18 region; resistance from CS.',
  ev='Results: QStb.ipk-7D2 in the region of Xgwm1220')
# ---- Stadlmeier et al. 2019 G3 9:1745 (BMWpop MAGIC), Table 3 Mbp
for q, c, p, pv in [('QStb.lfl-1A', '1A', 9, '6.1'), ('QStb.lfl-1B', '1B', 690, '5.9'), ('QStb.lfl-2B.1', '2B', 57, '5.3'), ('QStb.lfl-2B.2', '2B', 440, '4.7'), ('QStb.lfl-2B.3', '2B', 700, '6.3'),
                    ('QStb.lfl-2D', '2D', 10, '8.5'), ('QStb.lfl-4B', '4B', 24, '4.3')]:
    r('Stadlmeier2019G3', q, c, [], fixed=dict(peak=float(p), grade='A'), pve=pv, note='BMWpop 8-founder MAGIC; field STB; peak position (Mbp, rounded) printed in Table 3.', ev='Table 3: %s %s %s Mbp' % (q, c, p))
for q, c, p, pv in [('QTs.lfl-1A', '1A', 9, '6.8'), ('QTs.lfl-2A', '2A', 71, '6.6'), ('QTs.lfl-2B', '2B', 130, '4.6'), ('QTs.lfl-2D', '2D', 640, '6.7'), ('QTs.lfl-3D', '3D', 210, '7.1'),
                    ('QTs.lfl-4B', '4B', 600, '5.2'), ('QTs.lfl-5B', '5B', 41, '5.0'), ('QTs.lfl-7A.1', '7A', 130, '10.6'), ('QTs.lfl-7A.2', '7A', 710, '4.5')]:
    r('Stadlmeier2019G3', q, c, [], fixed=dict(peak=float(p), grade='A'), pve=pv, disease='Tan spot', note='BMWpop 8-founder MAGIC; field tan spot; peak position (Mbp, rounded) printed in Table 3.', ev='Table 3: %s %s %s Mbp' % (q, c, p))
# ---- Zwart et al. 2010 Mol Breed 26:107 (CPI133872 x Janz)
r('X5_1871', 'QStb.wai-3A', '3A', ['wPt-0714', 'wPt-1589'], pve='11.1', note='CPI133872 x Janz DH; resistance from Janz (Stb6 region).', ev='Table: QStb.wai-3A wPt-0714 (4.3) wPt-1589')
r('X5_1871', 'QStb.wai-3D', '3D', [], fixed=dict(start=590.7, end=603.5, grade='C'), pve='21.2-45.1',
  note='CPI133872 x Janz DH; 3DL between AFLP P37M54e and Lr24; interval as projected by He et al. (2021) and Odilbekov et al. (2019) (590.7-603.5 Mb). The atlas had it at Xgwm341 on 3DS.',
  ev='Table: QStb.wai-3D P37M54e (0.4) Lr24 (1.0)')
r('X5_1871', 'QStb.wai-7D', '7D', ['wPt-3727', 'Xwmc405'], pve='10.8', note='CPI133872 x Janz DH; resistance from CPI133872 (Stb5 region).', ev='Table: QStb.wai-7D wPt-3727 (15.0) Xwmc405b')

# ---- Riaz et al. 2020 Genes 11:887 (MAGIC): built from riaz.json (Table 5 Mbp) in stb44.py
