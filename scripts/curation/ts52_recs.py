"""Batch 5 (Nov 2026): tan spot papers from Archive_1.zip, re-read. Same tuple layout as stb44_recs.py.
Tsc1 is placed at the fine-mapped 184 kb interval (Running et al. 2022, 1A 4.60-4.76 Mb); Tsc2 at the cosegregating EST BE444541
(TraesCS2B02G048700, 2B 24.10 Mb; Corsi et al. 2020); Tsn1 at the gene (5B 546.70-546.84 Mb)."""
TS = 'Tan spot'
SNB = 'Septoria nodorum blotch'
R = []
def r(src, name, ch, markers=(), peak=None, fixed=None, pve='', stat='', note='', ev='', **extra):
    extra.setdefault('disease', TS)
    R.append(dict(source=src, name=name, chr=ch, markers=list(markers), peak_marker=peak, fixed=fixed or {}, pve=pve, stat=stat,
                  note=note, evidence=ev, extra=extra))
TSC1 = dict(start=4.60, end=4.76, grade='C')
TSC2 = dict(peak=24.10, grade='C')
TSN1 = ['Tsn1']

# ---- Abeysekara et al. 2010 Plant Genome 3:179 (Tsc2)
r('Abeysekara2010PG', 'Tsc2 (Abeysekara 2010)', '2B', ['XTC339813', 'XBE517745'], cls='DG', fixed=dict(peak=24.10, grade='B'), pve='54',
  note='Salamouni x Katepwa RIL; Ptr ToxB sensitivity gene Tsc2 on 2BS, cosegregating with EST-STS XBE444541 (TraesCS2B02G048700, 24.10 Mb) '
       'in a 3.3 cM interval flanked by Xmag681/XTC339813 (distal) and XBE517745 (proximal); compatible Tsc2-ToxB interaction explained 54% of race 5 disease.',
  ev='Abstract and results: Tsc2 cosegregated with XBE444541; 3.3 cM interval Xmag681/XTC339813-XBE517745; 54% of variation')
# ---- Chu et al. 2008 TAG 117:873 (TA4152-60 x ND495 DH), Table 3
CH8 = 'TA4152-60 x ND495 DH; seedling resistance to races 1, 2, 3 and 5 from TA4152-60; CIM.'
r('A_Tan_spot_024', 'QTs.fcu-2AS', '2A', ['Xgwm515', 'Xfcp526'], pve='14-22', stat='LOD 4.95-8.55', note=CH8, ev='Table 3: QTs.fcu-2AS Xgwm515-Xfcp526')
r('A_Tan_spot_024', 'QTs.fcu-4AL', '4A', ['Xbarc236', 'Xgwm644'], pve='10', stat='LOD 3.56', note=CH8 + ' Race 3 only.', ev='Table 3: QTs.fcu-4AL Xbarc236-Xgwm644')
r('A_Tan_spot_024', 'QTs.fcu-5AL', '5A', ['Xbarc1061', 'Xcfa2185'], pve='9-14', stat='LOD 3.47-4.77', note=CH8, ev='Table 3: QTs.fcu-5AL Xbarc1061-Xcfa2185')
r('A_Tan_spot_024', 'QTs.fcu-5BL.1', '5B', ['Xbarc138', 'Xgwm260'], pve='14-26', stat='LOD 4.88-9.02', note=CH8, ev='Table 3: QTs.fcu-5BL.1 Xbarc138-Xgwm260')
r('A_Tan_spot_024', 'QTs.fcu-5BL.2', '5B', ['Xfcp615', 'Xbarc142'], pve='14-17', stat='LOD 5.49-5.75', note=CH8 + ' Races 1 and 2 only.', ev='Table 3: QTs.fcu-5BL.2 Xfcp615-Xbarc142')
# ---- Chu et al. 2010 Mol Breed 25:327 (Lebsock x PI 94749 durum DH), Table 3
CH10 = 'Lebsock x PI 94749 (T. turgidum subsp. carthlicum) durum DH; seedling resistance to races 1 and 2; CIM.'
r('Chu2010MBts', 'QTs.fcu-3A', '3A', ['Xbarc321', 'Xwmc11'], peak='Xbarc321', pve='8-11', stat='LOD 3.10-4.62', note=CH10 + ' Resistance from Lebsock.', ev='Table 3: Xbarc321-Xwmc11, peak Xbarc321')
r('Chu2010MBts', 'QTs.fcu-3B', '3B', ['Xwmc231', 'Xwmc762'], peak='Xwmc529', pve='5-8', stat='LOD 2.36-3.88', note=CH10 + ' Resistance from Lebsock.', ev='Table 3: Xwmc231-Xwmc762, peak Xwmc529')
r('Chu2010MBts', 'QTs.fcu-5A.1', '5A', ['Xgwm425', 'Xgwm6'], peak='Xbarc360', pve='15-22', stat='LOD 4.05-5.14', note=CH10 + ' Resistance from Lebsock.', ev='Table 3: Xgwm425-Xgwm6.1, peak Xbarc360')
r('Chu2010MBts', 'QTs.fcu-5A.2', '5A', ['Xwmc110', 'Xgwm595'], pve='8-13', stat='LOD 3.32-4.11', note=CH10 + ' Resistance from PI 94749; peak at Xgwm6.2 (second locus, not placed).', ev='Table 3: Xwmc110-Xgwm595, peak Xgwm6.2')
r('Chu2010MBts', 'QTs.fcu-7B', '7B', ['Xwmc276', 'Xbarc182'], peak='Xwmc273', pve='6-8', stat='LOD 3.01-3.13', note=CH10 + ' Resistance from PI 94749.', ev='Table 3: Xwmc276-Xbarc182, peak Xwmc273')
# ---- Cheong et al. 2004 Aust J Agric Res 55:315 (Brookton, Cranbrook)
r('Cheong2004AJAR', 'QTs-5BL (Krichauff/Brookton)', '5B', ['Xgwm499', 'Xgwm408'], pve='36', stat='LOD 7.32',
  note='Krichauff/Brookton DH; seedling yellow spot resistance from Brookton on 5BL near tsn1; peak interval between AFLPs XAGA/CAA161 and XACC/CTA113, flanked by Xgwm499 and Xgwm408.',
  ev='Table 1: XAGA/CAA161-XACC/CTA113 LOD 7.32, 36%; SSRs Xgwm067, Xgwm499, Xgwm408, Xgwm604')
r('Cheong2004AJAR', 'QTs-5BL (Cranbrook/Halberd)', '5B', ['Xgwm408'], pve='60', stat='LOD 29.93',
  note='Cranbrook/Halberd DH; same 5BL locus; peak interval Xcdo400-Xgwm408 (RFLP Xcdo400 not placed).', ev='Table 1: Xcdo400-Xgwm408 LOD 29.93, 60%')
# ---- Corsi et al. 2020 TAG 133:935 (ToxB sensitivity, NIAB Elite MAGIC and association panel), Table 2
r('Corsi2020TAG', 'Tsc2 (QTL, NIAB Elite MAGIC)', '2B', [], fixed=dict(start=24.09, end=25.02, grade='A'), pve='41.0', stat='-log10P 109.19',
  note='NIAB Elite MAGIC; ToxB sensitivity (Tsc2); flanks BS00070050_51 and Kukri_c63748_1453 with RefSeq v1.0 positions printed in Table 2. The association panel gave 23.106-25.027 Mb.',
  ev='Table 2: 2B BS00070050_51 (24.09)-Kukri_c63748_1453 (25.02), PVE 41%')
r('Corsi2020TAG', 'QTL-2A ToxB sensitivity (Corsi 2020)', '2A', [], fixed=dict(start=613.00, end=639.99, grade='A'), pve='1.4', stat='-log10P 5.70',
  note='NIAB Elite MAGIC; minor ToxB sensitivity QTL on 2AL (CIM-cov10 only); flank positions printed in Table 2.',
  ev='Table 2: 2A CAP12_rep_c7918_56 (613.00)-RAC875_c38018_278 (639.99)')
# ---- Dinglasan et al. 2021 TAG 134:2823 (WLA-036 x Banks RIL313), Table 1 and text
DG = 'RIL313 (WLA-036 x Banks); tan spot adult-plant and seedling tests in WA and QLD.'
r('Z_1097', 'Qts.313-1A', '1A', [], fixed=dict(peak=4.662, grade='B'), pve='30.7', stat='LOD 3.38',
  note=DG + ' All-stage resistance from WLA-036; peak DArTseq 3033029 lies 4 kb from TraesCS1A02G008353 (1A 4.66 Mb, Tsc1 region).', ev='Table 1 and text: Qts.313-1A peak 3033029')
r('Z_1097', 'Qts.313-2B', '2B', [], fixed=dict(start=137.81, end=215.93, grade='A'), pve='24.6', stat='LOD 12.62',
  note=DG + ' Adult-plant resistance from WLA-036; physical span printed in the text (137.81-215.93 Mb).', ev='Text: Qts.313-2B spanning 137.81-215.93 Mbp')
r('Z_1097', 'Qts.313-5B', '5B', [], fixed=dict(peak=339.925, grade='B'), pve='27.4', stat='LOD 9.95',
  note=DG + ' Adult-plant resistance from Banks; peak DArTseq 7352232 within TraesCS5B02G187900 (339.92 Mb).', ev='Table 1 and text: Qts.313-5B peak 7352232 at TraesCS5B02G187900')
r('Z_1097', 'Qts.313-5D', '5D', [], fixed=dict(peak=152.749, grade='B'), pve='12.5-16.7', stat='LOD 8.57',
  note=DG + ' Resistance from Banks; peak DArTseq marker 6 kb from TraesCS5D02G117900 (152.75 Mb).', ev='Table 1 and text: Qts.313-5D near TraesCS5D02G117900')
r('Z_1097', 'Qts.313-6A', '6A', [], fixed=dict(peak=609.255, grade='B'), pve='24.9', stat='LOD 3.94',
  note=DG + ' Seedling QTL for the ToxA-deleted strain; peak DArTseq 5972240 3.4 kb from TraesCS6A02G399300 (609.25 Mb).', ev='Text: 5972240 located 3396 bp to TraesCS6A02G399300')
for q, c, m, pv, ld in [('Qts.313-3D.a', '3D', '995750', '9.6', '5.76'), ('Qts.313-3D.b', '3D', '4538913', '14.1', '3.67'),
                        ('Qts.313-4B', '4B', '1252679', '13.6', '6.64'), ('Qts.313-5A', '5A', '1089517', '14.2', '3.84')]:
    r('Z_1097', q, c, [], pve=pv, stat='LOD ' + ld, note=DG + ' Seedling QTL; peak DArTseq %s has no physical position in the paper.' % m, ev='Table 1: %s peak %s' % (q, m))
# ---- Effertz et al. 2002 Phytopathology 92:527 (Tsc1, ITMI)
r('Effertz2002Phyto', 'Tsc1 (Effertz 2002)', '1A', [], cls='DG', fixed=TSC1, pve='49.4',
  note='W-7984 x Opata 85 (ITMI) RIL; Ptr ToxC insensitivity locus on 1AS, 5.7 cM from RFLP XGli1 (R2 0.49). Placed at the fine-mapped Tsc1 interval (Running et al. 2022).',
  ev='Abstract and Table 2: XGli1 1AS P<0.0001 R2 0.494')
r('Effertz2002Phyto', 'ToxC culture filtrate QTL 3A (Effertz 2002)', '3A', [], pve='6.3', note='Minor locus at RFLP Xcdo395 (3A), allele from Opata 85; RFLP not placed.',
  ev='Table 2: Xcdo395 3A R2 0.063')
# ---- Faris and Friesen 2005 TAG 111:386 (BR34 x Grandin RIL), Tables 2 and 3
FF = 'Grandin x BR34 RIL; race-nonspecific resistance from BR34 to races 1, 2, 3 and 5.'
r('A_Tan_spot_022', 'QTs.fcu-1BS', '1B', ['Xgdm33', 'Xgdm125'], pve='13-29', stat='LOD 5.47-8.61', note=FF, ev='Table 3: Xgdm33-Xgdm125')
r('A_Tan_spot_022', 'QTs.fcu-3BL', '3B', ['Xbarc248', 'Xfcp83'], fixed=dict(start=466.2649, end=470.047, grade='C'), pve='13-41', stat='LOD 5.03-13.51', note=FF + ' Same locus as Tsr7: Faris et al. (2020) re-analysed this population and found the QTL peak at the seven Tsr7 SNPs, so it is placed at Tsr7 (3B 466.26-470.05 Mb). Flank Xbarc248 (3B 162.47 Mb) disagrees and is not used.', ev='Table 3: Xbarc248-Xfcp83; Faris et al. 2020 results')
r('A_Tan_spot_022', 'QTs.fcu-3BS', '3B', ['Xfcp311', 'Xfcp114'], pve='12', stat='LOD 4.60', note=FF + ' Race 3 only.', ev='Table 3: Xfcp311-Xfcp114')
for c, m, pv, rc in [('2D', 'Xgwm614', '13', 'race 2'), ('4A', 'Xbarc190', '7', 'race 5'), ('5A', 'Xfcp231', '7-11', 'races 1 and 5'), ('7D', 'Xbarc126', '9', 'race 3')]:
    r('A_Tan_spot_022', 'QTs-%s (Faris 2005, minor)' % c, c, [m], pve=pv, note=FF + ' Minor locus found by single-marker regression only (%s).' % rc, ev='Table 2: %s %s' % (m, c))
# ---- Faris et al. 2012 Mol Breed 30:1669 (Salamouni x Katepwa), Table 2
SK = 'Salamouni x Katepwa RIL; seedling tan spot; resistance alleles from Salamouni.'
r('Z_1110', 'QTs.fcu-5B (Tsn1)', '5B', TSN1, pve='25-32', stat='LOD 10.0-14.5', note=SK + ' QTL peak at Tsn1 (races 1 and 2).', ev='Table 2: QTs.fcu-5B Tsn1')
r('Z_1110', 'QTs.fcu-5D', '5D', ['Xwmc215'], pve='13', stat='LOD 4.68', note=SK + ' Flank Xcfd183 amplifies only at 5D 5.6 Mb (5DS) in silico, which contradicts the 5DL map position, so the locus is set at Xwmc215.', ev='Table 2: QTs.fcu-5D Xwmc215-Xcfd183')
r('Z_1110', 'QTs.fcu-7B', '7B', ['Xgwm537', 'Xbarc126'], pve='5-8', stat='LOD 3.15-4.71', note=SK + ' Flank Xbarc126.1 is the 7B locus of barc126.', ev='Table 2: QTs.fcu-7B Xgwm537-Xbarc126.1')
r('Z_1110', 'QTs.fcu-7D', '7D', ['Xcfd21', 'Xgwm44'], pve='7', stat='LOD 4.54', note=SK, ev='Table 2: QTs.fcu-7D Xcfd21-Xgwm44.2')
# ---- Faris et al. 2020 TAG 133:829 (Tsr7)
r('Faris2020TAG', 'Tsr7', '3B', ['IWB7782', 'IWB370', 'IWB40103', 'IWB80300', 'IWB45103', 'IWB80906', 'IWB80299', 'Xwmc69'], cls='DG', pve='17.5-55.1',
  note='Langdon-wild emmer (IsraelA) 3B substitution recombinants (LD3B); single dominant race-nonspecific gene Tsr7 on 3BL, cosegregating with wmc69 and seven SNPs '
       '(IWB7782, IWB370, IWB40103, IWB80300, IWB45103, IWB80906, IWB80299; STARP markers fcp735-fcp738); same locus as the BR34 and Penawawa QTL.',
  ev='Results: seven cosegregating SNPs and wmc69 define Tsr7')
# ---- Kalia et al. 2018 Plant Breed 137:313 (TA4161-L1 x TAM105 F2:3), Table 2
KA = 'TA4161-L1 (synthetic) x TAM105 F2:3; race 1; resistance from TA4161-L1.'
r('V_Kalia2018', 'QTs.ksu-1A', '1A', ['wPt-730618', 'wPt-669484'], pve='22', stat='LOD 11.2', note=KA, ev='Table 2: 1A wPt-730618-wPt-669484')
r('V_Kalia2018', 'QTs.ksu-6A', '6A', ['wPt-4047', 'wPt-664733'], pve='20', stat='LOD 10.7', note=KA, ev='Table 2: 6A wPt-4047-wPt-664733')
r('V_Kalia2018', 'QTs.ksu-7A.1', '7A', ['wPt-0827', 'wPt-3403'], pve='12', stat='LOD 7.4', note=KA, ev='Table 2: 7A wPt-0827-wPt-3403')
r('V_Kalia2018', 'QTs.ksu-7A.2', '7A', ['wPt-671669', 'wPt-740675'], pve='14', stat='LOD 6.4', note=KA, ev='Table 2: 7A wPt-671669-wPt-740675')
# ---- Kariyawasam et al. 2016 TAG 129:897 (Louise x Penawawa RIL), Table 2
LP = 'Louise x Penawawa RIL; races 1, 2, 3, 5 and AR CrossB10; resistance from Penawawa.'
r('Z_1198', 'QTs.zhl-1A', '1A', ['IWA6644', 'Xpsp2999'], pve='9-22', stat='LOD 14.4-23.1', note=LP + ' Ptr ToxC-Tsc1 interaction.', ev='Table 2: Xiwa6644-Xpsp2999')
r('Z_1198', 'QTs.zhl-2D', '2D', ['Xwmc41', 'Xgwm608'], pve='3-9', stat='LOD 4.7-8.6', note=LP, ev='Table 2: Xwmc41-Xgwm608')
r('Z_1198', 'QTs.zhl-3B', '3B', ['IWA1383', 'IWA4613'], pve='22-53', stat='LOD 13.6-44.0', note=LP + ' Race-nonspecific; same locus as Tsr7 (Faris et al. 2020).', ev='Table 2: Xiwa1383-Xiwa4613')
r('Z_1198', 'QTs.zhl-5A', '5A', ['IWA7025', 'IWA5173'], pve='6-14', stat='LOD 5.2-18.2', note=LP, ev='Table 2: Xiwa7025-Xiwa5173')
# ---- Kariyawasam et al. 2018 Mol Breed 38:140 (Harry x Wesley RIL)
HW = 'Harry x Wesley RIL (Nebraska winter wheat); flanking markers are GBS tags without published positions.'
r('V_Kariyawasam2018', 'QTs.zhl-1A (Tsc1)', '1A', [], fixed=TSC1, pve='10-64', stat='LOD 9.9-46.9', note=HW + ' Ptr ToxC-Tsc1 interaction (Pti2, 331-9, AR CrossB10); placed at Tsc1.', ev='Text and Table 2: QTs.zhl-1A HWGBS60-HWGBS150, Tsc1')
r('V_Kariyawasam2018', 'QTs.zhl-5B (Tsn1)', '5B', TSN1, note=HW + ' Ptr ToxA-Tsn1 interaction (Pti2, 86-124); placed at Tsn1.', ev='Text: QTs.zhl-5B mapped at the Tsn1 locus')
for q, c, f, pv, x in [('QTs.zhl-7A', '7A', 'HWGBS5420-HWGBS5422', '12', '86-124 delta ToxA'), ('QTs.zhl-7B', '7B', 'HWGBS5696-HWGBS5992', '8', 'DW5, LOD 6.6'),
                       ('QTs.zhl-7D', '7D', 'HWGBS6047-HWGBS6066', '13', '86-124 delta ToxA')]:
    r('V_Kariyawasam2018', q, c, [], pve=pv, note=HW + ' Minor QTL (%s).' % x, ev='Table 2: %s %s' % (q, f))
# ---- Li et al. 2011 TAG 122:395 (Batavia x Ernie DH), Table 3
LE = 'Batavia x Ernie DH; three seedling trials (TS-1 to TS-3); closest DArT marker.'
for c, m, pv, ld, s in [('1A', 'wPt-5776', '15.8-18.4', '2.72-3.02', 'Batavia'), ('2B', 'wPt-0289', '29.8-38.2', '6.11-7.36', 'Ernie'),
                        ('3B', 'wPt-0571', '25.3', '5.37', 'Ernie'), ('3D', 'wPt-7705', '23.2', '3.74', 'Batavia'), ('5B', 'wPt-3076', '11.9-24.4', '2.94-3.92', 'Ernie'),
                        ('6A', 'rPt-9065', '13.2-16.3', '2.74', 'Batavia'), ('7A', 'wPt-6447', '16.9', '3.63', 'Ernie'), ('7D', 'wPt-6769', '21.1', '', 'Batavia')]:
    r('Li2011TAG', 'QTs-%s (Batavia/Ernie)' % c, c, [m] + (['Xbarc35'] if c == '2B' else []), pve=pv, stat=('LOD ' + ld) if ld else '',
      note=LE + ' Resistance allele from %s.%s' % (s, ' Major QTL on 2BS, Xbarc35 2.7 cM from the peak.' if c == '2B' else ''), ev='Table 3: %s %s' % (c, m))
# ---- Liu et al. 2006 Genome 49:1265 (Grandin x BR34; SNB and Tsn1-ToxA), Table 2
BG = 'Grandin x BR34 RIL; SNB seedling resistance (isolate Sn2000).'
r('Liu2006Genome', 'QSnb.fcu-5BL.1 (Tsn1)', '5B', TSN1, pve='62', stat='LOD 24.62', note=BG + ' Major QTL at the Tsn1 locus (SnToxA-Tsn1).', ev='Table 2: QSnb.fcu-5BL.1 tsn1', disease=SNB)
r('Liu2006Genome', 'QSnb.fcu-5BL.2', '5B', ['Xbarc234', 'Xfcp273'], pve='6', stat='LOD 4.25', note=BG + ' 63 cM distal to Tsn1.', ev='Table 2: Xbarc234-Xfcp273', disease=SNB)
r('Liu2006Genome', 'QSnb.fcu-1BS', '1B', ['Xgdm33', 'Xgdm125'], pve='10', stat='LOD 7.92', note=BG, ev='Table 2: Xgdm33-Xgdm125', disease=SNB)
# ---- Liu et al. 2017 TAG 130:1267 (LMPG-6 x PI 626573 RIL), Table 2
LM = 'LMPG-6 x PI 626573 RIL; isolates Asc1, Pti2, 86-124, 331-9 and DW5; closest marker.'
r('V_Liu2017', 'QTs.zhl-1A (Tsc1)', '1A', ['IWA1367'], fixed=TSC1, pve='23-27', stat='LOD 14-16', note=LM + ' Ptr ToxC-Tsc1; resistance from PI 626573. IWA1367 has no RefSeq position, so the locus is placed at Tsc1.', ev='Table 2: QTs.zhl-1A IWA1367')
r('V_Liu2017', 'QTs.zhl-2A', '2A', ['IWA6514', 'IWA5893'], pve='9-10', stat='LOD 4.7-5.5', note=LM + ' Resistance from LMPG-6.', ev='Table 2: QTs.zhl-2A IWA6514/IWA5893')
r('V_Liu2017', 'QTs.zhl-2B', '2B', ['IWA3594'], pve='7', stat='LOD 3.8', note=LM + ' Asc1 only.', ev='Table 2: QTs.zhl-2B IWA3594')
r('V_Liu2017', 'QTs.zhl-5B.1 (Tsn1)', '5B', TSN1, pve='16-36', stat='LOD 8.9-23.0', note=LM + ' Ptr ToxA-Tsn1.', ev='Table 2: QTs.zhl-5B.1 Tsn1')
r('V_Liu2017', 'QTs.zhl-5B.2', '5B', ['IWA6211'], pve='9', stat='LOD 4.8', note=LM + ' DW5 only; resistance from PI 626573.', ev='Table 2: QTs.zhl-5B.2 IWA6211')
# ---- Liu et al. 2020 TAG 133:433 (Rusty-derived durum RIL populations, race 2), Table 4 (RefSeq v1.0 Mbp)
for pop, q, c, s, e, ld, pv in [('RP336', 'RP336_86-124_2B', '2B', 2.8, 30.5, 4.3, 6.9), ('RP336', 'RP336_86-124_5A.1', '5A', 46.6, 475.6, 4.9, 7.9),
        ('RP336', 'RP336_86-124_5A.2', '5A', 665.9, 679.6, 8.2, 13.8), ('RP336', 'RP336_86-124deltaToxA_1B', '1B', 445.6, 598.7, 4.2, 5.2),
        ('RP336', 'RP336_86-124deltaToxA_3A', '3A', 7.6, 20.2, 6.5, 8.3), ('RP336', 'RP336_86-124deltaToxA_5A', '5A', 665.9, 681.5, 23.1, 36.4),
        ('RP336', 'RP336_86-124deltaToxA_7A', '7A', 1.9, 12.2, 3.6, 4.4), ('RP696', 'RP696_86-124_2B', '2B', 58.1, 118.6, 4.0, 7.1),
        ('RP696', 'RP696_86-124_3A', '3A', 7.4, 17.4, 7.3, 13.6), ('RP696', 'RP696_86-124_5A.1', '5A', 45.0, 452.9, 3.9, 7.0),
        ('RP696', 'RP696_86-124_5A.2', '5A', 663.4, 688.3, 5.2, 9.6), ('RIum', 'RIum_86-124_3A', '3A', 10.8, 26.2, 4.2, 9.8)]:
    r('Liu2020TAGts', q, c, [], fixed=dict(start=s, end=e, grade='A'), pve=str(pv), stat='LOD %s' % ld,
      note='%s durum RIL population (Rusty-derived), race 2 isolate 86-124 or its ToxA knockout; 1.5-LOD physical interval printed in Table 4.' % pop,
      ev='Table 4: %s %s %s-%s Mbp' % (q, c, s, e))
# ---- Shankar et al. 2017 TAG 130:2637, Tables 2 (Calingiri/Wyalkatchem) and 3 (IGW2574/Annuello)
CW = 'Calingiri/Wyalkatchem DH; field and controlled tan spot severity.'
IA = 'IGW2574/Annuello DH; field and controlled tan spot severity.'
SH = [
 (CW, '1A', 'C/W', ['wPt-7296', 'wPt-6280'], '3.6', '1.7', 'Calingiri', 'wPt-7296-wPt-6280'),
 (CW, '2A', 'C/W', ['wPt-9320', 'IWA5893'], '10.3-33.5', '2.7-11.4', 'Wyalkatchem',
  'five intervals between wPt-9320 (22.7 cM) and wri64 (49.5 cM); most often near wri75 and wri79'),
 (CW, '2B', 'C/W, distal', ['wPt-7004', 'wPt-1394'], '3.7', '1.7', 'Calingiri', 'wPt-7004-wPt-1394'),
 (CW, '2B', 'C/W, proximal', ['wPt-6471', 'wPt-8548'], '14.2', '2.9', 'Calingiri', 'wPt-6471-wPt-8548'),
 (CW, '2D', 'C/W', ['wPt-4242', 'wPt-4413'], '3.8', '1.6', 'Wyalkatchem', 'wPt-4242-wPt-4413'),
 (CW, '3A', 'C/W', ['tPt-1143', 'wPt-3839'], '5.6', '2.1', 'Wyalkatchem', 'tPt-1143-wPt-3839'),
 (CW, '3D', 'C/W', ['wPt-10018', 'wPt-733412'], '6.4', '2.4', 'Wyalkatchem', 'wPt-10018-wPt-733412'),
 (CW, '6B', 'C/W', ['wPt-6208', 'wPt-729979', 'tPt-1723', 'wPt-9195', 'wPt-0446', 'wPt-8183'], '4.2-10', '1.9-3.8', 'Wyalkatchem', 'three intervals wPt-6208 (83.4 cM)-wPt-8183 (99.2 cM)'),
 (CW, '7B', 'C/W', ['wPt-7046', 'wPt-744987'], '9.9', '3.5', 'Calingiri', 'wPt-7046-wPt-744987'),
 (IA, '1A', 'I/A, distal', ['wPt-2527', 'wPt-744902'], '25.4', '8.9', 'Annuello', 'wPt-2527-wPt-744902'),
 (IA, '1A', 'I/A, proximal', ['wPt-1418', 'wPt-667155'], '16.6', '4.7', 'Annuello', 'wPt-1418-wPt-667155'),
 (IA, '1B', 'I/A', ['wPt-5347', 'wPt-6592'], '9.0', '2.2', 'IGW2574', 'wPt-5347-wPt-6592'),
 (IA, '2B', 'I/A', ['wPt-733929', 'tPt-9767'], '9.6', '3.5', 'IGW2574', 'wPt-733929-tPt-9767'),
 (IA, '2D', 'I/A, 21-28 cM', ['wPt-666656', 'wPt-3728'], '6.8', '1.7', 'Annuello', 'wPt-666656-wPt-3728'),
 (IA, '2D', 'I/A, 39-44 cM', ['wPt-668239', 'wPt-3692'], '19.7', '7.3', 'Annuello', 'wPt-668239-wPt-3692'),
 (IA, '4A', 'I/A', ['wPt-7612', 'wPt-9445', 'wPt-7524'], '6.0-12.9', '1.6-4.1', 'IGW2574', 'wPt-7612-wPt-9445-wPt-7524'),
 (IA, '4B', 'I/A', ['wPt-5559', 'wPt-731635'], '8.0', '3.0', 'Annuello', 'wPt-5559-wPt-731635'),
 (IA, '5B', 'I/A', ['wPt-3049', 'wPt-5896'], '21.0', '6.5', 'IGW2574', 'wPt-3049-wPt-5896'),
 (IA, '7A', 'I/A, distal', ['wPt-2266', 'wPt-3393'], '6.8', '1.4', 'Annuello', 'wPt-2266-wPt-3393'),
 (IA, '7A', 'I/A, proximal', ['wPt-744897', 'wPt-7122'], '13.7', '2.8', 'Annuello', 'wPt-744897-wPt-7122'),
 (IA, '7B', 'I/A, 126-137 cM', ['wPt-0635', 'wPt-2356'], '12.2', '1.6', 'IGW2574', 'wPt-0635-wPt-2356'),
 (IA, '7B', 'I/A, 155-160 cM', ['wPt-3190', 'wPt-1957', 'wPt-7720', 'wPt-8561'], '10.3-18.6', '1.5-7.2', 'IGW2574', 'wPt-3190-wPt-8561'),
 (IA, '7D', 'I/A', ['wPt-744843', 'wPt-744219'], '9.4', '2.2', 'IGW2574', 'wPt-744843-wPt-744219'),
]
XN = {'2A': ' The other interval markers place inconsistently (Xcfa2263 426 Mb, tPt-8937 630 Mb; wri markers not placed); the locus is set by wPt-9320 and IWA5893 on 2AS.',
      '4A': ' wPt-9445 and wPt-7524 are 2.7 cM apart but 291 Mb apart physically (low-recombination region), so the interval is wide.'}
for pop, c, tag, mk, pv, ld, s, ev in SH:
    r('V_Shankar2017', 'QTs-%s (%s; Shankar 2017)' % (c, tag), c, mk, pve=pv, stat='LOD ' + ld, note=pop + ' Resistance allele from %s.' % s + (XN.get(c, '') if tag in ('C/W', 'I/A') else ''),
      ev=('Table 2: ' if pop is CW else 'Table 3: ') + ev)
# ---- Singh et al. 2008 Mol Breed 21:293 (T. turgidum #283 x Coulter RIL)
r('Singh2008MB', 'tsn5 (Tsr5)', '3B', ['Xgwm285'], cls='DG', note='T. turgidum #283 (PI 352519) x Coulter RIL; recessive race 5 necrosis resistance gene tsn5 on 3BL, 8.3 cM proximal to tsn2; Xgwm285 11 cM. Xgwm285 lies at 3B 273 Mb (RefSeq v1.0, proximal 3BS), so the gene may lie on the long-arm side of this position.',
  ev='Abstract and discussion: Xgwm285 2.7 cM from tsn2 and 11 cM from tsn5')
r('Singh2008MB', 'tsn2 (Tsr2; Singh 2008)', '3B', ['Xgwm285'], cls='DG', note='Same population; race 3 necrosis resistance gene tsn2 on 3BL, Xgwm285 2.7 cM (Xgwm285 at 3B 273 Mb, RefSeq v1.0).',
  ev='Discussion: Xgwm285 2.7 cM from tsn2')
# ---- Sun et al. 2010 Phytopathology 100:468 (Wangshuibai x Ning7840 RIL)
r('Sun2010Phyto', 'QTs.ksu-1AS', '1A', ['Xcfa2153', 'Xgwm33'], pve='39', note='Wangshuibai x Ning7840 RIL; race 1; resistance from Wangshuibai.', ev='Abstract and results: QTs.ksu-1AS Xcfa2153-Xgwm33, 39%')
r('Sun2010Phyto', 'QTs.ksu-2BS', '2B', ['Xbarc7', 'Xbarc55'], pve='10-13', note='Wangshuibai x Ning7840 RIL; race 1; resistance from Wangshuibai.', ev='Results: QTs.ksu-2BS Xbarc7-Xbarc55, 13% (IM) and 10% (CIM)')
# ---- Zou et al. 2017 Mol Breed 37:144 (Attila x CDC Go RIL), Table 3
AC = 'Attila x CDC Go RIL; field tan spot; resistance from CDC Go.'
r('V_Zou2017', 'QTs.dms-2B', '2B', ['GENE-1343_556', 'Kukri_c148_1346'], pve='7.0', stat='LOD 3.0', note=AC, ev='Table 3: QTs.dms-2B GENE-1343_556-Kukri_c148_1346')
r('V_Zou2017', 'QTs.dms-2D', '2D', ['Ppd-D1', 'GENE-0787_85'], pve='10.0', stat='LOD 3.6', note=AC, ev='Table 3: QTs.dms-2D Ppd-D1-GENE-0787_85')
r('V_Zou2017', 'QTs.dms-6B', '6B', ['wsnp_Ex_c19082_27999258', 'Excalibur_c35713_106'], pve='6.9', stat='LOD 3.2', note=AC, ev='Table 3: QTs.dms-6B wsnp_Ex_c19082_27999258-Excalibur_c35713_106')
# ---- Zwart et al. 2010 Mol Breed 26:107 (CPI133872 x Janz DH), yellow leaf spot rows of the QTL table
ZW = 'CPI133872 x Janz DH; yellow leaf spot (tan spot) field trials 2001, 2004, 2007; resistance from CPI133872.'
for q, c, mk, pv, x in [('QYls.lrc-2B', '2B', ['wPt-5027', 'wPt-6148', 'Xwmc25', 'Xcfd238'], '14.7', ' Tsc2 region.'),
                        ('QYls.lrc-2D', '2D', ['Xbarc168', 'Xgwm484'], '8.6', ''), ('QYls.lrc-3A', '3A', ['Xwmc322', 'wPt-2698'], '6.5', ''),
                        ('QYls.lrc-3D', '3D', ['Xwmc169', 'Xgwm114'], '5.4-12.9', ''), ('QYls.lrc-4B', '4B', ['Xgwm66', 'Xwmc47'], '12.0', ''),
                        ('QYls.lrc-5A', '5A', ['wPt-1370', 'Xgdm132', 'Xgwm443'], '10.7-12.9', ' The paper relates it to QTs.fcu-5AL, but the only placed marker, Xgwm443, amplifies at 5A 11.3 Mb (5AS); the position is uncertain.'),
                        ('QYls.lrc-5B', '5B', ['Xgwm499', 'Xfcp2', 'Xgwm408'], '28.9-40.0', ' Tsn1 region.')]:
    r('X5_1871', q, c, mk, pve=pv, note=ZW + x, ev='QTL table: %s %s' % (q, ', '.join(mk)))
