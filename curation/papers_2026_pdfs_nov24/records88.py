# Loci read from the full-text PDFs the owner uploaded on 24 Nov 2026 (Paper_2026.zip).
V1,V2='RefSeq v1.0','RefSeq v2.1'
R=[]
def add(doi,name,cls,dis,ch,s,e,p,ref,mk,pve,stat,g,ev,note): R.append((doi,name,cls,dis,ch,s,e,p,ref,mk,pve,stat,g,ev,note))
# Zhang Z et al. 2026 TAG, Gansu landraces, Table 3 (RefSeq v2.1)
D='10.1007/s00122-026-05323-5'
for n,mk,ch,s,e,p,P,pve,env,typ in [
 ('Yr.hzau.2A','2A_434967820','2A',None,None,439.026532,'0.0009','2.68','CYR34','ASR'),
 ('Yr.hzau.2D','2D_56059640','2D',None,None,58.540354,'0.0009','6.12','CYR34','ASR'),
 ('Yr.hzau.3B.1','3B_81093919','3B',None,None,92.856437,'0.0002','10.91','CYR34','ASR'),
 ('Yr.hzau.3B.2','3B_111554295','3B',None,None,123.301914,'0.0001','9.8','CYR34','ASR'),
 ('Yr.hzau.3D','3D_99439594','3D',None,None,99.983133,'0.0004','2.92','CYR34','ASR'),
 ('Yr.hzau.4B','4B_454694596; 4B_453257700','4B',452.200895,453.860571,None,'0.000015','9.46-11.58','CYR34 (ASR); EZ2025, GS2025, mean (APR)','ASR and APR'),
 ('Yr.hzau.7A.1','7A_43256889','7A',None,None,45.690988,'0.0003','8.22','CYR34','ASR'),
 ('Yr.hzau.4A','4A_709619583','4A',714.36157,714.361645,None,'0.000254','3.15','EZ2023, GS2024, mean','APR'),
 ('Yr.hzau.6A','6A_454397082','6A',None,None,456.815077,'0.000014','3.45','EZ2023, EZ2025, mean','APR'),
 ('Yr.hzau.7A.2','7A_516104014','7A',None,None,520.196895,'0.0000104','6.43','EZ2025, GS2025, mean','APR')]:
    add(D,n,'QTL','Stripe rust',ch,s,e,p,V2,mk,pve,'P '+P,'A' if p else 'B',f'Table 3: {n}, peak {mk}, {ch}, '+(f'{p*1e6:,.0f} bp' if p else f'{s*1e6:,.0f}-{e*1e6:,.0f} bp')+f' (RefSeq v2.1), {typ}, {env}',f'GWAS (MLM and MLMM) in Gansu landraces; {typ} locus.')
# Bian R et al. 2026 TAG, native FHB QTL, Table 2 (RefSeq v1.0)
D='10.1007/s00122-026-05200-1'
for n,ch,mk,p,r2,ev in [('QPss.hwwg-1AS','1A','BS00107852_51',20.980029,'5.1-9.7','2011GH-PSS, 2011FD-PSS'),('QPss.hwwg-1DS','1D','D_F1BEJMU01CYX2V_133',89.59756,'4.5-9.7','GH-PSS BLUP, 2010FD, 2011FD, FD-PSS BLUP'),
 ('QPss.hwwg-2BL','2B','Excalibur_c18989_1026',629.022659,'5.2-9.4','2009GH, 2010GH, 2011GH, GH BLUP'),('QPss.hwwg-3AL','3A','Tdurum_contig63460_144',725.841228,'5.7-7.5','2010GH, 2011GH, GH BLUP'),
 ('QPss.hwwg-3BL','3B','Tdurum_contig9738_170',798.348368,'5.5-6.2','2010GH, 2011GH, GH BLUP'),('QPss.hwwg-5BL','5B','Excalibur_c24511_1282',696.023038,'5.1-8.2','2010GH, 2011GH, GH BLUP')]:
    add(D,n,'MTA','Fusarium head blight',ch,None,None,p,V1,mk,r2,'',"A",f'Table 2A: {n}, {mk}, {ch}, {p*1e6:,.0f} bp; experiments {ev}','Percentage of symptomatic spikelets (type II resistance) in US hard winter wheat; positions on RefSeq v1.0 (stated in the paper).')
for n,ch,s,e,r2 in [('QFDon.hwwg-4AL','4A',622.237212,622.237212,'5.1-9.4'),('QFDon.hwwg-5BL','5B',689.851808,690.331053,'4.9-9.9'),('QFDon.hwwg-5DL','5D',546.086597,547.273861,'4.5-12')]:
    add(D,n,'QTL','Fusarium head blight',ch,s if s!=e else None,e if s!=e else None,s if s==e else None,V1,'',r2,'',"A" if s==e else "B",f'Table 2B: {n}, {ch}, {s*1e6:,.0f}-{e*1e6:,.0f} bp; FDK and DON in 2010FD, 2011FD and BLUP','Low Fusarium-damaged kernels and DON content (RefSeq v1.0).'+(' Candidate TraesCS4A02G342500 (serine/threonine kinase).' if '4AL' in n else ''))
# Wang J et al. 2026 TAG, yellow rust MEQ (RefSeq v2.1)
D='10.1007/s00122-026-05312-8'
for n,ch,s,e,nt in [('Q1A_3.92_10.73','1A',3.92,10.73,''),('Q2A_0.47_39.79','2A',0.47,39.79,' Detected in Exp II and Exp III.'),('Q2B_112.61_207.39','2B',112.61,207.39,''),('Q5D_504.52_508.21','5D',504.52,508.21,''),('Q6A_607.21_617.66','6A',607.21,617.66,'')]:
    add(D,n,'QTL','Stripe rust',ch,s,e,None,V2,'','','',"B",f'Multi-environment QTL {n} named in the Results (chromosome_start_end in Mb on Chinese Spring RefSeq v2.1)','Multi-environment QTL from GWAS of 5,840 hybrids and elite lines (central Europe); the full list of 38 MEQ is in Supplemental Table S6.'+nt)
# Luo Z et al. 2026 BMC Plant Biol, Table 2 (RefSeq v2.1, genes TraesCS..03G)
D='10.1186/s12870-026-10229-x'
for mk,ch,p,P,pve,env in [('AX-94506771','1A',156.07707,'0.000033','7.1-10.4','22YL, 23YL, BLUP'),('AX-108773928','6B',227.641302,'0.000489','7.1-10.3','22YL, 23YL, BLUP'),('AX-109414509','6B',227.071968,'0.000647','6.4-11.0','22YL, 23YL, BLUP'),
 ('AX-109428896','6B',200.387641,'0.000535','11.1-12.9','22YL, 23YL, BLUP'),('AX-110380107; AX-94516403','6B',228.977167,'0.000173','8.7-10.8','22YL, 23YL, BLUP'),('AX-110948069','6B',229.669634,'0.000723','7.9-9.8','22YL, 23YL, BLUP'),
 ('AX-111123557','6B',201.993149,'0.000197','10.6-13.3','22YL, 23YL, BLUP'),('AX-89488392','6B',335.170825,'0.000213','10.1-12.3','22YL, 23YL, BLUP'),('AX-179560800','7B',675.149692,'0.000214','1.3-2.7','22YL, 23YL, 24TS, BLUP'),
 ('AX-111091126','7B',678.124448,'0.000070','3.2-6.0','23YL, 24TS, BLUP'),('AX-94537587','7B',680.853828,'0.000070','5.2-8.6','22YL, 23YL, BLUP')]:
    add(D,f'MTA Yr {mk.split(";")[0]}','MTA','Stripe rust',ch,None,None,p,V2,mk,pve,'BLUP P '+P,'A',f'Table 2: {mk}, {ch}, {p*1e6:,.0f} bp; {env}','Recurrent association with adult-plant stripe rust severity in 120 accessions (50K liquid SNP chip); positions on RefSeq v2.1.'+(' KASP assay developed.' if mk=='AX-111091126' else ''))
# Jan F et al. 2026 Funct Integr Genomics, Table 1 (RefSeq v1.0 genes)
D='10.1007/s10142-025-01799-4'
for mk,ch,p,env,gn in [('AX-92370339','1A',14.30,'E1, E3, E7','NB-ARC genes TraesCS1A02G028900 and others'),('AX-92746236','1D',12.42,'E1, E7','NB-ARC genes TraesCS1D02G029100 and others'),('AX-92397119','3A',714.77,'E1, E7','TraesCS3A02G334200 (kinase), TraesCS3A02G333900 (NB-ARC)'),
 ('AX-92597993','4B',612.03,'E1, E3, E4, E7','TraesCS4B02G322100 (Mlo-related)'),('AX-92334415','1D',318.50,'E7','TraesCS1D02G230100 and others'),('AX-92293824','3B',145.52,'E3, E7','TraesCS3B02G153700 and others')]:
    add(D,f'MTA Yr {mk}','MTA','Stripe rust',ch,None,None,p,V1,mk,'','FDR-corrected','A',f'Table 1: {mk}, {ch}, {p} Mb; environments {env}',f'Stable MTA for stripe rust in western Himalayan germplasm; candidate genes {gn}. Positions with RefSeq v1.0 (v1.1 gene IDs).')
add(D,'MTA Yr AX-92621629','MTA','Stripe rust','3B',None,None,None,'','AX-92621629','','FDR-corrected','D','Table 1: AX-92621629, 3B, "1498.52" Mb; environments E4, E7; candidate TraesCS3B02G157500','Printed position (1,498.52 Mb) exceeds the length of chromosome 3B, so no position was recorded.')
# Thapa S et al. 2026 Plant Disease, BLS, Table 3 (RefSeq v1.0)
D='10.1094/pdis-05-25-1021-re'
for n,mk,ch,p,P,lod in [('q1A.1','scaffold22435_16886199','1A',312.563843,'1.4E-05','4.85'),('q1D.1','scaffold68063_4254674','1D',467.861913,'1.3E-06','5.86'),('q3B.1','scaffold126088_1398217','3B',13.521717,'3.1E-05','4.51'),
 ('q4A.1','scaffold60431_17341658','4A',119.692981,'3.0E-07','6.52'),('q5A.1','scaffold15531_2782724','5A',376.097499,'4.1E-10','9.39'),('q5A.2','scaffold12747_4957938','5A',481.697821,'9.1E-06','5.04')]:
    add(D,f'{n} (BLS)','MTA','Bacterial leaf streak',ch,None,None,p,V1,mk,'',f'P {P}; -log10P {lod}','A',f'Table 3: {n}, {mk}, {ch}, {p*1e6:,.0f} bp (RefSeq v1.0)','GWAS of 355 diverse accessions for bacterial leaf streak (302,524 SNPs).'+(' Validates the 3B QTL of cultivar Boost.' if n=='q3B.1' else '')+(' Four more SNPs within 1 Mb.' if n=='q5A.2' else ''))
# Xiang M et al. 2026 TAG, Yr52 and QTL, Table 3 (RefSeq v2.1)
D='10.1007/s00122-026-05269-8'
add(D,'Yr52 (fine map)','DG','Stripe rust','7B',728.05,730.63,None,V2,'7B_716635041-7B_718449045','13.6-40.2','LOD 6.2-22.7','B','"Yr52 ... physical interval 728.05-730.63 Mb" in CS RefSeq v2.1 (Results); QTL intervals 718.1-741.6 Mb (Table 3)','High-temperature adult-plant resistance gene from PI 660057; 72 genes (29 high-confidence) in the interval.')
for n,ch,s,e,pve,lod,nt in [('QYr.nwafu-6BS.6','6B',110.1,158.2,'8.7-17.7','5.4-13.6',' Probably Yr78; interval spans the detections in all environments.'),('QYr.nwafu-1AL','1A',466.8,471.2,'4.0','3.8',' Resistance allele from KN9204.'),
 ('QYr.nwafu-2AL','2A',644.4,646.9,'5.2-5.6','4.3-5.1',' Resistance allele from PI 660057.'),('QYr.nwafu-2BL','2B',696.4,698.6,'21.6-21.8','6.4-7.0',' Resistance allele from KN9204.'),
 ('QYr.nwafu-3DL','3D',335.0,346.5,'4.9','4.8',' Resistance allele from KN9204.'),('QYr.nwafu-6AL','6A',282.3,526.3,'9.2-11.8','4.5-7.2',' Wide interval (pericentromeric).')]:
    add(D,n,'QTL','Stripe rust',ch,s,e,None,V2,'',pve,'LOD '+lod,'B',f'Table 3: {n}, physical interval {s}-{e} Mb (CS RefSeq v2.1)','PI 660057 x KN9204 RIL.'+nt)
# Zhuo S et al. 2026 BMC Plant Biol, Pm4 allele in YZ16 (RefSeq v2.1)
D='10.1186/s12870-026-08903-1'
add(D,'Pm4 allele (YZ16)','DG','Powdery mildew','2A',766.9,784.6,None,V2,'JS717xJS718','','','B','"BSR-Seq localized the resistance locus to a distal interval on chromosome 2A (766.9-784.6 Mb), and linkage mapping revealed complete co-segregation between resistance and JS717xJS718, a functional marker of the cloned Pm4 gene" (abstract)','Pm4 allele conferring stable powdery mildew resistance in breeding line YZ16.')
# Thielmann J et al. 2026 Curr Plant Biol, Table 4 (powdery mildew colonies), RefSeq v2.1
D='10.1016/j.cpb.2025.100568'
for n,mk,ch,p,lod,r2,tr in [('QTL.Bgt_Control','BS00040600_51','7A',105.8,'3.21','13.87','control'),('QTL1.Bgt_S.indica','IACX5764','1B',435.6,'3.71','9.91','S. indica'),('QTL2.Bgt_S.indica','Kukri_c9133_2211','1B',577.7,'3.04','8.15','S. indica'),
 ('QTL3.Bgt_S.indica','IACX5912','3B',601.9,'3.19','8.53','S. indica'),('QTL4.Bgt_S.indica','Excalibur_rep_c112888_602','4A',720.3,'4.86','13.04','S. indica'),('QTL5.Bgt_S.indica','Kukri_c5033_1815','4A',737.5,'3.06','8.20','S. indica'),
 ('QTL6.Bgt_S.indica','Tdurum_contig46583_2203','4A',747.7,'3.12','8.35','S. indica'),('QTL7.Bgt_S.indica','TA005831-0880','5A',43.2,'3.03','8.10','S. indica'),('QTL8.Bgt_S.indica','GENE-3726_860','6B',605.5,'3.29','8.80','S. indica'),
 ('QTL9.Bgt_S.indica','BS00094895_51','6D',420.4,'3.29','8.80','S. indica'),('QTL10.Bgt_S.indica','tplb0057f21_914','7A',1.3,'5.02','13.51','S. indica'),('QTL11.Bgt_S.indica','Excalibur_c33542_113','7A',19.8,'3.52','9.40','S. indica'),
 ('QTL12.Bgt_S.indica','RAC875_c7988_1588','7A',101.3,'3.48','9.29','S. indica'),('QTL13.Bgt_S.indica','wsnp_JG_c227_167774','7A',105.8,'3.38','9.02','S. indica'),('QTL14.Bgt_S.indica','tplb0032k20_90','7D',16.7,'3.22','8.60','S. indica'),('QTL15.Bgt_S.indica','TA002294-0887','7D',133.8,'3.05','8.17','S. indica')]:
    add(D,n,'MTA','Powdery mildew',ch,None,None,p,V2,mk,r2,'LOD '+lod,'A',f'Table 4: Bgt colonies, {tr}, {n}, {mk}, {ch}, {p} Mb (RefSeq v2.1)','Association with powdery mildew colony number '+('in untreated plants' if tr=='control' else 'after priming with the root endophyte Serendipita indica (priming-induced resistance)')+'; GWAS of a European winter wheat panel.')
