import json,copy,re
d=json.load(open('../all_now.json'));A={r['id']:r for r in d['reports']}
L={'IntegratedMQTLmap':'integrated wheat consensus map (Supplementary Tables 2-8 supplied by the owner)','Wen2017':'Wen et al. 2017 consensus map (Front Plant Sci 8:1389)',
   'Maccaferri2014':'Maccaferri et al. 2014 durum consensus map (BMC Genomics)','Marone2012':'Marone et al. 2012 durum consensus map (Theor Appl Genet 125:1619)','KOMUGIcomposite':'KOMUGI composite wheat map'}
def lab(s):
    for k,v in L.items(): s=s.replace('interpolated on the '+k,'interpolated on the '+v)
    return s
def setpos(r,lo,hi):
    if hi-lo<0.05: r['start']=r['end']=None; r['peak']=round((lo+hi)/2,4)
    else: r['start'],r['end'],r['peak']=round(lo,4),round(hi,4),None
    r['ref']='RefSeq v1.0'; r['grade']='C'; r['curation']='auto (in-silico from markers, consensus maps Nov 2026)'
    a={k:v for k,v in (r.get('asm') or {}).items() if isinstance(v,dict) and (v.get('start') is not None or v.get('peak') is not None)}
    r['asm']={}   # reconvert
CHK=' Interpolated positions come from published consensus genetic maps anchored by markers of known physical position; check by hand.'
out={}
# 1. general consensus-map placements of the remaining records
for o in json.load(open('proposals2.json')):
    if o['status']!='placed' or o['id'] in ('R0660','R0661','R0662'): continue
    r=copy.deepcopy(A[o['id']]); setpos(r,o['lo'],o['hi'])
    how='; '.join(f"{m} {a:.2f}"+(f"-{b:.2f}" if b-a>=0.05 else '')+f" Mb ({s})" for m,a,b,s in o['P'])
    r['note']=re.sub(r'\s{2,}',' ',((r.get('note') or '')+f" Placed in silico on RefSeq v1.0 on 25 Nov 2026 from the markers named by the source: {lab(how)}."+CHK).strip())
    out[r['id']]=r
# 2. Aguilar et al. 2005 (Forno x Oberkulmer), Table 5
AG='Aguilar et al. 2005 Table 5 (Forno x Oberkulmer RILs, CIM, LOD > 3.7)'
def ag(i,lo=None,hi=None,mk='',lod=None,pve=None,how='',extra=''):
    r=copy.deepcopy(A[i])
    if lo is not None: setpos(r,lo,hi)
    if mk: r['marker']=mk
    if lod: r['stat']=f'LOD {lod}'
    if pve: r['pve']=str(pve)
    r['note']=re.sub(r'\s{2,}',' ',((r.get('note') or '')+' '+(f'{AG}: flanking markers {mk}, LOD {lod}, R2 {pve}%.' if mk else '')+(' '+how if how else '')+extra).strip())
    out[i]=r
ag('V0075',507.85,527.27,'Xpsr932, Xpsr331a',8.6,20.8,lab('Placed in silico on RefSeq v1.0 on 25 Nov 2026: Xpsr331a at 167.2 cM on the IntegratedMQTLmap between anchors at 507.85 and 527.27 Mb (relaxed rule, anchors within 3 cM; interval wider than 10 Mb); Xpsr932 has no position.')+CHK)
ag('V0081',42.56,100.97,'Xglk538, Xpsr902',5.63,9.9,lab('Placed in silico on RefSeq v1.0 on 25 Nov 2026: Xglk538 42.56-43.74 Mb and Xpsr902 95.16-100.97 Mb (interpolated on the IntegratedMQTLmap, anchors within 3 cM); the flanking interval is wide because the support interval spans the distal 3BS.')+CHK)
ag('V0080',604.07,605.48,'Xglk163, Xpsr426',5.82,1.6,lab('Placed in silico on RefSeq v1.0 on 25 Nov 2026: Xpsr426 604.07-605.48 Mb (interpolated on the IntegratedMQTLmap, anchors within 1 cM; the KOMUGI composite map gives 577-598 Mb); Xglk163 has no position.')+CHK,' Atlas name QSnl.eth-5B is QSnl.eth-5B2 of the paper.')
ag('K0149',654.09,663.94,'Xpsr1194, Xpsr918b',27.5,35.8,lab('Placed in silico on RefSeq v1.0 on 25 Nov 2026: Xpsr1194 654.09-657.06 Mb and Xpsr918b 659.73-663.94 Mb (interpolated on the IntegratedMQTLmap); the q locus lies just proximal.')+CHK,' Same QTL as QSng.eth-5A2.')
ag('K0152',700.63,711.23,'Xmwg710a, Xglk576',6.72,12.8,lab('Placed in silico on RefSeq v1.0 on 25 Nov 2026: Xmwg710a 700.63 Mb and Xglk576 700.63-711.23 Mb (interpolated on the KOMUGIcomposite map, anchors within 3 cM).')+CHK,' Same QTL as QSnl.eth-7B3.')
ag('V0076',mk='Xglk348, Xpsr921',lod=14.3,pve=17,how='Neither flanking RFLP has a known physical position or a GenBank probe sequence; arm 4BS distal (4 cM).')
ag('K0151',mk='Xglk348, Xpsr921',lod=14.3,pve=17,how='Same QTL as V0076 (QSnl.eth-4B); flanking RFLPs have no known physical position.')
ag('V0082',mk='Xpsr131, Xglk687',lod=3.96,pve=1.8,how='Flanking RFLPs have no usable physical position. Atlas name QSnl.eth-2A is QSnl.eth-2A2 of the paper.')
ag('V0077',mk='Xglk301, Xglk317',pve=7.4,how='Flanking RFLPs have no known physical position.')
ag('V0078',mk='Xpsr961, Xglk687',pve=3.8,how='Flanking RFLPs have no known physical position.')
ag('V0079',mk='Xglk165, pwir232a',pve=5.2,how='Flanking markers have no known physical position.')
ag('V0083',mk='Xgwm18, Xglk483',lod=3.79,pve=7.3,how='Atlas name QSnl.eth-1BS is QSnl.eth-1BS2 of the paper.')
ag('Z0255',mk='Xpsr952, Xgwm46',lod=4.55,pve=0.9)
base={'cls':'QTL','disease':'Septoria nodorum blotch','source':'Z_1048','pve':'','stat':'','grade':'D','ref':'RefSeq v1.0','start':None,'end':None,'peak':None,'asm':{},'curation':'auto (paper full text)'}
NEW=[('QSng.eth-3A1','3A','Xpsr304, Xpsr598',3.95,0.3,'50 (36-52)',(56.68,57.42,'Xpsr598 56.68-57.42 Mb (interpolated on the IntegratedMQTLmap, anchors within 3 cM); Xpsr304 has no position on 3A')),
     ('QSng.eth-4A1','4A','CD16.2, Xglk331',5.18,1.4,'48 (38-52)',(617.95,617.95,'Xglk331 617.95 Mb (interpolated on the IntegratedMQTLmap; KOMUGI composite map 601-617 Mb)')),
     ('QSng.eth-4A2','4A','Xglk128, Xgwm710b',5.51,4.3,'182 (178-188)',(625.86,627.82,'Xglk128 625.86-627.82 Mb (interpolated on the IntegratedMQTLmap, anchors within 1 cM)')),
     ('QSnl.eth-2B1','2B','Xglk407, Xpsr924',3.87,8.1,'86 (82-94)',None),('QSnl.eth-2B2','2B','Xpsr644b, Xpsr956a',6.59,4.7,'164 (158-168)',None),
     ('QSng.eth-3A2','3A','Xglk118, Xglk577',3.79,6.4,'94 (82-96)',None),('QSng.eth-4DL2','4D','Xglk302, Xpsr1101a',4.11,1.3,'98 (90-100)',None),
     ('QSng.eth-5B1','5B','Xpsr128, Xglk614',9.72,6.3,'22 (16-32)',None)]
n=0
for name,ch,mk,lod,pve,ci,pl in NEW:
    n+=1;r=dict(base,id=f'PF6-{n:03d}',name=name,chr=ch,marker=mk,stat=f'LOD {lod}',pve=str(pve),evidence=f'Table 5: {ch} {name} support interval {ci} cM, flanking {mk}, LOD {lod}, R2 {pve}%')
    tr='glume blotch (SNG)' if 'Sng' in name else 'leaf blotch (SNL)'
    r['note']=f'Forno x Oberkulmer RILs (204), field resistance to S. nodorum {tr}, CIM across 2 years.'
    if pl: setpos(r,pl[0],pl[1]); r['note']+=lab(f' Placed in silico on RefSeq v1.0 on 25 Nov 2026: {pl[2]}.')+CHK
    elif name=='QSng.eth-5B1': r['note']+=' Xpsr128 gives conflicting positions on two consensus maps (435-438 Mb vs 304-331 Mb); left unplaced.'
    else: r['note']+=' Flanking RFLPs have no known physical position; left unplaced.'
    r['note']+=' Read from the full-text PDF supplied by the owner on 25 Nov 2026; not yet checked by hand.'
    out[r['id']]=r
# 3. Adhikari et al. 2012 (BLS and SB association mapping, DArT), Table 2
AD='Adhikari et al. 2012 Table 2'
def adset(i,lo=None,hi=None,dis=None,stat=None,pve=None,how=''):
    r=copy.deepcopy(A[i])
    if lo is not None: setpos(r,lo,hi)
    if dis: r['disease']=dis
    if stat: r['stat']=stat
    if pve: r['pve']=str(pve)
    r['note']=re.sub(r'\s{2,}',' ',((r.get('note') or '')+' '+how).strip()); out[i]=r
adset('V1164',dis='Bacterial leaf streak',stat='-log10(P) 2.59',pve=1.6,how=f'{AD} lists wPt-663764 (6B, 31.7 cM) under bacterial leaf streak, not spot blotch; disease corrected 25 Nov 2026.')
adset('V1165',619.0769,619.0769,dis='Bacterial leaf streak',stat='-log10(P) 3.82',pve=2.6,how=f'{AD} lists wPt-5674 (7D, 174.2 cM) under bacterial leaf streak, not spot blotch; disease corrected. Placed on RefSeq v1.0 at 619.08 Mb from the Crop J 2021 consensus marker table, 25 Nov 2026.')
adset('V1168',741.0,741.0,stat='-log10(P) 3.21',pve=2.2,how=lab(f'{AD}: 7B 223.5 cM. Placed in silico on RefSeq v1.0 on 25 Nov 2026 at 741.0 Mb (interpolated on the Marone2012 map, single anchor within 3 cM; relaxed rule).')+CHK)
adset('V1166',stat='-log10(P) 3.08',pve=2.1); adset('V1167',stat='-log10(P) 3.12',pve=2.2)
NA=[('wPt-4886','Bacterial leaf streak','1A','60.7',2.29,1.4,(133.425,'DArT tag BLAST (this atlas)')),('wPt-665259','Bacterial leaf streak','1A','75.2',2.51,1.6,(483.4889,'DArT tag BLAST (this atlas)')),
    ('wPt-671483','Bacterial leaf streak','1A','75.2',2.56,1.6,(483.4889,'DArT tag BLAST (this atlas)')),('wPt-2780','Bacterial leaf streak','4A','90.6',2.37,1.5,((715.94,720.43),'interpolated on the Maccaferri2014 map, anchors within 3 cM')),
    ('wPt-8292','Bacterial leaf streak','4B','110.8',2.39,1.5,(664.5267,'Crop J 2021 consensus marker table')),('wPt-1159','Spot blotch','3B','53.2',2.40,1.6,(82.449,'DArT tag BLAST (this atlas)')),
    ('wPt-664459','Spot blotch','7D','175.1',2.37,1.5,None)]
n=0
for m,dis,ch,cm,lp,pve,pl in NA:
    n+=1;r=dict(base,id=f'PG6-{n:03d}',name=m,cls='MTA',disease=dis,chr=ch,marker=m,source='V_Adhikari2012',stat=f'-log10(P) {lp}',pve=str(pve),evidence=f'Table 2: {dis} {m} {ch} {cm} cM, p-value {lp}, R2 {pve}%')
    r['note']=f'Association mapping in 566 spring wheat landraces (USDA NSGC), DArT markers; chromosome and cM from WheatDArTmaps (Table 2).'
    if pl:
        lo,hi=(pl[0] if isinstance(pl[0],tuple) else (pl[0],pl[0])); setpos(r,lo,hi)
        if 'interpolated' not in pl[1]: r['grade']='C'; r['curation']='auto (in-silico from markers, Nov 2026)'
        r['note']+=lab(f' Placed on RefSeq v1.0 on 25 Nov 2026 from {m}: {lo:.2f}'+(f'-{hi:.2f}' if hi-lo>=0.05 else '')+f' Mb ({pl[1]}).')+(CHK if 'interpolated' in pl[1] else '')
    else: r['note']+=' Marker has no known physical position; left unplaced.'
    r['note']+=' Read from the full-text PDF supplied by the owner on 25 Nov 2026; not yet checked by hand.'
    out[r['id']]=r
json.dump(list(out.values()),open('rec96.json','w'))
need=[{'rec':r} for r in out.values() if r.get('start') is not None or r.get('peak') is not None]
json.dump(need,open('/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad/pan/need98.json','w'))
print(len(out),len(need))
