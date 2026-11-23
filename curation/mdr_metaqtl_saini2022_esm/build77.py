"""GWAS MTAs (Online Resource 6) and two known genes (Online Resource 7) from the supplementary file of
Saini et al. 2022 (Planta, multi-disease meta-QTL; 11032_2022_1282_MOESM1_ESM.xls), credited to the original papers.
Meta-QTLs themselves are not added.
Position rule: published physical position of the marker from our marker table if available, else the bp value in
the supplement (dropped where it is inconsistent with the paper's other rows).
usage: python3 build77.py esm6_rows.json cr77.jsonl out.json"""
import json,sys,re
rows=json.load(open(sys.argv[1]));CR={json.loads(l)[0].lower():json.loads(l) for l in open(sys.argv[2])}
PAP={'Gerard et al. (2017)':('10.1007/s10681-016-1820-1','Septoria tritici blotch'),'Muqaddasi et al. (2018)':('10.3835/plantgenome2018.05.0029','Septoria tritici blotch'),
     'Gupta et al. (2019)':('10.3390/ijms20133124','Karnal bunt'),'Singh et al. (2020)':('10.1038/s41598-020-62711-7','Karnal bunt'),'Hu et al. (2020)':('10.1186/s12870-019-2177-0','Fusarium head blight')}
AB={'Septoria tritici blotch':'STB','Karnal bunt':'Karnal bunt','Fusarium head blight':'FHB'}
DROP={('Hu et al. (2020)','Excalibur_c22724_85')}   # 62.22 Mb in the supplement vs 622 Mb for all its neighbours (typo)
def cite(doi):
    c=CR[doi.lower()];return {'citation':f"{c[1]} ({c[2]}) {c[3]}. {c[4]} {c[5]}:{c[6]}",'url':'https://doi.org/'+doi,'short':c[1].split(',')[0].split(' ')[0]+(' et al.' if ',' in c[1] else '')+f" ({c[2]})"}
src={};out=[];seen=set();n=0;cnt={}
for ref,(doi,dis) in PAP.items():
    sid='E7_'+re.sub(r'\W','',doi)[-10:];src[sid]=cite(doi)
    for x in rows:
        if x['ref']!=ref: continue
        for mk in [m.strip().replace('.0','') for m in x['mk'].split(',')][:1]:
            if (ref,mk) in DROP or not mk or mk=='-': continue
            if x['v1'] is not None: p,rf,how=x['v1'],'RefSeq v1.0','published physical position of the marker'
            elif x['v2'] is not None: p,rf,how=x['v2'],'RefSeq v2.1','published physical position of the marker'
            elif x['tab']: p,rf,how=x['tab'],'RefSeq v1.0','bp position given in the supplement'
            else: continue
            k=(sid,mk,x['chr'])
            if k in seen: continue
            seen.add(k);n+=1;cnt[ref]=cnt.get(ref,0)+1
            out.append(dict(id='E7-%03d'%n,name=f"{mk} ({AB[dis]})",cls='MTA',disease=dis,chr=x['chr'],start=None,end=None,peak=round(p,4),ref=rf,
                marker=mk,pve='',stat='',source=sid,grade='C',curation='auto (supplementary compilation)',
                evidence=f"Online Resource 6 of Saini et al. 2022: {mk} {x['chr']}",
                note=f"Significant GWAS marker from the original paper, as compiled in Online Resource 6 of Saini et al. (2022, multi-disease resistance meta-analysis). Position: {how} ({rf}). Not yet checked by hand against the original paper."))
G=[('Stb8','Septoria tritici blotch','7B',740.7206,'10.1094/phyto.2003.93.9.1158','linked SSR marker (primers F CCAAAAAAACTGCCTGCATG / R CTCTGGCATTGCTCCTTGG)'),
   ('Snn6','Septoria nodorum blotch','6A',554.8337,'10.1094/mpmi-12-14-0396-r','EST marker BE403326')]
for nm,dis,ch,p,doi,mk in G:
    sid='E7_'+re.sub(r'\W','',doi)[-10:];src[sid]=cite(doi);n+=1
    out.append(dict(id='E7-%03d'%n,name=nm,cls='DG',disease=dis,chr=ch,start=None,end=None,peak=p,ref='RefSeq v1.0',marker=mk.split(' (')[0],pve='',stat='',source=sid,grade='C',
        curation='auto (supplementary compilation)',evidence=f"Online Resource 7 of Saini et al. 2022: {nm} {ch}",
        note=f"Position of the {mk} as given in Online Resource 7 of Saini et al. (2022) (RefSeq v1.0). Not yet checked by hand against the original paper."))
json.dump(dict(reports=out,sources=src,counts=cnt),open(sys.argv[3],'w'),indent=0);print(n,cnt)
