"""Spot blotch / black point MTAs from Supplementary Table S4 of the Sci Rep 2024 spot blotch meta-QTL study
(doi 10.1038/s41598-024-63924-w), for the original GWAS papers not yet in the atlas.
Only the original GWAS MTAs are used (no meta-QTL). Positions: our own lookup of the marker's published physical
position where available; otherwise the bp position in Table S4, used only for papers whose looked-up markers agree
with Table S4 (for papers where Table S4 gives genetic-map values, only looked-up markers are kept).
usage: python3 build76.py m4rows.json cr.jsonl out.json"""
import json,sys,re
rows=json.load(open(sys.argv[1]));CR=[json.loads(l) for l in open(sys.argv[2])]
PAPERS={'Adhikari et al. 2012':('P004','10.3835/plantgenome2011.12.0032','lookup'),'Ahirwar et al. 2018':('P011','10.1371/journal.pone.0208196','lookup'),
 'Gurung et al. 2014':('P007','10.1371/journal.pone.0108179','both'),'Jamil et al. 2018':('P013','10.1094/phyto-02-18-0047-r','both'),
 'Kumar et al. 2022':('P023','10.1186/s12870-022-04013-w','table'),'Lozano-Ramirez et al. 2022':('P022','10.3390/genes13081387','both'),
 'Liu et al. 2017':('P033','10.1186/s12870-017-1167-3','both')}
TR={'SB':'spot blotch','AUDPC':'spot blotch AUDPC','Severity':'spot blotch severity','LN':'lesion number','IP':'incubation period','BPR':'black point'}
src={};out=[];seen=set();n=0;cnt={}
for ref,(pid,doi,mode) in PAPERS.items():
    c=[x for x in CR if x[0].lower()==doi.lower()][0]
    pg=c[6];title=re.sub(r'</?i>','',c[3])
    sid='S4_'+pid
    src[sid]={'citation':f"{c[1]} ({c[2]}) {title}. {c[4].replace('®','')} {c[5]}:{pg}",'url':'https://doi.org/'+doi,'short':ref.replace(' 20',' (20')+')'}
    for x in rows:
        if x['ref']!=ref or x['mk'] in('NA','') : continue
        if re.fullmatch(r'\d+',x['mk']) and x['tab'] and abs(x['tab']-int(x['mk'])/1e6)<1e-3: continue   # position equals the marker ID: not a position
        if x['v1'] is not None: p,rf,how=x['v1'],'RefSeq v1.0','published physical position of the marker'
        elif x['v2'] is not None: p,rf,how=x['v2'],'RefSeq v2.1','published physical position of the marker'
        elif mode in('both','table') and x['tab']: p,rf,how=x['tab'],'RefSeq v1.0','bp position given in Table S4'
        else: continue
        k=(sid,x['mk'],x['chr'])
        if k in seen: continue
        seen.add(k);n+=1;cnt[ref]=cnt.get(ref,0)+1
        tr=TR.get(x['trait'],x['trait']);bp=ref=='Liu et al. 2017'
        out.append(dict(id='S4-%03d'%n,name=f"{x['mk']} ({'black point' if bp else 'spot blotch'})",cls='MTA',disease='Spot blotch',chr=x['chr'],start=None,end=None,peak=round(p,4),ref=rf,
            marker=x['mk'],pve='',stat='',source=sid,grade='C',curation='auto (supplementary compilation)',
            evidence=f"Table S4 of the 2024 Sci Rep spot blotch meta-QTL study: {x['mk']} {x['chr']}, trait {x['trait']}",
            note=f"Significant marker ({tr}) from the original GWAS, as compiled in Supplementary Table S4 of the 2024 Scientific Reports spot blotch meta-QTL study (doi 10.1038/s41598-024-63924-w). Position: {how} ({rf}). Not yet checked by hand against the original paper."))
json.dump(dict(reports=out,sources=src,counts=cnt),open(sys.argv[3],'w'),indent=0)
print(n,cnt)
