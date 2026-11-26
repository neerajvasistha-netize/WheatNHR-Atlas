import json,collections,glob
from openpyxl import Workbook
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.utils import get_column_letter
G=json.load(open('gg_list.json'));d=json.load(open('all_now.json'));R={r['id']:r for r in d['reports']};S=d['sources']
WHEAT={f'{i}{g}' for i in range(1,8) for g in 'ABD'}
# results of the 26 Nov in-silico PCR of owner-supplied primers
W='../wt1/';P=json.load(open(W+'primers.json'));O=json.load(open(W+'post.json'));SQ=json.load(open(W+'seqs.json'))
hits=collections.defaultdict(set)
for f in glob.glob(W+'res_*.json'):
    ch=f.split('res_')[1][:-5]
    for m,hs in json.load(open(f))['epcr'].items():
        if any(h[3]+h[4]<=4 for h in hs): hits[m].add(ch)
def tried(name):
    for m in list(P)+list(SQ):
        if m.lower().lstrip('x')==name.lower().lstrip('x'):
            if m in P:
                o=sorted(hits.get(m,set())-set(P[m][2]))
                return 'primers tested 26 Nov: '+(f'product only on {", ".join(o)}' if o else 'no product in Chinese Spring')
            return 'sequence tested 26 Nov: '+O.get('SEQ:'+m,{}).get('status','no alignment')
    return ''
pri={'high':0,'medium':1,'low':2}
rows=sorted(G.values(),key=lambda g:(pri[g['t'][2]],g['t'][0],g['m']))
M=[];PAP=collections.OrderedDict()
for g in rows:
    links=[];
    for i in dict.fromkeys(g['ids']):
        s=R.get(i,{}).get('source');v=S.get(s,{});u=v.get('url') or ''
        if u and u not in links: links.append(u)
        PAP.setdefault(s,{'short':v.get('short',s),'cit':v.get('citation',''),'url':u,'markers':set(),'genes':set()})['markers'].add(g['m'])
    M.append([g['m'],g['t'][0],g['t'][2],g['t'][1],', '.join(g['ch']),len(set(g['ids'])),'; '.join(dict.fromkeys(g['ids'])),'; '.join(list(dict.fromkeys(g['loci']))[:6]),'; '.join(g['refs']),'\n'.join(links) or 'no link in atlas',tried(g['m'])])
gn=[]
for r in d['reports']:
    if r.get('start') is not None or r.get('peak') is not None: continue
    if r.get('cls') not in ('DG','SG') or r.get('chr') not in WHEAT: continue
    s=S.get(r.get('source'),{})
    gn.append([r.get('name',''),r.get('cls'),r.get('disease',''),r.get('chr'),r['id'],r.get('marker',''),s.get('short',''),s.get('url',''),'GrainGenes gene report: linked markers and map position; a sequence only if the gene is cloned'])
    PAP.setdefault(r.get('source'),{'short':s.get('short',''),'cit':s.get('citation',''),'url':s.get('url',''),'markers':set(),'genes':set()})['genes'].add(r.get('name',''))
gn.sort(key=lambda x:(x[2],x[0]))
wb=Workbook();ws0=wb.active;ws0.title='Summary'
H=Font(name='Arial',size=10,bold=True,color='FFFFFF');FL=PatternFill('solid',fgColor='2F5D3A');F=Font(name='Arial',size=10)
def sheet(t,head,data,w):
    ws=wb.create_sheet(t);ws.append(head)
    for x in data: ws.append(x)
    for c in ws[1]: c.font=H;c.fill=FL;c.alignment=Alignment(wrap_text=True,vertical='top')
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font=F;c.alignment=Alignment(wrap_text=True,vertical='top')
            if isinstance(c.value,str) and c.value.startswith('http') and '\n' not in c.value: c.hyperlink=c.value;c.font=Font(name='Arial',size=10,color='1F4E9E',underline='single')
    for i,x in enumerate(w,1): ws.column_dimensions[get_column_letter(i)].width=x
    ws.freeze_panes='B2';ws.auto_filter.ref=ws.dimensions
sheet('Unpositioned markers',['Marker','Type','Chance GrainGenes has the sequence','What to take from GrainGenes','Chromosome(s) needed','Atlas records','Atlas IDs','Loci','References','Paper links','Already tested'],M,[14,20,12,44,12,9,30,36,36,36,30])
sheet('Unpositioned genes',['Gene','Class','Disease','Chromosome','Atlas ID','Markers in atlas','Reference','DOI / link','What GrainGenes has'],gn,[22,7,22,9,10,24,22,30,50])
sheet('Papers',['Reference','Full citation','Link','Unpositioned markers','Unpositioned genes'],[[x['short'],x['cit'],x['url'] or 'no link in atlas (search the title)','; '.join(sorted(m for m in x['markers'] if m)),'; '.join(sorted(m for m in x['genes'] if m))] for x in sorted(PAP.values(),key=lambda x:x['short'] or '')],[24,70,36,40,40])
c=collections.Counter((g['t'][0],g['t'][2]) for g in rows)
info=[['Unpositioned markers and genes to look up in GrainGenes (updated 26 Nov 2026)'],[],
['Scope','Only markers and genes of atlas records that still have no physical position (wheat chromosomes; meta-QTL excluded). Markers that already have a RefSeq v1.0 position are left out.'],
['Markers',f'{len(M)} markers from {len({i for g in rows for i in g["ids"]})} atlas records'],]+[[f'  {k[0]}',f'{v} (chance the sequence is in GrainGenes: {k[1]})'] for k,v in sorted(c.items(),key=lambda t:pri[t[0][1]])]+[
['Genes',f'{len(gn)} named genes / effector-sensitivity genes without a position'],
['Papers',f'{len(PAP)} papers ({sum(1 for x in PAP.values() if x["url"])} with a link)'],
['Already tested','Column K of the markers sheet: primers or sequences you supplied on 26 Nov that were tested and gave no usable position (no need to look them up again).'],[],
['How to search','GrainGenes (https://wheat.pw.usda.gov or https://graingenes.org): type the marker name in the search box; the marker report gives primer sequences (SSR, STS) or the probe report links the clone sequence (RFLP). Save them in this workbook and send it back for in-silico PCR or BLAST on RefSeq v1.0.']]
for x in info: ws0.append(x)
ws0['A1'].font=Font(name='Arial',size=12,bold=True)
for row in ws0.iter_rows(min_row=2):
    for c_ in row: c_.font=Font(name='Arial',size=10,bold=(c_.column==1 and bool(c_.value) and not str(c_.value).startswith('  ')));c_.alignment=Alignment(wrap_text=True,vertical='top')
ws0.column_dimensions['A'].width=28;ws0.column_dimensions['B'].width=110
out='/home/claude/wheatnhr-atlas/curation/WheatTraitMap_unpositioned_markers_genes_for_GrainGenes.xlsx'
wb.save(out);print(len(M),len(gn),len(PAP),sum(1 for x in M if x[-1]))
