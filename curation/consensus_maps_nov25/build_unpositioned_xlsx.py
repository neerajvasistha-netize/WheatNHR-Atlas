import json,re,collections,sys
sys.path.insert(0,'w'); from nm import split_markers
from openpyxl import Workbook
from openpyxl.styles import Font,PatternFill,Alignment
from openpyxl.utils import get_column_letter
d=json.load(open('all_now.json'));R=d['reports'];S=d['sources']
WHEAT={f'{i}{g}' for i in range(1,8) for g in 'ABD'}
CYTO=re.compile(r'monosomic|telosomic|nullisomic|ditelo|substitution line|c-banding|\bgish\b|\bfish\b|aneuploid|deletion bin|bin map',re.I)
NOPOS=re.compile(r'\b(AFLP|RAPD|isozyme|E-?[ACGT]{3}/?M-?[ACGT]{3}|P\d+/M\d+|OP[A-Z]\d+)\b',re.I)
MK=re.compile(r'\b(X?(?:gwm|wmc|barc|cfd|cfa|gdm|wms|psp|cfe|swes|ksum|wmm|cwem|ncw|fcp)\d+[a-z]?|X?(?:psr|bcd|cdo|wg|mwg|abg|glk|ksu[a-z]?|fba|fbb|abc|rz|tam|bcd|mta|whs)\d+[a-z]?|IW[AB]\d+|AX-\d+|wsnp_\w+|BS\d{8}_\d+|(?:Kukri|Excalibur|RAC875|Tdurum_contig|BobWhite|Ku|Ra|RFL_Contig|IAAV|CAP\d+|GENE|BS)_c?\w+_\d+|wPt-?\d+|tPt-?\d+|rPt-?\d+|DArT\w*\d+|S\d[ABD]_\d{5,})\b',re.I)
def mtype(m):
    k=m.lower().lstrip('x')
    if re.match(r'(gwm|wmc|barc|cfd|cfa|gdm|wms|psp|cfe|swes|ksum|wmm|cwem|ncw)\d',k): return 'SSR'
    if re.match(r'(psr|bcd|cdo|wg|mwg|abg|glk|ksu|fba|fbb|abc|rz|tam|mta|whs)\d',k): return 'RFLP'
    if re.match(r'(wpt|tpt|rpt|dart)',k) or re.match(r'^\d{6,8}$',k): return 'DArT'
    if re.match(r'(iw[ab]|ax-|wsnp|bs\d|kukri|excalibur|rac875|tdurum|bobwhite|ku_|ra_|rfl|iaav|cap|gene|s\d[abd]_)',k): return 'SNP'
    if re.search(r'kasp|snp',k): return 'SNP'
    return 'other'
HINT={'SSR':'SSR primer sequences (GrainGenes marker report or the original SSR paper, e.g. Röder 1998, Somers 2004, Song 2005 for barc, Sourdille/Guyomarc\'h for cfa/cfd/gpw); run in-silico PCR on RefSeq v1.0',
 'RFLP':'RFLP probe/clone sequence in GenBank (search the clone name, e.g. "PSR546" or "BCD1434"); BLAST it against RefSeq v1.0 and take the copy on this chromosome',
 'DArT':'DArT clone sequence (paper supplement, or the Diversity Arrays marker sequence file); BLAST against RefSeq v1.0',
 'SNP':'SNP flanking sequence from the paper\'s supplementary table or the array manifest (90K iSelect, 35K/660K Axiom, DArTseq tag); BLAST against RefSeq v1.0 or read the position if given',
 'other':'marker sequence (primer or probe) from the paper or its supplementary file; run in-silico PCR / BLAST on RefSeq v1.0'}
rows=[];skip=collections.Counter()
for r in R:
    if r.get('start') is not None or r.get('peak') is not None: continue
    if r.get('cls')=='MQTL': skip['meta-QTL']+=1; continue
    ch=r.get('chr')
    if ch not in WHEAT: skip['alien or unknown chromosome']+=1; continue
    ev=r.get('evidence') or ''; note=r.get('note') or ''; txt=ev+' '+note
    if 'reject' in note.lower() or 'paralog' in note.lower(): skip['earlier placement rejected (paralogue / wrong arm)']+=1; continue
    mks=split_markers(r.get('marker'))
    found=[m for m in dict.fromkeys(x.group(0) for x in MK.finditer(r.get('name','')+' '+ev+' '+note)) if m.lower() not in {x.lower() for x in mks}]
    allm=mks+found
    types={mtype(m) for m in allm}
    if not allm:
        if CYTO.search(txt): skip['cytogenetic stock only (monosomic, substitution, GISH etc.)']+=1; continue
        if NOPOS.search(txt): skip['AFLP / RAPD / isozyme markers only']+=1; continue
        cm=re.search(r'(\d+(?:\.\d+)?)\s*cM',txt)
        if r.get('cls') in ('MTA','QTN'): hint='GWAS: the supplementary MTA table normally gives the SNP name, its sequence or its physical position; take the peak SNP'
        elif r.get('cls')=='DG': hint='Named gene: look for a later molecular-mapping or fine-mapping paper (search "<gene> fine mapping" or "<gene> molecular markers"); flanking markers there can be placed in silico'
        elif cm: hint=f'Paper gives a genetic position ({cm.group(0)}): read the flanking markers at that position from the linkage map figure or table, then place them in silico'
        elif re.search(r'catalogue of gene symbols',ev,re.I): hint='Source is the Catalogue of Gene Symbols: find the original paper it cites (number in braces), then take the flanking markers from its QTL table'
        else: hint='Read the QTL table (main text or supplementary file) for the flanking or peak markers and their cM; then place them in silico'
        how='no marker in atlas'
    else:
        if types<= {'other'} and NOPOS.search(' '.join(allm)): skip['AFLP / RAPD / isozyme markers only']+=1; continue
        hint=' | '.join(HINT[t] for t in sorted(types))
        how='markers known but no physical position'
    s=S.get(r.get('source'),{})
    rows.append([r['id'],r.get('name',''),r.get('cls',''),r.get('disease',''),ch,'; '.join(mks),'; '.join(found[:8]),how,hint,s.get('short',r.get('source','')),s.get('citation',''),s.get('url',''),ev[:300],note[:300]])
rows.sort(key=lambda x:(x[9],x[4],x[1]))
# LOD sheet: linkage-mapping QTL without LOD
lod=[]
for r in R:
    if r.get('cls')!='QTL': continue
    st=(r.get('stat') or '')
    if re.search(r'\bLOD\b',st,re.I): continue
    s=S.get(r.get('source'),{}); pos=r.get('peak') if r.get('peak') is not None else (f"{r['start']}-{r['end']}" if r.get('start') is not None else '')
    gw=bool(re.search(r'log10|p-?value|\bP\b',st))
    lod.append([r['id'],r.get('name',''),r.get('disease',''),r.get('chr',''),pos,r.get('pve',''),st,
       'Statistic given is -log10(P) (association test), not a LOD; check whether the paper also reports LOD from interval mapping' if gw else 'QTL table of the paper (LOD column), or its supplementary file; for ICIM/CIM papers the LOD is usually beside PVE',
       s.get('short',r.get('source','')),s.get('citation',''),s.get('url','')])
lod.sort(key=lambda x:(x[8],x[3],x[1]))
wb=Workbook();F=Font(name='Arial',size=10);H=Font(name='Arial',size=10,bold=True,color='FFFFFF');FILL=PatternFill('solid',fgColor='2F5D3A')
def sheet(ws,head,data,widths):
    ws.append(head)
    from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE as IL
    for row in data: ws.append([IL.sub(' ',v).replace('—','-') if isinstance(v,str) else v for v in row])
    for c in ws[1]: c.font=H;c.fill=FILL;c.alignment=Alignment(wrap_text=True,vertical='top')
    for row in ws.iter_rows(min_row=2):
        for c in row: c.font=F;c.alignment=Alignment(wrap_text=True,vertical='top')
    for i,w in enumerate(widths,1): ws.column_dimensions[get_column_letter(i)].width=w
    ws.freeze_panes='C2';ws.auto_filter.ref=ws.dimensions
    for row in ws.iter_rows(min_row=2):
        for c in row:
            if isinstance(c.value,str) and c.value.startswith('http'): c.hyperlink=c.value;c.font=Font(name='Arial',size=10,color='1F4E9E',underline='single')
ws=wb.active;ws.title='Summary'
ws1=wb.create_sheet('Unpositioned records');ws2=wb.create_sheet('QTL without LOD')
sheet(ws1,['Atlas ID','Locus','Class','Disease','Chromosome','Markers in atlas record','Other marker names in the evidence or note','Status','Where to get the position (hint)','Reference','Full citation','DOI / link','Evidence (from the paper)','Atlas note'],rows,[10,28,7,18,8,22,22,16,60,20,50,28,45,45])
sheet(ws2,['Atlas ID','QTL','Disease','Chromosome','Position (Mb, RefSeq v1.0)','PVE (%)','Statistic in atlas','Where to get the LOD (hint)','Reference','Full citation','DOI / link'],lod,[10,28,18,8,14,8,20,55,20,50,28])
by=collections.Counter(x[7] for x in rows)
info=[['WheatTraitMap: records still needing a position or a LOD value (25 Nov 2026)'],[],
 ['Sheet "Unpositioned records"',f'{len(rows)} records with a wheat chromosome but no physical position, for which a position can still be found'],
 ['  markers known but no physical position',by.get('markers known but no physical position',0)],['  no marker in atlas (paper must be read)',by.get('no marker in atlas',0)],[],
 ['Left out (position cannot be found this way)']]+[[f'  {k}',v] for k,v in skip.most_common()]+[[],
 ['Sheet "QTL without LOD"',f'{len(lod)} linkage-mapping QTL records whose LOD value is not yet in the atlas (needed for meta-QTL projection). Meta-QTL records are not included.'],[],
 ['How to use','Filter by Reference to work paper by paper. The hint column says where the missing value is usually found. Positions found can be sent back for in-silico placement.']]
for x in info: ws.append(x)
ws['A1'].font=Font(name='Arial',size=12,bold=True)
for row in ws.iter_rows(min_row=2):
    for c in row: c.font=Font(name='Arial',size=10,bold=(c.column==1 and c.value and not str(c.value).startswith('  ')))
ws.column_dimensions['A'].width=48;ws.column_dimensions['B'].width=110
out='/home/claude/wheatnhr-atlas/curation/WheatTraitMap_unpositioned_records_and_missing_LOD.xlsx'
wb.save(out);print(len(rows),by,dict(skip),len(lod),len(set(x[9] for x in rows)),'papers')
