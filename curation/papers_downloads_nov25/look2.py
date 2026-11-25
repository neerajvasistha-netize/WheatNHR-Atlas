import sys,re,csv
sys.path.insert(0,'/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad/pp')
from mp import look,norm
C={}
for r in csv.DictReader(open('/home/claude/wheatnhr-atlas/curation/fhb_metaqtl_cropj2021_esm/consensus_marker_positions_refseq_v1.tsv'),delimiter='\t'):
    C.setdefault(norm(r[list(r)[0]]),[]).append(r)
def pos(m,ch):
    L=[x for x in look(m,ch) if x[2] in ('RefSeq v1.0','RefSeq v2.1')]
    v1=[x for x in L if x[2]=='RefSeq v1.0']
    if v1 or L: return tuple((v1 or L)[0])
    for r in C.get(norm(m),[]):
        if r['chr']==ch and r['pos_bp']:
            return (ch,round(int(r['pos_bp'])/1e6,4),'RefSeq v1.0','consensus marker table of the Crop J 2021 FHB meta-QTL study (Table S3)')
    return None
