import json,re,sys,collections
sys.path.insert(0,'w');sys.path.insert(0,'../p26')
from nm import n2,split_markers; from look2 import pos
exec(open('build_xl.py').read().split('rows=[];skip')[0])   # MK regex, mtype, CYTO, NOPOS, WHEAT
d=json.load(open('all_now.json'));R=d['reports'];S=d['sources']
def gtype(m):
    k=n2(m)
    if re.match(r'(gwm|wmc|barc|cfa|cfd|gdm|psp|wms|cfe|ksum|gpw|swes|wmm|cwem|ncw|fcp|cnl|dupw|edm|mag|bcd|cdo)\d',k) and re.match(r'(gwm|wmc|barc|cfa|cfd|gdm|psp|wms|cfe|ksum|gpw|swes|wmm|cwem|ncw|fcp|cnl|dupw|edm|mag)\d',k): return 'SSR','primer sequences in the GrainGenes marker report (or the original SSR paper)','high'
    if re.match(r'(bcd|cdo|wg|mwg|abg|ksu[a-z]?|fba|fbb|abc|rz|tam|umc|bcd|cmwg|mta|whs|ksud|ksue|ksuf|ksug|ksuh)\d',k): return 'RFLP probe','probe/clone report with a GenBank accession for the clone sequence (BLAST it to RefSeq v1.0)','medium'
    if re.match(r'(psr|psb|pwir|wir)\d',k): return 'RFLP probe (PSR / pWIR)','probe report; few PSR clones have a public sequence','low'
    if re.match(r'glk\d',k): return 'RFLP probe (GLK, ETH Zurich)','probe report; GLK clone sequences are mostly not public','low'
    if re.match(r'(sts|scar|sc|bcd|mag|ubc|opa|opb|opc|opd|ope|opf|opg|oph|opi|opj|opk|opl|opm|opn|opo|opp|opq|opr|ops|opt|opu|opv|opw|opx|opy|opz)',k): return 'STS / SCAR','primer sequences if the marker was converted to an STS/SCAR','medium'
    return None
G=collections.OrderedDict()
for r in R:
    if r.get('start') is not None or r.get('peak') is not None: continue
    ch=r.get('chr')
    if ch not in WHEAT or r.get('cls')=='MQTL': continue
    ev=r.get('evidence') or '';note=r.get('note') or ''
    ms=split_markers(r.get('marker'))+[x.group(0) for x in MK.finditer(r.get('name','')+' '+ev+' '+note)]
    for m in dict.fromkeys(ms):
        t=gtype(m)
        if not t or pos(m,ch): continue
        k=n2(m); g=G.setdefault(k,{'m':m if m[0] in 'Xx' else 'X'+m,'t':t,'ch':set(),'ids':[],'loci':[],'refs':set()})
        g['ch'].add(ch); g['ids'].append(r['id']); g['loci'].append(r.get('name','')); g['refs'].add(S.get(r.get('source'),{}).get('short',r.get('source','')))
print(len(G),collections.Counter(g['t'][0] for g in G.values()))
json.dump({k:dict(v,ch=sorted(v['ch']),refs=sorted(v['refs'])) for k,v in G.items()},open('gg_list.json','w'))
