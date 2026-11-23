import csv,json,math,re,collections,sys
# usage: python3 build_tissue.py E-MTAB-4484-tpms.tsv web/genes_all.json web/tissue_expr.json
O={'Root':['g11','g4','g9'],'Leaf':['g10','g1','g18'],'Stem':['g5','g7','g3'],'Spike':['g6','g8','g2'],'Grain':['g17','g19','g20']}
G=json.load(open(sys.argv[2]));have={c:{g[1] for g in v} for c,v in G['g'].items()};out=collections.defaultdict(list)
with open(sys.argv[1]) as f:
    r=csv.reader(f,delimiter='\t');h=next(r);ix={k:i for i,k in enumerate(h)}
    for row in r:
        m=re.fullmatch(r'TraesCS([1-7][ABD])02G(\d{6})',row[0])
        if not m or int(m.group(2)) not in have.get(m.group(1),()): continue
        vals=[sum(float(row[ix[g]].split(',')[2]) if row[ix[g]] else 0.0 for g in gs)/len(gs) for gs in O.values()]
        if max(vals)<0.5: continue
        out[m.group(1)].append(str(int(m.group(2)))+','+','.join(f'{math.log2(v+1):.1f}'.rstrip('0').rstrip('.') or '0' for v in vals))
D={'src':'Pingault et al. 2015, Genome Biology 16:29 (doi 10.1186/s13059-015-0601-9): Chinese Spring root, leaf, stem, spike and grain, three developmental stages each; Expression Atlas E-MTAB-4484. Value: median TPM of each stage, averaged over the three stages of the organ, as log2(TPM+1).',
   'organs':list(O),'unit':'log2(TPM+1)','g':{c:';'.join(v) for c,v in out.items()}}
open(sys.argv[3],'w').write(json.dumps(D,separators=(',',':')))
