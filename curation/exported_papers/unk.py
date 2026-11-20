import json,re,collections,sys
sys.path.insert(0,'../pp')
from mp import look
need=json.load(open('need_fig.json'))
done={(r['h'],r['name']) for f in ['../fig/fig_res.json','../supp/supp_res.json','gbs_res.json'] for r in json.load(open(f))}
MK=re.compile(r'\b(?:X?(?:gwm|wmc|barc|cfd|cfa|gdm|wms|psp|cfb|wpt|wPt|gpw|swes|cwem|wmm|cau|xue|hbe|hbg|hbd|bwm|ksm|stars|ucw|wem)[-_]?\d+[A-Za-z]?|IW[AB]\d{2,6}|BS\d{8}(?:_\d+)?|AX-\d{8,9}|(?:Kukri|Excalibur|RAC875|BobWhite|Tdurum|wsnp|Ku|Ex|BS|Jagger|CAP|GENE|IAAV|tplb|TA|RFL|D_GDS7LZN|D_contig|Ra|Tdurum_contig)[A-Za-z0-9_]*_\d+[A-Za-z0-9_]*|KASP[-_][A-Za-z0-9_\-]+|K[-_][A-Za-z0-9_\-]{4,}|Xbcd\d+|Xcdo\d+|Xpsr\d+|Xglk\d+|Xmwg\d+)\b')
unk=collections.defaultdict(set);cnt=collections.Counter()
for h,n,ch in need:
    n=n.replace('‐','-').replace('–','-')
    if (h,n) in done: continue
    t=re.sub(r'\s+',' ',open('txr/'+h+'.txt',errors='ignore').read())
    for m in re.finditer(re.escape(n),t):
        w=t[max(0,m.start()-150):m.end()+250]
        for tok in MK.findall(w):
            if not ch: continue
            if look(tok,ch) or (tok[:1]=='X' and look(tok[1:],ch)): continue
            if look(tok) or (tok[:1]=='X' and look(tok[1:])): continue  # known but on other chromosome
            unk[tok].add((h,n,ch));cnt[(h,n)]+=1
json.dump({k:sorted(v) for k,v in unk.items()},open('unknown_markers.json','w'))
print(len(unk),'unknown markers for',len(cnt),'loci')
print(collections.Counter(re.sub(r'[\d_]+.*','',k) for k in unk).most_common(25))
print(list(unk)[:40])
