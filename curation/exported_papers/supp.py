import json,re,collections
SUP=re.compile(r'((?:Supplementa(?:l|ry)|Suppl\.?)\s+(?:Tables?|Data(?:\s*Sheet)?|Files?|Figures?|Figs?\.?|Materials?|Information)\s*S?\d+[A-Za-z]?|(?:Tables?|Figs?\.?|Figures?|Files?|Data(?:sets?)?)\s+S\d+[A-Za-z]?|Additional\s+[Ff]iles?\s+\d+|Online\s+Resources?\s+\d+|ESM\s*\d+|Supplementary\s+Data\s+\d+|Data\s+Sheet\s+\d+)')
GOOD=re.compile(r'marker|SNP|QTL|MTA|position|physical|linkage|map|chromosom|flank|interval|associat|loc(us|i)|cM|Mb|DArT|KASP|SSR|genotyp',re.I)
def supp(h,loci):
    t=re.sub(r'\s+',' ',open('txr/'+h+'.txt',errors='ignore').read())
    # cut references section
    i=max(t.rfind(' References '),t.rfind(' REFERENCES '),t.rfind('Literature Cited'))
    if i>len(t)*0.5: t=t[:i]
    sc=collections.Counter();ctx={}
    for m in SUP.finditer(t):
        lab=re.sub(r'\s+',' ',m.group(1)).strip().rstrip('.')
        lab=re.sub(r'^(Tables?)\s+S','Table S',lab);lab=re.sub(r'^Supplementa(l|ry)\s+Tables\b','Supplementary Table',lab)
        w=t[max(0,m.start()-200):m.end()+200]
        s=1+len(GOOD.findall(w))+5*sum(1 for n in loci if n in w)
        if re.search(r'Fig',lab): s*=0.5
        sc[lab]+=s
        if lab not in ctx or s>ctx[lab][0]: ctx[lab]=(s,t[max(0,m.start()-160):m.end()+120])
    return sc,ctx
if __name__=='__main__':
    import sys
    sc,ctx=supp(sys.argv[1],sys.argv[2:])
    for k,v in sc.most_common(8): print(round(v),k,'|',ctx[k][1][:260])
