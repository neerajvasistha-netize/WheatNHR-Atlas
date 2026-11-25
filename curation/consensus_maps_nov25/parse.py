import pandas as pd,re,json
I='../in/';rows=[]
def add(mp,ch,mk,cm):
    if not isinstance(mk,str) or pd.isna(cm): return
    try: cm=float(cm)
    except: return
    ch=str(ch).strip().upper()[:2]
    if re.match(r'^[1-7][ABD]$',ch): rows.append((mp,ch,mk.strip(),cm))
# Marone 2012 durum
d=pd.read_excel(I+'122_2012_1939_MOESM2_ESM.xls',header=None).iloc[9:]
for _,r in d.iterrows(): add('Marone2012',str(r[2]).split('-')[0],r[3],r[4]/100 if pd.notna(r[4]) else None)
# Maccaferri 2014 BMC Genomics durum framework
d=pd.read_excel(I+'12864_2014_6782_MOESM2_ESM.xlsx',header=None).iloc[6:]
for _,r in d.iterrows():
    m=r[0]
    if isinstance(m,str) and re.search(r'-[1-7][AB]',m):
        mk,ch=re.match(r'(.+)-([1-7][AB])',m.rstrip('+#*')).groups(); add('Maccaferri2014',ch,mk,r[1])
# Quraishi-type integrated map (Supplementary Table 2-8)
for i in range(2,9):
    x=pd.ExcelFile(I+f'Supplementary Table {i}.XLSX')
    for s in x.sheet_names:
        if not s.startswith('ch'): continue
        d=x.parse(s,header=None); c=[c for c in d.columns if (d[c].astype(str)==s[2:]).sum()>100][0]; d=d.iloc[6:]
        for _,r in d.iterrows(): add('IntegratedMQTLmap',r[c],r[c-2],r[c-1])
# Wen 2017 integrated consensus map (Table S6)
d=pd.read_excel(I+'Table_6.XLSX',header=None).iloc[3:]
for _,r in d.iterrows(): add('Wen2017',r[2],r[1],r[3])
# KOMUGI (NBRP Japan) Composite Wheat Map, ChromosomeMap.xls
x=pd.ExcelFile('../gg/komugi_zip/x/ChromosomeMap.xls')
for sname in x.sheet_names[1:]:
    ch=sname[-2:]; d=x.parse(sname,header=None).iloc[1:]
    for _,r in d.iterrows():
        m=str(r[1]); m=re.sub(r'-'+ch+r'(\.\d+)?$','',m); add('KOMUGIcomposite',ch,m,r[2])
df=pd.DataFrame(rows,columns=['map','chr','marker','cM']);df.to_csv('consensus_maps.tsv',sep='\t',index=False)
print(df.groupby('map').size())
