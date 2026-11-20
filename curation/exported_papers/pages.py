import subprocess,json,glob,re,collections
out={}
for f in glob.glob('/mnt/user-data/uploads/Downloads/Exported /**/*.pdf',recursive=True):
    L=subprocess.run(['pdfimages','-list',f],capture_output=True,text=True).stdout.splitlines()[2:]
    pg=collections.Counter()
    for l in L:
        c=l.split()
        try:
            if int(c[3])>=400 and int(c[4])>=250: pg[int(c[0])]+=1
        except: pass
    # vector figure pages: caption on page
    n=int(re.search(r'Pages:\s+(\d+)',subprocess.run(['pdfinfo',f],capture_output=True,text=True).stdout).group(1))
    vec=set()
    for p in range(1,n+1):
        t=subprocess.run(['pdftotext','-f',str(p),'-l',str(p),f,'-'],capture_output=True,text=True).stdout
        if re.search(r'\n\s*(Fig\.?|Figure)\s*\d+[\s.:|]',t) : vec.add(p)
    out[f]=sorted(set(pg)|vec)
json.dump(out,open('pages.json','w'))
print(sum(len(v) for v in out.values()))
