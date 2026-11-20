import subprocess,json,os,hashlib,sys
from concurrent.futures import ThreadPoolExecutor
P=json.load(open('pages.json'));os.makedirs('ocr',exist_ok=True)
jobs=[(f,p) for f,ps in P.items() for p in ps]
def run(j):
    f,p=j;k=hashlib.md5(f.encode()).hexdigest()[:10]+'_%03d'%p
    if os.path.exists('ocr/'+k+'.tsv'): return
    subprocess.run(['pdftoppm','-r','300','-gray','-f',str(p),'-l',str(p),'-png',f,'ocr/'+k],check=False)
    png=[x for x in os.listdir('ocr') if x.startswith(k) and x.endswith('.png')]
    if not png: return
    subprocess.run(['tesseract','ocr/'+png[0],'ocr/'+k,'--psm','11','tsv'],capture_output=True)
    os.remove('ocr/'+png[0])
with ThreadPoolExecutor(2) as ex: list(ex.map(run,jobs))
json.dump({hashlib.md5(f.encode()).hexdigest()[:10]:f for f in P},open('keys.json','w'))
print('done')
