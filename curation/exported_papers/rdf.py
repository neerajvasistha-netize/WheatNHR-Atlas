import xml.etree.ElementTree as ET,json,re,glob
NS={'rdf':'http://www.w3.org/1999/02/22-rdf-syntax-ns#','z':'http://www.zotero.org/namespaces/export#','dc':'http://purl.org/dc/elements/1.1/','bib':'http://purl.org/net/biblio#','foaf':'http://xmlns.com/foaf/0.1/','link':'http://purl.org/rss/1.0/modules/link/','dcterms':'http://purl.org/dc/terms/','prism':'http://prismstandard.org/namespaces/1.2/basic/'}
R='{%s}'%NS['rdf']
meta={}
for f,pre in [('/mnt/user-data/uploads/Downloads/Exported /Exported Items/Exported Items.rdf','./Exported Items/'),('/mnt/user-data/uploads/Downloads/Exported /Exported Items_2/Exported Items_2.rdf','./Exported Items_2/'),('/mnt/user-data/uploads/Downloads/Exported /Exported .rdf','./')]:
    root=ET.parse(f).getroot()
    att={}
    for a in root.findall('z:Attachment',NS):
        p=a.find('z:path',NS)
        if p is not None: att['#'+a.get(R+'about').split('#')[-1]]=p.get(R+'resource')
    for it in root:
        if it.tag.endswith('Attachment') or it.tag.endswith('Memo'): continue
        ti=it.findtext('dc:title',default='',namespaces=NS)
        au=[(p.findtext('foaf:surname',default='',namespaces=NS),p.findtext('foaf:givenName',default='',namespaces=NS)) for p in it.findall('bib:authors/rdf:Seq/rdf:li/foaf:Person',NS)]
        j=it.find('dcterms:isPartOf/bib:Journal',NS)
        jt=j.findtext('dc:title',default='',namespaces=NS) if j is not None else ''
        vol=j.findtext('prism:volume',default='',namespaces=NS) if j is not None else ''
        pg=it.findtext('bib:pages',default='',namespaces=NS)
        dt=it.findtext('dc:date',default='',namespaces=NS)
        doi=''
        for v in it.iter(R+'value'):
            m=re.search(r'10\.\d{4,9}/\S+',v.text or '')
            if m: doi=m.group(0);break
        if not doi and j is not None:
            m=re.search(r'10\.\d{4,9}/\S+',j.findtext('dc:identifier',default='',namespaces=NS));doi=m.group(0) if m else ''
        for l in it.findall('link:link',NS):
            p=att.get(l.get(R+'resource'))
            if p: meta[pre+p if not p.startswith('./') else p]=dict(ti=ti,au=au,j=jt,vol=vol,pg=pg,yr=(re.search(r'\d{4}',dt) or [''])[0] if dt else '',doi=doi)
json.dump(meta,open('rdfmeta.json','w'))
o=json.load(open('s1.json'));hit=sum(1 for x in o if x['path'] in meta);print(len(meta),'matched',hit,'of',len(o))
print([x['path'] for x in o if x['path'] not in meta][:5])
