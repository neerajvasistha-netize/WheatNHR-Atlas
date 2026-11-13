"""Position fix (Nov 2026): Pundir et al. 2023 BMC Plant Biology (cereal cyst nematode), from the paper and Additional file 1.
- Main-effect QTL (Table 2): the "Position" column gives the RefSeq v1.1 position of the marker closest to the QTL (Table S3
  calls it the peak position). Stored as the peak (grade A); the interval of the flanking markers (Table S3) is kept where it spans more
  than one SNP. For QCcnr.ccsu-2A.1 and -2A.2 the CIM and ICIM positions differ; both are given in the note.
- Epistatic QCcnr.ccsu-2A.5 (WSNP6862-WSNP4934): both flanking markers have positions in Table S3 (WSNP6862 31.045 Mb,
  WSNP4934 50.039 Mb), so it is now placed (grade B). The other epistatic loci have only study-internal WSNP codes without
  positions and stay unplaced.
Writes data/merged_v13f.json, web/index.html (seed) and /home/claude/work/pu23_db.json."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift

T = ' Table 2 / Table S3 (RefSeq v1.1).'
def Q(ids, peak, start=None, end=None, ev='', extra=''):
    f = dict(peak=peak, grade='A', curation='hand', evidence=ev,
             note='Main-effect QTL for Heterodera avenae resistance, HUW 468 x C 306 RIL population. Peak = RefSeq v1.1 position of the marker closest to the QTL (Table 2; called peak position in Table S3).' + extra +
                  ' WSNP#### are the study internal labels for its 90K SNPs.')
    if start is not None: f.update(start=start, end=end)
    return (ids, f)
FIX = [
 Q(['V3162'], 49.940952, ev='Table 2 (CIM, R2): QCcnr.ccsu-1B.1 WSNP4047-WSNP3401 79.5 cM, 49,940,952; LOD 2.80, PVE 7%', extra=' Closest marker WSNP4047.'),
 Q(['V3163'], 49.927166, ev='Table 2 (CIM, CD): QCcnr.ccsu-1B.2 WSNP4047-WSNP3401 81.4 cM, 49,927,166; LOD 2.76, PVE 7%', extra=' Closest marker WSNP3401.'),
 Q(['V3164'], 70.663028, 59.193613, 70.663028, 'Table 2 (ICIM, CD): QCcnr.ccsu-1B.3 WSNP1974-WSNP5998 83.0 cM, 70,663,028; Table S3 1B:59193613-70663028',
   ' Interval of the flanking markers WSNP1974-WSNP5998 from Table S3.'),
 Q(['V3165'], 31.045433, 31.045433, 50.038871, 'Table 2: QCcnr.ccsu-2A.1 WSNP6862-WSNP4934, CIM 18.0 cM 31,045,433 (CD, LOD 3.62, PVE 14%); ICIM 19.0 cM 42,996,081 (CD)',
   ' CIM position 31.045 Mb (stored); ICIM gave 42.996 Mb. Interval of the flanking markers WSNP6862 (31.045 Mb) and WSNP4934 (50.039 Mb), Table S3. Also significant in the permutation test.'),
 Q(['V3166'], 50.038871, 31.045433, 50.038871, 'Table 2: QCcnr.ccsu-2A.2 WSNP6862-WSNP4934, CIM 22.8 cM 50,038,871 (R1, LOD 2.59); ICIM 22.0 cM 43,263,254 (R1, LOD 3.37)',
   ' CIM position 50.039 Mb (stored); ICIM gave 43.263 Mb. Interval of the flanking markers WSNP6862 (31.045 Mb) and WSNP4934 (50.039 Mb), Table S3.'),
 Q(['V3167'], 43.476992, 42.194565, 43.476992, 'Table 2 (CIM, R2): QCcnr.ccsu-2A.3 WSNP208-WSNP5183 27.5 cM, 43,476,992; Table S3 2A:42194565-43476992',
   ' Interval of the flanking markers WSNP5183-WSNP208 from Table S3.'),
 Q(['V3168'], 747.233694, 746.982564, 747.26645, 'Table 2 (ICIM, R2): QCcnr.ccsu-2A.4 WSNP181-WSNP6519 210.0 cM, 747,233,694; Table S3 2A:746982564-747266450',
   ' Interval of the flanking markers WSNP6519-WSNP181 from Table S3. Also significant in the permutation test.'),
 Q(['V3169'], 498.529974, 432.953612, 498.529974, 'Table 2 (CIM, R1): QCcnr.ccsu-3A.1 WSNP1702-WSNP4151 46.6 cM, 498,529,974; Table S3 3A:432953612-498529974',
   ' Interval of the flanking markers WSNP1702-WSNP4151 from Table S3.'),
 (['V3332'], dict(start=31.045433, end=50.038871, grade='B', curation='hand',
     note='QTL for Heterodera avenae resistance from the digenic epistatic analysis of the HUW 468 x C 306 RIL population (combined data; 7 cM on the study map), interacting with QCcnr.ccsu-4A.2 (LOD 3.40, PVE 5.76%).'
          ' Placed between its flanking markers WSNP6862 (31.045 Mb) and WSNP4934 (50.039 Mb), whose RefSeq v1.1 positions are given in Table S3.',
     evidence='Table 3: QCcnr.ccsu-2A.5 (WSNP6862-WSNP4934) x QCcnr.ccsu-4A.2; Table S3 WSNP6862 2A:31045433, WSNP4934 2A:50038871')),
]
SID = lambda i: 'S81-%03d' % i
def build(reports):
    byid = {r['id']: r for r in reports}
    old = {i for ids, _ in FIX for i in ids}
    new = []
    for n, (ids, f) in enumerate(FIX, 1):
        r = dict(byid[ids[0]]); r.update(f); r['id'] = SID(n); r['ref'] = f.get('ref', 'RefSeq v1.0')
        for k in ('start', 'end', 'peak'):
            if k not in f: r[k] = None
        r['note'] += ' Corrected from the paper, Nov 2026 (pu23.py).'
        r['asm'] = lift.lift_record(r)
        new.append(r)
    keep = [r for r in reports if r['id'] not in old and not r['id'].startswith('S81-')]
    return keep + new, new, sorted(old)

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], new, rem = build(M['reports'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</'))
    Sd['reports'], new2, rem2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': rem2}, open('/home/claude/work/pu23_db.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['chr'], r['start'], r['end'], r['peak'], r['grade'], r['asm'].get('RefSeq v2.1'))
    print('removed', rem2)
