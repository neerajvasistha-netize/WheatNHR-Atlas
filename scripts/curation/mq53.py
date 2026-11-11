"""Remove all meta-QTL records (user request, Nov 2026: no meta-QTL papers in the atlas).
Whole sources that are meta-analyses: Liu et al. 2020 (tan spot MQTL), Sharma et al. 2024 (powdery mildew), Mao et al. 2010 (FHB),
and the Goudemand et al. 2013 meta-analysis locus. In primary-data papers only the MQTL records are removed (Pundir 2023, Lhamo 2024,
Arruda 2016). Writes data/merged_v13f.json, web/index.html (seed) and /home/claude/work/mq53_del.json."""
import json, re
WHOLE = {'Liu2020MQTL_TS', 'Sharma2024FPS', 'V_Mao2010'}
PART = {'Pundir2023BMC'}   # meta-QTL study; its MQTL records only (corrected: Lhamo 2024, Arruda 2016 and Goudemand 2013 are not meta-QTL studies and are kept)
def is_meta(r):
    if r['source'] in WHOLE: return True
    # Pundir et al. 2023: drop only its meta-QTL; keep its own QTL and the earlier CCN loci it anchored physically (user request)
    return r['source'] in PART and (r.get('cls') == 'MQTL' or r['name'].startswith('MQTL'))
def strip(reports):
    gone = [r['id'] for r in reports if is_meta(r)]
    return [r for r in reports if not is_meta(r)], gone
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], g1 = strip(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    S = json.loads(m.group(2).replace('<\\/', '</')); S['reports'], g2 = strip(S['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(S) + H[m.end(2):])
    json.dump(sorted(set(g1) | set(g2)), open('/home/claude/work/mq53_del.json', 'w'))
    print(len(g1), len(g2), len(set(g1) | set(g2)))
