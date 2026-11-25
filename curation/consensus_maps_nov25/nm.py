import re
def n2(m):
    m=str(m).strip().lower()
    m=re.sub(r'[\s\-_\.()]','',m)
    m=re.sub(r'^x(?=[a-z])','',m).replace('wms','gwm')
    return re.sub(r'^(gwm|wmc|barc|cfd|cfa|gdm|psp|psr|bcd|cdo|wg|mwg|abg|glk|ksu|fba|fbb|abc|rz|tam|ksm|ksum|fcp)0+(\d)',r'\1\2',m)
def base(k): return re.sub(r'(d\d|[a-e]|\d?s|\d?l)$','',k) if re.search(r'\d(d\d|[a-e])$',k) else k
def split_markers(s):
    return [x.strip() for x in re.split(r'[;,/]| and |–|—| - ',s or '') if x.strip() and re.search(r'\d',x)]
