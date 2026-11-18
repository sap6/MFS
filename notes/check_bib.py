"""Compare the titles and years in refs.bib with the Crossref record of each DOI.

    python3 check_bib.py > check_bib_report.txt     (about 10 s per entry; api.crossref.org is rate-limited)
"""
import json
import re
import time
import urllib.request

entries = re.split(r"\n(?=@)", open("refs.bib").read())
for e in entries:
    key = re.match(r"@\w+\{([^,]+),", e)
    doi = re.search(r"doi\s*=\s*\{([^}]+)\}", e)
    title = re.search(r"\btitle\s*=\s*\{(.+?)\},?\s*\n", e, re.S)
    if not key:
        continue
    if not doi:
        print(f"{key.group(1):18s} no DOI in the entry")
        continue
    try:
        req = urllib.request.Request("https://api.crossref.org/works/" + doi.group(1), headers={"User-Agent": "bibcheck/0.1"})
        rec = json.load(urllib.request.urlopen(req, timeout=40))["message"]
        t = (rec.get("title") or ["?"])[0]
        yr = (rec.get("issued", {}).get("date-parts") or [["?"]])[0][0]
        fam = [a.get("family", "?") for a in rec.get("author", rec.get("editor", []))]
        mine = re.sub(r"[{}\\$']", "", title.group(1)).lower() if title else ""
        norm = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())
        same = norm(mine) == norm(re.sub(r"<[^>]+>", "", t))
        print(f"{key.group(1):18s} {'OK   ' if same else 'DIFF '} {yr} {'; '.join(fam)[:60]} | {t}")
    except Exception as ex:
        print(f"{key.group(1):18s} lookup failed: {ex}")
    time.sleep(10)
