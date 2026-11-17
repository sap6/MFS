"""Titles and authors of the arXiv entries cited in refs.bib (slow: the API is rate-limited)."""
import re
import time
import urllib.request

for aid in ("2310.01665", "2610.26540", "2605.15527", "2603.19113"):
    for attempt in range(6):
        try:
            req = urllib.request.Request("https://export.arxiv.org/api/query?id_list=" + aid, headers={"User-Agent": "bibcheck/0.1"})
            x = urllib.request.urlopen(req, timeout=60).read().decode()
            e = x[x.index("<entry>"):]
            t = re.sub(r"\s+", " ", re.search(r"<title>(.*?)</title>", e, re.S).group(1))
            print(aid, "|", t, "|", "; ".join(re.findall(r"<name>(.*?)</name>", e)), "|",
                  re.search(r"<published>(.*?)</published>", e).group(1)[:10], flush=True)
            break
        except Exception as ex:
            print(aid, "attempt", attempt, "failed:", ex, flush=True)
            time.sleep(60)
    time.sleep(20)
