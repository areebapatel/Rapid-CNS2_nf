#!/usr/bin/env python3
"""Fetch MNP-Flex classification results for an uploaded sample.

Credentials come from EPIGNOSTIX_USER / EPIGNOSTIX_PASSWORD, as for the upload.
Exits 0 with a message if the analysis is still running, so it can be re-run.
"""
import argparse
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

API_ARTEFACTS = {
    "bundle_summary": "bundle_summary",
    "bundle": "bundle",
    "qc_coverage_plot": "qc_coverage_plot",
    "qc_methylation_density_plot": "qc_methylation_density_plot",
    "mgmt_region_plot": "mgmt_region_plot",
}


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def call(url, headers, data=None, method=None, raw=False):
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            body = r.read()
            if raw:
                return r.status, body, r.headers.get("Content-Type", "")
            return r.status, (json.loads(body.decode()) if body.strip() else {}), ""
    except urllib.error.HTTPError as e:
        return e.code, e.read()[:300], ""
    except urllib.error.URLError as e:
        die(f"cannot reach {url}: {e.reason}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample-id", type=int, required=True, help="platform sample id")
    p.add_argument("--api", default="https://app.epignostix.com/api")
    p.add_argument("--outdir", default=".")
    p.add_argument("--prefix", default="mnpflex")
    p.add_argument("--wait", type=int, default=0, help="seconds to poll for the run to finish")
    a = p.parse_args()

    user, password = os.environ.get("EPIGNOSTIX_USER"), os.environ.get("EPIGNOSTIX_PASSWORD")
    if not user or not password:
        die("EPIGNOSTIX_USER and EPIGNOSTIX_PASSWORD must be set")

    api = a.api.rstrip("/")
    form = urllib.parse.urlencode({"grant_type": "password",
                                   "username": user, "password": password}).encode()
    s, body, _ = call(f"{api}/v1/auth/token", {"Content-Type": "application/x-www-form-urlencoded"},
                      data=form, method="POST")
    if s != 200:
        die(f"login failed: HTTP {s}")
    H = {"Authorization": f"Bearer {body['access_token']}"}

    s, runs, _ = call(f"{api}/v1/workflow_runs/by_entity/{a.sample_id}", H)
    if s != 200 or not runs:
        die(f"no workflow run for sample {a.sample_id} (HTTP {s})")
    run = runs[0]
    run_id = run["id"]

    import time
    waited = 0
    while True:
        s, detail, _ = call(f"{api}/v1/workflow_runs/{run_id}", H)
        tasks = detail.get("task_runs", []) if s == 200 else []
        done = [t for t in tasks if t.get("task_result_id")]
        if done or waited >= a.wait:
            break
        time.sleep(15)
        waited += 15
    if not done:
        states = ", ".join(sorted({str(t.get("status")) for t in tasks})) or "unknown"
        print(f"analysis not finished for sample {a.sample_id} "
              f"(run {run_id}, status {detail.get('status')}, tasks: {states}) - try again later")
        return

    os.makedirs(a.outdir, exist_ok=True)
    result_id = done[-1]["task_result_id"]
    for name, ep in API_ARTEFACTS.items():
        s, blob, ctype = call(f"{api}/v1/mnpflex_sample/analysis/{ep}/{run_id}/{result_id}",
                              H, raw=True)
        if s != 200:
            print(f"  {name}: HTTP {s} (skipped)")
            continue
        ext = ".png" if "image" in ctype else (".zip" if "zip" in ctype else ".json")
        path = os.path.join(a.outdir, f"{a.prefix}_{name}{ext}")
        with open(path, "wb") as fh:
            fh.write(blob)
        print(f"  {name}: {len(blob)/1024:.0f} KB -> {path}")

    # compact table for the report: top predictions plus the QC line
    summary_path = os.path.join(a.outdir, f"{a.prefix}_bundle_summary.json")
    if os.path.exists(summary_path):
        d = json.load(open(summary_path))
        cs = d.get("classifier_summary", {})
        rows = sorted(cs.get("scores", []), key=lambda x: -float(x.get("score", 0)))
        out = os.path.join(a.outdir, f"{a.prefix}_mnpflex_predictions.tsv")
        with open(out, "w") as fh:
            clf = cs.get("classifier", {})
            fh.write(f"# classifier\t{clf.get('name','')} {clf.get('version','')}\n")
            fh.write(f"# qc\t{d.get('qc', {}).get('status','')}\n")
            fh.write("methylation_class\tscore\n")
            for r in rows[:10]:
                rg = r.get("reference_group") or {}
                name = rg.get("name") or f"group {r.get('reference_group_id')}"
                fh.write(f"{name}\t{float(r['score']):.4f}\n")
        print(f"  predictions table -> {out}")


if __name__ == "__main__":
    main()
