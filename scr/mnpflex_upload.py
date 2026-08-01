#!/usr/bin/env python3
"""Upload an MNP-Flex bed file to the Epignostix research platform.

Credentials come from the environment (EPIGNOSTIX_USER, EPIGNOSTIX_PASSWORD)
and are never passed as command-line arguments, so they do not appear in `ps`
or in the Nextflow .command.* files.

This sends patient-derived methylation data to a third-party service. It runs
only when --mnpFlexUpload is set.
"""
import argparse
import json
import mimetypes
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
import uuid


def die(msg, code=1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def request(url, data=None, headers=None, method=None):
    req = urllib.request.Request(url, data=data, headers=headers or {}, method=method)
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            body = r.read().decode()
            return r.status, json.loads(body) if body.strip() else {}
    except urllib.error.HTTPError as e:
        detail = e.read().decode()[:400]
        die(f"{method or 'GET'} {url} -> HTTP {e.code}: {detail}")
    except urllib.error.URLError as e:
        die(f"cannot reach {url}: {e.reason}")


def get_token(api, user, password):
    form = urllib.parse.urlencode({
        "grant_type": "password", "username": user, "password": password,
    }).encode()
    _, body = request(f"{api}/v1/auth/token", data=form,
                      headers={"Content-Type": "application/x-www-form-urlencoded"},
                      method="POST")
    token = body.get("access_token")
    if not token:
        die("no access_token in the login response")
    return token


def multipart(fields, filename, filepath):
    """Build a multipart/form-data body; the file is the 'bed_file' part."""
    boundary = uuid.uuid4().hex
    out = []
    for k, v in fields.items():
        if v is None:
            continue
        out.append(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n".encode())
    ctype = mimetypes.guess_type(filename)[0] or "application/octet-stream"
    out.append(
        f"--{boundary}\r\nContent-Disposition: form-data; name=\"bed_file\"; "
        f"filename=\"{filename}\"\r\nContent-Type: {ctype}\r\n\r\n".encode())
    with open(filepath, "rb") as fh:
        out.append(fh.read())
    out.append(f"\r\n--{boundary}--\r\n".encode())
    return b"".join(out), f"multipart/form-data; boundary={boundary}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bed", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--api", default="https://app.epignostix.com/api")
    p.add_argument("--workflow-id", type=int, required=True)
    p.add_argument("--technology", required=True)
    p.add_argument("--target-coverage", required=True)
    p.add_argument("--extraction-type", required=True)
    p.add_argument("--sex", required=True)
    p.add_argument("--localisation")
    p.add_argument("--diagnosis")
    p.add_argument("--comment")
    p.add_argument("--out", default="mnpflex_upload.json")
    a = p.parse_args()

    user = os.environ.get("EPIGNOSTIX_USER")
    password = os.environ.get("EPIGNOSTIX_PASSWORD")
    if not user or not password:
        die("EPIGNOSTIX_USER and EPIGNOSTIX_PASSWORD must be set in the environment")
    if not os.path.isfile(a.bed):
        die(f"bed file not found: {a.bed}")

    api = a.api.rstrip("/")
    token = get_token(api, user, password)
    auth = {"Authorization": f"Bearer {token}"}

    query = urllib.parse.urlencode({
        "sample_identifier": a.sample,
        "used_technology": a.technology,
        "target_coverage": a.target_coverage,
        "extraction_type": a.extraction_type,
        "sex": a.sex,
        "workflow_id": a.workflow_id,
        "keep_filename": "true",
    })
    body, ctype = multipart(
        {"localisation": a.localisation, "diagnosis": a.diagnosis, "comment": a.comment},
        os.path.basename(a.bed), a.bed)
    status, sample = request(f"{api}/v1/mnpflex_sample?{query}", data=body,
                             headers={**auth, "Content-Type": ctype}, method="PUT")

    sample_id = sample.get("id") or sample.get("sample_id")
    print(f"uploaded {a.bed} as sample {a.sample} (id={sample_id}, HTTP {status})")

    run = {}
    if sample_id:
        _, run = request(
            f"{api}/v1/mnpflex_sample/{sample_id}/execute_workflow?workflow_id={a.workflow_id}",
            data=b"", headers=auth, method="POST")
        print(f"started workflow {a.workflow_id} -> run {run.get('id', '?')}")

    with open(a.out, "w") as fh:
        json.dump({"sample": sample, "workflow_run": run}, fh, indent=2)


if __name__ == "__main__":
    main()
