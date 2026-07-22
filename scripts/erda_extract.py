#!/usr/bin/env python3
"""ERDA byte-range FASTA + duckdb/httpfs parquet extraction helpers.

Self-contained copy of the proven primitives from
``FDZ004.1_aws/scripts/build_subset.py`` so the RNA-seq pipeline deploys
standalone (no cross-subproject import, no dependency on an untracked sibling).
Kept deliberately small: just what reference-building needs — pull a few contig
records from ``all_contigs.fna`` and subset ``all_meta_ncbi_gff.parquet``.

The ERDA base URL is a BEARER credential: callers pass it in; it is never
written to any output file or log here.
"""
from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

FETCH_WORKERS = 8
MERGE_GAP = 262_144           # merge records into one GET if the gap <= 256 KB
MAX_CHUNK = 8 * 1024 * 1024   # cap a merged chunk at 8 MB (a single bigger record stands alone)

_tls = threading.local()


def _session():
    import requests  # lazy: only needed when actually fetching from ERDA
    s = getattr(_tls, "session", None)
    if s is None:
        s = requests.Session()
        _tls.session = s
    return s


def compute_read_bytes(length: int, line_bases: int, line_width: int) -> int:
    """Bytes to read for a FASTA record given its residue length and wrapping."""
    breaks_per_line = line_width - line_bases
    newline_count = (length - 1) // line_bases if line_bases > 0 else 0
    return length + newline_count * breaks_per_line


def http_range(url: str, start: int, nbytes: int) -> bytes:
    end = start + nbytes - 1
    r = _session().get(url, headers={"Range": f"bytes={start}-{end}"}, timeout=180)
    if r.status_code != 206:
        raise RuntimeError(f"expected HTTP 206, got {r.status_code} (range {start}-{end})")
    if len(r.content) != nbytes:
        raise RuntimeError(f"expected {nbytes} bytes, got {len(r.content)}")
    return r.content


def fai_rows(con, fai_url: str, id_key: str, len_key: str, ids: list[str]) -> list[dict]:
    """Fetch .fai.parquet index rows (offset + wrapping) for the requested ids."""
    con.execute("CREATE OR REPLACE TEMP TABLE want_fai(id VARCHAR)")
    con.executemany("INSERT INTO want_fai VALUES (?)", [(i,) for i in ids])
    q = (
        f'SELECT f."{id_key}", f."{len_key}", f.fasta_offset, f.line_bases, f.line_width '
        f"FROM read_parquet('{fai_url}') f SEMI JOIN want_fai w ON f.\"{id_key}\" = w.id"
    )
    return [
        {id_key: r[0], len_key: r[1], "fasta_offset": r[2], "line_bases": r[3], "line_width": r[4]}
        for r in con.execute(q).fetchall()
    ]


def _plan_records(rows: list[dict], id_key: str, len_key: str) -> list[dict]:
    recs = []
    for r in rows:
        length, lb, lw, off = int(r[len_key]), int(r["line_bases"]), int(r["line_width"]), int(r["fasta_offset"])
        if length < 1 or lb < 1:
            raise RuntimeError(f"bad length/line_bases for {r[id_key]}: {length}/{lb}")
        recs.append({"id": str(r[id_key]), "length": length, "offset": off,
                     "read_bytes": compute_read_bytes(length, lb, lw)})
    recs.sort(key=lambda x: x["offset"])
    for rec in recs:
        rec["end"] = rec["offset"] + rec["read_bytes"]
    return recs


def _plan_chunks(recs: list[dict]) -> list[dict]:
    chunks: list[dict] = []
    for rec in recs:
        if (chunks and rec["offset"] - chunks[-1]["end"] <= MERGE_GAP
                and (rec["end"] - chunks[-1]["start"]) <= MAX_CHUNK):
            chunks[-1]["end"] = max(chunks[-1]["end"], rec["end"])
            chunks[-1]["members"].append(rec)
        else:
            chunks.append({"start": rec["offset"], "end": rec["end"], "members": [rec]})
    return chunks


def _fetch_chunk(fasta_url: str, chunk: dict) -> dict:
    buf = http_range(fasta_url, chunk["start"], chunk["end"] - chunk["start"])
    out = {}
    for rec in chunk["members"]:
        s = rec["offset"] - chunk["start"]
        seq = buf[s:s + rec["read_bytes"]].replace(b"\n", b"").replace(b"\r", b"")
        if len(seq) != rec["length"]:
            raise RuntimeError(f"{rec['id']}: extracted {len(seq)} residues != indexed {rec['length']}")
        if b">" in seq:
            raise RuntimeError(f"{rec['id']}: '>' found in sequence bytes (offset misaligned)")
        out[rec["id"]] = seq
    return out


def extract_contigs(con, fasta_url: str, fai_url: str, contig_ids: list[str], out_fasta) -> dict[str, int]:
    """Extract the given contigs from ERDA into a single-line multi-FASTA.

    Returns {contig_id: length}. Raises if any requested contig is absent from
    the FAI or a byte-range slice fails its length/'>' integrity checks.
    """
    rows = fai_rows(con, fai_url, "contig_id", "contig_length", contig_ids)
    found = {str(r["contig_id"]) for r in rows}
    missing = sorted(set(contig_ids) - found)
    if missing:
        raise RuntimeError(f"{len(missing)} contigs absent from FAI: {missing}")
    recs = _plan_records(rows, "contig_id", "contig_length")
    chunks = _plan_chunks(recs)
    seqs: dict[str, bytes] = {}
    if chunks:
        with ThreadPoolExecutor(max_workers=FETCH_WORKERS) as ex:
            futures = [ex.submit(_fetch_chunk, fasta_url, c) for c in chunks]
            for fut in as_completed(futures):
                seqs.update(fut.result())
    lengths: dict[str, int] = {}
    with open(out_fasta, "wb") as fh:
        for rec in recs:  # deterministic offset order
            fh.write(f">{rec['id']}\n".encode())
            fh.write(seqs[rec["id"]])
            fh.write(b"\n")
            lengths[rec["id"]] = rec["length"]
    return lengths


def gff_rows_for_contig(con, gff_url: str, contig_id: str) -> list[dict]:
    """All GFF gene rows for one contig, ordered by start. Fails loud if the
    required columns are absent, so a schema drift can't silently return junk."""
    cols = {r[0] for r in con.execute(
        f"DESCRIBE SELECT * FROM read_parquet('{gff_url}')").fetchall()}
    required = {"contig", "start", "end", "strand", "type", "protein_id"}
    missing = required - cols
    if missing:
        raise RuntimeError(f"GFF parquet missing required columns {sorted(missing)}; has {sorted(cols)}")
    # quote start/end — both are duckdb reserved words. `type` carried so downstream can filter CDS.
    q = (f'SELECT contig, "start", "end", strand, type, protein_id FROM read_parquet(\'{gff_url}\') '
         f'WHERE contig = ? ORDER BY "start"')
    return [
        {"contig": r[0], "start": int(r[1]), "end": int(r[2]), "strand": r[3], "type": r[4], "protein_id": r[5]}
        for r in con.execute(q, [contig_id]).fetchall()
    ]
