#!/usr/bin/env python3
"""Restore the four original files from a Parquet bundle written by ossuary.py:

    <name>.state          nucleotide ASR file
    <name>_binary.state   binary (gap/indel) ASR file
    <name>.fasta          MSA file
    <name>.treefile       tree file

The FASTA and treefile come back from their xz blobs byte for byte. The two
.state files are regenerated from the packed columns with the original newline
style and a fresh "Node Site State p_X..." header line; the comment block IQ-TREE
writes above that header is not restored.

Files land in outdir under the original PHACT result-directory layout:
    outdir/1_preprocessed/<name>.fasta
    outdir/2_iqtree_ancestral/<name>.state, <name>.treefile
    outdir/3_binary_iqtree_ancestral/<name>_binary.state

Usage:
    python reborn.py bundle.parquet [outdir]
    python reborn.py bundle.parquet --only state   # state, fasta, tree or all
    python reborn.py bundle.parquet --list         # names and sizes, write nothing
"""

import argparse
import json
import lzma
import os
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

SCALE = 100000


class Bundle:
    """Bundle opened for restoring: metadata up front, columns on demand."""

    def __init__(self, path):
        self._pf = pq.ParquetFile(path)
        md = self._pf.schema_arrow.metadata
        if md is None or b"bundle_node_ids" not in md:
            raise SystemExit(f"error: not an ossuary bundle: {path}")

        self.node_ids = json.loads(md[b"bundle_node_ids"])
        self.n_sites = int(md[b"bundle_n_sites"])
        self.n_nodes = len(self.node_ids)
        self.alphabet = {"nt": tuple(json.loads(md[b"bundle_alphabet_nt"])),
                         "gap": tuple(json.loads(md[b"bundle_alphabet_gap"]))}
        self.exceptions = {
            w: {int(k): v for k, v in
                json.loads(md[f"bundle_exceptions_{w}".encode()]).items()}
            for w in ("nt", "gap")}
        self.has_adj = {w: md[f"bundle_has_sum_adj_{w}".encode()] == b"1"
                        for w in ("nt", "gap")}
        self.newline = {w: ("\r\n" if md.get(f"bundle_newline_{w}".encode()) == b"crlf"
                            else "\n") for w in ("nt", "gap")}
        self.scale = int(md.get(b"bundle_scale", str(SCALE).encode()))
        self.filenames = json.loads(md.get(b"bundle_filenames", b"{}"))
        self._blob_rows = json.loads(md[b"bundle_blob_rows"])
        self._blobs = None

    def blob(self, name):
        """Decompressed blob by name, or None when the bundle predates it."""
        if self._blobs is None:
            col = self._pf.read(columns=["blob"]).column("blob")
            self._blobs = {k: lzma.decompress(col[r].as_py())
                           for k, r in self._blob_rows.items()}
        return self._blobs.get(name)

    def name(self, role, default):
        return self.filenames.get(role) or default

    def values(self, which):
        """Rebuild the (n_rows, k) probability matrix and its argmax column."""
        k = len(self.alphabet[which])
        t = self._pf.read(columns=[c for c in self._pf.schema_arrow.names
                                   if c.endswith(f"_{which}")])
        n = t.num_rows
        dropped = t.column(f"dropped_{which}").to_numpy(zero_copy_only=False)
        rest = np.empty((n, k - 1), dtype=np.int32)
        for c in range(k - 1):
            rest[:, c] = t.column(f"v{c}_{which}").to_numpy(zero_copy_only=False)
        adj = (t.column(f"sum_adj_{which}").to_numpy(zero_copy_only=False).astype(np.int32)
               if self.has_adj[which] else np.zeros(n, dtype=np.int32))
        vals = np.empty((n, k), dtype=np.int32)
        missing = (self.scale + adj) - rest.sum(axis=1)
        for c in range(k):
            m = dropped == c
            if not m.any():
                continue
            vals[np.ix_(m, [j for j in range(k) if j != c])] = rest[m]
            vals[m, c] = missing[m]
        return vals, dropped


def _decimal_table(hi, scale):
    """Map a scaled integer to 5-decimal text, whatever the bundle's scale is."""
    mul = SCALE // scale
    return np.array([f"{(v * mul) // SCALE}.{(v * mul) % SCALE:05d}"
                     for v in range(hi + 1)], dtype=object)


def _state_labels(b, which, vals_dropped):
    """Per-row State letters, with the stored non-alphabet exceptions applied."""
    _, dropped = vals_dropped
    alphabet = np.array(b.alphabet[which], dtype=object)
    labels = alphabet[dropped]
    for row, lab in b.exceptions[which].items():
        labels[row] = lab
    return labels


def write_state(b, which, out_path):
    """Write one .state file in the original node-major order."""
    vals, dropped = b.values(which)
    labels = _state_labels(b, which, (vals, dropped))
    fmt = _decimal_table(int(vals.max()), b.scale)
    nl = b.newline[which]
    header = "Node\tSite\tState\t" + "\t".join(f"p_{a}" for a in b.alphabet[which])

    site_col = np.array([str(s) for s in range(1, b.n_sites + 1)], dtype=object)
    k = vals.shape[1]

    with open(out_path, "wb") as fout:
        fout.write((header + nl).encode())
        for pos, node_id in enumerate(b.node_ids):
            sl = slice(pos, None, b.n_nodes)
            node_col = np.full(b.n_sites, f"Node{node_id}", dtype=object)
            cols = [node_col, site_col, labels[sl]]
            cols += [fmt[vals[sl, c]] for c in range(k)]
            fout.write((nl.join(map("\t".join, zip(*cols))) + nl).encode())


SUBDIRS = {"fasta": "1_preprocessed",
          "nt": "2_iqtree_ancestral", "tree": "2_iqtree_ancestral",
          "gap": "3_binary_iqtree_ancestral"}


def restore(bundle_path, outdir, only="all"):
    """Restore the requested files under outdir, mirroring the PHACT layout."""
    b = Bundle(bundle_path)
    stem = Path(bundle_path).stem
    outdir = Path(outdir)
    dirs = {role: outdir / sub for role, sub in SUBDIRS.items()}
    for d in set(dirs.values()):
        d.mkdir(parents=True, exist_ok=True)
    written = []

    if only in ("all", "state"):
        for which, default in (("nt", f"{stem}.state"),
                               ("gap", f"{stem}_binary.state")):
            path = dirs[which] / b.name(which, default)
            write_state(b, which, path)
            written.append(path)

    if only in ("all", "fasta"):
        path = dirs["fasta"] / b.name("fasta", f"{stem}.fasta")
        path.write_bytes(b.blob("fasta_xz"))
        written.append(path)

    if only in ("all", "tree"):
        path = dirs["tree"] / b.name("tree", f"{stem}.treefile")
        path.write_bytes(b.blob("treefile_xz"))
        written.append(path)

    return [(p, p.stat().st_size) for p in written]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Restore the four original files from an ossuary bundle.")
    ap.add_argument("bundle", metavar="BUNDLE.parquet")
    ap.add_argument("outdir", nargs="?", default=".",
                    help="output directory; the PHACT subfolders are created "
                         "under it (default: current)")
    ap.add_argument("--only", default="all", choices=("all", "state", "fasta", "tree"),
                    help="restore just one kind of file (default: all)")
    ap.add_argument("--list", action="store_true",
                    help="show what the bundle holds and exit")
    args = ap.parse_args(argv)

    if args.list:
        b = Bundle(args.bundle)
        stem = Path(args.bundle).stem
        print(f"nodes x sites : {b.n_nodes} x {b.n_sites}")
        print(f"alphabets     : nt {b.alphabet['nt']}  gap {b.alphabet['gap']}")
        for role, default, w in (("nt", f"{stem}.state", "nt"),
                                 ("gap", f"{stem}_binary.state", "gap")):
            style = "CRLF" if b.newline[w] == "\r\n" else "LF"
            print(f"{role:5s}         : {b.name(role, default)} ({style})")
        for role, default in (("fasta", f"{stem}.fasta"), ("tree", f"{stem}.treefile")):
            key = "fasta_xz" if role == "fasta" else "treefile_xz"
            print(f"{role:5s}         : {b.name(role, default)} "
                  f"({len(b.blob(key))/1e6:.2f} MB)")
        return

    t0 = time.time()
    written = restore(args.bundle, args.outdir, args.only)
    for path, size in written:
        print(f"{path}  ({size/1e6:.2f} MB)", file=sys.stderr)
    print(f"done: {len(written)} files in {time.time() - t0:.1f} s", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        sys.exit(0)
