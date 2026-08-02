#!/usr/bin/env python3
"""Read a Parquet bundle produced by ossuary.py: nucleotide and gap ancestral-state
probabilities, the alignment (FASTA) and the tree (Newick) from one interface.

Library usage:
    from ossuary_reader import BundleReader
    b = BundleReader("bundle.parquet")

    b.node_ids                  # internal nodes, in file order
    b.sites_per_node            # 10000
    b.species                   # leaf names, in FASTA order

    b.get(node_id=2, site=11)   # -> {'nt': ('C', 0.00083, 0.99751, ...),
                                #     'gap': ('1', 0.0, 1.0)}
    b.tree()                    # Newick text (str)
    b.fasta()                   # {'hg38': 'ACTAAG...', ...}
    b.column(site=11)           # the leaf letters at that site

CLI usage:
    python ossuary_reader.py bundle.parquet                  # summary
    python ossuary_reader.py bundle.parquet <node> <site>    # one record
    python ossuary_reader.py bundle.parquet --tree           # print Newick
    python ossuary_reader.py bundle.parquet --fasta          # print FASTA
"""

import json
import lzma
import os
import sys

import numpy as np
import pyarrow.parquet as pq

SCALE = 100000


class BundleReader:
    """Lazy reader for an ossuary bundle; probability columns load on first use."""

    def __init__(self, path: str, preload: bool = False):
        self._pf = pq.ParquetFile(path)
        md = self._pf.schema_arrow.metadata

        self.node_ids = json.loads(md[b"bundle_node_ids"])
        self.sites_per_node = int(md[b"bundle_n_sites"])
        self.alphabet_nt = tuple(json.loads(md[b"bundle_alphabet_nt"]))
        self.alphabet_gap = tuple(json.loads(md[b"bundle_alphabet_gap"]))
        self._exc_nt = {int(k): v for k, v in
                        json.loads(md[b"bundle_exceptions_nt"]).items()}
        self._exc_gap = {int(k): v for k, v in
                         json.loads(md[b"bundle_exceptions_gap"]).items()}
        self._has_adj_nt = md[b"bundle_has_sum_adj_nt"] == b"1"
        self._has_adj_gap = md[b"bundle_has_sum_adj_gap"] == b"1"
        self.scale = int(md.get(b"bundle_scale", str(SCALE).encode()))
        self._blob_rows = json.loads(md[b"bundle_blob_rows"])

        self._n_nodes = len(self.node_ids)
        self._node_pos = {n: i for i, n in enumerate(self.node_ids)}
        self._probs = None
        self._fasta = None

        if preload:
            self._load()

    def _load(self):
        """Read the probability columns once and restore the full matrices."""
        if self._probs is not None:
            return
        t = self._pf.read(columns=[c for c in self._pf.schema_arrow.names
                                   if c != "blob"])
        self._probs = {
            "nt": self._rebuild(t, "nt", len(self.alphabet_nt), self._has_adj_nt,
                                self.scale),
            "gap": self._rebuild(t, "gap", len(self.alphabet_gap), self._has_adj_gap,
                                 self.scale),
        }

    @staticmethod
    def _rebuild(t, suffix, k, has_adj, scale):
        """Invert the packing: reinsert the dropped column from the row sum."""
        n = t.num_rows
        dropped = t.column(f"dropped_{suffix}").to_numpy(zero_copy_only=False)
        rest = np.empty((n, k - 1), dtype=np.int32)
        for c in range(k - 1):
            rest[:, c] = t.column(f"v{c}_{suffix}").to_numpy(zero_copy_only=False)
        adj = (t.column(f"sum_adj_{suffix}").to_numpy(zero_copy_only=False).astype(np.int32)
               if has_adj else np.zeros(n, dtype=np.int32))
        vals = np.empty((n, k), dtype=np.int32)
        missing = (scale + adj) - rest.sum(axis=1)
        for c in range(k):
            m = dropped == c
            if not m.any():
                continue
            vals[np.ix_(m, [j for j in range(k) if j != c])] = rest[m]
            vals[m, c] = missing[m]
        return vals, dropped.astype(np.int64)

    def _row_index(self, node_id, site):
        """Map (node, site) to a site-major row index, or None if out of range."""
        pos = self._node_pos.get(node_id)
        if pos is None or not (1 <= site <= self.sites_per_node):
            return None
        return (site - 1) * self._n_nodes + pos

    def _record(self, which, r):
        """Return (state, p0, p1, ...) for one row of the nt or gap table."""
        self._load()
        vals, dropped = self._probs[which]
        alphabet = self.alphabet_nt if which == "nt" else self.alphabet_gap
        exc = self._exc_nt if which == "nt" else self._exc_gap
        state = exc.get(r) or alphabet[dropped[r]]
        return (state, *(float(v) / self.scale for v in vals[r]))

    def _blob(self, name):
        """Read and decompress one xz blob by name."""
        row = self._blob_rows[name]
        t = self._pf.read(columns=["blob"])
        return lzma.decompress(t.column("blob")[row].as_py())

    def get_nt(self, node_id: int, site: int):
        """Nucleotide record for (node, site), or None."""
        r = self._row_index(node_id, site)
        return None if r is None else self._record("nt", r)

    def get_gap(self, node_id: int, site: int):
        """Gap/indel record for (node, site), or None."""
        r = self._row_index(node_id, site)
        return None if r is None else self._record("gap", r)

    def get(self, node_id: int, site: int):
        """Both records for (node, site) as {"nt": ..., "gap": ...}, or None."""
        r = self._row_index(node_id, site)
        if r is None:
            return None
        return {"nt": self._record("nt", r), "gap": self._record("gap", r)}

    def tree(self) -> str:
        """The Newick tree as text."""
        return self._blob("treefile_xz").decode()

    def fasta(self) -> dict:
        """{species_name: aligned_sequence}, decompressed once and cached."""
        if self._fasta is None:
            seqs, name, cur = {}, None, []
            for line in self._blob("fasta_xz").decode().splitlines():
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(cur)
                    name, cur = line[1:].strip(), []
                else:
                    cur.append(line.strip())
            if name:
                seqs[name] = "".join(cur)
            self._fasta = seqs
        return self._fasta

    @property
    def species(self):
        """Leaf names in FASTA order."""
        return list(self.fasta().keys())

    def column(self, site: int) -> dict:
        """One column of the alignment: {species_name: letter}."""
        return {n: s[site - 1] for n, s in self.fasta().items()}


def main():
    if len(sys.argv) < 2:
        print(f"usage: python {sys.argv[0]} <bundle.parquet> "
              f"[<node> <site> | --tree | --fasta]", file=sys.stderr)
        sys.exit(1)

    b = BundleReader(sys.argv[1])
    args = sys.argv[2:]

    if not args:
        print(f"file            : {sys.argv[1]} "
              f"({os.path.getsize(sys.argv[1])/1e6:.2f} MB)")
        print(f"internal nodes  : {len(b.node_ids)}")
        print(f"sites per node  : {b.sites_per_node}")
        print(f"leaves (species): {len(b.species)}")
        print(f"nt alphabet     : {b.alphabet_nt}")
        print(f"gap alphabet    : {b.alphabet_gap}")
        print(f"tree            : {len(b.tree())} characters of Newick")
    elif args[0] == "--tree":
        print(b.tree())
    elif args[0] == "--fasta":
        for name, seq in b.fasta().items():
            print(f">{name}")
            for i in range(0, len(seq), 60):
                print(seq[i:i + 60])
    else:
        node_id, site = int(args[0]), int(args[1])
        rec = b.get(node_id, site)
        if rec is None:
            print(f"node {node_id} / site {site} not found", file=sys.stderr)
            sys.exit(1)
        for which in ("nt", "gap"):
            state, *vals = rec[which]
            print(f"{which:4s} {node_id}\t{site}\t{state}\t" +
                  "\t".join(str(v) for v in vals))


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        sys.exit(0)
