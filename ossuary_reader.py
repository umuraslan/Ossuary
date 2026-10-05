#!/usr/bin/env python3
"""Read a Parquet bundle produced by ossuary.py: nucleotide and gap ancestral-state
probabilities, the alignment (FASTA) and the tree (Newick) from one interface.

Three things the four raw files carried are not stored, because they can be
derived, and are reconstructed here from the row index alone: the Node and Site
columns (rows are site-major, so row `r` is site `r // n_nodes + 1` and node
`node_ids[r % n_nodes]`), the State column (the argmax of the row, with the rows
IQ-TREE labelled otherwise listed as exceptions), and the argmax probability
itself (dropped from storage and recovered as whatever is left of the row's
total, which is what keeps every remaining value inside a uint16).

Probabilities are integers scaled by `bundle_scale`, so the bundle is exact only
to that many decimals. See README.md for the BundleReader API and the command
line; the reading/decoding internals here are kept in sync with the backend's
copy in ePHACT (backend/app/ossuary_reader.py).
"""

import json
import lzma
import os
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

# The positional order results are reported in, checked against the bundle's
# own alphabet at open time rather than assumed to match.
NT_ORDER = ("A", "C", "G", "T")

# The gap table's "state 0" is the gap probability.
GAP_LABEL = "0"


class BundleReader:
    """Lazy reader for one bundle, holding its probability columns in the compact form they are stored in until a site is asked for."""

    def __init__(self, path: str, preload: bool = False):
        self._path = str(path)
        # Errors name the file, never the full path, for tidier CLI/API messages.
        self._label = Path(self._path).name
        self._pf = pq.ParquetFile(self._path)
        md = self._pf.schema_arrow.metadata or {}
        try:
            self.node_ids = json.loads(md[b"bundle_node_ids"])
            self.sites_per_node = int(md[b"bundle_n_sites"])
            self.alphabet_nt = tuple(json.loads(md[b"bundle_alphabet_nt"]))
            self.alphabet_gap = tuple(json.loads(md[b"bundle_alphabet_gap"]))
            self.scale = int(md[b"bundle_scale"])
            self._blob_rows = json.loads(md[b"bundle_blob_rows"])
            self._has_adj = {"nt": md[b"bundle_has_sum_adj_nt"] == b"1",
                             "gap": md[b"bundle_has_sum_adj_gap"] == b"1"}
            # Rows whose State was not the row's argmax, keyed by row index.
            self._exceptions = {
                "nt": {int(k): v for k, v in
                       json.loads(md[b"bundle_exceptions_nt"]).items()},
                "gap": {int(k): v for k, v in
                        json.loads(md[b"bundle_exceptions_gap"]).items()},
            }
        except KeyError as e:
            raise ValueError(f"{self._label}: not a PHACTn Parquet bundle "
                             f"(missing {e.args[0].decode()}).") from None

        if set(self.alphabet_nt) != set(NT_ORDER):
            raise ValueError(f"{self._label}: nucleotide alphabet is "
                             f"{self.alphabet_nt}, expected {NT_ORDER}.")
        if GAP_LABEL not in self.alphabet_gap:
            raise ValueError(f"{self._label}: gap alphabet {self.alphabet_gap} "
                             f"has no '{GAP_LABEL}' state.")

        # The .state files spelled internal nodes "Node2" and the bundle stores
        # the bare number, so the prefix is put back to stay closest to the
        # original files.
        self.node_names = [f"Node{n}" for n in self.node_ids]
        self._n_nodes = len(self.node_ids)
        self._node_pos = {n: i for i, n in enumerate(self.node_ids)}
        self._nt_order = [self.alphabet_nt.index(nt) for nt in NT_ORDER]
        self._gap_col = self.alphabet_gap.index(GAP_LABEL)
        self._compact = None
        self._blobs = {}
        self._sequences = None

        if preload:
            self._load()

    # --- decoding ------------------------------------------------------------

    def _column(self, name):
        """One whole column as a numpy array, read on its own because the bundle is a single row group with 64 MB data pages and decoding every column at once peaks near 500 MB against ~50 MB of actual data."""
        return self._pf.read(columns=[name]).column(name).to_numpy(
            zero_copy_only=False)

    def _load(self):
        """Read every probability column once, keeping it in the compact form it is stored in rather than expanding all rows."""
        if self._compact is None:
            self._compact = {
                which: self._compact_columns(which, len(alphabet))
                for which, alphabet in (("nt", self.alphabet_nt),
                                        ("gap", self.alphabet_gap))
            }
        return self._compact

    def _compact_columns(self, which, k):
        """The (dropped, kept values, sum drift) triple one table is stored as."""
        dropped = self._column(f"dropped_{which}")
        rest = np.empty((len(dropped), k - 1), dtype=np.uint16)
        for c in range(k - 1):
            rest[:, c] = self._column(f"v{c}_{which}")
        adj = self._column(f"sum_adj_{which}") if self._has_adj[which] else None
        return dropped, rest, adj

    def _check_site(self, site):
        """The site as an int, or a ValueError if it is out of range."""
        site = int(site)
        if not 1 <= site <= self.sites_per_node:
            raise ValueError(f"site {site} outside 1..{self.sites_per_node}")
        return site

    def _block(self, which, site):
        """One site's probabilities for every node, as a (n_nodes, k) float array in the bundle's own alphabet order."""
        dropped, rest, adj = self._load()[which]
        k = rest.shape[1] + 1
        lo = (self._check_site(site) - 1) * self._n_nodes
        hi = lo + self._n_nodes

        d = dropped[lo:hi].astype(np.int64)
        r = rest[lo:hi].astype(np.int32)
        total = self.scale + (adj[lo:hi].astype(np.int32)
                              if adj is not None else 0)

        # Put the kept values back in their original columns and recover the
        # dropped one as whatever is left of the row's total.
        vals = np.empty((self._n_nodes, k), dtype=np.int32)
        missing = total - r.sum(axis=1)
        for c in range(k):
            m = d == c
            if not m.any():
                continue
            vals[np.ix_(m, [j for j in range(k) if j != c])] = r[m]
            vals[m, c] = missing[m]
        return vals / self.scale

    def state_at(self, which, site, node_id):
        """The State column's value for one row: the exception label if IQ-TREE recorded one, otherwise the alphabet letter at the dropped column."""
        pos = self._node_pos.get(node_id)
        if pos is None:
            raise ValueError(f"no node {node_id} in this bundle")
        r = (self._check_site(site) - 1) * self._n_nodes + pos
        exception = self._exceptions[which].get(r)
        if exception is not None:
            return exception
        alphabet = self.alphabet_nt if which == "nt" else self.alphabet_gap
        dropped, _rest, _adj = self._load()[which]
        return alphabet[int(dropped[r])]

    # --- the four inputs -------------------------------------------------------

    def _blob(self, name):
        """One of the two embedded files, decompressed once and kept, since both share a column that would otherwise be re-read on every request."""
        if not self._blobs:
            column = self._pf.read(columns=["blob"]).column("blob")
            self._blobs = {k: lzma.decompress(column[row].as_py())
                           for k, row in self._blob_rows.items()}
        return self._blobs[name]

    @property
    def treefile(self) -> bytes:
        """The Newick tree, byte for byte as the .treefile held it."""
        return self._blob("treefile_xz")

    @property
    def fasta(self) -> bytes:
        """The alignment, byte for byte as the _noGapped.fasta held it."""
        return self._blob("fasta_xz")

    def tree(self) -> str:
        """The Newick tree as text."""
        return self.treefile.decode()

    def sequences(self) -> dict:
        """{species_name: aligned_sequence}, parsed from the embedded FASTA and cached."""
        if self._sequences is None:
            seqs, name, cur = {}, None, []
            for line in self.fasta.decode().splitlines():
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(cur)
                    name, cur = line[1:].strip(), []
                else:
                    cur.append(line.strip())
            if name:
                seqs[name] = "".join(cur)
            self._sequences = seqs
        return self._sequences

    @property
    def species(self):
        """Leaf names in FASTA order."""
        return list(self.sequences().keys())

    def column(self, site: int) -> dict:
        """One column of the alignment: {species_name: letter}."""
        return {n: s[site - 1] for n, s in self.sequences().items()}

    # --- per-node/site lookups ---------------------------------------------

    def _record(self, which, node_id, site):
        """(state, p0, p1, ...) for one (node, site) row, in the bundle's own alphabet order."""
        pos = self._node_pos[node_id]
        row = self._block(which, site)[pos]
        state = self.state_at(which, site, node_id)
        return (state, *(float(v) for v in row))

    def get_nt(self, node_id: int, site: int):
        """Nucleotide record for (node, site), or None."""
        if node_id not in self._node_pos or not (1 <= site <= self.sites_per_node):
            return None
        return self._record("nt", node_id, site)

    def get_gap(self, node_id: int, site: int):
        """Gap/indel record for (node, site), or None."""
        if node_id not in self._node_pos or not (1 <= site <= self.sites_per_node):
            return None
        return self._record("gap", node_id, site)

    def get(self, node_id: int, site: int):
        """Both records for (node, site) as {"nt": ..., "gap": ...}, or None."""
        if node_id not in self._node_pos or not (1 <= site <= self.sites_per_node):
            return None
        return {"nt": self._record("nt", node_id, site),
                "gap": self._record("gap", node_id, site)}


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(f"usage: python {sys.argv[0]} <bundle.parquet> "
              f"[<node> <site> | --tree | --fasta]",
              file=sys.stderr if len(sys.argv) < 2 else sys.stdout)
        sys.exit(1 if len(sys.argv) < 2 else 0)

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
        for name, seq in b.sequences().items():
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
