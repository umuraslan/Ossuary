#!/usr/bin/env python3
"""Pack an ePHACTn input set into a single Parquet bundle.

The four input files may be given in any order; roles are detected from the file
names, the same way an input directory is scanned. See README.md for usage, the
command-line options and the packing method.
"""

import argparse
import bz2
import gzip
import io
import json
import lzma
import os
import sys
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path, PurePosixPath

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pv
import pyarrow.parquet as pq

SCALE = 100000

FASTA_ROW = 0
TREE_ROW = 1

JOBS = 3

STATE_SUFFIXES = (".state", ".state.gz", ".state.xz", ".state.bz2")
FASTA_SUFFIXES = (".fasta", ".fa", ".fasta.gz", ".fa.gz")
TAR_SUFFIXES = (".tar", ".tar.gz", ".tgz", ".tar.bz2", ".tbz2", ".tar.xz", ".txz")

_STREAM_DECOMPRESSORS = {".gz": lambda f: gzip.GzipFile(fileobj=f),
                         ".xz": lzma.LZMAFile, ".bz2": bz2.BZ2File}
_ARROW_CODECS = {".gz": "gzip", ".bz2": "bz2"}


class _Source:
    """An input to read: a filesystem path or a member of a tar archive."""

    __slots__ = ("name", "tar_path", "member")

    def __init__(self, name, tar_path=None, member=None):
        self.name = str(name)
        self.tar_path = tar_path
        self.member = member

    def __str__(self):
        return self.name if self.tar_path is None else f"{self.tar_path}::{self.member}"

    @property
    def suffix(self):
        return PurePosixPath(self.name).suffix

    def open_raw(self):
        """Raw, possibly compressed byte stream; tar members get their own handle."""
        if self.tar_path is None:
            return open(self.name, "rb")
        tar = tarfile.open(self.tar_path, "r:*")
        stream = tar.extractfile(self.member)
        if stream is None:
            raise SystemExit(f"error: cannot read tar member: {self}")
        stream._tar = tar  # GC of the TarFile would close the stream
        return stream

    def open_binary(self):
        """Decompressed byte stream."""
        raw = self.open_raw()
        dec = _STREAM_DECOMPRESSORS.get(self.suffix)
        return dec(raw) if dec else raw

    def open_text(self):
        return io.TextIOWrapper(self.open_binary(), encoding="utf-8", newline="")

    def read_bytes(self) -> bytes:
        with self.open_binary() as fin:
            return fin.read()

    def csv_input(self):
        """Input for pyarrow.csv: a path when possible, else a stream Arrow can decode."""
        if self.tar_path is None:
            return self.name
        codec = _ARROW_CODECS.get(self.suffix)
        return (pa.CompressedInputStream(self.open_raw(), codec) if codec
                else self.open_binary())


def _as_source(obj) -> _Source:
    return obj if isinstance(obj, _Source) else _Source(obj)


def read_bytes(path) -> bytes:
    """Read a file as raw bytes, decompressing by suffix when needed."""
    return _as_source(path).read_bytes()


def open_text(path):
    """Open a file as a text stream; .gz/.xz/.bz2 are handled transparently."""
    return _as_source(path).open_text()


def _detect_newline(src):
    """Detect a .state file's newline style from its first line."""
    with src.open_binary() as fin:
        first = fin.readline()
    return b"\r\n" if first.endswith(b"\r\n") else b"\n"


def parse_exact_5dec(raw: str) -> int:
    """Parse a 5-decimal fixed-point value exactly: "0.99999" -> 99999."""
    if raw == "":
        return 0
    if "." in raw:
        intpart, decpart = raw.split(".")
    else:
        intpart, decpart = raw, ""
    decpart = (decpart + "00000")[:5]
    return int(intpart) * 100000 + int(decpart)


def strip_node_prefix(raw: str) -> str:
    return raw[len("Node"):] if raw.startswith("Node") else raw


def _detect_layout(cols, prob_columns, has_site_col):
    """Derive (prob_columns, has_site_col, is_header) from a header or first data line."""
    if cols[0] == "Node" and cols[1] in ("Site", "State"):
        has_site_col = cols[1] == "Site"
        p_cols = cols[3:] if has_site_col else cols[2:]
        prob_columns = tuple(c[2:] if c.startswith("p_") else c for c in p_cols)
        return prob_columns, has_site_col, True
    if prob_columns is None:
        has_site_col = False
        num_prob_cols = len(cols) - 2
        prob_columns = (
            tuple(str(i) for i in range(num_prob_cols))
            if num_prob_cols != 4
            else ("A", "C", "G", "T")
        )
    return prob_columns, has_site_col, False


def _iter_data_lines(in_path):
    """Yield (cols, prob_columns, has_site_col) for each data line, skipping comments."""
    with _as_source(in_path).open_text() as fin:
        prob_columns = None
        has_site_col = None
        for line in fin:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            prob_columns, has_site_col, is_header = _detect_layout(
                cols, prob_columns, has_site_col)
            if is_header:
                continue
            yield cols, prob_columns, has_site_col


def _encode_labels(labels):
    """Turn a label list into (int16 codes, vocabulary)."""
    vocab = {}
    codes = np.empty(len(labels), dtype=np.int16)
    for i, lab in enumerate(labels):
        code = vocab.get(lab)
        if code is None:
            code = vocab[lab] = len(vocab)
        codes[i] = code
    return codes, tuple(vocab)


def _read_state_slow(in_path):
    """Pure-Python fallback parser; handles both raw and cleaned .state layouts."""
    values = []
    state_labels = []
    node_ids = []
    prob_columns = None
    prev_node = None

    for cols, prob_columns, has_site_col in _iter_data_lines(in_path):
        k = len(prob_columns)
        if has_site_col:
            node_raw, state, vals = cols[0], cols[2], cols[3:3 + k]
        else:
            node_raw, state, vals = cols[0], cols[1], cols[2:2 + k]

        node = strip_node_prefix(node_raw) if node_raw != "" else prev_node
        if node != prev_node:
            node_ids.append(int(node))
            prev_node = node

        state_labels.append(state)
        values.append([parse_exact_5dec(v) for v in vals])

    if prob_columns is None:
        raise ValueError("no data rows found in file")

    codes, vocab = _encode_labels(state_labels)
    return np.array(values, dtype=np.int32), codes, vocab, node_ids, tuple(prob_columns)


def _read_state_fast(in_path):
    """Vectorised parser for the raw IQ-TREE layout; None if it does not apply."""
    src = _as_source(in_path)
    n_comment = 0
    with src.open_text() as fin:
        for line in fin:
            if line.startswith("#"):
                n_comment += 1
                continue
            header = line.rstrip("\r\n").split("\t")
            break
        else:
            return None

    if len(header) < 4 or header[:3] != ["Node", "Site", "State"]:
        return None
    prob_columns = tuple(c[2:] if c.startswith("p_") else c for c in header[3:])

    table = pv.read_csv(
        src.csv_input(),
        read_options=pv.ReadOptions(skip_rows=n_comment, block_size=16 << 20),
        parse_options=pv.ParseOptions(delimiter="\t"),
        convert_options=pv.ConvertOptions(
            column_types={"Node": pa.string(), "Site": pa.int32(),
                          "State": pa.string(),
                          **{c: pa.float64() for c in header[3:]}}),
    )
    n_rows = table.num_rows

    sites = table.column("Site").combine_chunks().to_numpy(zero_copy_only=False)
    n_sites = int(sites.max())
    if n_rows % n_sites:
        raise ValueError(f"{n_rows} rows do not divide into {n_sites} sites")
    n_nodes = n_rows // n_sites
    if not np.array_equal(sites, np.tile(np.arange(1, n_sites + 1, dtype=sites.dtype),
                                         n_nodes)):
        return None

    nodes_col = table.column("Node").combine_chunks()
    heads = pc.take(nodes_col, pa.array(np.arange(n_nodes) * n_sites)).to_pylist()
    node_ids = [int(strip_node_prefix(n)) for n in heads]

    states = pc.dictionary_encode(table.column("State").combine_chunks())
    codes = states.indices.to_numpy(zero_copy_only=False).astype(np.int16)
    vocab = tuple(states.dictionary.to_pylist())

    values = np.empty((n_rows, len(prob_columns)), dtype=np.int32)
    scratch = np.empty(n_rows, dtype=np.float64)
    for i, col in enumerate(header[3:]):
        f = table.column(col).combine_chunks().to_numpy(zero_copy_only=False)
        np.multiply(f, SCALE, out=scratch)
        np.rint(scratch, out=scratch)
        values[:, i] = scratch
    return values, codes, vocab, node_ids, prob_columns


def _read_state_file(in_path):
    """Read a .state file into (values, codes, vocab, node_ids, prob_columns)."""
    fast = _read_state_fast(in_path)
    return fast if fast is not None else _read_state_slow(in_path)


def quantize(values, out_scale):
    """Round to out_scale by largest remainder, keeping each row's sum exact."""
    scaled = values / (SCALE / out_scale)
    floor = np.floor(scaled).astype(np.int32)
    deficit = out_scale - floor.sum(axis=1)
    rank = np.argsort(np.argsort(floor - scaled, axis=1), axis=1)
    return floor + (rank < deficit[:, None]).astype(np.int32)


def build_arrays(in_path: str, out_scale: int = SCALE):
    """Pack a .state file into (arrays, metadata, n_nodes, n_rows)."""
    newline = _detect_newline(_as_source(in_path))
    values, codes, vocab, node_ids, prob_columns = _read_state_file(in_path)
    if out_scale != SCALE:
        values = quantize(values, out_scale)

    n_rows, k = values.shape
    n_nodes = len(node_ids)
    if n_rows % n_nodes:
        raise ValueError(f"{n_rows} rows do not divide into {n_nodes} nodes")
    n_sites = n_rows // n_nodes

    sm = values.reshape(n_nodes, n_sites, k).transpose(1, 0, 2).reshape(-1, k).copy()
    codes_sm = codes.reshape(n_nodes, n_sites).T.reshape(-1)

    dropped = np.argmax(sm, axis=1).astype(np.int8)
    sum_adj = (sm.sum(axis=1) - out_scale).astype(np.int8)

    rest = np.empty((n_rows, k - 1), dtype=np.int32)
    for c in range(k):
        m = dropped == c
        if m.any():
            rest[m] = sm[np.ix_(m, [j for j in range(k) if j != c])]
    if rest.size and rest.max() > np.iinfo(np.uint16).max:
        raise ValueError("remaining value does not fit in uint16")
    rest = rest.astype(np.uint16)

    lookup = {lab: i for i, lab in enumerate(vocab)}
    expected = np.array([lookup.get(p, -1) for p in prob_columns], dtype=np.int16)
    exceptions = {
        str(i): vocab[codes_sm[i]]
        for i in np.flatnonzero(codes_sm != expected[dropped]).tolist()
    }

    arrays = {"dropped": pa.array(dropped, type=pa.int8())}
    for c in range(k - 1):
        arrays[f"v{c}"] = pa.array(np.ascontiguousarray(rest[:, c]), type=pa.uint16())
    has_sum_adj = bool((sum_adj != 0).any())
    if has_sum_adj:
        arrays["sum_adj"] = pa.array(sum_adj, type=pa.int8())

    metadata = {
        b"stb_alphabet": json.dumps(list(prob_columns)).encode(),
        b"stb_node_ids": json.dumps(node_ids).encode(),
        b"stb_n_sites": str(n_sites).encode(),
        b"stb_ordering": b"site-major",
        b"stb_state_exceptions": json.dumps(exceptions).encode(),
        b"stb_has_sum_adj": (b"1" if has_sum_adj else b"0"),
        b"stb_newline": (b"crlf" if newline == b"\r\n" else b"lf"),
    }
    return arrays, metadata, n_nodes, n_rows


def _one(role, hits):
    """Return the single match for a role, or exit with a readable error."""
    if not hits:
        raise SystemExit(f"error: no {role} file found")
    if len(hits) > 1:
        listed = "\n  ".join(str(src) for _, src in hits)
        raise SystemExit(f"error: multiple {role} files found, "
                         f"pass paths explicitly:\n  {listed}")
    return hits[0]


def _select(entries):
    """Pick (nt, gap, fasta, tree) from [(relative_path, _Source)] by path shape."""
    entries = sorted(entries, key=lambda e: e[0])
    states = [e for e in entries if e[0].endswith(STATE_SUFFIXES)]
    nt_rel, nt = _one("nucleotide .state",
                      [e for e in states if "binary" not in e[0].lower()])
    _, gap = _one("binary .state", [e for e in states if "binary" in e[0].lower()])
    _, fasta = _one("alignment .fasta",
                    [e for e in entries if e[0].endswith(FASTA_SUFFIXES)])

    trees = [e for e in entries if e[0].endswith(".treefile")]
    parent = PurePosixPath(nt_rel).parent
    beside = [e for e in trees if PurePosixPath(e[0]).parent == parent]
    _, tree = _one("tree .treefile", beside or trees)

    return nt, gap, fasta, tree


def discover(root):
    """Locate the four inputs in a PHACT result directory."""
    root = Path(root)
    if not root.is_dir():
        raise SystemExit(f"error: not a directory: {root}")
    entries = [(p.relative_to(root).as_posix(), _Source(p))
               for p in root.rglob("*") if p.is_file()]
    return _select(entries)


def _strip_common_root(names):
    """Drop the shared top-level directory when every member sits under it."""
    first = names[0].split("/")[0]
    if all("/" in n and n.split("/")[0] == first for n in names):
        return [n.split("/", 1)[1] for n in names]
    return list(names)


def discover_tar(tar_path):
    """Locate the four inputs in a tar archive, streaming members individually."""
    with tarfile.open(tar_path, "r:*") as tar:
        names = [m.name for m in tar.getmembers() if m.isfile()]
    if not names:
        raise SystemExit(f"error: tar archive contains no files: {tar_path}")
    entries = [(rel, _Source(name, tar_path=tar_path, member=name))
               for name, rel in zip(names, _strip_common_root(names))]
    return _select(entries)


def _write_settings(names, compression, compression_level):
    """Encoding and compression settings for pq.write_table (see OZET.md 2.6)."""
    wide = [n for n in names if n.startswith("v")]
    small = [n for n in names if not n.startswith("v")]
    settings = {
        "use_dictionary": small,
        "column_encoding": {n: "BYTE_STREAM_SPLIT" for n in wide},
        "compression": {**{n: compression for n in names}, "blob": "none"},
        "data_page_size": 64 * 1024 * 1024,
    }
    if compression_level is not None:
        settings["compression_level"] = {n: compression_level for n in names}
    return settings


def build(nt_path, gap_path, fasta_path, tree_path, out_path,
          compression="brotli", compression_level=9, out_scale=SCALE):
    """Write the four inputs to out_path as a single Parquet bundle."""
    nt_path, gap_path = _as_source(nt_path), _as_source(gap_path)
    fasta_path, tree_path = _as_source(fasta_path), _as_source(tree_path)

    def _blobs():
        return (lzma.compress(fasta_path.read_bytes(), preset=9 | lzma.PRESET_EXTREME),
                lzma.compress(tree_path.read_bytes(), preset=9))

    # Three independent jobs; numpy, pyarrow.csv and lzma all release the GIL.
    with ThreadPoolExecutor(max_workers=JOBS) as pool:
        fut_nt = pool.submit(build_arrays, nt_path, out_scale)
        fut_gap = pool.submit(build_arrays, gap_path, out_scale)
        fut_blob = pool.submit(_blobs)
        nt_arr, nt_md, n_nodes, n_rows = fut_nt.result()
        gap_arr, gap_md, g_nodes, g_rows = fut_gap.result()
        fasta_xz, tree_xz = fut_blob.result()

    if nt_md[b"stb_node_ids"] != gap_md[b"stb_node_ids"]:
        raise ValueError("the two .state files have different node order")
    if nt_md[b"stb_n_sites"] != gap_md[b"stb_n_sites"]:
        raise ValueError("the two .state files have different site counts")
    if n_rows != g_rows:
        raise ValueError(f"row counts differ: {n_rows} vs {g_rows}")

    cols = {f"{k}_nt": v for k, v in nt_arr.items()}
    cols.update({f"{k}_gap": v for k, v in gap_arr.items()})
    names = list(cols)

    blob = [None] * n_rows
    blob[FASTA_ROW] = fasta_xz
    blob[TREE_ROW] = tree_xz
    cols["blob"] = pa.array(blob, type=pa.binary())

    filenames = {"nt": PurePosixPath(nt_path.name).name,
                 "gap": PurePosixPath(gap_path.name).name,
                 "fasta": PurePosixPath(fasta_path.name).name,
                 "tree": PurePosixPath(tree_path.name).name}

    metadata = {
        b"bundle_node_ids": nt_md[b"stb_node_ids"],
        b"bundle_n_sites": nt_md[b"stb_n_sites"],
        b"bundle_ordering": b"site-major",
        b"bundle_alphabet_nt": nt_md[b"stb_alphabet"],
        b"bundle_alphabet_gap": gap_md[b"stb_alphabet"],
        b"bundle_exceptions_nt": nt_md[b"stb_state_exceptions"],
        b"bundle_exceptions_gap": gap_md[b"stb_state_exceptions"],
        b"bundle_has_sum_adj_nt": nt_md[b"stb_has_sum_adj"],
        b"bundle_has_sum_adj_gap": gap_md[b"stb_has_sum_adj"],
        b"bundle_newline_nt": nt_md[b"stb_newline"],
        b"bundle_newline_gap": gap_md[b"stb_newline"],
        b"bundle_scale": str(out_scale).encode(),
        b"bundle_filenames": json.dumps(filenames).encode(),
        b"bundle_blob_rows": json.dumps(
            {"fasta_xz": FASTA_ROW, "treefile_xz": TREE_ROW}).encode(),
    }

    table = pa.table(cols).replace_schema_metadata(metadata)
    pq.write_table(
        table, out_path,
        row_group_size=n_rows,
        **_write_settings(names, compression, compression_level),
    )


DEFAULT_LEVEL = {"brotli": 9, "zstd": 13, "gzip": 9, "snappy": None, "lz4": None}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Pack PHACT results into a single Parquet bundle.")
    ap.add_argument("inputs", nargs="+", metavar="INPUT",
                    help="either <result_dir|result.tar> [out.parquet], or the "
                         "four input files in any order followed by <out.parquet>")
    ap.add_argument("--outdir", default=None,
                    help="directory to write the output file into "
                         "(default: alongside the given output path, or the "
                         "current directory)")
    ap.add_argument("--codec", default="brotli",
                    choices=sorted(DEFAULT_LEVEL), help="default: brotli")
    ap.add_argument("--level", type=int, default=None,
                    help="compression level (default 9 for brotli, 13 for zstd). "
                         "brotli 11 is 6%% smaller but ~25x slower.")
    ap.add_argument("--precision", type=int, default=5, choices=(5, 4, 3, 2),
                    help="decimal places kept for probabilities (default 5 = "
                         "lossless). 3 gives a ~64%% smaller file, rounding each "
                         "row so its probabilities still sum to exactly 1.")
    args = ap.parse_args(argv)

    first = args.inputs[0]
    is_tar = first.endswith(TAR_SUFFIXES) and os.path.isfile(first)
    if os.path.isdir(first) or is_tar:
        if len(args.inputs) > 2:
            ap.error("with a directory or tar, at most one output path is allowed")
        if is_tar:
            nt, gap, fasta, tree = discover_tar(first)
            base = os.path.basename(first)
            for suf in TAR_SUFFIXES:
                if base.endswith(suf):
                    base = base[:-len(suf)]
                    break
        else:
            root = first.rstrip("/\\")
            nt, gap, fasta, tree = discover(root)
            base = os.path.basename(root)
        out = args.inputs[1] if len(args.inputs) == 2 else f"{base}.parquet"
    elif len(args.inputs) == 5:
        *paths, out = args.inputs
        missing = [p for p in paths if not os.path.isfile(p)]
        if missing:
            ap.error("no such file: " + ", ".join(missing))
        nt, gap, fasta, tree = _select([(os.path.basename(p), _Source(p))
                                        for p in paths])
    else:
        ap.error("pass either a result directory/tar, or exactly 5 paths")

    if args.outdir is not None:
        os.makedirs(args.outdir, exist_ok=True)
        out = os.path.join(args.outdir, os.path.basename(out))

    print(f"nt     : {nt}\ngap    : {gap}\nfasta  : {fasta}\n"
          f"tree   : {tree}\noutput : {out}", file=sys.stderr)

    level = args.level if args.level is not None else DEFAULT_LEVEL[args.codec]
    codec = args.codec + (str(level) if level is not None else "")

    slow = ((args.codec == "brotli" and level >= 10)
            or (args.codec == "zstd" and level >= 18))
    print(f"compression: {codec}" + (" (~25 s for a 10k-site block; the default "
          "level is ~25x faster and 6% larger)" if slow else ""),
          file=sys.stderr, flush=True)
    if args.precision != 5:
        print(f"precision  : {args.precision} decimals (lossy)", file=sys.stderr)

    t0 = time.time()
    build(nt, gap, fasta, tree, out, compression=args.codec,
          compression_level=level, out_scale=10 ** args.precision)
    print(f"done: {out} ({os.path.getsize(out)/1e6:.2f} MB, "
          f"{time.time() - t0:.1f} s)", file=sys.stderr)


if __name__ == "__main__":
    main()
