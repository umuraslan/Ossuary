# Ossuary

Ossuary is a tool that packs the processed ePHACTn inputs into a single Parquet file.

A ePHACTn input set consists of four files:

| File | Role |
| --- | --- |
| `<name>.state` | Nucleotide ancestral-state reconstruction (ASR) |
| `<name>_binary.state` | Binary (gap/indel) ASR |
| `<name>.fasta` | Multiple sequence alignment (MSA) |
| `<name>.treefile` | Newick tree |

`ossuary.py` packs them into one `.parquet` bundle, `ossuary_reader.py` queries that
bundle, and `reborn.py` restores the four original files from it.

## Requirements

- Python 3.10+ (tested on 3.13). The code itself uses nothing newer than 3.8; the
  floor comes from numpy 2.1.3, which requires 3.10.
- [numpy](https://numpy.org/) (tested with 2.1.3)
- [pyarrow](https://arrow.apache.org/docs/python/) (tested with 19.0.0)

```bash
pip install numpy pyarrow
```

These are the only third-party packages. Everything else the scripts import —
`argparse`, `bz2`, `concurrent.futures`, `gzip`, `io`, `json`, `lzma`, `os`,
`pathlib`, `sys`, `tarfile`, `time` — is in the standard library.

Inputs may be `.gz` / `.xz` / `.bz2` compressed, and the result directory may be
given as a tar archive, which is streamed rather than extracted.

## Usage

Pack a input set:

```bash
python ossuary.py chr1_100001-110000.tar                   # -> chr1_100001-110000.parquet
python ossuary.py chr1_100001-110000/ out.parquet                  # a directory works too
python ossuary.py <nt.state> <binary.state> <fasta> <treefile> out.parquet

python ossuary.py chr1_100001-110000.tar --precision 3
# 3 decimals for probabilities (default: 5 = lossless)
```

Options:

| Option | Description |
| --- | --- |
| `--codec {brotli,gzip,lz4,snappy,zstd}` | Compression codec (default: `brotli`) |
| `--level LEVEL` | Compression level (default 9 for brotli, 13 for zstd). brotli 11 is 6% smaller but ~25x slower. |
| `--precision {5,4,3,2}` | Decimal places kept for probabilities (default 5 = lossless). 3 gives a ~64% smaller file, rounding each row so its probabilities still sum to exactly 1. |

Read your parquet file:

```bash
python ossuary_reader.py chr1_100001-110000.parquet            # summary of the bundle
python ossuary_reader.py chr1_100001-110000.parquet <node> <site>
python ossuary_reader.py chr1_100001-110000.parquet --tree     # print the Newick tree
python ossuary_reader.py chr1_100001-110000.parquet --fasta     # print the alignment
```

Reacquire your processed inputs:

```bash
python reborn.py chr1_100001-110000.parquet out/
python reborn.py chr1_100001-110000.parquet --list             # see what is inside
python reborn.py chr1_100001-110000.parquet out/ --only state  # state, fasta, tree or all
```

Restored files land under the original ePHACTn input-directory layout:

```
out/1_preprocessed/<name>.fasta
out/2_iqtree_ancestral/<name>.state, <name>.treefile
out/3_binary_iqtree_ancestral/<name>_binary.state
```
## As a python module

```python
from ossuary_reader import BundleReader
b = BundleReader("chr1_100001-110000.parquet")

b.node_ids          # [2, 22, 25, 26, 24, ...] internal nodes, in file order
b.sites_per_node    # 10000
b.species           # ['hg38', 'panPan3', 'panTro6', ...] leaf names, in FASTA order

b.get(node_id=2, site=11)
# {'nt': ('C', 0.001, 0.997, 0.001, 0.001), 'gap': ('0', 0.99, 0.01)}
# each tuple is (state, p_0, p_1, ...) over b.alphabet_nt / b.alphabet_gap

b.get_nt(node_id=2, site=11)   # ('C', 0.001, 0.997, 0.001, 0.001)
b.get_gap(node_id=2, site=11)  # ('0', 0.99, 0.01)

b.tree()            # Newick text: '(hg38:0.0106668292,(panPan3:0.0069829796,...'
b.fasta()           # {'hg38': 'ACTAAGCACACAGAGAATAATGTCTAGAATCTGAGTGCCA...', ...}
b.column(site=11)   # {'hg38': 'C', 'panPan3': '-', 'panTro6': '-', ...}
```

`get`, `get_nt` and `get_gap` return `None` when the node is unknown or the site is
out of range. Probability columns are read lazily on the first lookup; pass
`BundleReader(path, preload=True)` to load them up front.

The probabilities above come from a bundle written with `--precision 3`; a default
`--precision 5` bundle returns values such as `('C', 0.00083, 0.99751, ...)`.

## Method

Almost all of the bulk of a ePHACTn input set is the two `.state` files, so that is
where the packing happens. Each is read into one integer matrix and then reduced:

1. **Fixed point instead of text.** IQ-TREE writes probabilities with 5 decimals, so
   each one is stored as an integer scaled by `10^precision` (`0.99751` -> `99751`).
   At the default precision this is exact, not lossy. `--precision` rescales with a
   largest-remainder rounding that keeps every row summing to exactly 1.
2. **One column is never stored.** The probabilities of a row sum to 1, so the
   largest of the `k` letters is dropped and recomputed on read as
   `scale - sum(rest)`. Only `k - 1` `uint16` columns survive, plus one `int8`
   column recording which letter was dropped.
3. **The State column comes for free.** The letter dropped in step 2 is the argmax,
   which is exactly what the source file's `State` column holds — so `State` is read
   back as `alphabet[dropped]` instead of being stored. The only disagreements are
   ties, where two letters share the top probability and IQ-TREE picked the one
   `argmax` did not. Those rows go into a small map in the metadata: 90 out of
   4,460,000 in `chr1_100001-110000`.
4. **Rounding is repaired, not assumed away.** `scale - sum(rest)` is only correct if
   the row truly sums to 1, but IQ-TREE rounds each probability on its own, so a row
   can land on `0.99999`. An `int8` `sum_adj` column carries that difference where it
   occurs, and is left out of the file when no row needs it — as in both bundles here.
5. **Site-major reordering.** Rows are transposed from node-major to site-major, so
   all nodes for site 1 come first. Neighbouring nodes at the same site hold nearly
   identical probabilities, which is what makes the columns compress well.
6. **The other two files are kept verbatim** as xz blobs in a `blob` column, and
   everything needed to undo the above (node ids, alphabets, site count, scale,
   original file names and newline style) lives in the Parquet schema metadata.

The nucleotide file, the binary file and the two blobs are processed on three
threads. For `chr1_100001-110000` at `--precision 3`, 349.32 MB of uncompressed input
(shipped as a 43.36 MB tar of gzipped `.state` files) becomes a 1.67 MB bundle.

Restoration is therefore not symmetric. The FASTA and treefile come back from their
xz blobs byte for byte. The two `.state` files are rebuilt from the packed columns
with the original newline style and a fresh `Node Site State p_X...` header line;
the comment block IQ-TREE writes above that header is not kept. File names come from
the bundle metadata, so a source file named `<name>.state.gz` is restored under that
same name even though its contents are written as plain text.