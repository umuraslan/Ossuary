#!/usr/bin/env python3
"""
Bir PHACTn/IQ-TREE analizinin DORT ciktisini tek bir kayipsiz Parquet
dosyasinda birlestirir:

    <ad>.state          nukleotid ata-durumu olasiliklari  (ic node x site)
    <ad>.state (binary) gap/indel ata-durumu olasiliklari  (ayni node x site)
    <ad>.fasta          yapraklarin hizalanmis dizileri    (ayni site sayisi)
    <ad>.treefile       Newick agaci

Bunlar tek bir veri kumesinin yuzleri: FASTA basliklari agacin yapraklariyla
birebir ayni, iki .state dosyasi ayni node sirasini ve ayni site sayisini
paylasiyor. Bu yuzden iki olasilik tablosu (node, site) anahtari uzerinde
DOGAL bir birlestirmeyle yan yana konur; node listesi/site duzeni tek kez
saklanir.

FASTA ve treefile ise TABLOYA SUTUN OLARAK ACILMAZ, xz'lenmis birer blob
olarak `blob` adli BINARY sutunda tasinir. Sebebi olculdu:
  - FASTA'yi site-major sutunlu duzene acmak sikismayi 149 KB'dan 300 KB'a
    cikariyor (turler arasi uzun LZ eslesmeleri dagiliyor).
  - Ayni blob'u parquet key-value metadata'sina koymak ise UTF-8 genislemesi
    yuzunden 156 KB'lik yuku 365 KB'a sisiriyor; BINARY sutunda sisme yok.

Sema (satirlar site-major, n_nodes*n_sites satir):
    dropped_nt, v0_nt, v1_nt, v2_nt, sum_adj_nt   nukleotid (bkz. build_arrays)
    dropped_gap, v0_gap                            gap/indel
    blob   BINARY  satir0 = fasta.xz, satir1 = treefile.xz, digerleri NULL

Bu dosya kendi kendine yeter: .state okuma ve kompakt sutunlari kurma isi
(state_to_parquet.py / state_to_binary.py'dekiyle ayni mantik) asagida
bulunuyor, disaridan sadece numpy + pyarrow gerekiyor.

Girdiler .gz/.xz/.bz2 sikistirilmis olabilir; acmaya gerek yok. Klasorun
kendisi bir tar arsivi olarak da verilebilir; arsiv DISARI ACILMAZ, uyeler
dogrudan akis olarak okunur.

Kullanim:
    # PHACT sonuc klasoru (1_preprocessed / 2_iqtree_ancestral /
    # 3_binary_iqtree_ancestral) verilir, dosyalar otomatik bulunur:
    python make_bundle.py chr2_10001-20000/ [cikti.parquet]

    # ayni klasorun tar'i (.tar/.tar.gz/.tgz/.tar.bz2/.tar.xz):
    python make_bundle.py chr1_100001-110000.tar [cikti.parquet]

    # ya da dort dosya tek tek:
    python make_bundle.py <nt.state> <binary.state> <fasta> <treefile> <cikti.parquet>

    # sikistirma secenekleri (hiz/boyut dengesi icin):
    python make_bundle.py chr2_10001-20000/ --codec zstd --level 19
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

STATE_SUFFIXES = (".state", ".state.gz", ".state.xz", ".state.bz2")
FASTA_SUFFIXES = (".fasta", ".fa", ".fasta.gz", ".fa.gz")
TAR_SUFFIXES = (".tar", ".tar.gz", ".tgz", ".tar.bz2", ".tbz2", ".tar.xz", ".txz")

_STREAM_DECOMPRESSORS = {".gz": lambda f: gzip.GzipFile(fileobj=f),
                         ".xz": lzma.LZMAFile, ".bz2": bz2.BZ2File}
_ARROW_CODECS = {".gz": "gzip", ".bz2": "bz2"}  # Arrow'da xz codec'i yok


class _Source:
    """Okunacak tek bir girdi: dosya yolu ya da bir tar arsivinin uyesi.

    Ikisi de ayni arayuzden okunur. Tek fark tar uyesinin her acilista KENDI
    TarFile tutamagini almasi: tarfile is parcaciklari arasinda paylasilamaz,
    build() ise uc okumayi paralel yapiyor.
    """

    __slots__ = ("name", "tar_path", "member")

    def __init__(self, name, tar_path=None, member=None):
        self.name = str(name)  # uzanti kontrolleri hep bunun uzerinden
        self.tar_path = tar_path
        self.member = member

    def __str__(self):
        return self.name if self.tar_path is None else f"{self.tar_path}::{self.member}"

    @property
    def suffix(self):
        return PurePosixPath(self.name).suffix

    def open_raw(self):
        """Ham (belki sikistirilmis) bayt akisi."""
        if self.tar_path is None:
            return open(self.name, "rb")
        tar = tarfile.open(self.tar_path, "r:*")
        stream = tar.extractfile(self.member)
        if stream is None:
            raise SystemExit(f"hata: tar uyesi okunamadi: {self}")
        stream._tar = tar  # tar cop toplanirsa akis kapanir; referansi canli tut
        return stream

    def open_binary(self):
        """Cozulmus bayt akisi."""
        raw = self.open_raw()
        dec = _STREAM_DECOMPRESSORS.get(self.suffix)
        return dec(raw) if dec else raw

    def open_text(self):
        return io.TextIOWrapper(self.open_binary(), encoding="utf-8", newline="")

    def read_bytes(self) -> bytes:
        with self.open_binary() as fin:
            return fin.read()

    def csv_input(self):
        """pyarrow.csv'ye verilecek girdi.

        Dosya yolunda YOLUN KENDISI verilir: Arrow dosyayi kendi acip bloklari
        paralel isliyor. Tar uyesinde bu mumkun degil, akis verilir; cozme isini
        yine Arrow'un C++ codec'i yapsin diye CompressedInputStream'e sarilir
        (Python'un gzip modulunden belirgin hizli).
        """
        if self.tar_path is None:
            return self.name
        codec = _ARROW_CODECS.get(self.suffix)
        return (pa.CompressedInputStream(self.open_raw(), codec) if codec
                else self.open_binary())


def _as_source(obj) -> _Source:
    return obj if isinstance(obj, _Source) else _Source(obj)


def read_bytes(path) -> bytes:
    """Dosyayi ham bayt olarak okur; uzantisi sikistirmaliysa acar."""
    return _as_source(path).read_bytes()


def open_text(path):
    """Dosyayi metin akisi olarak acar; .gz/.xz/.bz2 seffaf sekilde cozulur."""
    return _as_source(path).open_text()


# --------------------------------------------------------------------------
# .state okuma
#
# Iki girdi sekli desteklenir:
#   - Ham (IQ-TREE ciktisi): "# ..." yorumlari, bir header satiri
#     (Node  Site  State  p_X ...), Node her satirda yazili, Site sutunu var.
#   - Temiz (clean_state.py ciktisi): yorum/header yok, Site sutunu yok, Node
#     sadece blok basinda yazili, olasilik 0 ise bos, 1 ise sadece "1".
# --------------------------------------------------------------------------


def parse_exact_5dec(raw: str) -> int:
    # "0.99999" -> 99999 ; "0.5" -> 50000 ; hicbir float() kullanmadan, tam ondalikli
    if raw == "":
        return 0
    if "." in raw:
        intpart, decpart = raw.split(".")
    else:
        intpart, decpart = raw, ""
    decpart = (decpart + "00000")[:5]
    return int(intpart) * 100000 + int(decpart)


def strip_node_prefix(raw: str) -> str:
    if raw.startswith("Node"):
        return raw[len("Node"):]
    return raw


def _detect_layout(cols, prob_columns, has_site_col):
    """Header satirindan ya da (yoksa) ilk veri satirindan p_X sutunlarini cikarir."""
    if cols[0] == "Node" and cols[1] in ("Site", "State"):
        has_site_col = cols[1] == "Site"
        p_cols = cols[3:] if has_site_col else cols[2:]
        prob_columns = tuple(c[2:] if c.startswith("p_") else c for c in p_cols)
        return prob_columns, has_site_col, True  # True = bu satir header, atlanmali
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
    """Etiket listesini (kodlar, sozluk) haline getirir."""
    vocab = {}
    codes = np.empty(len(labels), dtype=np.int16)
    for i, lab in enumerate(labels):
        code = vocab.get(lab)
        if code is None:
            code = vocab[lab] = len(vocab)
        codes[i] = code
    return codes, tuple(vocab)


def _read_state_slow(in_path):
    """Satir satir, saf Python okuma. Her girdi seklini kaldirir."""
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
        raise ValueError("Dosyada veri satiri bulunamadi")

    codes, vocab = _encode_labels(state_labels)
    return np.array(values, dtype=np.int32), codes, vocab, node_ids, tuple(prob_columns)


def _read_state_fast(in_path):
    """Ham IQ-TREE duzenini pyarrow.csv ile vektorize okur.

    Ham .state dosyasi CSV motorunun sevdigi bicimde: sabit yorum blogu, tek
    header satiri, sekme ayrilmis, her satirda ayni sutunlar. Bu yolda 2.5M
    satir saf Python'a gore ~10 kat hizli cozulur; ayrica pyarrow blok blok
    coklu is parcacigi kullanir.

    Beklenen duzen tutmuyorsa (orn. clean_state.py ciktisi: header yok, Site
    sutunu yok) None doner ve cagiran yavas yola duser.

    Olasiliklar float64 okunup 100000 ile carpilarak yuvarlaniyor; 5 ondalikli
    ve <= 1 olan degerlerde bu KAYIPSIZ: double'in bagil hatasi 1e-16
    seviyesinde, yani carpim gercek tam sayidan 1e-10'dan fazla sapamaz.
    """
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

    # Site sutunu, dosyanin node-major ve her node'da esit site'li oldugunu
    # kanitlar; asagidaki stride hilesi ve build_arrays'teki reshape buna dayanir.
    sites = table.column("Site").combine_chunks().to_numpy(zero_copy_only=False)
    n_sites = int(sites.max())
    if n_rows % n_sites:
        raise ValueError(f"{n_rows} satir {n_sites} site'a bolunmuyor")
    n_nodes = n_rows // n_sites
    if not np.array_equal(sites, np.tile(np.arange(1, n_sites + 1, dtype=sites.dtype),
                                         n_nodes)):
        return None  # beklenmeyen satir sirasi -> guvenli yola dus

    # Node her blokta sabit; sadece blok baslarindaki degeri okumak yeterli.
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
    """Dosyayi okuyup (values, codes, vocab, node_ids, prob_columns) dondurur.

    values    : (n_rows, k) int32, dosyadaki sira (node-major)
    codes     : n_rows uzunlugunda int16, State sutununun vocab'daki indeksi
    vocab     : codes'un cozuldugu etiket demeti
    node_ids  : node bloklarinin dosyadaki gorulme sirasi
    """
    fast = _read_state_fast(in_path)
    return fast if fast is not None else _read_state_slow(in_path)


def build_arrays(in_path: str):
    """Bir .state dosyasini okuyup kompakt sutunlari ve metadata'sini kurar.

    Doner: (arrays, metadata, n_nodes, n_rows)
      arrays   : {"dropped": .., "v0": .., ..., "sum_adj": ..} pyarrow dizileri
      metadata : stb_* anahtarli bytes sozlugu

    Kucultmeyi saglayan 4 yapisal gercek (gercek veri uzerinde olculdu):
      1. State sutunu = argmax; ayrica saklanmaz, argmax'tan turetilir.
         (Nadiren IQ-TREE State'e alfabede olmayan bir deger yazar, orn. "-";
         bunlar metadata'da istisna olarak tutulur.)
      2. Olasiliklarin toplami sabite yakindir (100000 +/- 1); tek baytlik
         sum_adj ile bir olasilik sutunu tamamen turetilebilir hale gelir.
      3. Dusurulecek sutun argmax olmali: kalanlar hep kucuk kalir ve uint16'ya
         sigar.
      4. Satirlar site-major dizilir; ardisik satirlar ayni site'in farkli
         node'lari olunca filogenetik korelasyon sikismayi belirgin iyilestirir.

    Node ve Site sutunlari hic saklanmaz; satir indeksinden turetilir:
        site = r // n_nodes,  node = node_ids[r % n_nodes]
    """
    values, codes, vocab, node_ids, prob_columns = _read_state_file(in_path)

    n_rows, k = values.shape
    n_nodes = len(node_ids)
    if n_rows % n_nodes:
        raise ValueError(f"{n_rows} satir {n_nodes} node'a bolunmuyor")
    n_sites = n_rows // n_nodes

    # --- node-major -> site-major ---
    sm = values.reshape(n_nodes, n_sites, k).transpose(1, 0, 2).reshape(-1, k).copy()
    codes_sm = codes.reshape(n_nodes, n_sites).T.reshape(-1)

    dropped = np.argmax(sm, axis=1).astype(np.int8)
    sum_adj = (sm.sum(axis=1) - SCALE).astype(np.int8)

    # argmax sutununu cikar, kalanlari orijinal sutun sirasinda topla
    rest = np.empty((n_rows, k - 1), dtype=np.int32)
    for c in range(k):
        m = dropped == c
        if m.any():
            rest[m] = sm[np.ix_(m, [j for j in range(k) if j != c])]
    if rest.size and rest.max() > np.iinfo(np.uint16).max:
        raise ValueError("kalan deger uint16'ya sigmiyor")
    rest = rest.astype(np.uint16)

    # State sutunu argmax'tan turetilemeyen satirlar (orn. "-") istisna listesi.
    # Her satirin "beklenen" etiketi prob_columns[dropped]; kod uzayinda
    # karsilastirmak 2.5M satirda Python dongusunden kurtarir.
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
    }
    return arrays, metadata, n_nodes, n_rows


def _one(role, hits):
    if not hits:
        raise SystemExit(f"hata: {role} dosyasi bulunamadi")
    if len(hits) > 1:
        listed = "\n  ".join(str(src) for _, src in hits)
        raise SystemExit(f"hata: birden fazla {role} dosyasi var, "
                         f"tek tek yol verin:\n  {listed}")
    return hits[0]


def _select(entries):
    """[(goreli_yol, _Source)] icinden dort girdiyi secer -> (nt, gap, fasta, tree).

    Klasor duzeni sabit degil; ayrimi goreli yoldaki "binary" kelimesi yapar
    (3_binary_iqtree_ancestral/... -> gap/indel, 2_iqtree_ancestral/... ->
    nukleotid). Dosya adinda da arandigi icin duz duzenler (raw.state +
    raw_binary.state) da calisir.
    """
    entries = sorted(entries, key=lambda e: e[0])
    states = [e for e in entries if e[0].endswith(STATE_SUFFIXES)]
    nt_rel, nt = _one("nukleotid .state",
                      [e for e in states if "binary" not in e[0].lower()])
    _, gap = _one("binary .state", [e for e in states if "binary" in e[0].lower()])
    _, fasta = _one("hizalama .fasta",
                    [e for e in entries if e[0].endswith(FASTA_SUFFIXES)])

    # Agac tercihen nukleotid .state'in yanindaki olsun: bundle'in sakladigi
    # node adlari o agactan geliyor.
    trees = [e for e in entries if e[0].endswith(".treefile")]
    parent = PurePosixPath(nt_rel).parent
    beside = [e for e in trees if PurePosixPath(e[0]).parent == parent]
    _, tree = _one("agac .treefile", beside or trees)

    return nt, gap, fasta, tree


def discover(root):
    """PHACT sonuc klasorunden dort girdiyi bulur -> (nt, gap, fasta, tree)."""
    root = Path(root)
    if not root.is_dir():
        raise SystemExit(f"hata: klasor degil: {root}")
    entries = [(p.relative_to(root).as_posix(), _Source(p))
               for p in root.rglob("*") if p.is_file()]
    return _select(entries)


def _strip_common_root(names):
    """Tum uyeler ayni ust klasordeyse o klasoru atar.

    Tar'lar genelde "chr1_100001-110000/..." diye tek bir kok klasor tasiyor.
    O kok, klasor yolundan gelen goreli yola karsilik gelmiyor; ustelik adinda
    "binary" gecerse nt/gap ayrimini bozardi.
    """
    first = names[0].split("/")[0]
    if all("/" in n and n.split("/")[0] == first for n in names):
        return [n.split("/", 1)[1] for n in names]
    return list(names)


def discover_tar(tar_path):
    """Tar arsivinden dort girdiyi bulur -> (nt, gap, fasta, tree).

    Arsiv diske ACILMAZ; her uye okunacagi zaman kendi akisi olarak acilir.
    """
    with tarfile.open(tar_path, "r:*") as tar:
        names = [m.name for m in tar.getmembers() if m.isfile()]
    if not names:
        raise SystemExit(f"hata: tar arsivinde dosya yok: {tar_path}")
    entries = [(rel, _Source(name, tar_path=tar_path, member=name))
               for name, rel in zip(names, _strip_common_root(names))]
    return _select(entries)


def _write_settings(names, compression, compression_level):
    """pq.write_table icin kodlama/sikistirma ayarlari.

    Isin tamami burada donuyor: Parquet'in INT32'den dar bir fiziksel tipi yok,
    yani int8/uint16 sutunlar PLAIN kodlamayla 4 bayt/deger yaziliyor ve
    sikistirici mantiksal verinin ~2 katini, cogu sifir dolgu, ciginemek zorunda
    kaliyordu. Cozum sikistirmayi zayiflatmak degil, ona giden bayt hacmini
    Parquet kodlamalariyla dusurmek; sonuc hem cok hizli hem daha kucuk.

    chr2_10001-20000 (2.51M satir) uzerinde olculdu:
        PLAIN     + brotli11   74.4 s   2.962 MB   <- eski varsayilan
        BSS       + brotli11   34.4 s   2.681 MB
        BSS+dict  + brotli11   23.7 s   2.581 MB
        BSS+dict  + brotli9     0.9 s   2.753 MB   <- yeni varsayilan
        BSS+dict  + zstd13      0.8 s   2.842 MB
        BSS+dict  + zstd19      4.4 s   2.670 MB

      - BYTE_STREAM_SPLIT (v* sutunlari): INT32'yi 4 bayt duzlemine ayirir;
        uint16 degerlerin ust iki duzlemi tamamen sifir oldugu icin hem hacim
        duser hem oran iyilesir.
      - sozluk (dropped_*, sum_adj_*): kardinalite 2-4, RLE_DICTIONARY bu
        sutunlari neredeyse sifira indirir. v* sutunlarinda ise sozluk ZARARLI
        (binlerce farkli deger; indeksler ham degerlerden kotu sikisiyor).
      - blob sutunu sikistirilmaz: icerik zaten xz'li, ikinci tur bosa zaman.
      - buyuk sayfa (64 MB) sikistiriciya genis pencere verir; varsayilan 1 MB'lik
        sayfalar sikistirmayi surekli bastan baslatiyordu.
    """
    wide = [n for n in names if n.startswith("v")]           # uint16 olasiliklar
    small = [n for n in names if not n.startswith("v")]      # int8 dropped/sum_adj
    settings = {
        "use_dictionary": small,
        "column_encoding": {n: "BYTE_STREAM_SPLIT" for n in wide},
        "compression": {**{n: compression for n in names}, "blob": "none"},
        "data_page_size": 64 * 1024 * 1024,
    }
    if compression_level is not None:
        # blob'un codec'i "none"; ona seviye verilirse pyarrow hata veriyor.
        settings["compression_level"] = {n: compression_level for n in names}
    return settings


def build(nt_path, gap_path, fasta_path, tree_path, out_path,
          compression="brotli", compression_level=9, jobs=3):
    # Girdiler ya duz yol ya da tar uyesi olabilir; ikisi de _Source'a sarilir.
    nt_path, gap_path = _as_source(nt_path), _as_source(gap_path)
    fasta_path, tree_path = _as_source(fasta_path), _as_source(tree_path)

    def _blobs():
        return (lzma.compress(fasta_path.read_bytes(), preset=9 | lzma.PRESET_EXTREME),
                lzma.compress(tree_path.read_bytes(), preset=9))

    # Uc is de birbirinden bagimsiz: iki .state cozumu ve blob'larin xz'lenmesi.
    # pyarrow.csv, numpy ve lzma GIL'i biraktigi icin is parcaciklari gercek
    # paralellik veriyor. Ucten fazlasinin anlami yok.
    if jobs > 1:
        with ThreadPoolExecutor(max_workers=min(jobs, 3)) as pool:
            fut_nt = pool.submit(build_arrays, nt_path)
            fut_gap = pool.submit(build_arrays, gap_path)
            fut_blob = pool.submit(_blobs)
            nt_arr, nt_md, n_nodes, n_rows = fut_nt.result()
            gap_arr, gap_md, g_nodes, g_rows = fut_gap.result()
            fasta_xz, tree_xz = fut_blob.result()
    else:
        nt_arr, nt_md, n_nodes, n_rows = build_arrays(nt_path)
        gap_arr, gap_md, g_nodes, g_rows = build_arrays(gap_path)
        fasta_xz, tree_xz = _blobs()

    # Iki tablonun ayni (node, site) izgarasini paylastigini dogrula --
    # paylasmiyorlarsa yan yana koymak sessizce yanlis veri uretirdi.
    if nt_md[b"stb_node_ids"] != gap_md[b"stb_node_ids"]:
        raise ValueError("iki .state dosyasinin node sirasi ayni degil")
    if nt_md[b"stb_n_sites"] != gap_md[b"stb_n_sites"]:
        raise ValueError("iki .state dosyasinin site sayisi ayni degil")
    if n_rows != g_rows:
        raise ValueError(f"satir sayilari farkli: {n_rows} vs {g_rows}")

    cols = {f"{k}_nt": v for k, v in nt_arr.items()}
    cols.update({f"{k}_gap": v for k, v in gap_arr.items()})
    names = list(cols)

    blob = [None] * n_rows
    blob[FASTA_ROW] = fasta_xz
    blob[TREE_ROW] = tree_xz
    cols["blob"] = pa.array(blob, type=pa.binary())

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
        description="PHACT sonuclarini tek Parquet bundle'a paketler.")
    ap.add_argument("inputs", nargs="+", metavar="GIRDI",
                    help="ya <sonuc_klasoru|sonuc.tar> [cikti.parquet], ya da "
                         "<nt.state> <binary.state> <fasta> <treefile> <cikti.parquet>")
    ap.add_argument("--codec", default="brotli",
                    choices=sorted(DEFAULT_LEVEL), help="varsayilan: brotli")
    ap.add_argument("--level", type=int, default=None,
                    help="sikistirma seviyesi (brotli 9, zstd 13 varsayilan). "
                         "brotli 11 %6 daha kucuk dosya verir ama ~25 kat yavas.")
    ap.add_argument("--jobs", type=int, default=3,
                    help="hazirlik is parcacigi sayisi: iki .state cozumu + "
                         "blob xz'i. 1 = seri (~150 MB az bellek, ~2 s yavas); "
                         "3'ten fazlasinin etkisi yok. Sikistirmayi etkilemez.")
    args = ap.parse_args(argv)

    first = args.inputs[0]
    is_tar = first.endswith(TAR_SUFFIXES) and os.path.isfile(first)
    if os.path.isdir(first) or is_tar:
        if len(args.inputs) > 2:
            ap.error("klasor/tar verildiginde en fazla bir cikti yolu alinir")
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
        print(f"nt    : {nt}\ngap   : {gap}\nfasta : {fasta}\n"
              f"tree  : {tree}\ncikti : {out}", file=sys.stderr)
    elif len(args.inputs) == 5:
        nt, gap, fasta, tree, out = args.inputs
    else:
        ap.error("ya bir sonuc klasoru/tar'i, ya da tam 5 yol verin")

    level = args.level if args.level is not None else DEFAULT_LEVEL[args.codec]
    codec = args.codec + (str(level) if level is not None else "")

    # Yuksek seviyeler tek bir C++ cagrisinda geciyor ve o sirada ekrana hicbir
    # sey akmiyor; Ctrl+C bile cagri donene kadar islenmiyor. En azindan neyin
    # beklendigi yazili olsun.
    slow = ((args.codec == "brotli" and level >= 10)
            or (args.codec == "zstd" and level >= 18))
    print(f"sikistirma: {codec}" + (" (10k site'lik blokta ~25 s surer; "
          "acele varsa varsayilan seviye ~25x hizli, %6 daha buyuk)" if slow else ""),
          file=sys.stderr, flush=True)

    t0 = time.time()
    build(nt, gap, fasta, tree, out,
          compression=args.codec, compression_level=level, jobs=args.jobs)
    print(f"bitti: {out} ({os.path.getsize(out)/1e6:.2f} MB, "
          f"{time.time() - t0:.1f} s)", file=sys.stderr)


if __name__ == "__main__":
    main()
