#!/usr/bin/env python3
"""
Temizlenmis .state dosyasini (clean_state.py ciktisi) kayipsiz, sikistirma
kullanmadan kompakt bir ikili formata cevirir.

Kazanc kaynaklari (hepsi kayipsiz):
  - Node ve Site sutunlari hic saklanmiyor; her node'un site'lari 1..N
    ardisik oldugu icin pozisyondan geri turetiliyor.
  - State harfi (A/C/G/T) 2 bit'e sikisiyor (8 bit yerine).
  - p_A/p_C/p_G/p_T: deger 0 ise hic bit harcanmiyor (sadece mask'ta 1 bit),
    deger tam 1.00000 ise (butun site "kesin") sadece 1 bit ile isaretleniyor,
    diger her deger tam 5 ondalik hassasiyetle (0-99999 arasi tamsayi, 17 bit)
    aynen saklaniyor.

Format:
  magic:            4 byte  b"STB1"
  num_nodes:        4 byte  uint32 LE
  sites_per_node:   4 byte  uint32 LE
  node_ids:         num_nodes x 2 byte uint16 LE  (dosyadaki sira ile)
  bitstream:        her satir icin:
      [2 bit state_code: A=0,C=1,G=2,T=3]
      [1 bit certain_flag]
      certain_flag==1 ise:  (baska bit yok; state kolonu=100000, digerleri=0)
      certain_flag==0 ise:
          [4 bit mask: A,C,G,T sirasiyla, o kolon sifirdan farkli mi]
          mask'ta 1 olan her kolon icin: [17 bit deger, 0-99999]
      (bitstream'in sonu 0 bit ile bayt sinirina tamamlanir)

Kullanim:
    python state_to_binary.py <temiz.state> <cikti.stb>
"""

import sys

STATE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
COLUMNS = ("A", "C", "G", "T")  # p_A, p_C, p_G, p_T sirasi


class BitWriter:
    __slots__ = ("buf", "acc", "nbits")

    def __init__(self):
        self.buf = bytearray()
        self.acc = 0
        self.nbits = 0

    def write(self, value, n):
        self.acc = (self.acc << n) | value
        self.nbits += n
        buf = self.buf
        acc = self.acc
        nbits = self.nbits
        while nbits >= 8:
            nbits -= 8
            buf.append((acc >> nbits) & 0xFF)
        self.acc = acc & ((1 << nbits) - 1) if nbits else 0
        self.nbits = nbits

    def flush(self):
        if self.nbits:
            self.buf.append((self.acc << (8 - self.nbits)) & 0xFF)
            self.acc = 0
            self.nbits = 0
        return self.buf


def parse_exact_5dec(raw: str) -> int:
    # "0.99999" -> 99999 ; "0.5" -> 50000 ; hicbir float() kullanmadan, tam ondalikli
    if "." in raw:
        intpart, decpart = raw.split(".")
    else:
        intpart, decpart = raw, ""
    decpart = (decpart + "00000")[:5]
    return int(intpart) * 100000 + int(decpart)


def encode(in_path: str, out_path: str) -> None:
    node_ids = []
    prev_node = None
    sites_per_node = None
    site_count_this_node = 0

    bw = BitWriter()

    with open(in_path, "r", encoding="utf-8") as fin:
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            node, site, state = cols[0], cols[1], cols[2]
            vals = cols[3:7]

            if node != prev_node:
                if prev_node is not None:
                    if sites_per_node is None:
                        sites_per_node = site_count_this_node
                    elif site_count_this_node != sites_per_node:
                        raise ValueError(f"Node {prev_node}: {site_count_this_node} site (beklenen {sites_per_node})")
                node_ids.append(int(node))
                prev_node = node
                site_count_this_node = 0

            site_count_this_node += 1

            state_code = STATE_CODE[state]
            bw.write(state_code, 2)

            state_col = COLUMNS.index(state)
            if vals[state_col] == "1" and all(vals[i] == "" for i in range(4) if i != state_col):
                bw.write(1, 1)  # certain
                continue

            bw.write(0, 1)
            mask = 0
            entries = []
            for i in range(4):
                if vals[i] != "":
                    mask |= (1 << (3 - i))
                    entries.append(vals[i])
            bw.write(mask, 4)
            for raw in entries:
                bw.write(parse_exact_5dec(raw), 17)

        # son node'un site sayisini da dogrula
        if sites_per_node is not None and site_count_this_node != sites_per_node:
            raise ValueError(f"Node {prev_node}: {site_count_this_node} site (beklenen {sites_per_node})")

    payload = bw.flush()

    with open(out_path, "wb") as fout:
        fout.write(b"STB1")
        fout.write(len(node_ids).to_bytes(4, "little"))
        fout.write((sites_per_node or 0).to_bytes(4, "little"))
        for nid in node_ids:
            fout.write(nid.to_bytes(2, "little"))
        fout.write(payload)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Kullanim: python {sys.argv[0]} <temiz.state> <cikti.stb>", file=sys.stderr)
        sys.exit(1)
    encode(sys.argv[1], sys.argv[2])
