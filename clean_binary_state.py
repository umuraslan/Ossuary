#!/usr/bin/env python3
"""
IQ-TREE .state dosyasini temizler:
  - Yorum satirlarini kaldirir (# ile baslayan satirlar)
  - Header satirlarini kaldirir (Node  Site  State  p_A  p_C  p_G  p_T)
  - 4,5,6,7. sutunlarda (p_A, p_C, p_G, p_T):
      * deger 0 ise -> hucre bos birakilir
      * deger 1 ise -> "1" yazilir (1.00000 degil)
      * diger degerler oldugu gibi birakilir
  - 1. sutundaki "NodeXXX" degerinden "Node" yazisi kaldirilir (sadece "XXX" kalir)

Kullanim:
    python clean_state.py <girdi.state> <cikti.state>
"""

import sys


def clean_value(raw: str) -> str:
    try:
        f = float(raw)
    except ValueError:
        return raw
    if f == 0:
        return ""
    if f == 1:
        return "1"
    return raw


def clean_node(raw: str) -> str:
    if raw.startswith("Node"):
        return raw[len("Node"):]
    return raw


def drop_site_column(cols: list) -> list:
    # Site (2. sutun) her node icinde 1'den baslayip kesintisiz arttigi
    # icin satir sirasindan geri turetilebilir; dosyada tutmaya gerek yok.
    return [cols[0]] + cols[2:]


def sparse_node(cols: list, prev_node: str) -> tuple:
    # Node bir onceki satirla ayniysa tekrar yazma, bos birak;
    # sadece node degistiginde (blogun ilk satirinda) yaz.
    node = cols[0]
    if node == prev_node:
        cols = [""] + cols[1:]
    return cols, node


def process(in_path: str, out_path: str) -> None:
    value_cols = (3, 4)  # 0-indexli: 4,5. sutunlar (p_0, p_1 - bu dosyada sadece 2 olasilik sutunu var)
    prev_node = None

    with open(in_path, "r", encoding="utf-8", newline="") as fin, \
         open(out_path, "w", encoding="utf-8", newline="") as fout:
        for line in fin:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                continue

            cols = line.split("\t")

            # Header satirini atla
            if cols[0] == "Node":
                continue

            cols[0] = clean_node(cols[0])

            for i in value_cols:
                if i < len(cols):
                    cols[i] = clean_value(cols[i])

            cols = drop_site_column(cols)
            cols, prev_node = sparse_node(cols, prev_node)
            fout.write("\t".join(cols) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Kullanim: python {sys.argv[0]} <girdi.state> <cikti.state>", file=sys.stderr)
        sys.exit(1)

    process(sys.argv[1], sys.argv[2])
