#!/usr/bin/env python3
"""
state_to_binary.py ile uretilen .stb dosyasini, clean_state.py'nin ciktisiyla
BAYT BAYT AYNI olacak sekilde tekrar tab-separated .state dosyasina cevirir.

Kullanim:
    python binary_to_state.py <girdi.stb> <cikti.state>
"""

import sys

COLUMNS = ("A", "C", "G", "T")


class BitReader:
    __slots__ = ("data", "pos", "acc", "nbits")

    def __init__(self, data: bytes):
        self.data = data
        self.pos = 0
        self.acc = 0
        self.nbits = 0

    def read(self, n):
        while self.nbits < n:
            self.acc = (self.acc << 8) | self.data[self.pos]
            self.pos += 1
            self.nbits += 8
        self.nbits -= n
        value = (self.acc >> self.nbits) & ((1 << n) - 1)
        self.acc &= (1 << self.nbits) - 1
        return value


def format_5dec(intval: int) -> str:
    whole, frac = divmod(intval, 100000)
    return f"{whole}.{frac:05d}"


def decode(in_path: str, out_path: str) -> None:
    with open(in_path, "rb") as f:
        data = f.read()

    if data[:4] != b"STB1":
        raise ValueError("Gecersiz dosya (magic uyusmuyor)")

    num_nodes = int.from_bytes(data[4:8], "little")
    sites_per_node = int.from_bytes(data[8:12], "little")

    offset = 12
    node_ids = []
    for _ in range(num_nodes):
        node_ids.append(int.from_bytes(data[offset:offset + 2], "little"))
        offset += 2

    br = BitReader(data[offset:])

    with open(out_path, "w", encoding="utf-8", newline="\n") as fout:
        for node_id in node_ids:
            node_str = str(node_id)
            for site in range(1, sites_per_node + 1):
                state_code = br.read(2)
                certain = br.read(1)
                state = COLUMNS[state_code]

                if certain:
                    vals = ["", "", "", ""]
                    vals[state_code] = "1"
                else:
                    mask = br.read(4)
                    vals = ["", "", "", ""]
                    for i in range(4):
                        if mask & (1 << (3 - i)):
                            vals[i] = format_5dec(br.read(17))

                fout.write(f"{node_str}\t{site}\t{state}\t{vals[0]}\t{vals[1]}\t{vals[2]}\t{vals[3]}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Kullanim: python {sys.argv[0]} <girdi.stb> <cikti.state>", file=sys.stderr)
        sys.exit(1)
    decode(sys.argv[1], sys.argv[2])
