# Ossuary

Bir PHACT/IQ-TREE ata-durumu analizinin dort ciktisini -- nukleotid olasiliklari,
gap/indel olasiliklari, hizalama (FASTA) ve agac (Newick) -- tek bir **kayipsiz**
Parquet dosyasinda toplar. Ad buradan geliyor: atalarin hepsi tek bir kapta.

```bash
# paketle (klasor, tar ya da dort dosya tek tek)
python ossuary.py chr2_10001-20000/                # -> chr2_10001-20000.parquet
python ossuary.py chr1_100001-110000.tar

# kucultmek icin hassasiyeti dusur (varsayilan 5 = kayipsiz)
python ossuary.py chr1_100001-110000.tar --precision 3   # 4.58 -> 1.67 MB

# oku (sorgula)
python ossuary_reader.py bundle.parquet            # ozet
python ossuary_reader.py bundle.parquet <node> <site>

# geri ac (dort orijinal dosyayi byte-byte geri uretir)
python reborn.py bundle.parquet cikti_klasoru/
python reborn.py bundle.parquet --list             # ne var icinde
```

Kutuphane olarak:

```python
from ossuary_reader import BundleReader
b = BundleReader("bundle.parquet")
b.get(node_id=2, site=11)   # {'nt': ('C', 0.00083, ...), 'gap': ('1', 0.0, 1.0)}
b.tree()                    # Newick metni
b.fasta()                   # {tur_adi: hizalanmis_dizi}
```

Ayrinti (sema, sikistirma olcumleri, girdi bicimleri) `ossuary.py`'nin
basindaki aciklamada. Satir satir Turkce anlatim icin: `OZET.md`.
