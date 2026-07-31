# Ossuary

Bir PHACT/IQ-TREE ata-durumu analizinin dort ciktisini -- nukleotid olasiliklari,
gap/indel olasiliklari, hizalama (FASTA) ve agac (Newick) -- tek bir **kayipsiz**
Parquet dosyasinda toplar. Ad buradan geliyor: atalarin hepsi tek bir kapta.

```bash
# paketle (klasor, tar ya da dort dosya tek tek)
python make_bundle.py chr2_10001-20000/            # -> chr2_10001-20000.parquet
python make_bundle.py chr1_100001-110000.tar

# oku
python bundle_reader.py bundle.parquet             # ozet
python bundle_reader.py bundle.parquet <node> <site>
```

Kutuphane olarak:

```python
from bundle_reader import BundleReader
b = BundleReader("bundle.parquet")
b.get(node_id=2, site=11)   # {'nt': ('C', 0.00083, ...), 'gap': ('1', 0.0, 1.0)}
b.tree()                    # Newick metni
b.fasta()                   # {tur_adi: hizalanmis_dizi}
```

Ayrinti (sema, sikistirma olcumleri, girdi bicimleri) `make_bundle.py`'nin
basindaki aciklamada.
