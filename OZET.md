# Ossuary — Satir Satir Turkce Anlatim

Iki dosya var: `ossuary.py` (paketleyici) ve `ossuary_reader.py` (okuyucu).
Asagida her ikisinin de ne yaptigi, bolum bolum ve satir numaralariyla anlatiliyor.

> Kod tarafinda yorumlar bilincli olarak minimumda tutuldu (fonksiyon basina en
> fazla bir satir). Tasarim gerekceleri ve olcum sonuclari **bu dosyada** duruyor;
> `ossuary.py`'de gormeyeceginiz her aciklama burada.

---

## 1. Genel Fikir

Bir PHACT/IQ-TREE ata-durumu analizi dort dosya uretir:

| Dosya | Icerik |
|---|---|
| `<ad>.state` | Nukleotid ata-durumu olasiliklari (ic node x site) |
| `<ad>.state` (binary) | Gap/indel ata-durumu olasiliklari (ayni node x site izgarasi) |
| `<ad>.fasta` | Yapraklarin hizalanmis dizileri (ayni site sayisi) |
| `<ad>.treefile` | Newick agaci |

Bunlar ayri veri kumeleri degil, **tek bir veri kumesinin dort yuzu**: FASTA
basliklari agacin yapraklariyla birebir ayni, iki `.state` dosyasi ayni node
sirasini ve ayni site sayisini paylasiyor.

`ossuary.py` bu dordunu **tek bir kayipsiz Parquet dosyasina** paketler.
`ossuary_reader.py` o dosyayi geri okur. Kayipsizlik gercek: dosyadan geri
okunan degerler orijinal `.state` dosyasindaki degerlerin aynisidir.

**Kazanc:** 212 MB `.state` + 137 MB binary `.state` + 4.5 MB FASTA + 20 KB
agac ≈ 354 MB girdi → **4.58 MB** tek Parquet dosyasi (chr1_100001-110000
ornegi, ~5 saniyede).

---

## 2. `ossuary.py` — Paketleyici

### 2.1 Modul basligi ve sabitler (satir 2-62)

- **Satir 3-26 (docstring):** Sema, tasarim gerekceleri ve kullanim ornekleri.
  Iki olculmus karar burada aciklaniyor:
  - FASTA'yi tabloya sutun olarak acmak sikismayi **149 KB'dan 300 KB'a**
    cikariyor (turler arasi uzun LZ eslesmeleri dagiliyor). Bu yuzden FASTA ve
    agac, xz'lenmis birer blob olarak `blob` adli BINARY sutunda tasiniyor.
  - Ayni blob'u Parquet key-value metadata'sina koymak UTF-8 genislemesi
    yuzunden **156 KB'lik yuku 365 KB'a** sisiriyor; BINARY sutunda sisme yok.
- **Satir 49 `SCALE = 100000`:** Olasiliklar `.state` dosyasinda 5 ondalikli
  yaziliyor (`0.99751`). Float yerine 100000 ile carpilmis **tam sayi** olarak
  saklaniyor — hem kayipsiz hem cok daha iyi sikisiyor.
- **Satir 51-52 `FASTA_ROW`, `TREE_ROW`:** `blob` sutununda FASTA 0. satirda,
  agac 1. satirda; gerisi NULL.
- **Satir 54 `JOBS = 3`:** hazirlik is parcacigi sayisi; bkz. bolum 2.9.
- **Satir 56-58:** Girdi dosyalarini uzantisindan tanimak icin son ek listeleri.
- **Satir 60-62:** `.gz`, `.xz`, `.bz2` icin cozucu esleme tablolari. Arrow'un
  xz codec'i olmadigi icin `_ARROW_CODECS` sadece gzip ve bz2 iceriyor.

### 2.2 `_Source` sinifi (satir 65-116)

Okunacak bir girdiyi temsil eder: ya **diskteki bir dosya yolu**, ya da **bir
tar arsivinin uyesi**. Ikisi de ayni arayuzden okunur, boylece geri kalan tum
kod "bu dosya mi tar icinde mi" diye sormak zorunda kalmaz.

- **Satir 68 `__slots__`:** Sinifa `__dict__` verilmez; hafif nesne.
- **Satir 82-91 `open_raw()`:** Ham (belki sikistirilmis) bayt akisi.
  - **Satir 90 `stream._tar = tar`** — kritik incelik: `tarfile` nesnesi cop
    toplanirsa acik akis da kapanir. Referansi akisin uzerinde canli tutuyoruz.
  - Her acilista **yeni bir TarFile** aciliyor, cunku `tarfile` is parcaciklari
    arasinda paylasilamaz ve `build()` uc okumayi paralel yapiyor.
- **Satir 93-97 `open_binary()`:** Uzantiya bakip gerekiyorsa cozucu sarar.
- **Satir 99-100 `open_text()`:** UTF-8 metin akisi.
- **Satir 106-113 `csv_input()`:** pyarrow.csv'ye verilecek girdiyi hazirlar.
  - Normal dosyada **yolun kendisi** verilir → Arrow dosyayi kendi acip
    bloklari paralel isler (en hizli yol).
  - Tar uyesinde bu mumkun degil, akis verilir; cozmeyi yine Arrow'un C++
    codec'i yapsin diye `CompressedInputStream`'e sarilir (Python'un `gzip`
    modulunden belirgin hizli).

### 2.3 `.state` okuma (satir 129-289)

Iki farkli girdi sekli desteklenir:
- **Ham IQ-TREE ciktisi:** `# ...` yorum satirlari, bir header satiri
  (`Node  Site  State  p_A  p_C  p_G  p_T`), Node her satirda yazili.
- **Temizlenmis cikti:** yorum/header yok, Site sutunu yok, Node sadece blok
  basinda yazili, olasilik 0 ise bos, 1 ise sadece `1`.

**`parse_exact_5dec` (satir 136-146).** `"0.99999"` → `99999`, `"0.5"` →
`50000`. Hicbir yerde `float()` cagrilmaz: metin nokta yerinden ikiye bolunur,
ondalik kismi 5 haneye tamamlanir, tam sayi aritmetigiyle birlestirilir. Bu
sayede yuvarlama hatasi teorik olarak bile imkansiz.

**`strip_node_prefix` (satir 148-149).** `"Node12"` → `"12"`.

**`_detect_layout` (satir 152-167).** Header satirindan p_X sutunlarini cikarir.
Header yoksa ilk veri satirinin sutun sayisindan tahmin eder: 4 olasilik sutunu
varsa alfabe `("A","C","G","T")`, degilse `("0","1",...)` gibi sayisal.
Ucuncu donus degeri `is_header`, satirin atlanmasi gerektigini soyler.

**`_iter_data_lines` (satir 170-184).** Bos satirlari ve `#` yorumlarini atlar,
her veri satirini `(cols, prob_columns, has_site_col)` olarak uretir.

**`_encode_labels` (satir 187-196).** Etiket listesini `(kodlar, sozluk)`
haline getirir — `"A"`,`"C"`,... yerine `int16` indeksler.

**`_read_state_slow` (satir 199-226).** Saf Python, satir satir okuma.
**Her** girdi seklini kaldirir, bu yuzden guvenli yedek yoldur.
- **Satir 214:** Node sutunu bossa (temizlenmis format) bir onceki node
  kullanilir; degistiginde `node_ids` listesine eklenir.

**`_read_state_fast` (satir 229-282).** Ham IQ-TREE duzenini `pyarrow.csv` ile
**vektorize** okur. Bu yolda 2.5M satir saf Python'a gore **~10 kat hizli**
cozulur; ayrica pyarrow blok blok coklu is parcacigi kullanir.
- **Satir 232-241:** Yorum satirlari sayilir, header okunur.
- **Satir 243-244:** Beklenen duzen (`Node  Site  State  ...`) tutmuyorsa
  **`None` doner** ve cagiran yavas yola duser. Guvenlik supabi budur.
- **Satir 247-256:** Sutun tipleri acikca verilerek CSV okunur (16 MB'lik
  bloklar halinde).
- **Satir 258-265:** `Site` sutunu, dosyanin gercekten node-major oldugunu ve
  her node'da esit sayida site bulundugunu **kanitlar**. Beklenen desen
  (`1,2,...,n_sites` n_nodes kere tekrar) tutmuyorsa yine `None` donulur.
  Asagidaki tum reshape hilesi bu garantiye dayaniyor.
- **Satir 267-269:** Node her blokta sabit oldugu icin sadece blok
  baslarindaki deger okunur — 2.5M satir yerine 446 deger.
- **Satir 271-273:** `State` sutunu Arrow'un sozluk kodlamasiyla `int16`
  kodlara cevrilir.
- **Satir 275-281:** Olasiliklar `float64` okunup 100000 ile carpilarak
  yuvarlanir. **Bu kayipsizdir:** 5 ondalikli ve `<= 1` olan degerlerde
  double'in bagil hatasi 1e-16 seviyesinde, yani carpim gercek tam sayidan
  1e-10'dan fazla sapamaz — `np.rint` her zaman dogru tam sayiyi verir.

**`_read_state_file` (satir 285-289).** Once hizli yolu dener, olmazsa yavas
yola duser. Donus: `(values, codes, vocab, node_ids, prob_columns)`.

### 2.4 `build_arrays` — Isin Kalbi (satir 300-351)

Bir `.state` dosyasini okuyup **kompakt sutunlara** cevirir. Kucultmeyi
saglayan **dort yapisal gercek** (hepsi gercek veri uzerinde olculdu):

1. **`State` sutunu = argmax.** En yuksek olasilikli harf hangisiyse `State`
   odur. Ayrica saklanmaz, argmax'tan turetilir. (Nadiren IQ-TREE alfabede
   olmayan bir deger yazar, orn. `-`; bunlar metadata'da **istisna** olarak
   tutulur.)
2. **Olasiliklarin toplami sabite yakindir** (100000 ± 1). Tek baytlik
   `sum_adj` ile bir olasilik sutunu **tamamen turetilebilir** hale gelir —
   yani 4 sutundan sadece 3'u saklanir.
3. **Dusurulecek sutun argmax olmali.** Cunku o zaman kalan degerler hep
   kucuk kalir ve `int32` yerine `uint16`'ya sigar (yari yer).
4. **Satirlar site-major dizilir.** Ardisik satirlar ayni site'in farkli
   node'lari olunca **filogenetik korelasyon** sikismayi belirgin iyilestirir
   (yakin akraba node'lar benzer olasiliklara sahip).

Ayrica **`Node` ve `Site` sutunlari hic saklanmaz**; satir indeksinden
turetilir:

```
site = r // n_nodes,   node = node_ids[r % n_nodes]
```

Satir satir:
- **Satir 302-311:** Satir sonu uslubu tespit edilir, dosya okunur; `out_scale`
  varsayilandan farkliysa `quantize()` uygulanir (bkz. bolum 2.10).
- **Satir 313-314:** `reshape → transpose → reshape` ile **node-major duzenden
  site-major duzene** gecilir. (4. gercek.)
- **Satir 316:** `dropped` = her satirin argmax'i, `int8`.
- **Satir 317:** `sum_adj` = satir toplami − 100000, `int8` (2. gercek).
- **Satir 319-323:** argmax sutunu cikarilir, kalanlar **orijinal sutun
  sirasinda** toplanir. Dongu k (=4) kez doner, her seferinde argmax'i o sutun
  olan tum satirlar maskeyle toplu islenir.
- **Satir 324-326:** `uint16` tasma kontrolu (3. gercege guvenmiyoruz, dogruluyoruz).
- **Satir 328-333:** `State`'i argmax'tan turetilemeyen satirlarin **istisna
  listesi**. Her satirin "beklenen" etiketi `prob_columns[dropped]`;
  karsilastirma kod uzayinda vektorize yapilir — 2.5M satirda Python
  dongusunden kurtarir.
- **Satir 335-340:** pyarrow dizileri kurulur: `dropped` (int8), `v0..v{k-2}`
  (uint16), ve **sadece gerekiyorsa** `sum_adj`.
- **Satir 342-351:** Metadata: alfabe, node listesi, site sayisi, siralama,
  istisnalar ve `sum_adj` var mi bilgisi.

### 2.5 Girdi bulma (satir 354-410)

**`_one` (satir 354-362).** Bir rolde tam bir dosya bekler; hic yoksa ya da
birden fazlaysa **anlasilir bir hata** verir ve bulunanlari listeler.

**`_select` (satir 365-380).** Dort girdiyi secer. Klasor duzeni sabit
olmadigi icin nukleotid/gap ayrimini **goreli yoldaki "binary" kelimesi**
yapar:
- `3_binary_iqtree_ancestral/...` → gap/indel
- `2_iqtree_ancestral/...` → nukleotid

Dosya adinda da arandigi icin duz duzenler de calisir (`raw.state` +
`raw_binary.state`).
- **Satir 375-378:** Agac tercihen **nukleotid `.state`'in yanindaki** olsun,
  cunku bundle'in sakladigi node adlari o agactan geliyor.

**`discover` (satir 383-390).** Klasoru `rglob` ile tarar.

**`_strip_common_root` (satir 393-398).** Tar'lar genelde
`chr1_100001-110000/...` diye tek bir kok klasor tasir. O kok, klasor
yolundan gelen goreli yola karsilik gelmiyor; ustelik adinda "binary" gecerse
nt/gap ayrimini bozardi. Bu yuzden atilir.

**`discover_tar` (satir 401-410).** Arsiv **diske acilmaz**; uye listesi
okunur, her uye okunacagi zaman kendi akisi olarak acilir.

### 2.6 `_write_settings` — Parquet Ayarlari (satir 412-424)

Isin performans tarafinin tamami burada donuyor.

**Problem:** Parquet'in INT32'den dar bir fiziksel tipi yok. Yani `int8`/
`uint16` sutunlar PLAIN kodlamayla **4 bayt/deger** yaziliyor ve sikistirici
mantiksal verinin ~2 katini, cogu sifir dolgu, ciginemek zorunda kaliyordu.

**Cozum:** Sikistirmayi zayiflatmak degil, **ona giden bayt hacmini Parquet
kodlamalariyla dusurmek**. Sonuc hem cok hizli hem daha kucuk.

chr2_10001-20000 (2.51M satir) uzerinde olculdu:

| Kodlama + Sikistirma | Sure | Boyut | |
|---|---|---|---|
| PLAIN + brotli11 | 74.4 s | 2.962 MB | ← eski varsayilan |
| BSS + brotli11 | 34.4 s | 2.681 MB | |
| BSS+dict + brotli11 | 23.7 s | 2.581 MB | |
| **BSS+dict + brotli9** | **0.9 s** | **2.753 MB** | ← **yeni varsayilan** |
| BSS+dict + zstd13 | 0.8 s | 2.842 MB | |
| BSS+dict + zstd19 | 4.4 s | 2.670 MB | |

Yani: **80 kat hizlanma, uustelik eski varsayilandan %7 daha kucuk.**

Dort ayar:
- **`BYTE_STREAM_SPLIT` (v* sutunlari, satir 418):** INT32'yi 4 bayt duzlemine
  ayirir. `uint16` degerlerin **ust iki duzlemi tamamen sifir** oldugu icin hem
  hacim duser hem oran iyilesir.
- **Sozluk (`dropped_*`, `sum_adj_*`, satir 417):** kardinalite 2-4;
  RLE_DICTIONARY bu sutunlari neredeyse sifira indirir. `v*` sutunlarinda ise
  sozluk **zararli** (binlerce farkli deger; indeksler ham degerlerden kotu
  sikisiyor) — bu yuzden acikca sadece `small` listesine veriliyor.
- **`blob` sutunu sikistirilmaz (satir 419):** icerik zaten xz'li, ikinci tur
  bosa zaman.
- **Buyuk sayfa, 64 MB (satir 420):** sikistiriciya genis pencere verir;
  varsayilan 1 MB'lik sayfalar sikistirmayi surekli bastan baslatiyordu.
- **Satir 422-423:** `blob`'un codec'i `"none"`; ona seviye verilirse pyarrow
  hata veriyor, o yuzden seviye sadece gercek sutunlara veriliyor.

### 2.7 `build` (satir 427-491)

- **Satir 433-435 `_blobs()`:** FASTA `preset=9|PRESET_EXTREME` ile, agac
  `preset=9` ile xz'lenir.
- **Satir 438-444:** **Uc is birbirinden bagimsiz** — iki `.state` cozumu ve
  blob'larin xz'i. `pyarrow.csv`, `numpy` ve `lzma` GIL'i biraktigi icin is
  parcaciklari **gercek paralellik** verir. Is parcacigi sayisi `JOBS = 3`
  sabitiyle **sabitlenmis** durumda (ayrintisi icin bkz. bolum 2.9).
- **Satir 446-451:** Iki tablonun **ayni (node, site) izgarasini** paylastigi
  dogrulanir. Paylasmiyorlarsa yan yana koymak **sessizce yanlis veri**
  uretirdi — bu yuzden node sirasi, site sayisi ve satir sayisi tek tek
  kontrol ediliyor.
- **Satir 453-455:** Sutunlar `_nt` / `_gap` son ekleriyle birlestirilir.
- **Satir 457-460:** `blob` sutunu kurulur: 0. satir FASTA, 1. satir agac,
  gerisi NULL (artik sadece bu ikisi — bkz. bolum 4.1).
- **Satir 462-483:** Okuyucunun ve `reborn.py`'nin ihtiyaci olan her sey
  `bundle_*` anahtarlariyla sema metadata'sina yazilir (yeni eklenenler:
  `bundle_newline_*` ve `bundle_filenames`).
- **Satir 485-490:** Tablo tek bir row group olarak yazilir
  (`row_group_size=n_rows`) — sikistirici tum veriyi tek pencerede gorur.

### 2.8 `main` (satir 493-559)

- **Satir 493:** Codec basina varsayilan seviye tablosu.
- **Satir 513-538:** Girdi cozumleme — klasor mu, tar mi, yoksa bes acik yol mu.
  Cikti adi verilmemisse klasor/tar adindan turetilir.
- **Satir 530-536:** Bes-yol modunda **sira onemsiz**. Dort yol once
  `os.path.isfile` ile kontrol ediliyor, sonra roller `_select()` ile dosya
  adindan tespit ediliyor — klasor taramasinda kullanilan mantigin aynisi.
  Gerekcesi bir sonraki bolumde.
- **Satir 540-541:** Hangi dosyanin hangi role atandigi **her zaman** yazdiriliyor
  (eskiden sadece klasor/tar modunda yazdiriliyordu), boylece otomatik tespitin
  dogru karar verdigi gozle gorulebiliyor.
- **Satir 546-550:** Yuksek sikistirma seviyeleri **tek bir C++ cagrisinda**
  geciyor ve o sirada ekrana hicbir sey akmiyor; Ctrl+C bile cagri donene kadar
  islenmiyor. En azindan neyin beklendigi yaziliyor.
- **Satir 554-558:** Paketleme calistirilir, sonunda boyut ve sure yazilir.

### 2.8.1 Neden sira onemsiz hale getirildi

Bes yol elle verildiginde eskiden sira **kati** idi
(`nt.state`, `binary.state`, `fasta`, `treefile`, `cikti`). Yanlis sirada
verilince ortaya cok kotu bir davranis cikiyordu. Ornek:

```bash
python ossuary.py raw.fasta raw.state raw_binary.state raw.treefile raw.parquet
```

Burada `raw.fasta` nukleotid `.state` sanildi, `raw_binary.state` (137 MB) de
FASTA sanildi. Sonuc:

1. `build_arrays(raw.fasta)` **aninda** `IndexError` firlatti.
2. Ama ayni anda calisan ucuncu is `raw_binary.state`'i `xz -9e` ile
   sikistiriyordu — 137 MB icin **~190 saniye** (olculdu: 20 MB → 27.7 s).
3. `with ThreadPoolExecutor(...)` blogundan cikarken `shutdown(wait=True)`
   devreye girdigi icin Python, **sonucu cope gidecek** o isin bitmesini
   bekledi.

Yani kullanici acisindan: ekranda hicbir sey yokken uc dakika donan, sonra
alakasiz bir `IndexError` veren bir program. Hata mesajinin kendisi de yanlis
yeri isaret ediyordu.

Cozum agir isten **once** karar vermek: roller dosya adindan tespit ediliyor,
eksik/yanlis girdi `_select()` icindeki `_one()` tarafindan anlasilir bir
mesajla reddediliyor. Olculen fark:

| Durum | Once | Sonra |
|---|---|---|
| Yanlis sira | ~190 s sonra `IndexError` | 4.9 s'de dogru sonuc |
| Olmayan dosya | ~190 s sonra `FileNotFoundError` | 0.25 s'de `no such file: ...` |
| Rol eksik | ~190 s sonra belirsiz hata | 0.25 s'de `no binary .state file found` |

### 2.9 `JOBS = 3` — neden sabit

Eskiden `--jobs` diye bir CLI parametresiydi; **kaldirildi**, artik satir 55'de
sabit. Sebebi: ayarlanacak bir sey yok. Havuza atilan is **her zaman tam ucun**:
nt `.state` cozumu, gap `.state` cozumu, bir de FASTA+agac xz'i. Ucten fazla is
parcacigi olusturmanin karsiligi yok, azaltmak da sadece yavaslatiyor.

Bu makinede olculdu (chr1, 446 node × 10000 site):

| Is parcacigi | Sure | Tepe bellek |
|---|---|---|
| 1 (seri) | 9.34 s | 969 MB |
| 2 | 7.24 s | 1326 MB |
| **3 (sabit)** | **5.78 s** | **1333 MB** |
| 6 | 5.64 s | 1310 MB |

Uretilen dosya dort durumda da **byte-byte ayni**; bu yalnizca bir hiz/bellek
ayariydi. Bellek farki `build_arrays`'in gecici ara dizilerinden geliyor (CSV'nin
float64 tablosu + `values` + `sm` kopyasi): seri calisinca ayni anda bunlardan
sadece biri hayatta oluyor.

**Neyi etkilemedigi:** Parquet sikistirmasini/yazmasini degil (`pq.write_table`
tek bir C++ cagrisi) — onun icin `--codec` / `--level` var. Ayrica `pyarrow.csv`
kendi icinde zaten coklu is parcacigi kullaniyor; bu sabit ona dokunmuyor.

### 2.10 `--precision` ve `quantize()` — opsiyonel hassasiyet dusurme

Varsayilan `--precision 5` **kayipsizdir ve davranis hic degismez**. Daha dusuk
bir deger verilince olasiliklar o kadar ondaliga yuvarlanir (`0.99877` → `0.999`).

**`quantize()` (satir 291-297)** — en-buyuk-kalan (largest remainder / Hamilton)
yontemi. Neden duz `np.round` degil: 4 olasilik tek tek yuvarlanirsa toplamlari
1.000 tutmayabiliyor, bu veride **satirlarin %27.1'inde** oluyor. `build_arrays`
en buyuk olasiligi saklamayip `out_scale - sum(rest)` diye geri hesapladigi icin
sapmayi `sum_adj` sutununda tutmak zorunda kalirdi — 3 ondalikta **276 KB**.

`Node2 / site 11` uzerinden:

```
gercek           : 0.00083  0.99751  0.00083  0.00083   toplam 1.00000
1) tabana yuvarla: 0.000    0.997    0.000    0.000     toplam 0.997 -> 3 eksik
   kalanlar      : 0.83     0.51     0.83     0.83
2) 3 birimi en buyuk kalanlara dagit:
                   0.001    0.997    0.001    0.001     toplam 1.000  <- tam
   (duz yuvarlama: 0.001 0.998 0.001 0.001 = 1.001, sapma +1)
```

Toplam her satirda tam `out_scale` oldugu icin `sum_adj` hep sifir kalir ve mevcut
`has_sum_adj` kontrolu (satir 338) sutunu **hic yazmaz** — ek koda gerek yok.
Ayni sekilde argmax'in kaydigi 120 satiri (nt 90, gap 30) mevcut istisna mantigi
(satir 328-333) otomatik yakalayip gercek `State` harfini metadata'ya yaziyor;
maliyeti 1.9 KB.

Olculen sonuclar (chr1_100001-110000, gercek dosyalar uretilerek):

| `--precision` | Boyut | Maks hata | `State` degisen | Toplam 1.000 mu |
|---|---|---|---|---|
| **5 (varsayilan)** | **4.577 MB** | **0 (kayipsiz)** | 0 | 1.00000 ± 0.00001 |
| 3 | 1.670 MB (−%64) | 0.00075 | 120 | **tam 1.0** |
| 2 | 1.070 MB (−%77) | 0.005 | ~1120 | **tam 1.0** |

Duz yuvarlama ile 3 ondalik 1.801 MB verirdi; toplam-koruyan yontem `sum_adj`'i
yok ederek **142 KB** daha kazandiriyor.

`out_scale` sema metadata'sina `bundle_scale` olarak yaziliyor (satir 479);
`ossuary_reader.py` ve `reborn.py` bunu okuyup dogru olcekle bolüyor. Eski
bundle'larda bu anahtar yok, ikisi de 100000 varsayiyor.

**Elenen alternatifler** (3 ondalikta olculdu): `RLE_DICTIONARY` 2.396 MB,
delta kodlama 1.680 MB — ikisi de BSS+brotli9'un 1.659 MB'indan kotu. 1.5 MB'a
inmek icin ya brotli11 (1.485 MB ama **33 sn**) ya da olasiliklari xz blob yapmak
(1.483 MB, ama disaridan pandas/DuckDB ile okunamaz) gerekiyordu; ikisi de
reddedildi.

---

## 3. `ossuary_reader.py` — Okuyucu

Paketlemenin **tam tersi**: kompakt sutunlardan orijinal olasilik tablosunu
geri kurar.

### 3.1 `BundleReader.__init__` (satir 40-64)

Parquet dosyasi **acilir ama okunmaz** — sadece sema metadata'si alinir
(`bundle_node_ids`, `bundle_n_sites`, alfabeler, istisnalar, `blob` satir
haritasi). Agir olasilik sutunlari ilk ihtiyac aninda okunur (**tembel
yukleme**); `preload=True` verilirse hemen okunur.
- **Satir 57 `_node_pos`:** node id → dosyadaki sira; satir indeksi hesabi
  icin.

### 3.2 `_load` (satir 66-77)

`blob` **haric** tum sutunlari bir kere okur ve `nt`/`gap` matrislerini kurar.
Ikinci cagride hicbir sey yapmaz.

### 3.3 `_rebuild` — Geri Kurma (satir 79-97)

Paketlemenin tersini yapar: **dusurulen sutunu satir toplamindan geri koyar.**

```
missing = (100000 + sum_adj) − (kalan sutunlarin toplami)
```

- **Satir 83-86:** `dropped` ve `v*` sutunlari okunur.
- **Satir 87-88:** `sum_adj` yoksa sifir kabul edilir.
- **Satir 90:** Eksik degerin **tek hamlede tum satirlar icin** hesabi.
- **Satir 91-96:** Her olasi argmax sutunu icin maskeyle: kalanlar yerlerine,
  eksik olan da argmax pozisyonuna yazilir.

Sonuc, orijinal `.state` degerlerinin **birebir aynisi** (100000 ile carpilmis
tam sayi olarak).

### 3.4 Konumlandirma ve kayit (satir 99-119)

- **`_row_index` (satir 99-104):** `(site − 1) * n_nodes + node_pos`.
  Paketlemedeki site-major dizilimin tersi. Gecersiz node/site icin `None`.
- **`_record` (satir 106-113):** Kaydi `(state, p0, p1, ...)` demeti olarak
  dondurur. `state` once **istisna listesinden** aranir (orn. `-`), yoksa
  argmax'tan alfabe uzerinden turetilir. Olasiliklar 100000'e bolunup float'a
  cevrilir.
- **`_blob` (satir 115-119):** `blob` sutununun ilgili satirini okuyup xz'i
  cozer.

### 3.5 Genel API (satir 121-165)

| Metot | Ne yapar |
|---|---|
| `get_nt(node, site)` | Sadece nukleotid kaydi |
| `get_gap(node, site)` | Sadece gap/indel kaydi |
| `get(node, site)` | Ikisi birden, `{"nt": ..., "gap": ...}` |
| `tree()` | Newick metni (str) |
| `fasta()` | `{tur_adi: hizalanmis_dizi}` — bir kere acilip saklanir |
| `species` | Yaprak adlari listesi (FASTA sirasi) |
| `column(site)` | Hizalamanin bir sutunu: `{tur_adi: harf}` |

- **Satir 142-157 `fasta()`:** xz cozulur, FASTA satir satir ayristirilir
  (`>` ile baslayan satir yeni tur), sonuc **onbellege alinir**.

### 3.6 `main` — CLI (satir 167-202)

Dort mod:
- **Argumansiz:** ozet (dosya boyutu, node sayisi, site/node, tur sayisi,
  alfabeler, agac uzunlugu).
- **`--tree`:** Newick'i basar.
- **`--fasta`:** FASTA'yi 60 karakterlik satirlar halinde basar.
- **`<node> <site>`:** tek kaydi nt ve gap olarak sekme ayrilmis basar.

- **Satir 205-210:** `BrokenPipeError` yakalanir. `python ossuary_reader.py x
  --fasta | head` gibi bir kullanimda `head` cikinca Python cirkin bir yigin
  izi basardi; burada stdout `/dev/null`'a yonlendirilip sessizce cikiliyor.

---

## 4. `reborn.py` — Geri Dirilten

`ossuary_reader.py` bundle'i **sorgulamak** icin (tek tek kayit, agac, hizalama).
`reborn.py` ise bundle'i **dort orijinal dosyaya geri acmak** icin. Ciktisi
byte-byte girdinin aynisi.

```bash
python reborn.py bundle.parquet [cikti_klasoru]
python reborn.py bundle.parquet --only state    # state | fasta | tree | all
python reborn.py bundle.parquet --list          # icerigi yaz, dosya yazma
```

Verilen `cikti_klasoru` icine PHACT'in orijinal sonuc-klasoru duzeni kuruluyor —
`ossuary.py`'nin `discover()` fonksiyonunun taradigi duzenin ayni:

```
cikti_klasoru/
├── 1_preprocessed/<ad>.fasta
├── 2_iqtree_ancestral/<ad>.state, <ad>.treefile
└── 3_binary_iqtree_ancestral/<ad>_binary.state
```

Bu sayede `python ossuary.py cikti_klasoru/ yeni.parquet` calistirildiginda
girdileri **tekrar otomatik bulabiliyor** — dogrulamada da bu yapildi:
bundle → `reborn.py` ile diske ac → `ossuary.py` ile tekrar paketle → sonuc
ilk bundle ile **byte-byte ayni**.

### 4.1 Neyin saklandigi, neyin saklanmadigi

FASTA ve treefile zaten xz blob olarak duruyor, onlar bedava geri geliyor.
Iki `.state` dosyasi ise sutunlardan **yeniden uretiliyor**. IQ-TREE bu
dosyalarin basina iki farkli blok koyuyor:

1. **`#` ile baslayan yorum blogu** — dosya yolu, R/Excel'de nasil okunacagi,
   sutunlarin anlami gibi bilgiler. Dosyadan dosyaya degisiyor; bu veride
   `raw.state`'te 1 satir, `raw_binary.state`'te 8 satir.
2. **Header satiri** — `Node\tSite\tState\tp_A\tp_C\tp_G\tp_T` gibi sabit bir
   satir.

Bundle **yalnizca ikinciyi** geri uretiyor; yorum blogu hic saklanmiyor. Sebep:
header'in tum icerigi zaten bundle'da duran `stb_alphabet` / `bundle_alphabet_*`
metadata'sindan (`("A","C","G","T")`, `("0","1")`) **turetilebiliyor** —
`reborn.py` onu orijinal dosyadan bayt bayt kopyalamak yerine **yeniden
yaziyor**. Yorum blogunun ise turetilecek bir kaynagi yok (serbest metin),
saklamak sadece boyut ve karmasiklik ekler; PHACT/IQ-TREE'yi tekrar
calistirmak icin de zaten gerekmiyor.

Satir sonu uslubu (`raw.state` **CRLF**, `raw_binary.state` **LF** — ayni
klasordeki iki dosya, iki farkli uslup) hala saklaniyor, cunku veri satirlarinin
kendisini de etkiliyor:

- **`_detect_newline()` (ossuary.py satir 129-133):** Dosyanin ilk satirina
  bakip `\r\n` mi `\n` mi oldugunu tespit eder.
- Sonuc `bundle_newline_nt` / `bundle_newline_gap` metadata'sinda (`crlf` /
  `lf`); orijinal dosya adlari da `bundle_filenames`'te.

**Sonuc:** `.state` dosyalari orijinalle **veri satirlari ve header'da**
birebir ayni, ama IQ-TREE'nin yorum blogu **yeniden uretilmiyor** —
byte-byte tam eslesme artik yok, ama pratikte hicbir seyi kaybettirmiyor.

### 4.2 Yapisi (satir satir)

**`Bundle` sinifi (satir 39-97).** Bundle'i acar; metadata hemen, agir sutunlar
istendiginde okunur.
- **Satir 45-47:** Dosyanin gercekten bir ossuary bundle'i oldugu dogrulanir.
- **Satir 59-63:** `bundle_newline_*`, `bundle_filenames` ve `bundle_scale`
  **`md.get()` ile** okunuyor — bu alanlar eklenmeden once uretilmis bundle'lar da
  acilabilsin diye (`bundle_scale` yoksa 100000 varsayilir).
- **`blob()` (satir 65-71):** `blob` sutununu **tek seferde** okuyup dordunu de
  cozer ve saklar. Bilinmeyen ad icin `None` doner (eski bundle'lar).
- **`values()` (satir 77-97):** Paketlemenin tersi. Dusurulen sutun satir
  toplamindan geri konur:
  `missing = (100000 + sum_adj) − (kalan sutunlarin toplami)`.
  Sadece istenen tarafin (`_nt` ya da `_gap`) sutunlarini okur.

**`_decimal_table()` (satir 100-104).** Olasiliklari metne cevirmenin hizli yolu.
`float` uzerinden gitmek hem yavas hem yuvarlama riski; onun yerine
`0..max` arasindaki **her tam sayi icin** metin karsiligi bir kere uretilip
(`f"{v//100000}.{v%100000:05d}"`) diziye konuyor. Sonrasi saf indeksleme.
Boylece `0.99751` gibi degerler **kesinlikle** dogru basiliyor.

**`_state_labels()` (satir 107-114).** `State` sutunu argmax'tan turetilir;
alfabede olmayan istisnalar (orn. `-`) uzerine yazilir.

**`write_state()` (satir 117-136).** Bir `.state` dosyasini yazar.
- **Satir 123:** Header satiri `bundle_alphabet_*`'ten yeniden kuruluyor —
  orijinal dosyadan hicbir bayt kopyalanmiyor (bkz. bolum 4.1).
- **Satir 130-136:** Depolama **site-major**, cikti ise **node-major** olmali.
  Cozum tek satirlik: `vals[pos::n_nodes]` dilimi o node'un tum site'larini
  dogru sirada veriyor. Her node icin 10 000 satir tek `zip` + `"\t".join`
  ile kuruluyor, tek `write` ile diske gidiyor.

**`SUBDIRS` (satir 142-144).** Rol → alt klasor eslemesi:
`fasta → 1_preprocessed`, `nt/tree → 2_iqtree_ancestral`,
`gap → 3_binary_iqtree_ancestral`. `ossuary.py`'nin `_select()` fonksiyonundaki
"binary gecen .state → gap" kuralinin ters yonu.

**`restore()` (satir 147-174).** Istenen dosyalari `SUBDIRS`'e gore dogru alt
klasore yazar (satir 151-154: kullanilacak klasorler once olusturuluyor). Ad
once `bundle_filenames`'ten alinir, yoksa bundle dosya adindan turetilir.

**`main()` (satir 173-207).** `--list` bundle icerigini ozetler (node×site,
alfabeler, dosya adlari, satir sonu uslubu); arguman yoksa dort dosyayi da
yazar.

### 4.3 Dogrulama

`chr1_100001-110000` verisinde tam tur atildi — paketle, sonra geri ac:

| Dosya | Karsilastirma | Sonuc |
|---|---|---|
| `raw.fasta` | tum dosya | byte-byte ayni |
| `raw.treefile` | tum dosya | byte-byte ayni |
| `raw.state` | header + veri satirlari (yorum blogu haric) | byte-byte ayni |
| `raw_binary.state` | header + veri satirlari (yorum blogu haric) | byte-byte ayni |

Uretilen header satirlari da orijinaliyle birebir ayni cikti:
`Node\tSite\tState\tp_A\tp_C\tp_G\tp_T` (nt) ve `Node\tSite\tState\tp_0\tp_1`
(gap) — cunku IQ-TREE'nin kullandigi kalip zaten sabit.

Geri acma suresi ~4 saniye.

---

## 5. Ozet Akis

**Paketleme (`ossuary.py`):**
```
.state (nt) ─┐
.state (gap)─┼→ oku → site-major → argmax dusur → uint16 ─┐
             │   └→ satir sonu tespiti ────────────────┐   ├→ Parquet
.fasta ──────┼→ xz ─────────────────────────────────┐ │   │  (tek dosya)
.treefile ───┘                                      └─┴───┴→ blob sutunu
```

**Sorgulama (`ossuary_reader.py`):**
```
Parquet → metadata (node/site/alfabe) ─┐
       → v* + dropped + sum_adj ────────┼→ olasilik matrisi → get(node, site)
       → blob ─────────────────────────→ xz coz → FASTA / Newick
```

**Geri acma (`reborn.py`):**
```
Parquet → v* + dropped + sum_adj → olasilik matrisi ─┐
       → metadata → alfabe → yeni header satiri ─────┼→ .state x2 (node-major)
       → blob[0],[1] ────────────────────────────────→ .fasta + .treefile
```
