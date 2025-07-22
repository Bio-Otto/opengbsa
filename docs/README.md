# OpenGBSA Documentation

Hoş geldiniz! Bu dökümantasyon, OpenGBSA'nın tüm özelliklerini, kurulumunu, kullanımını ve geliştirici rehberini profesyonel ve modüler şekilde sunar.

---

## 📚 İçindekiler
- [Genel Bakış](#genel-bakış)
- [Kurulum](#kurulum)
- [Temel Kullanım](#temel-kullanim)
- [Gelişmiş Analizler](#gelismis-analizler)
- [YAML Konfigürasyon Rehberi](#yaml-konfigürasyon-rehberi)
- [Çıktı Dosyaları ve Formatları](#çıktı-dosyaları-ve-formatları)
- [Sık Karşılaşılan Sorunlar (Troubleshooting)](#sik-karsilasilan-sorunlar)
- [Katkı ve Geliştirici Rehberi](#katki-ve-gelistirici-rehberi)

---

## Genel Bakış
OpenGBSA, protein-ligand bağlanma serbest enerjisi hesaplamaları için gelişmiş MM/GBSA analizleri sunan, YAML tabanlı konfigürasyon ve kapsamlı raporlama özelliklerine sahip bir Python paketidir.

- Çoklu GB modelleri (OBC2, OBC1, HCT, GBn, GBn2)
- Entropi analizi (Normal Mode Analysis)
- Per-residue decomposition (ayrıntılı kalıntı katkısı)
- CUDA/CPU desteği, paralel analiz
- Otomatik validasyon ve raporlama

---

## Kurulum
**En güvenli ve önerilen yol:**
```bash
conda create -n opengbsa -c bio-otto opengbsa
conda activate opengbsa
```
> Python >=3.8,<3.13, numpy <2.0, openmm >=8.0.0,<8.3 ile tam uyumludur.

Ortamı temizlemek için:
```bash
conda deactivate
conda remove -n opengbsa --all
```

---

## Temel Kullanım
1. **Konfigürasyon dosyası oluşturun:**
   ```bash
   opengbsa --create-config
   # veya
   python mmgbsa_cli.py --create-config
   ```
2. **Analizi başlatın:**
   ```bash
   opengbsa mmgbsa_config.yaml
   # veya
   python mmgbsa_cli.py mmgbsa_config.yaml
   ```
3. **Sonuçları görüntüleyin:**
   - `mmgbsa_results/analysis_YYYYMMDD_HHMMSS/` klasöründe tüm rapor ve çıktılar bulunur.

---

## Gelişmiş Analizler
- **Frame aralığı seçimi:**
  ```yaml
  analysis_settings:
    frame_start: 200
    frame_end: 1000
    frame_stride: 1
    # max_frames ve decomp_frames otomatik hesaplanır
  ```
- **Per-residue decomposition ve frame-by-frame:**
  ```yaml
  analysis_settings:
    run_per_residue_decomposition: true
    save_frame_by_frame_csv: true
    frame_output_format: "csv"
  ```
- **Entropi analizi:**
  ```yaml
  analysis_settings:
    run_entropy_analysis: true
  ```

---

## YAML Konfigürasyon Rehberi
- Tüm parametreler ve örnekler için [../config/COMPLETE_CONFIG_GUIDE.md](../config/COMPLETE_CONFIG_GUIDE.md) dosyasına bakınız.
- Sık kullanılan parametreler:
  - `input_files`, `analysis_settings`, `output_settings`, `advanced_settings`, `platform_settings`

---

## Çıktı Dosyaları ve Formatları
- `final_report.txt`: Kapsamlı analiz raporu
- `results_summary.yaml`: Sonuç özeti
- `per_residue_detailed.csv`: Kalıntı bazında enerji katkıları
- `frame_by_frame_decomposition.csv`: (isteğe bağlı) Her frame için decomposition
- `binding_hot_spots.csv`: En güçlü bağlanma katkısı yapan kalıntılar
- `decomposition_summary.csv`: Özet istatistikler

---

## Sık Karşılaşılan Sorunlar
- **No module named 'numpy.compat'**: Ortamda numpy 2.x varsa, ortamı silip yeniden oluşturun.
- **openmm/cudatoolkit uyumsuzluğu**: openmm 8.3+ yüklenirse, ortamı silip yeniden oluşturun.
- **pip ile kurulumda dependency hatası**: conda ile kurulum yapın.
- **CUDA bulunamıyor**: CUDA destekli ortamda çalıştığınızdan emin olun veya CPU platformunu seçin.

---

## Katkı ve Geliştirici Rehberi
- Pull request ve issue açabilirsiniz.
- Kod stili, testler ve dökümantasyon için [CONTRIBUTING.md](../CONTRIBUTING.md) dosyasına bakınız.
- Geliştirici ortamı için:
  ```bash
  git clone https://github.com/Bio-Otto/opengbsa.git
  cd opengbsa
  conda env create -f environment.yml  # (varsa)
  pip install -e .
  ```

---

**Not:** Daha fazla örnek, detay ve güncel bilgi için ana [README.md](../README.md) dosyasına ve config/ klasöründeki rehberlere bakınız. 