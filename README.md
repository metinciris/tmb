# Tumor Mutational Burden (TMB) Bulk Analyzer

Bu Python aracı, iç içe klasörlerdeki “passing_filters.vcf” dosyalarını otomatik olarak bulur, her birini analiz eder ve **TMB** sonucunu hem renkli grafik hem de Excel tablosu olarak sunar.

## Özellikler

- Görsel klasör seçimi
- Otomatik alt klasör tarama
- Filtrelenmiş VCF dosyası analizleri
- Her hasta/vaka için TMB hesabı
- Renkli, cutoff çizgili yatay barplot (Uzun vaka listeleri için ideal)
- Excel (.xlsx) formatında otomatik çıktı

## Kullanım

1. Tüm `.vcf` dosyalarını ve scripti uygun klasöre yerleştirin.
2. Scripti çalıştırın:
    ```bash
    python tmb_bulk_analyzer.py
    ```
3. Kök klasörü görsel olarak seçin ve işlem tamamlandığında:
    - Excel dosyası otomatik olarak kaydedilir (tmb_results.xlsx)
    - Tüm hastaların TMB sonuçları renkli olarak görselleştirilir.

## Sonuç Ekranı

Aşağıda örnek analiz sonuç ekranı:
![Örnek analiz grafiği](https://raw.githubusercontent.com/metinciris/tmb/refs/heads/main/Figure_1.png)


## Gereksinimler

- Python 3.x
- pandas
- matplotlib
- tkinter

```bash
pip install matplotlib pandas
```

## Lisans

MIT

***

**Not:** “Figure_1.jpg” görselini `README.md` ile aynı klasöre koy ve yukarıdaki görsel referansını kullanarak doğrudan README’de paylaşabilirsin. Görselin ismini ve yolunu kendi tercihine göre değiştirebilirsin.
