import os
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import pandas as pd

def count_variants(filepath):
    count = 0
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) >= 5:
                ref = fields[3]
                alt = fields[4]
                if (len(ref) == 1 and len(alt) == 1) or (len(ref) > 1 or len(alt) > 1):
                    count += 1
    return count

def estimate_tmb(variant_count, sequenced_mb=35):
    return variant_count / sequenced_mb

# 1. Klasör seçim
root = tk.Tk()
root.withdraw()
selected_folder = filedialog.askdirectory(title="Kök klasörü seçin")

# 2. Alt klasörlerde arama ve analiz
results = []
for subdir, _, files in os.walk(selected_folder):
    for file in files:
        if 'passing_filters' in file and file.endswith('.vcf'):
            filepath = os.path.join(subdir, file)
            variant_count = count_variants(filepath)
            tmb = estimate_tmb(variant_count)
            basename = file.split('_')[0]  # MP255 gibi örnek ismi dosya başından alıyor
            results.append({'Sample': basename, 'TMB': tmb})

# 3. Excel çıktısı (tmb_results.xlsx)
df = pd.DataFrame(results)
excel_path = os.path.join(selected_folder, "tmb_results.xlsx")
df.to_excel(excel_path, index=False)
print(f"Excel çıktı: {excel_path}")

# 4. Renkli yatay barplot (long list için ideal)
colors = ['royalblue' if r['TMB'] >= 10 else 'gray' for r in results]
plt.figure(figsize=(7, len(results) * 0.5 + 2))  # Liste uzunluğuna göre yükseklik
plt.barh(df['Sample'], df['TMB'], color=colors, edgecolor='black')
plt.axvline(10, color='red', linestyle='dashed', label='Klinik TMB Sınırı (10)')
plt.xlabel('TMB (mut/Mb)', fontsize=12)
plt.ylabel('Vaka', fontsize=12)
plt.title('Tümör Mutasyon Yükü (TMB) - Vakalar', fontsize=14)
plt.legend()
plt.tight_layout()
plt.show()
