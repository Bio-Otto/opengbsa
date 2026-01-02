#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap


def plot_ss_data(dat_file, output_file):
    # Dosyayı oku
    with open(dat_file, 'r') as f:
        content = f.read()

    # Her satırı işle
    lines = [line.strip() for line in content.split('\n') if line.strip()]

    if not lines:
        print("Error: No data found in the file.")
        return

    print(f"Found {len(lines)} non-empty lines.")

    # Her satır bir kalıntı dizisini temsil eder
    ss_matrix = []
    for line in lines:
        # Her karakter bir zaman noktası için sekonder yapıyı temsil eder
        ss_matrix.append(list(line))

    # Numpy dizisine dönüştür
    ss_matrix = np.array(ss_matrix, dtype=str)

    print(f"Matrix shape: {ss_matrix.shape}")

    # SS kodlarını sayısal değerlere dönüştür
    num_matrix = np.zeros(ss_matrix.shape, dtype=int)

    ss_code_map = {
        '~': 0,  # Loop
        'H': 1,  # Alpha helix
        'E': 2,  # Beta strand
        'B': 3,  # Beta bridge
        'G': 4,  # 3-10 helix
        'I': 5,  # Pi helix
        'T': 6,  # Turn
        'S': 7,  # Bend
        'P': 8,  # Poly-proline II helix
        '=': 9,  # Break
        'C': 0,  # Coil (alternatif)
        ' ': 0  # Boşluk
    }

    # Matris değerlerini doldur
    for i in range(ss_matrix.shape[0]):
        for j in range(ss_matrix.shape[1]):
            num_matrix[i, j] = ss_code_map.get(ss_matrix[i, j], 0)

    # Renk haritası oluştur
    colors = [
        '#FFFFFF',  # Loop - beyaz
        '#FF0000',  # Alpha helix - kırmızı
        '#0000FF',  # Beta strand - mavi
        '#A0A0FF',  # Beta bridge - açık mavi
        '#FF00FF',  # 3-10 helix - magenta
        '#FFFF00',  # Pi helix - sarı
        '#00FFFF',  # Turn - cyan
        '#FFA000',  # Bend - turuncu
        '#A0FFA0',  # Poly-proline - açık yeşil
        '#808080'  # Break - gri
    ]
    cmap = ListedColormap(colors)

    # Görselleştirme
    plt.figure(figsize=(12, 8))
    plt.imshow(num_matrix, aspect='auto', cmap=cmap, interpolation='none', vmin=0, vmax=9)

    # Renk çubuğu
    cbar = plt.colorbar(ticks=np.arange(0, 10))
    cbar.set_ticklabels(['Loop (~)', 'Alpha helix (H)', 'Beta strand (E)', 'Beta bridge (B)',
                         '3-10 helix (G)', 'Pi helix (I)', 'Turn (T)', 'Bend (S)', 'PP-helix (P)', 'Break (=)'])

    # Eksen etiketleri
    plt.xlabel('Time Point')
    plt.ylabel('Residue')

    # X ekseni etiketlerini ayarla
    if num_matrix.shape[1] > 10:
        tick_indices = np.linspace(0, num_matrix.shape[1] - 1, 10, dtype=int)
        plt.xticks(tick_indices)
    else:
        plt.xticks(range(num_matrix.shape[1]))

    # Y ekseni etiketlerini ayarla
    if num_matrix.shape[0] > 20:
        tick_indices = np.linspace(0, num_matrix.shape[0] - 1, 20, dtype=int)
        plt.yticks(tick_indices)
    else:
        plt.yticks(range(num_matrix.shape[0]))

    plt.title('Protein Secondary Structure')
    plt.tight_layout()

    # Kaydet
    plt.savefig(output_file, dpi=300)
    print(f"Secondary structure plot saved to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 plot_ss.py input.dat output.png")
        sys.exit(1)

    dat_file = sys.argv[1]
    output_file = sys.argv[2]
    plot_ss_data(dat_file, output_file)