import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

def lorentzian(x, x0, fwhm, intensity):
    return intensity * (fwhm**2 / 4) / ((x - x0)**2 + (fwhm**2 / 4))

def generate_multiplet(center, j_constants, nuclei_counts, fwhm, total_integration, resolution=0.001):
    # Start with a single peak at the center
    peaks = [(center, 1.0)]
    
    for j, n in zip(j_constants, nuclei_counts):
        new_peaks = []
        j_ppm = j / 400.0 # Assuming 400 MHz as per image
        for _ in range(n):
            current_peaks = []
            for pos, amp in peaks:
                current_peaks.append((pos + j_ppm/2, amp * 0.5))
                current_peaks.append((pos - j_ppm/2, amp * 0.5))
            # Merge overlapping positions
            merged = {}
            for p, a in current_peaks:
                p_round = round(p, 6)
                merged[p_round] = merged.get(p_round, 0) + a
            peaks = list(merged.items())
            
    # Normalize peaks to match integration
    current_sum = sum(amp for pos, amp in peaks)
    peaks = [(pos, (amp / current_sum) * total_integration) for pos, amp in peaks]
    return peaks

x = np.linspace(0, 5, 10000)
y = np.zeros_like(x)

# Signal Data from image
# C2-H: 3.72, sextet, 1H, J ~ 6
c2_peaks = generate_multiplet(3.72, [6.0], [5], 0.01, 1.0)

# O-H: ~2.2, broad singlet, 1H
# Using a much wider FWHM for "broad singlet"
oh_peaks = [(2.2, 1.0)] # Lorentzian will handle width

# C3-H2: 1.46, multiplet, 2H (Diastereotopic)
# Simplified as a complex multiplet for visualization
c3_peaks = generate_multiplet(1.46, [7.4, 6.0], [3, 1], 0.01, 2.0)

# C1-H3: 1.16, doublet, 3H, J = 6.2
c1_peaks = generate_multiplet(1.16, [6.2], [1], 0.01, 3.0)

# C4-H3: 0.94, triplet, 3H, J = 7.4
c4_peaks = generate_multiplet(0.94, [7.4], [2], 0.01, 3.0)

all_signals = [
    (c2_peaks, 0.01),
    (oh_peaks, 0.15), # Broad OH
    (c3_peaks, 0.01),
    (c1_peaks, 0.01),
    (c4_peaks, 0.01)
]

for peaks, fwhm in all_signals:
    for pos, amp in peaks:
        y += lorentzian(x, pos, fwhm, amp)

plt.figure(figsize=(12, 6))
plt.plot(x, y, color='black', linewidth=1)
plt.xlim(5, 0) # Downfield to upfield
plt.ylim(-0.1, max(y)*1.1)
plt.xlabel('$\delta$ (ppm)', fontsize=12)
plt.ylabel('Intensity (arbitrary units)', fontsize=12)
plt.title('$^1$H NMR Spectrum of 2-butanol (400 MHz, $CDCl_3$)', fontsize=14)
plt.grid(True, which='both', linestyle='--', alpha=0.5)

# Labeling
plt.annotate('C2-H', xy=(3.72, lorentzian(3.72, 3.72, 0.01, c2_peaks[0][1])), xytext=(3.72, max(y)*0.3), arrowprops=dict(arrowstyle='->'))
plt.annotate('O-H', xy=(2.2, lorentzian(2.2, 2.2, 0.15, 1.0)), xytext=(2.2, max(y)*0.4), arrowprops=dict(arrowstyle='->'))
plt.annotate('C3-H$_2$', xy=(1.46, 5), xytext=(1.46, max(y)*0.6), arrowprops=dict(arrowstyle='->'))
plt.annotate('C1-H$_3$', xy=(1.16, 15), xytext=(1.16, max(y)*0.8), arrowprops=dict(arrowstyle='->'))
plt.annotate('C4-H$_3$', xy=(0.94, 10), xytext=(0.94, max(y)*0.9), arrowprops=dict(arrowstyle='->'))

plt.tight_layout()
plt.savefig('nmr_spectrum.png')
plt.show()
