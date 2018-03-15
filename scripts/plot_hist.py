import sys
import numpy as np
import matplotlib.pyplot as plt


x = []
for ln in sys.stdin:
    x.append(int(ln.strip()))
n_b = 50
# x.remove(max(x))
bins = np.arange(min(x),max(x),n_b)
# the histogram of the data
counts, bins = np.histogram(x, bins=bins)
sum_counts = float(sum(counts))
i = 0
while i < len(counts):
    print str(bins[i])+"-"+str(bins[i+1])+"\t"+str(counts[i] / sum_counts * 100)+"%"
    i = i + 1
width = 1 * (bins[1] - bins[0]) / 2
center = bins[:-1] + bins[1:] / 2
plt.figure(figsize=(18,12))
plt.bar(center, counts, align='center',width = width)
for i, v in enumerate(counts):
    if v == max(counts):
        plt.text(i, v-1, str(v), color='b',fontweight='bold')
plt.xlim(min(x),max(x))
plt.tight_layout()
plt.show()


# Tweak spacing to prevent clipping of ylabel
plt.show()
