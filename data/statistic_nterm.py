import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

seq_list = []
with open('HuOR_E60SecondE60Res_OR_sequence-lengfilter280.fas') as temf:
    text = temf.read()
    text = text.replace('\r', '')
    seqs = text.split('>')[1:]
    for seq in seqs:
        lines = seq.split('\n')
        name = lines[0]
        aads = "".join(lines[1:])
        seq_list.append((name, aads))

iaad = 0
tm1 = 0
for i in range(len(seq_list[0][1])):
    aad = seq_list[0][1][i]
    if aad != "-": iaad += 1
    if iaad == 23: tm1 = i

len_list = []
for name, seq in seq_list:
    nterm = seq[:tm1].replace('-', '')
    len_list.append(len(nterm))
count = Counter(len_list)
total = sum(count.values())
count = sorted(count.items(), key=lambda x: x[1], reverse=True)
count = [(x, round(y*100.0/total, 2)) for x, y in count]
with open("N-term_length.tab", 'w') as outf:
    lines = ['Nterm_length\tpercentage(%)\n']
    lines += [str(x)+'\t'+str(y)+'\n' for x, y in count]
    outf.writelines(lines)

plt.hist(len_list, bins=102, density=True)
plt.text(17.5, 0.08, "19")
plt.plot([19]*100, np.linspace(0, 0.07, 100), color='green')
plt.text(25, 0.08, "26")
plt.plot([26]*100, np.linspace(0, 0.07, 100), color='green')
plt.text(14.5, 0.06, "16")
plt.plot([16]*100, np.linspace(0, 0.05, 100), color='yellow')
plt.text(35, 0.06, "36")
plt.plot([36]*100, np.linspace(0, 0.05, 100), color='yellow')
plt.text(10.5, 0.04, "12")
plt.plot([12]*100, np.linspace(0, 0.03, 100), color='red')
plt.text(47, 0.04, "48")
plt.plot([48]*100, np.linspace(0, 0.03, 100), color='red')
plt.title("N-terminal sequence length statistics")
plt.xlabel("Sequence Length")
plt.ylabel("Frequency")
plt.xlim(0, 60)
plt.savefig("N-term_length.png", format='png')
plt.show()
