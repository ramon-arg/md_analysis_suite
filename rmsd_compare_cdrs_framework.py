import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

df1 = pd.read_csv("./cdrs_data/05_rmsd_romo_cluster_1_prod1_cdrs.dat", delim_whitespace=True, header=0)
df2 = pd.read_csv("./framework_data/05_rmsd_romo_cluster_1_prod1_framework.dat", delim_whitespace=True, header=0)

rmsd = pd.DataFrame({
    'time': df1['#Frame'],
    'r1': df1.iloc[:, 1],
}).dropna()
rmsd['time'] = rmsd['time'] / 100

rmsd2 = pd.DataFrame({
    'time': df2['#Frame'],
    'r1': df2.iloc[:, 1],
}).dropna()
rmsd2['time'] = rmsd2['time'] / 100

# Faz o plot
fig, ax = plt.subplots(figsize=(5, 2))

ax.plot(rmsd["time"], rmsd["r1"], color='#A0526A', linewidth=0.5, alpha=0.8, label="CDRs")
ax.plot(rmsd2["time"], rmsd2["r1"], color='#008B8B', linewidth=0.5, alpha=0.8, label="Framework")

#ax.xaxis.set_major_locator(MultipleLocator(100))
ax.set_xlim([0, 500])
ax.set_ylim([0, 4])
ax.grid(True, linestyle="--", alpha=0.6)
plt.title('Valores de RMSD de CDRs vs Framework para Romosozumab do complexo I', fontsize='8')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (Å)')
plt.legend(fontsize=8, handlelength=1.5, handleheight=0.7, labelspacing=0.3, borderpad=0.4)
plt.tight_layout()
plt.savefig('./05_rmsd_romo_cluster1_cdrsvsframework.png', dpi=300)
plt.show()
