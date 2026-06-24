import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

df1  = pd.read_csv("./05_rmsd_cluster_1_prod1_onlyromo.dat", delim_whitespace=True, header=0)
df2  = pd.read_csv("./05_rmsd_cluster_1_prod1_onlyromo.dat", delim_whitespace=True, header=0)
df3  = pd.read_csv("./05_rmsd_cluster_1_prod1_onlyromo.dat", delim_whitespace=True, header=0)
rmsd = pd.DataFrame({
    'time': df1['#Frame'],
    'r1': df1.iloc[:, 1],
    'r2': df2.iloc[:, 1],
    'r3': df3.iloc[:, 1],
}); rmsd = rmsd.dropna()

rmsd['time'] = rmsd['time'] / 100

# Calcula a média e desvio padrão
rmsd['mean'] = rmsd[['r1','r2','r3']].mean(axis=1)
rmsd['std']  = rmsd[['r1','r2','r3']].std(axis=1)

# Segundo conjunto de dados
da1 = pd.read_csv("../../../../3_dynamics/romosozumab/03_equ/analysis/analise_dissertação/full_data/05_rmsd_romo_full_equ2.dat", delim_whitespace=True, header=0)
da2 = pd.read_csv("../../../../3_dynamics/romosozumab/03_equ/analysis/analise_dissertação/full_data/05_rmsd_romo_full_equ2.dat", delim_whitespace=True, header=0)
da3 = pd.read_csv("../../../../3_dynamics/romosozumab/03_equ/analysis/analise_dissertação/full_data/05_rmsd_romo_full_equ2.dat", delim_whitespace=True, header=0)

rmsd2 = pd.DataFrame({
    'time': da1['#Frame'],
    'r1': da1.iloc[:, 1],
    'r2': da2.iloc[:, 1],
    'r3': da3.iloc[:, 1],
}).dropna()

rmsd2['time'] = rmsd2['time'] / 100

rmsd2['mean'] = rmsd2[['r1','r2','r3']].mean(axis=1)
rmsd2['std']  = rmsd2[['r1','r2','r3']].std(axis=1)

# Faz o plot
fig, ax = plt.subplots(figsize=(5, 2))

# Dataset 2
ax.plot(rmsd2["time"], rmsd2["mean"], c='black', linewidth=0.5, alpha=0.8, label="ROMO livre")
ax.fill_between(rmsd2["time"],
    rmsd2["mean"] - rmsd2["std"],
    rmsd2["mean"] + rmsd2["std"],
    color="gray", alpha=0.2)

# Dataset 1
ax.plot(rmsd["time"], rmsd["mean"], c='steelblue', linewidth=0.5, alpha=0.8, label="ROMO ligado")
ax.fill_between(rmsd["time"], 
    rmsd["mean"] - rmsd["std"],
    rmsd["mean"] + rmsd["std"],
    color="steelblue",
    alpha=0.3)

#ax.xaxis.set_major_locator(MultipleLocator(100))
ax.set_xlim([0, 500])
ax.set_ylim([0, 4])
ax.grid(True, linestyle="--", alpha=0.6)
plt.title('Valores de RMSD para Romosozumab do complexo I e para Romosozumab livre', fontsize='8')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (Å)')
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig('./05_rmsd_cluster_1_prod1_onlyromo_compared.png', dpi=300)
plt.show()
