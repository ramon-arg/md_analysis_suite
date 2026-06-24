import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

df1  = pd.read_csv("./05_rmsf_cluster_1_prod1_onlyromo.dat", delim_whitespace=True, header=0)
df2  = pd.read_csv("./05_rmsf_cluster_1_prod1_onlyromo.dat", delim_whitespace=True, header=0)
df3  = pd.read_csv("./05_rmsf_cluster_1_prod1_onlyromo.dat", delim_whitespace=True, header=0)
rmsf = pd.DataFrame({
    'res': df1['#Res'],
    'r1': df1.iloc[:, 1],
    'r2': df2.iloc[:, 1],
    'r3': df3.iloc[:, 1],
}); rmsf = rmsf.dropna()

# Calcula a média e desvio padrão
rmsf['mean'] = rmsf[['r1','r2','r3']].mean(axis=1)
rmsf['std']  = rmsf[['r1','r2','r3']].std(axis=1)

da1 = pd.read_csv("../../../../3_dynamics/romosozumab/03_equ/analysis/analise_dissertação/full_data/05_rmsf_romo_full_equ2.dat", delim_whitespace=True, header=0)
da2 = pd.read_csv("../../../../3_dynamics/romosozumab/03_equ/analysis/analise_dissertação/full_data/05_rmsf_romo_full_equ2.dat", delim_whitespace=True, header=0)
da3 = pd.read_csv("../../../../3_dynamics/romosozumab/03_equ/analysis/analise_dissertação/full_data/05_rmsf_romo_full_equ2.dat", delim_whitespace=True, header=0)
rmsf2 = pd.DataFrame({
    'res': da1['#Res'],
    'r1': da1.iloc[:, 1],
    'r2': da2.iloc[:, 1],
    'r3': da3.iloc[:, 1],
}).dropna()

# Calcula a média e desvio padrão
rmsf2['mean'] = rmsf2[['r1','r2','r3']].mean(axis=1)
rmsf2['std']  = rmsf2[['r1','r2','r3']].std(axis=1)

# Faz o plot
fig, ax = plt.subplots(figsize=(5, 2))

ax.plot(rmsf["res"], rmsf["mean"], c='steelblue', linewidth=0.5, alpha=0.8, label="Ligado")
ax.fill_between(rmsf["res"], 
    rmsf["mean"] - rmsf["std"],
    rmsf["mean"] + rmsf["std"],
    color="steelblue",
    alpha=0.3)

ax.plot(rmsf2["res"], rmsf2["mean"], c='black', linewidth=0.5, alpha=0.8, label="Livre")
ax.fill_between(rmsf2["res"], 
    rmsf2["mean"] - rmsf2["std"],
    rmsf2["mean"] + rmsf2["std"],
    color="gray",
    alpha=0.3)

for start, end in [(26,33), (51,58), (97,112)]:
    ax.axvspan(start, end, alpha=0.15, color='#CC6600')
for start, end in [(150,156), (173,175), (212,220)]:
    ax.axvspan(start, end, alpha=0.15, color='#4B0082')

handles, labels = ax.get_legend_handles_labels()
handles.append(Patch(color='#CC6600', alpha=0.15))
labels.append('CDR-H')
handles.append(Patch(color='#4B0082', alpha=0.15))
labels.append('CDR-L')
ax.legend(handles, labels, fontsize=6)

ax.xaxis.set_major_locator(MultipleLocator(20))
ax.set_xlim([1, 230])
ax.set_ylim([0, 3])
ax.grid(True, linestyle="--", alpha=0.6)

plt.title('Comparação de RMSF entre ROMO livre e ROMO Ligado do complexo I', fontsize=8)
plt.xlabel('Posição do Resíduo')
plt.ylabel('RMSF (Å)')
plt.legend(fontsize=8)

ax.legend(handles, labels, fontsize=6, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)

plt.tight_layout()
plt.savefig('./05_rmsf_cluster_1_prod1_onlyromo.png', dpi=300)
plt.show()
