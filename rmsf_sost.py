import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

df1  = pd.read_csv("./05_rmsf_cluster_1_prod1_onlysost.dat", delim_whitespace=True, header=0)
df2  = pd.read_csv("./05_rmsf_cluster_1_prod1_onlysost.dat", delim_whitespace=True, header=0)
df3  = pd.read_csv("./05_rmsf_cluster_1_prod1_onlysost.dat", delim_whitespace=True, header=0)
rmsf = pd.DataFrame({
    'res': df1['#Res'],
    'r1': df1.iloc[:, 1],
    'r2': df2.iloc[:, 1],
    'r3': df3.iloc[:, 1],
}).dropna()

# Renumera resíduos para começar em 1
rmsf['res'] = rmsf['res'] - rmsf['res'].min() + 1

# Calcula a média e desvio padrão
rmsf['mean'] = rmsf[['r1','r2','r3']].mean(axis=1)
rmsf['std']  = rmsf[['r1','r2','r3']].std(axis=1)

da1 = pd.read_csv("../../../../3_dynamics/sost/04_prod/analysis/analise_dissertação/full_data/05_rmsf_sost_full_equ2.dat", delim_whitespace=True, header=0)
da2 = pd.read_csv("../../../../3_dynamics/sost/04_prod/analysis/analise_dissertação/full_data/05_rmsf_sost_full_prod2.dat", delim_whitespace=True, header=0)
da3 = pd.read_csv("../../../../3_dynamics/sost/04_prod/analysis/analise_dissertação/full_data/05_rmsf_sost_full_prod3.dat", delim_whitespace=True, header=0)
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

ax.axvspan(2,  25, alpha=0.15, color='#333384', label='Loop 1')
ax.axvspan(31, 54, alpha=0.15, color='#aa2424', label='Loop 2')
ax.axvspan(56, 85, alpha=0.15, color='#267326', label='Loop 3')

ax.xaxis.set_major_locator(MultipleLocator(10))
ax.set_xlim([1, 89])
ax.set_ylim([0, 10])
ax.grid(True, linestyle="--", alpha=0.6)

plt.title('Comparação de RMSF entre SOST livre e SOST Ligada do complexo I', fontsize=8)
plt.xlabel('Posição do Resíduo')
plt.ylabel('RMSF (Å)')
plt.legend(fontsize=8)

ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)

plt.tight_layout()
plt.savefig('./05_rmsf_cluster_1_prod1_onlysost.png', dpi=300)
plt.show()
