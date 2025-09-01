import pandas as pd
import matplotlib.pyplot as plt

# Lê o CSV
df = pd.read_csv("mc_results.csv")

# Mostra as primeiras linhas
print(df.head())

# Exemplo: histograma da coluna 'uy'
plt.figure(figsize=(8,5))
plt.hist(df["uy"], bins=20, color="skyblue", edgecolor="black")
plt.xlabel("u_y")
plt.ylabel("Frequência")
plt.title("Histograma das amostras de u_y")
plt.grid(True)
plt.tight_layout()
plt.show()

