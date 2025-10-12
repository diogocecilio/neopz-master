# post.py
import pandas as pd
import matplotlib.pyplot as plt

def main():
    # Lê os dados do CSV
    df = pd.read_csv("loadsweep.csv")

    # Se houver coluna "ok", filtra apenas os convergidos
    if "ok" in df.columns:
        df = df[df["ok"] == 1]

    # Ordena por fator para curva mais suave
    df = df.sort_values("factor")

    # Plota curva Uy × fator
    plt.figure()
    plt.plot(df["factor"], df["uy"], marker="o")
    plt.xlabel("Uy no ponto")
    plt.ylabel("fator de carga")
    plt.title("Curva carga–deslocamento")
    plt.grid(True)
    plt.tight_layout()

    # Salva figura
    plt.savefig("loadsweep_plot.png", dpi=150)
    print("Figura salva em loadsweep_plot.png")

    # Mostra na tela
    plt.show()

if __name__ == "__main__":
    main()

