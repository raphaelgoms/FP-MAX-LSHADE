import os
import pandas as pd

def comparar_csvs_por_instancia(arquivo1, arquivo2):
    # Lê os arquivos com separador ;
    df1 = pd.read_csv(arquivo1, sep=';')
    df2 = pd.read_csv(arquivo2, sep=';')

    # Faz o merge com base nas colunas 'Func' e 'Dim'
    df_merged = pd.merge(df1, df2, on=['Func', 'Dim'], suffixes=('_alg1', '_alg2'))

    # Inicializa contadores de vitórias
    vitorias_alg1 = {'AvgFO': 0, 'BestFO': 0, 'Time': 0}
    vitorias_alg2 = {'AvgFO': 0, 'BestFO': 0, 'Time': 0}

    # Compara as métricas
    for index, row in df_merged.iterrows():
        for metrica in ['AvgFO', 'BestFO', 'Time']:
            v1 = row[f'{metrica}_alg1']
            v2 = row[f'{metrica}_alg2']

            if v1 < v2:
                vitorias_alg1[metrica] += 1
            elif v2 < v1:
                vitorias_alg2[metrica] += 1
            # empates são ignorados

    # Mostra resultados
    print("Vitórias do L-SHADE:")
    for metrica, v in vitorias_alg1.items():
        print(f"  {metrica}: {v}")

    print("\nVitórias do DM-L-SHADE:")
    for metrica, v in vitorias_alg2.items():
        print(f"  {metrica}: {v}")



from pathlib import Path

pasta = Path('results/uv=10')  # ou substitua por Path('caminho/da/pasta')
arquivos = [f for f in pasta.iterdir() if f.is_file()]

if arquivos:
    arquivo_mais_recente = max(arquivos, key=lambda f: f.stat().st_mtime)
else:
    print("Nenhum arquivo encontrado.")

LSHADE = "results/out(--lshade --sf ).csv"
DMLSHADE = arquivo_mais_recente

comparar_csvs_por_instancia(LSHADE, DMLSHADE)
