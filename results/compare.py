import pandas as pd
import argparse # Adicionado para lidar com argumentos de linha de comando
from pathlib import Path
# Adicionado para significância estatística
import numpy as np
from scipy.stats import mannwhitneyu # Alterado de wilcoxon para mannwhitneyu

# Helper function to parse Func and Dim from the first column
def parse_func_dim_from_first_column(df):
    """
    Tenta parsear as colunas 'Func' e 'Dim' a partir da primeira coluna do DataFrame
    se elas já não existirem. Espera um formato como 'F-D' (ex: '1-10').
    """
    if 'Func' in df.columns and 'Dim' in df.columns:
        try:
            df['Dim'] = pd.to_numeric(df['Dim']) # Garante que Dim seja numérica
            if df['Dim'].isnull().any():
                print(f"Aviso: Valores nulos encontrados na coluna 'Dim' existente após conversão para numérico. Removendo linhas com Dim nula.")
                df.dropna(subset=['Dim'], inplace=True)
            df['Dim'] = df['Dim'].astype(int)
        except ValueError:
            print(f"Erro: Não foi possível converter a coluna 'Dim' existente para numérico/inteiro. Verifique os dados.")
            return None # Indica falha
        return df

    if not df.columns.empty:
        first_col_name = df.columns[0]
        # print(f"Info: Tentando parsear Func/Dim da primeira coluna: '{first_col_name}'")
        try:
            split_data = df[first_col_name].astype(str).str.split('-', n=1, expand=True)
            if split_data.shape[1] == 2:
                df['Func'] = split_data[0]
                df['Dim'] = pd.to_numeric(split_data[1])
                df.dropna(subset=['Dim'], inplace=True) # Remove linhas onde Dim não pôde ser convertida
                df['Dim'] = df['Dim'].astype(int)
                df.drop(columns=[first_col_name], inplace=True)
            else:
                print(f"Erro: A primeira coluna '{first_col_name}' não pôde ser dividida em duas partes com '-'.")
                return None # Indica falha
        except Exception as e:
            print(f"Erro ao parsear Func/Dim da coluna '{first_col_name}': {e}")
            return None # Indica falha
    else:
        print("Erro: DataFrame vazio, não é possível parsear Func/Dim.")
        return None # Indica falha

    if 'Func' not in df.columns or 'Dim' not in df.columns:
        print(f"Alerta: Colunas 'Func' e/ou 'Dim' não puderam ser criadas. A comparação pode falhar.")
        return None # Indica falha
    return df

# Main comparison function
def comparar_csvs_por_instancia(arquivo1_path, arquivo2_path, alg1_name="Algoritmo 1", alg2_name="Algoritmo 2", dim_filter=None):
    # Lê os arquivos com separador ;
    try:
        df1_raw = pd.read_csv(arquivo1_path, sep=';')
        df2_raw = pd.read_csv(arquivo2_path, sep=';')
    except FileNotFoundError as e:
        print(f"Erro ao ler arquivo: {e}")
        return
    except pd.errors.EmptyDataError:
        print(f"Erro: Um dos arquivos CSV está vazio ou mal formatado ({arquivo1_path} ou {arquivo2_path}).")
        return
    except Exception as e:
        print(f"Erro inesperado ao ler CSVs: {e}")
        return

    # Parseia Func e Dim
    df1 = parse_func_dim_from_first_column(df1_raw.copy()) # Use .copy() to avoid modifying original raw df if parse fails mid-way
    df2 = parse_func_dim_from_first_column(df2_raw.copy())

    if df1 is None or df2 is None:
        print("Falha ao parsear Func/Dim de um ou ambos os arquivos. Abortando comparação.")
        return

    # Filtra por dimensão, se especificado
    if dim_filter is not None:
        if 'Dim' not in df1.columns or 'Dim' not in df2.columns: # Should not happen if parse_func_dim_from_first_column worked
            print(f"Erro: Coluna 'Dim' não encontrada em um dos arquivos para aplicar o filtro de dimensão.")
            return
        # 'Dim' já deve ser int devido ao parse_func_dim_from_first_column
        df1 = df1[df1['Dim'] == dim_filter]
        df2 = df2[df2['Dim'] == dim_filter]
    
    # Faz o merge com base nas colunas 'Func' e 'Dim', mantendo apenas instâncias comuns
    try:
        df_merged = pd.merge(df1, df2, on=['Func', 'Dim'], suffixes=('_alg1', '_alg2'), how='inner')
    except KeyError as e: # Should not happen if parse_func_dim_from_first_column worked
        print(f"Erro no merge: Coluna 'Func' ou 'Dim' não encontrada em um dos arquivos. Detalhe: {e}")
        print(f"  Colunas em {arquivo1_path} (após parse): {df1.columns.tolist()}")
        print(f"  Colunas em {arquivo2_path} (após parse): {df2.columns.tolist()}")
        return
    
    if df_merged.empty:
        if dim_filter is not None:
            print(f"Nenhuma instância (Func, Dim={dim_filter}) em comum encontrada entre os dois arquivos CSV após filtragem.")
        else:
            print("Nenhuma instância (Func, Dim) em comum encontrada entre os dois arquivos CSV.")
        return
    
    dim_info = f" para Dim={dim_filter}" if dim_filter is not None else ""
    print(f"Comparando {len(df_merged)} instâncias (Func, Dim{dim_info}) em comum.")

    vitorias_media_alg1 = 0
    vitorias_media_alg2 = 0
    vitorias_significativas_alg1 = 0
    vitorias_significativas_alg2 = 0

    # Identifica as colunas de "run" originais do df1 (L-SHADE base)
    # Assume que são todas as colunas exceto 'Func' e 'Dim' (que foram criadas/validadas)
    # Usamos as colunas de df1 ANTES do merge para identificar as colunas de run
    original_run_cols_from_df1 = [col for col in df1.columns if col not in ['Func', 'Dim']]


    if not original_run_cols_from_df1:
        print(f"Erro: Nenhuma coluna de 'run' (dados de execução) encontrada no arquivo base {arquivo1_path} (após parse).")
        print(f"  Colunas encontradas em df1 (após parse): {df1.columns.tolist()}")
        return

    # Constrói os nomes das colunas de run como apareceriam em df_merged
    # Se df1 e df2 tinham originalmente os mesmos nomes de colunas de run (ex: s1, s2...),
    # então no df_merged elas serão sufixadas com _alg1 e _alg2.
    # As colunas de run para df1 serão as originais com sufixo _alg1
    run_cols_alg1_merged = [f"{col}_alg1" for col in original_run_cols_from_df1]
    run_cols_alg2_merged = [f"{col}_alg2" for col in original_run_cols_from_df1] # Assume que df2 tem as mesmas colunas de run
    
    # Verifica se todas as colunas de run esperadas existem no df_merged
    missing_cols_alg1 = [col for col in run_cols_alg1_merged if col not in df_merged.columns]
    missing_cols_alg2 = [col for col in run_cols_alg2_merged if col not in df_merged.columns]

    if missing_cols_alg1:
        print(f"Erro: Colunas de run ausentes para Algoritmo 1 no DataFrame mesclado: {missing_cols_alg1}")
        print(f"  Verifique se os arquivos CSV originais têm o mesmo número e nomes de colunas de dados de execução.")
        return
    if missing_cols_alg2:
        print(f"Erro: Colunas de run ausentes para Algoritmo 2 no DataFrame mesclado: {missing_cols_alg2}")
        print(f"  Verifique se os arquivos CSV originais têm o mesmo número e nomes de colunas de dados de execução.")
        return


    for index, row in df_merged.iterrows():
        func_name = str(row['Func'])
        dim_val = str(row['Dim'])

        try:
            # Extrai os valores das execuções para cada algoritmo, convertendo para float
            runs_alg1_values = row[run_cols_alg1_merged].values.astype(float)
            runs_alg2_values = row[run_cols_alg2_merged].values.astype(float)
        except KeyError as e:
            print(f"  Aviso: Coluna de run esperada não encontrada para Func={func_name}, Dim={dim_val}. Detalhe: {e}. Ignorando esta instância.")
            continue
        except ValueError as e: # Erro na conversão para float
            print(f"  Aviso: Erro ao converter dados de run para numérico para Func={func_name}, Dim={dim_val}. Detalhe: {e}. Ignorando esta instância.")
            continue

        # Remove NaNs que podem ter surgido de colunas vazias ou erros de conversão
        runs_alg1_values = runs_alg1_values[~np.isnan(runs_alg1_values)]
        runs_alg2_values = runs_alg2_values[~np.isnan(runs_alg2_values)]

        if len(runs_alg1_values) == 0 or len(runs_alg2_values) == 0:
            print(f"  Aviso: Dados de execução insuficientes para Func={func_name}, Dim={dim_val} após remover NaNs. Ignorando.")
            continue
        
        media_alg1 = np.mean(runs_alg1_values)
        media_alg2 = np.mean(runs_alg2_values)

        # Contabiliza vitória pela média
        if media_alg1 < media_alg2:
            vitorias_media_alg1 += 1
        elif media_alg2 < media_alg1:
            vitorias_media_alg2 += 1
        
        # Teste de significância estatística (Mann-Whitney U) para esta instância
        try:
            # Evita erro no teste se ambos os arrays forem constantes e idênticos
            # ou se uma das amostras tiver variação zero (o que pode causar problemas com mannwhitneyu)
            if (len(np.unique(runs_alg1_values)) <= 1 and len(np.unique(runs_alg2_values)) <= 1 and 
                np.array_equal(np.unique(runs_alg1_values), np.unique(runs_alg2_values))):
                p_value = 1.0 # Nenhuma diferença ou dados insuficientes para diferenciar
            elif len(runs_alg1_values) < 1 or len(runs_alg2_values) < 1: # Mann-Whitney U precisa de pelo menos 1 em cada
                 p_value = 1.0
            else:
                # alternative='two-sided' é o padrão, mas pode ser explícito
                # Usar 'less' ou 'greater' se tiver uma hipótese direcional ANTES de ver os dados
                stat, p_value = mannwhitneyu(runs_alg1_values, runs_alg2_values, alternative='two-sided', use_continuity=True)
            
            if p_value < 0.05:
                if media_alg1 < media_alg2: # alg1 é significativamente melhor
                    vitorias_significativas_alg1 += 1
                elif media_alg2 < media_alg1: # alg2 é significativamente melhor
                    vitorias_significativas_alg2 += 1
        except ValueError as e:
            # print(f"  Aviso: Teste de Mann-Whitney U não pôde ser realizado para Func={func_name}, Dim={dim_val}. Erro: {e}")
            pass # Não contabiliza como vitória significativa

    print(f"\nVitórias de {alg1_name} (Total de instâncias (Instâncias com significância)):")
    print(f"  Performance: {vitorias_media_alg1} ({vitorias_significativas_alg1})")

    print(f"\nVitórias de {alg2_name} (Total de instâncias (Instâncias com significância)):")
    print(f"  Performance: {vitorias_media_alg2} ({vitorias_significativas_alg2})")


if __name__ == "__main__":
    DEFAULT_SOURCE_DIR = str(Path("/home/raphael/Code/FP-MAX-LSHADE/results/stats_test_data").resolve())
    LSHADE_BASE_FILENAME = "out(--lshade --sf ).csv" # Nome do arquivo L-SHADE base dentro do source_dir

    parser = argparse.ArgumentParser(
        description=f"Compara resultados de LSHADE (base) com FIM-LSHADE (variantes) a partir de arquivos CSV "
                    f"contendo múltiplas execuções por instância. O diretório padrão é '{DEFAULT_SOURCE_DIR}'. "
                    f"O arquivo L-SHADE base deve se chamar '{LSHADE_BASE_FILENAME}' e estar neste diretório."
    )
    parser.add_argument('--source_dir', type=str, default=DEFAULT_SOURCE_DIR,
                        help=f"Diretório contendo os arquivos CSV para comparação. Padrão: {DEFAULT_SOURCE_DIR}")
    parser.add_argument('--lshade_base_name', type=str, default=LSHADE_BASE_FILENAME,
                        help=f"Nome do arquivo L-SHADE base dentro do source_dir. Padrão: {LSHADE_BASE_FILENAME}")
    parser.add_argument('--d', '--dim', type=int, default=None,
                        help="Dimensão específica para filtrar os resultados (ex: 10, 30, 50). Se não especificado, todas as dimensões são consideradas.")

    args = parser.parse_args()

    source_dir_path = Path(args.source_dir)
    lshade_base_file_path = source_dir_path / args.lshade_base_name

    if not source_dir_path.is_dir():
        print(f"Erro: O diretório fonte '{source_dir_path}' não foi encontrado ou não é um diretório.")
        exit(1)

    if not lshade_base_file_path.is_file():
        print(f"Erro: Arquivo L-SHADE base '{args.lshade_base_name}' não encontrado em '{source_dir_path}'.")
        exit(1)

    print(f"Usando L-SHADE base: {lshade_base_file_path}")
    alg1_display_name = "L-SHADE (base)"

    # Encontra outros arquivos CSV no diretório para comparar com o base
    other_csv_files = [
        f for f in source_dir_path.iterdir()
        if f.is_file() and f.suffix == '.csv' and f.resolve() != lshade_base_file_path.resolve()
    ]

    if not other_csv_files:
        print(f"Nenhum outro arquivo CSV encontrado em '{source_dir_path}' para comparar com o L-SHADE base.")
        exit(0) 
    
    # Seleciona apenas o arquivo mais recente
    arquivo_mais_recente_para_comparar = max(other_csv_files, key=lambda f: f.stat().st_mtime)
    print(f"Arquivo mais recente selecionado para comparação: {arquivo_mais_recente_para_comparar.name}")

    print(f"\n--- Comparando com: {arquivo_mais_recente_para_comparar.name} ---")
    
    alg2_display_name = ""
    filename_str = arquivo_mais_recente_para_comparar.name
    # Tenta extrair configuração do nome do arquivo, ex: "out(uv=10,fim=0.5).csv"
    start_paren_idx = filename_str.find('(')
    end_paren_idx = filename_str.rfind(')')

    if start_paren_idx != -1 and end_paren_idx != -1 and start_paren_idx < end_paren_idx:
        config_params = filename_str[start_paren_idx + 1 : end_paren_idx]
        alg2_display_name = f"FIM-L-SHADE ({config_params})"
    else:
        # Fallback: usa o nome do arquivo (sem extensão) se o padrão não for encontrado
        alg2_display_name = f"FIM-L-SHADE ({arquivo_mais_recente_para_comparar.stem})"
    
    comparar_csvs_por_instancia(str(lshade_base_file_path), str(arquivo_mais_recente_para_comparar), 
                                alg1_display_name, alg2_display_name, dim_filter=args.d)
