import pandas as pd
import os
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from scipy.stats import chi2_contingency
import numpy as np


# --- Constantes ---
# NOVO GRUPO DE GENES: IL5, IL33, IL25, TSLP
GENES_INTERESSE = ['IL5', 'IL33', 'IL25', 'TSLP']
COL_TEMPO = 'OS_Time_nature2012'
COL_EVENTO = 'OS_event_nature2012'
COL_PAM50 = 'PAM50Call_RNAseq'

# --- Paleta de Cores Customizada para PAM50 ---
# Tons de Verde, Amarelo/Dourado, Vermelho e Marrom para melhor visualização
PAM50_COLORS = {
    "LumA": "#06B600",   # Verde
    "LumB": "#B9E90E",   # Amarelo Esverdeado
    "Her2": "#F59505",   # Laranja/Dourado
    "Basal": "#D60000",  # Vermelho
    "Normal": "#008080"  # Turquesa (se não for filtrado)
}

def formatar_pval(p_value):
    """Formata um p-valor para exibição no gráfico."""
    p_val_float = float(p_value) 
    if p_val_float < 0.001:
        return "p < 0.001"
    else:
        return f"p = {p_val_float:.3f}"

def definir_grupos_comparacao(df):
    """
    Cria as colunas de grupo específicas para cada comparação com base nos GENES_INTERESSE.
    """
    print("Definindo grupos de comparação...")
    
    # 1. Criar colunas booleanas para expressão (> 0)
    for gene in GENES_INTERESSE:
        if gene in df.columns:
            # Garante que a coluna de expressão seja booleana (True/False)
            df[f'{gene}_expresso'] = df[gene] > 0
        else:
            # Para a execução se algum dos genes essenciais estiver faltando
            print(f"Atenção: Gene {gene} não encontrado. Abortando a definição de grupos.")
            return None

    # 2. Definir o grupo complexo (4 Genes Expressos Simultaneamente)
    g_4_all = (
        df['IL5RA_expresso'] & 
        df['IL33_expresso'] & 
        df['IL25_expresso'] & 
        df['TSLP_expresso']
    )
    
    # 3. Criar as colunas de comparação para os gráficos
    
    # Comparação 1: (Todos os 4 Genes vs. Todos os Outros)
    df['grupo_4_genes_vs_outros'] = np.where(
        g_4_all, 
        'IL5RA_IL33_IL25_TSLP', # Todos expressos
        'Demais'                # Pelo menos um não expresso
    )
    
    # Comparações Individuais (Gene vs. Nao Expressa)
    df['grupo_IL5RA_vs_Nao'] = np.where(
        df['IL5RA_expresso'], 
        'Expressa IL5RA', 
        'Não Expressa IL5RA'
    )
    
    df['grupo_IL33_vs_Nao'] = np.where(
        df['IL33_expresso'], 
        'Expressa IL33', 
        'Não Expressa IL33'
    )
    
    df['grupo_IL25_vs_Nao'] = np.where(
        df['IL25_expresso'], 
        'Expressa IL25', 
        'Não Expressa IL25'
    )
    
    df['grupo_TSLP_vs_Nao'] = np.where(
        df['TSLP_expresso'], 
        'Expressa TSLP', 
        'Não Expressa TSLP'
    )
    
    return df

def plotar_sobrevida(df_dados, coluna_grupo, titulo, nome_arquivo):
    """
    Gera e salva um gráfico de sobrevida (Kaplan-Meier) para os grupos definidos.
    """
    print(f"Gerando gráfico de Sobrevida: {titulo}")
    
    df_plot = df_dados.dropna(subset=[COL_TEMPO, COL_EVENTO, coluna_grupo])
    # Assume-se que 'Ignorar' não existe, mas mantemos o filtro por segurança
    df_plot = df_plot[df_plot[coluna_grupo] != 'Ignorar']
    
    if df_plot.empty or df_plot[coluna_grupo].nunique() < 2:
        print(f"  -> Aviso: Dados insuficientes para '{titulo}' após filtragem.")
        return

    T = df_plot[COL_TEMPO] / 30.44 # Converte dias para meses
    E = df_plot[COL_EVENTO]
    grupos = df_plot[coluna_grupo]
    # np.sort garante ordem alfabética e facilita a reprodutibilidade do logrank_test
    grupos_unicos = np.sort(grupos.unique()) 

    # --- Análise Estatística (Log-Rank Test) ---
    p_valor = 0.99
    try:
        if len(grupos_unicos) == 2:
            # Teste Log-Rank padrão (logrank_test) para 2 grupos
            grupo_A = grupos_unicos[0]
            grupo_B = grupos_unicos[1]
            idx_A = (grupos == grupo_A)
            idx_B = (grupos == grupo_B)
            resultado_stats = logrank_test(T[idx_A], T[idx_B], E[idx_A], E[idx_B])
            p_valor = resultado_stats.p_value 
        
        elif len(grupos_unicos) > 2:
            # Teste Log-Rank multivariado (para 3+ grupos)
            resultado_stats = multivariate_logrank_test(T, grupos, E)
            p_valor = resultado_stats.p_value
            
    except Exception as e:
        print(f"  -> Erro no teste estatístico (logrank_test): {e}")

    # --- Plotagem (Kaplan-Meier) ---
    plt.figure(figsize=(10, 7))
    ax = plt.subplot(111)
    kmf = KaplanMeierFitter()

    for grupo in grupos_unicos:
        idx = (grupos == grupo)
        if idx.sum() > 0:
            kmf.fit(T[idx], E[idx], label=f'{grupo} (n={idx.sum()})')
            # Plotando a função de sobrevivência com intervalo de confiança
            kmf.plot_survival_function(ax=ax, ci_show=True)

    # Exibe o P-Valor no canto inferior esquerdo
    plt.text(0.05, 0.05, f"Log-Rank: {formatar_pval(p_valor)}", transform=ax.transAxes,
             fontsize=12, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))

    plt.title(titulo, fontsize=16)
    plt.xlabel('Tempo (Meses)', fontsize=12)
    plt.ylabel('Probabilidade de Sobrevida', fontsize=12)
    plt.legend(title="Grupo", loc='upper right')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    try:
        plt.savefig(nome_arquivo)
        print(f"  -> Imagem salva: {os.path.basename(nome_arquivo)}")
    except Exception as e:
        print(f"  -> ERRO ao salvar imagem: {e}")
        
    plt.close()

def plotar_pam50(df_dados, coluna_grupo, titulo, nome_arquivo):
    """
    Gera e salva um gráfico de barras empilhadas (PAM50) para os grupos definidos.
    """
    print(f"Gerando gráfico de PAM50: {titulo}")
    
    df_plot = df_dados.dropna(subset=[COL_PAM50, coluna_grupo])
    df_plot = df_plot[df_plot[coluna_grupo] != 'Ignorar']
    
    # --- Filtra o subtipo "Normal" ---
    if 'Normal' in df_plot[COL_PAM50].values:
        df_plot = df_plot[df_plot[COL_PAM50] != 'Normal'].copy()
        print("  -> Subtipo 'Normal' removido da análise PAM50.")

    if df_plot.empty or df_plot[COL_PAM50].nunique() < 2:
        print(f"  -> Aviso: Dados insuficientes para '{titulo}' após filtragem.")
        return

    # --- Preparação dos Dados (Crosstab) ---
    df_cross = pd.crosstab(df_plot[coluna_grupo], df_plot[COL_PAM50])
    df_perc = df_cross.div(df_cross.sum(axis=1), axis=0) * 100
    
    # Ordena as colunas do PAM50 para consistência (LumA, LumB, Her2, Basal)
    ordered_keys = ['LumA', 'LumB', 'Her2', 'Basal', 'Normal']
    col_order = [col for col in ordered_keys if col in df_perc.columns]
    df_perc = df_perc[col_order]

    # --- Análise Estatística (Qui-Quadrado / Chi-Square) ---
    p_valor = 0.99
    try:
        # Tenta calcular o teste de Qui-quadrado
        chi2, p_valor, dof, expected = chi2_contingency(df_cross)
    except ValueError as e:
        # Isso pode ocorrer se houver células com frequência zero ou muito pequenas
        print(f"  -> Erro no teste estatístico (chi2_contingency): {e}. P-valor definido como 0.99.")

    # --- Plotagem (Gráfico de Barras Empilhadas) ---
    
    # Pega as cores na ordem correta das colunas
    color_list = [PAM50_COLORS.get(col, '#999999') for col in df_perc.columns]

    ax = df_perc.plot(
        kind='bar', 
        stacked=True, 
        figsize=(12, 8),
        color=color_list # Usa a lista de cores customizada
    )

    # Exibe o P-Valor (abaixo da legenda)
    plt.text(1.02, 0.5, f"Teste Qui-quadrado:\n{formatar_pval(p_valor)}", 
             transform=ax.transAxes, fontsize=12, fontweight='bold',
             bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
    
    plt.title(titulo, fontsize=16)
    plt.xlabel('Grupo de Comparação', fontsize=12)
    plt.ylabel('Distribuição Percentual (%)', fontsize=12)
    plt.xticks(rotation=0)
    
    plt.legend(title='Subtipo PAM50', bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    # Ajusta o layout para caber a legenda e o texto do p-valor
    plt.tight_layout(rect=[0, 0, 0.82, 1]) 
    
    try:
        plt.savefig(nome_arquivo)
        print(f"  -> Imagem salva: {os.path.basename(nome_arquivo)}")
    except Exception as e:
        print(f"  -> ERRO ao salvar imagem: {e}")
        
    plt.close()


def main():
    """
    Função principal: Carrega dados, define grupos e chama as funções de plotagem.
    """
    print("--- Iniciando Script de Geração de Imagens (IL5RA, IL33, IL25, TSLP) ---")
    
    # Tenta criar o objeto Tkinter antes de usá-lo
    try:
        root = tk.Tk()
        root.withdraw() # Esconde a janela principal do Tkinter
    except Exception as e:
        print(f"Aviso: Não foi possível inicializar o Tkinter ({e}). A caixa de diálogo de seleção de arquivo pode não funcionar em todos os ambientes.")
        return

    print("\nSelecione o arquivo CSV processado (do Script 1)")
    
    try:
        caminho_csv = filedialog.askopenfilename(
            title="Selecione o arquivo CSV processado (com dados de expressão e clínicos)",
            filetypes=[("CSV files", "*.csv")]
        )
        if not caminho_csv:
            print("Seleção de arquivo cancelada. Encerrando.")
            return
            
        print(f"Carregando dados de: {caminho_csv}")
        df_merged = pd.read_csv(caminho_csv, index_col=0, low_memory=False)
        caminho_pasta = os.path.dirname(caminho_csv)
        print("Dados carregados com sucesso.")
        
    except Exception as e:
        print(f"ERRO ao carregar o arquivo CSV: {e}")
        return

    colunas_essenciais = GENES_INTERESSE + [COL_TEMPO, COL_EVENTO, COL_PAM50]
    # Verifica se as colunas de genes e clínicas existem
    colunas_faltando = [col for col in colunas_essenciais if col not in df_merged.columns]
    
    if colunas_faltando:
        print("ERRO: O arquivo CSV não contém todas as colunas necessárias.")
        print(f"Colunas ausentes: {', '.join(colunas_faltando)}")
        return

    # Define os grupos de comparação no DataFrame
    df_merged = definir_grupos_comparacao(df_merged)
    if df_merged is None:
        # Se definir_grupos_comparacao retornar None, ocorreu um erro de gene ausente
        print("ERRO: Falha ao definir os grupos. Encerrando.")
        return

    # Lista de colunas de grupo a serem plotadas
    grupos_para_plotar = [
        ('grupo_4_genes_vs_outros', 'IL5RA_IL33_IL25_TSLP vs. Demais', '4_genes'),
        ('grupo_IL5RA_vs_Nao', 'Expressa IL5RA vs. Não Expressa IL5RA', 'IL5RA'),
        ('grupo_IL33_vs_Nao', 'Expressa IL33 vs. Não Expressa IL33', 'IL33'),
        ('grupo_IL25_vs_Nao', 'Expressa IL25 vs. Não Expressa IL25', 'IL25'),
        ('grupo_TSLP_vs_Nao', 'Expressa TSLP vs. Não Expressa TSLP', 'TSLP'),
    ]

    # --- 4. Gerar Gráficos de Sobrevida ---
    print("\n--- Gerando Gráficos de Sobrevida ---")
    for col_grupo, titulo_base, nome_curto in grupos_para_plotar:
        plotar_sobrevida(
            df_merged, col_grupo,
            f'Sobrevida: {titulo_base}',
            os.path.join(caminho_pasta, f'sobrevida_{nome_curto}.png')
        )
    
    # --- 5. Gerar Gráficos de PAM50 ---
    print("\n--- Gerando Gráficos de PAM50 ---")
    for col_grupo, titulo_base, nome_curto in grupos_para_plotar:
        plotar_pam50(
            df_merged, col_grupo,
            f'Distribuição PAM50: {titulo_base}',
            os.path.join(caminho_pasta, f'pam50_{nome_curto}.png')
        )
    
    print("\n--- Processo de Geração de Imagens Concluído ---")

if __name__ == "__main__":
    # Configurações de estilo para os gráficos
    sns.set_theme(style='whitegrid', palette='deep')
    try:
        plt.rcParams['font.family'] = 'Times New Roman'
    except:
        print("Fonte Arial não encontrada, usando fonte padrão.")
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['figure.dpi'] = 100
    
    main()