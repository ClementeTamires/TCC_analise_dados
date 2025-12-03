# ---- Plotagem de gráficos de sobrevida e distribuição PAM50 para genes CCL ----

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
# Genes de interesse atualizados para o grupo CCL
GENES_INTERESSE = ['CCL11', 'CCL24', 'CCL26']
COL_TEMPO = 'OS_Time_nature2012'
COL_EVENTO = 'OS_event_nature2012'
COL_PAM50 = 'PAM50Call_RNAseq'

# --- Paleta de Cores Customizada para PAM50 ---
# Tons de Azul, Roxo, Rosa, Vermelho, Verde-azulado
PAM50_COLORS = {
    "LumA": "#06B600",    # Azul (tipo "Azul Marinho")
    "LumB": "#B9E90E",    # Roxo
    "Her2": "#F59505",    # Vermelho
    "Basal": "#D60000",   # Rosa/Pink
    # "Normal" é filtrado, mas mantemos a cor caso o filtro seja removido
    "Normal": "#008080" 
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
    Cria as colunas de grupo específicas para as comparações dos genes CCL.
    """
    print("Definindo grupos de comparação...")
    
    # 1. Criar colunas booleanas para expressão (> 0)
    for gene in GENES_INTERESSE:
        if gene in df.columns:
            df[f'{gene}_expresso'] = df[gene] > 0
        else:
            print(f"Atenção: Gene {gene} não encontrado para criar grupos.")
            return None

    # 2. Definir o grupo complexo (todos os 3 genes expressos)
    
    # Grupo "3 Genes Expressos"
    g_3_all = (
        df['CCL11_expresso'] & 
        df['CCL24_expresso'] & 
        df['CCL26_expresso']
    )
    
    # 3. Criar as colunas de comparação para os gráficos
    
    # Comparações 1 e 5: (3 Genes vs. Todos os Outros)
    df['grupo_3_vs_outros'] = np.where(
        g_3_all, 
        'CCL11_CCL24_CCL26_Expresso', # Nome do grupo com os 3 expressos
        'Demais' 
    )
    
    # Comparações 2 e 6: (CCL11 vs. Nao)
    df['grupo_CCL11_vs_Nao'] = np.where(
        df['CCL11_expresso'], 
        'Expressa CCL11', 
        'Não Expressa CCL11'
    )
    
    # Comparações 3 e 7: (CCL24 vs. Nao)
    df['grupo_CCL24_vs_Nao'] = np.where(
        df['CCL24_expresso'], 
        'Expressa CCL24', 
        'Não Expressa CCL24'
    )
    
    # Comparações 4 e 8: (CCL26 vs. Nao)
    df['grupo_CCL26_vs_Nao'] = np.where(
        df['CCL26_expresso'], 
        'Expressa CCL26', 
        'Não Expressa CCL26'
    )
    
    return df

def plotar_sobrevida(df_dados, coluna_grupo, titulo, nome_arquivo):
    """
    Gera e salva um gráfico de sobrevida (Kaplan-Meier) para os grupos definidos.
    """
    print(f"Gerando gráfico: {titulo}")
    
    df_plot = df_dados.dropna(subset=[COL_TEMPO, COL_EVENTO, coluna_grupo])
    df_plot = df_plot[df_plot[coluna_grupo] != 'Ignorar']
    
    if df_plot.empty or df_plot[coluna_grupo].nunique() < 2:
        print(f"  -> Aviso: Dados insuficientes para '{titulo}' após filtragem.")
        return

    T = df_plot[COL_TEMPO] / 30.44 # Converte dias para meses
    E = df_plot[COL_EVENTO]
    grupos = df_plot[coluna_grupo]
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
            kmf.plot_survival_function(ax=ax, ci_show=True)

    plt.text(0.05, 0.05, formatar_pval(p_valor), transform=ax.transAxes,
             fontsize=14, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))

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
    print(f"Gerando gráfico: {titulo}")
    
    df_plot = df_dados.dropna(subset=[COL_PAM50, coluna_grupo])
    df_plot = df_plot[df_plot[coluna_grupo] != 'Ignorar']
    
    # --- FILTRO: Filtra o subtipo "Normal" ---
    if 'Normal' in df_plot[COL_PAM50].values:
        df_plot = df_plot[df_plot[COL_PAM50] != 'Normal'].copy()
        print("  -> Subtipo 'Normal' removido da análise PAM50.")

    if df_plot.empty or df_plot[COL_PAM50].nunique() < 2:
        print(f"  -> Aviso: Dados insuficientes para '{titulo}' após filtragem.")
        return

    # --- Preparação dos Dados (Crosstab) ---
    df_cross = pd.crosstab(df_plot[coluna_grupo], df_plot[COL_PAM50])
    df_perc = df_cross.div(df_cross.sum(axis=1), axis=0) * 100
    
    # Ordena as colunas do PAM50 para consistência
    col_order = [col for col in PAM50_COLORS.keys() if col in df_perc.columns]
    df_perc = df_perc[col_order]

    # --- Análise Estatística (Qui-Quadrado / Chi-Square) ---
    p_valor = 0.99
    try:
        chi2, p_valor, dof, expected = chi2_contingency(df_cross)
    except ValueError as e:
        print(f"  -> Erro no teste estatístico (chi2_contingency): {e}")

    # --- Plotagem (Gráfico de Barras Empilhadas) ---
    
    # Pega as cores na ordem correta das colunas
    color_list = [PAM50_COLORS.get(col, '#999999') for col in df_perc.columns]

    ax = df_perc.plot(
        kind='bar', 
        stacked=True, 
        figsize=(12, 8),
        color=color_list # Usa a lista de cores customizada
    )

    # Posição do P-Valor (abaixo da legenda)
    plt.text(1.02, 0.5, f"Teste Chi-Square:\n{formatar_pval(p_valor)}", 
             transform=ax.transAxes, fontsize=14, fontweight='bold',
             bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
    
    plt.title(titulo, fontsize=16)
    plt.xlabel('Grupo de Comparação', fontsize=12)
    plt.ylabel('Distribuição Percentual (%)', fontsize=12)
    plt.xticks(rotation=0)
    
    plt.legend(title='PAM50', bbox_to_anchor=(1.02, 1), loc='upper left')
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
    print("--- Iniciando Script de Análise de Genes CCL ---")
    print("\nSelecione o arquivo CSV processado (contendo dados dos genes CCL)")
    
    try:
        root = tk.Tk()
        root.withdraw()
        caminho_csv = filedialog.askopenfilename(
            title="Selecione o arquivo CSV processado",
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

    df_merged = definir_grupos_comparacao(df_merged)
    if df_merged is None:
        print("ERRO: Falha ao definir os grupos (genes essenciais ausentes). Encerrando.")
        return

    # --- 4. Gerar Gráficos de Sobrevida ---
    print("\n--- Gerando Gráficos de Sobrevida (Kaplan-Meier) ---")
    
    # Imagem 1: 3 Genes CCL vs. Demais
    plotar_sobrevida(
        df_merged, 'grupo_3_vs_outros',
        'Sobrevida: CCL11, CCL24, CCL26 Expressos vs. Demais',
        os.path.join(caminho_pasta, 'sobrevida_1_3_CCL_vs_demais.png')
    )
    
    # Imagem 2: CCL11 vs Nao
    plotar_sobrevida(
        df_merged, 'grupo_CCL11_vs_Nao',
        'Sobrevida: Expressa CCL11 vs. Não Expressa CCL11',
        os.path.join(caminho_pasta, 'sobrevida_2_CCL11_vs_Nao.png')
    )
    
    # Imagem 3: CCL24 vs Nao
    plotar_sobrevida(
        df_merged, 'grupo_CCL24_vs_Nao',
        'Sobrevida: Expressa CCL24 vs. Não Expressa CCL24',
        os.path.join(caminho_pasta, 'sobrevida_3_CCL24_vs_Nao.png')
    )

    # Imagem 4: CCL26 vs Nao
    plotar_sobrevida(
        df_merged, 'grupo_CCL26_vs_Nao',
        'Sobrevida: Expressa CCL26 vs. Não Expressa CCL26',
        os.path.join(caminho_pasta, 'sobrevida_4_CCL26_vs_Nao.png')
    )
    
    # --- 5. Gerar Gráficos de PAM50 ---
    print("\n--- Gerando Gráficos de PAM50 ---")
    
    # Imagem 5: 3 Genes CCL vs. Demais
    plotar_pam50(
        df_merged, 'grupo_3_vs_outros',
        'Distribuição PAM50: CCL11, CCL24, CCL26 Expressos vs. Demais',
        os.path.join(caminho_pasta, 'pam50_1_3_CCL_vs_demais.png')
    )
    
    # Imagem 6: CCL11 vs Nao
    plotar_pam50(
        df_merged, 'grupo_CCL11_vs_Nao',
        'Distribuição PAM50: Expressa CCL11 vs. Não Expressa CCL11',
        os.path.join(caminho_pasta, 'pam50_2_CCL11_vs_Nao.png')
    )
    
    # Imagem 7: CCL24 vs Nao
    plotar_pam50(
        df_merged, 'grupo_CCL24_vs_Nao',
        'Distribuição PAM50: Expressa CCL24 vs. Não Expressa CCL24',
        os.path.join(caminho_pasta, 'pam50_3_CCL24_vs_Nao.png')
    )
    
    # Imagem 8: CCL26 vs Nao
    plotar_pam50(
        df_merged, 'grupo_CCL26_vs_Nao',
        'Distribuição PAM50: Expressa CCL26 vs. Não Expressa CCL26',
        os.path.join(caminho_pasta, 'pam50_4_CCL26_vs_Nao.png')
    )
    
    print("\n--- Processo de Geração de Imagens Concluído ---")

if __name__ == "__main__":
    sns.set_theme(style='whitegrid', palette='deep')
    # Tenta usar uma fonte mais profissional se disponível
    try:
        plt.rcParams['font.family'] = 'Arial'
    except:
        print("Fonte Arial não encontrada, usando fonte padrão.")
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['figure.dpi'] = 100
    
    main()