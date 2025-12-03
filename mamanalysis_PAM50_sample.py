import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import sys

# Dados de fallback (dados originais extraídos)
# O subtipo 'Normal' (118) foi removido das contagens para plotagem e análise.
DADOS_FALLBACK = {
    'LumA': 433,
    'LumB': 194,
    'Her2': 67,
    'Basal': 141,
    'Total_Classificado_PAM50': 953, # 433 + 194 + 67 + 141 + 118 (Normal)
    'Total_Dataset': 1215 # Total de pacientes listado no documento
}

def carregar_dados_pam50():
    """Abre uma caixa de diálogo para o usuário selecionar o arquivo CSV e extrai os dados PAM50 (excluindo 'Normal')."""
    
    # Esconde a janela principal do Tkinter
    root = tk.Tk()
    root.withdraw() 
    
    # Abre a caixa de diálogo para seleção de arquivo CSV
    caminho_arquivo = filedialog.askopenfilename(
        title="Selecione o arquivo CSV da Análise Estatística (Analise_Estatistica.csv)",
        filetypes=[("Arquivos CSV", "*.csv"), ("Todos os arquivos", "*.*")]
    )
    
    if not caminho_arquivo:
        # Retorna os dados de fallback se nenhum arquivo for selecionado
        print("Nenhum arquivo selecionado. Usando dados internos (hardcoded) para demonstração.")
        return DADOS_FALLBACK['LumA'], DADOS_FALLBACK['LumB'], DADOS_FALLBACK['Her2'], DADOS_FALLBACK['Basal'], DADOS_FALLBACK['Total_Classificado_PAM50'], DADOS_FALLBACK['Total_Dataset']

    try:
        # Carrega o arquivo inteiro.
        df_completo = pd.read_csv(caminho_arquivo, header=None)
        
        # 1. Tenta extrair o Total de Pacientes (Valor '1215' esperado na linha de índice 1, coluna 0)
        total_dataset = int(df_completo.iloc[1, 0])
        
        # 2. Localiza a tabela PAM50
        start_row_pam50 = df_completo[df_completo.iloc[:, 0] == 'PAM50 Subtipo'].index
        
        if start_row_pam50.empty:
            raise ValueError("Não foi possível encontrar o cabeçalho 'PAM50 Subtipo' no arquivo.")
            
        # O início dos dados é a linha seguinte ao cabeçalho (data frame tem 5 linhas: LumA, LumB, Basal, Normal, Her2)
        start_row_data = start_row_pam50[0] + 1
        
        # Seleciona as colunas de Subtipo e Contagem Absoluta (colunas 0 e 1)
        df_pam50 = df_completo.iloc[start_row_data : start_row_data + 5, :2].copy()
        df_pam50.columns = ['Subtipo', 'Contagem']
        df_pam50 = df_pam50.set_index('Subtipo')['Contagem'].to_dict()

        # Extrai os valores dos 4 subtipos clínicos
        lum_a = df_pam50.get('LumA', 0)
        lum_b = df_pam50.get('LumB', 0)
        her2 = df_pam50.get('Her2', 0)
        basal = df_pam50.get('Basal', 0)
        normal = df_pam50.get('Normal', 0) # Ainda extrai 'Normal' para calcular o total classificado
        
        total_classificado_pam50 = lum_a + lum_b + her2 + basal + normal

        print(f"Dados carregados com sucesso do arquivo: {caminho_arquivo}")
        return lum_a, lum_b, her2, basal, total_classificado_pam50, total_dataset

    except Exception as e:
        print(f"Erro ao carregar ou processar o arquivo: {e}. Usando dados internos (hardcoded) para demonstração.")
        # Retorna os dados de fallback em caso de qualquer erro de processamento
        return DADOS_FALLBACK['LumA'], DADOS_FALLBACK['LumB'], DADOS_FALLBACK['Her2'], DADOS_FALLBACK['Basal'], DADOS_FALLBACK['Total_Classificado_PAM50'], DADOS_FALLBACK['Total_Dataset']

# Chama a função de carregamento para obter os dados
lum_a_count, lum_b_count, her2_count, basal_count, TOTAL_CLASSIFICADO_PAM50, TOTAL_DATASET = carregar_dados_pam50()

# -----------------------------------------------------------
# Configuração dos Dados para Plotagem (Apenas 4 subtipos)
# -----------------------------------------------------------

# Organiza os dados para plotagem (apenas os 4 subtipos clínicos solicitados)
contagens = np.array([lum_a_count, lum_b_count, her2_count, basal_count])
subtipos = ['LumA', 'LumB', 'Her2', 'Basal']
total_amostras_plotadas = contagens.sum()

# Imprime o resumo da análise
print("==================================================================================")
print(f"Total de Amostras no Dataset (Total de Pacientes): {TOTAL_DATASET} amostras.")
print(f"Total de Amostras com Classificação PAM50 (5 subtipos): {TOTAL_CLASSIFICADO_PAM50} amostras.")
print(f"Total de Amostras ANALISADAS/PLOTADAS (4 Subtipos Clínicos): {total_amostras_plotadas} amostras.")
print("Amostras do subtipo 'Normal' foram excluídas da plotagem.")
print("==================================================================================")

# -----------------------------------------------------------
# 1. Preparação dos Dados para Plotagem
# -----------------------------------------------------------

# Calcula as porcentagens (Baseadas no TOTAL CLASSIFICADO, para refletir os valores da tabela original)
porcentagens = (contagens / TOTAL_CLASSIFICADO_PAM50) * 100

# Formata as porcentagens para serem exibidas como rótulos
rotulos_porcentagem = [f'{p:.1f}%' for p in porcentagens]

# -----------------------------------------------------------
# 2. Configuração e Criação do Gráfico
# -----------------------------------------------------------

# Configurações do Matplotlib para usar a fonte Times New Roman (se disponível)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']

# Cria a figura e os eixos
fig, ax = plt.subplots(figsize=(10, 6))

# Cria o gráfico de barras
barras = ax.bar(
    subtipos, 
    contagens, 
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'] 
)

# Adiciona o rótulo de porcentagem em cima de cada barra
for i, barra in enumerate(barras):
    altura = barra.get_height()
    # Adiciona o texto no topo da barra
    ax.text(
        barra.get_x() + barra.get_width() / 2, # Posição X (centro da barra)
        altura + 5,                             # Posição Y (um pouco acima da barra)
        rotulos_porcentagem[i],                 # O texto da porcentagem
        ha='center',                            # Alinhamento horizontal (centralizado)
        va='bottom',                            # Alinhamento vertical (abaixo da posição)
        fontsize=16,                            # AJUSTE: Aumentado para 16
        fontweight='bold'
    )

# -----------------------------------------------------------
# 3. Configuração dos Eixos e Títulos
# -----------------------------------------------------------

# Eixo Y (Contagem Absoluta)
ax.set_ylabel('Contagem Absoluta de Amostras', fontsize=16)

# AJUSTE: Define o tamanho da fonte dos rótulos do eixo X (os subtipos) para 18
ax.tick_params(axis='x', labelsize=18)

# Limite do Eixo Y (ajustado para que as porcentagens caibam)
ax.set_ylim(0, max(contagens) * 1.15) 

# Adiciona uma linha horizontal para o zero (opcional, para clareza)
ax.axhline(0, color='grey', linewidth=0.8)

# Remove as molduras superior e direita (melhora a estética acadêmica)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Adiciona grid ao eixo Y (opcional, para melhor leitura dos valores)
# ax.grid(axis='y', linestyle='--', alpha=0.7) 

# Exibe o gráfico
plt.tight_layout() # Ajusta o layout para evitar sobreposição
plt.show()