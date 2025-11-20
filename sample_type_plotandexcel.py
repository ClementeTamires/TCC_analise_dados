import pandas as pd
import os
import tkinter as tk
from tkinter import filedialog
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import math
import openpyxl

# --- Mapeamento dos Códigos de Tipo de Amostra (Sample-Type) do TCGA ---
# Estes códigos são os dois dígitos na 4ª posição do barcode (ex: TCGA-XX-YYYY-ZZ-A)
SAMPLE_TYPE_CODES = {
    '01': 'Primary Solid Tumor (Tumor Sólido Primário - Geralmente Ressecção Cirúrgica)',
    '02': 'Recurrent Solid Tumor (Tumor Recorrente)',
    '03': 'Primary Blood Derived Cancer (Câncer Sanguíneo Primário)',
    '04': 'Recurrent Blood Derived Cancer - Bone Marrow (Câncer Recorrente Medular)',
    '05': 'Additional - New Primary (Adicional - Novo Primário)',
    '06': 'Metastatic (Metástase - Frequentemente Biópsia)',
    '07': 'Additional Metastatic (Metástase Adicional)',
    '08': 'Human Tumor Original Cells (Células Originais de Tumor Humano)',
    '09': 'Primary Blood Derived Cancer - Bone Marrow (Câncer Primário Medular)',
    '10': 'Blood Derived Normal (Normal Derivado de Sangue)',
    '11': 'Solid Tissue Normal (Tecido Sólido Normal - Adjacente à Ressecção)',
    '12': 'Buccal Cell Normal (Célula Bucal Normal)',
    '13': 'EBV Immortalized Normal (Normal Imortalizado por EBV)',
    '14': 'Bone Marrow Normal (Medula Óssea Normal)',
    '20': 'Control Analyte (Analito de Controle)',
    '40': 'Recurrent Blood Derived Cancer - Peripheral Blood (Câncer Recorrente Sanguíneo)',
    '50': 'Cell Lines (Linhagens Celulares)',
    '60': 'Primary Xenograft Tissue (Tecido Xenográfico Primário)',
    '61': 'Cell Line Derived Xenograft Tissue (Tecido Xenográfico Derivado de Linhagem)',
}

def extract_sample_type_code(tcga_barcode):
    """Extrai o código de 2 dígitos do tipo de amostra do barcode TCGA."""
    try:
        parts = tcga_barcode.split('-')
        if len(parts) >= 4:
            return parts[3][:2]
        return 'UNKNOWN'
    except Exception:
        return 'UNKNOWN'

def tentar_ler_csv_ou_tsv(file_path):
    """Tenta ler o arquivo como TSV e, se falhar, tenta como CSV. Retorna o DataFrame e o separador usado."""
    print("Tentando ler o arquivo como TSV (separador '\\t')...")
    try:
        df = pd.read_csv(file_path, sep='\t')
        if not df.empty and len(df.columns) > 1:
            print("Sucesso na leitura como TSV.")
            return df, '\t'
        print("Falha na leitura TSV ou DataFrame vazio/mal-formado. Tentando CSV...")
    except Exception:
        print("Falha na leitura TSV. Tentando CSV...")
        
    print("Tentando ler o arquivo como CSV (separador ',')...")
    try:
        df = pd.read_csv(file_path, sep=',')
        if not df.empty and len(df.columns) > 1:
            print("Sucesso na leitura como CSV.")
            return df, ','
        print("Falha na leitura CSV ou DataFrame vazio/mal-formado.")
    except Exception:
        print("Falha na leitura CSV.")
        
    return None, None

def calcular_intervalo_confianca_wilson(k, n, nivel_confianca=0.95):
    """
    Calcula o Intervalo de Confiança (CI) usando o Wilson Score Interval para proporções.
    k: número de sucessos (count), n: número total de tentativas (total_samples).
    Retorna a margem de erro.
    """
    if n == 0:
        return 0.0 # Sem amostras
    
    # Nível de significância (alfa)
    alfa = 1 - nivel_confianca
    
    # Z-score (valor Z para a distribuição normal padrão)
    z = norm.ppf(1 - alfa / 2)
    
    # Proporção
    p = k / n
    
    # Cálculo do intervalo de Wilson
    termo1 = p + (z*z) / (2*n)
    termo2 = z * math.sqrt(p * (1 - p) / n + (z * z) / (4 * n * n))
    denominador = 1 + (z * z) / n
    
    lower_bound = (termo1 - termo2) / denominador
    
    # A margem de erro (ME) para plotagem é a diferença entre a proporção e o limite inferior.
    # Usamos o limite inferior pois ele é o termo mais conservador do erro.
    margin_of_error_proportion = p - lower_bound
    
    return margin_of_error_proportion

def gerar_grafico_barras_estatistico(summary_df, total_samples, output_path):
    """Gera um gráfico de barras em Python com Contagem (n), Porcentagem (%) e Barras de Erro (CI)."""
    
    # 1. Preparar os dados para plotagem
    codes = summary_df['Código do Tipo'].tolist()
    
    # ALTERAÇÃO SOLICITADA: Extrair APENAS o termo principal da tradução (antes do primeiro traço "-")
    # Ex: de "Tumor Sólido Primário - Geralmente Ressecção Cirúrgica" para "Tumor Sólido Primário"
    descriptions = summary_df['Descrição da Amostra'].apply(
        lambda x: x.split('(')[1].rstrip(')').split('-')[0].strip() if '(' in x else x
    ).tolist() 
    
    counts = summary_df['Contagem Absoluta'].tolist()
    percentages = summary_df['Porcentagem (%)'].tolist()
    
    # Rótulos combinados (ex: '01 - Tumor Sólido Primário')
    labels = [f"{c} - {d}" for c, d in zip(codes, descriptions)]

    # 2. Calcular as Barras de Erro (Margem de Erro do CI 95% do Wilson Score)
    error_bars_count = []
    for count in counts:
        # Erro é calculado na proporção
        me_proportion = calcular_intervalo_confianca_wilson(count, total_samples)
        # Convertemos para a escala da contagem (y-axis)
        error_count = me_proportion * total_samples
        error_bars_count.append(error_count)
    
    # 3. Criar o Gráfico
    # Ajuste o tamanho da figura para acomodar rótulos longos
    fig, ax = plt.subplots(figsize=(16, 9)) 
    
    # Plotar as barras com as barras de erro
    # Transpõe as barras de erro para serem simétricas em torno da média para visualização
    error_array = np.array(error_bars_count)
    bar_container = ax.bar(labels, counts, yerr=error_array, capsize=5, 
                           color='#059669', edgecolor='black', linewidth=0.7)
    
    # 4. Adicionar anotações (n e %) acima de cada barra
    y_lim_max = ax.get_ylim()[1] # Limite superior do eixo Y para referência
    
    for bar, count, percent, error in zip(bar_container, counts, percentages, error_array):
        height = bar.get_height()
        
        # Posição da anotação: no topo da barra + erro + uma pequena margem (2% do limite max do Y)
        text_pos_y = height + error + (y_lim_max * 0.02)

        # O texto é formatado como "n=X (Y.Y%)"
        label = f'n={count}\n({percent:.2f}%)'
        
        ax.text(bar.get_x() + bar.get_width() / 2., text_pos_y,
                label,
                ha='center', va='bottom', fontsize=10, 
                fontweight='bold', 
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', boxstyle="round,pad=0.3"))

    # 5. Configurações de Título e Eixos
    
    # Aumentar tamanho da fonte dos TÍTULOS dos eixos (mantido em 14)
    # Título do Eixo X em Português
    ax.set_xlabel('Tipo de Amostra (Código e Descrição em Português)', fontsize=14) 
    # Título do Eixo Y em Português
    ax.set_ylabel('Contagem Absoluta (n)', fontsize=14)
    
    # Ajustar o limite Y para garantir que as anotações caibam
    # Aumentamos o multiplicador para garantir mais espaço no topo
    ax.set_ylim(top=max(counts) * 1.35 if counts else 10) 
    
    # Aumentar tamanho da fonte dos RÓTULOS (tick labels) dos eixos (mantido em 11)
    plt.xticks(rotation=45, ha='right', fontsize=11)
    plt.yticks(fontsize=11)
    
    # Adicionar a informação de Confidence Interval na legenda
    legend_text = f"Barra de Erro: Margem de Erro do Intervalo de Confiança de 95% para a Proporção (Wilson Score)."
    
    # Adicionar uma caixa de texto no gráfico para a legenda 
    ax.text(0.01, 0.98, legend_text, 
            transform=ax.transAxes, 
            bbox=dict(boxstyle="round,pad=0.5", fc="lightgray", alpha=0.6, ec="gray"), 
            fontsize=9, va='top', ha='left')
            
    plt.grid(axis='y', linestyle='--', alpha=0.7) # Linhas de grade no eixo Y
    plt.tight_layout() # Ajusta o layout automaticamente
    
    # 6. Salvar o gráfico
    caminho_grafico_saida = output_path.replace('.xlsx', '_Grafico.png')
    plt.savefig(caminho_grafico_saida, dpi=300)
    print(f"\nGráfico de barras gerado com sucesso em: {caminho_grafico_saida}")
    plt.close(fig)
    
    return caminho_grafico_saida

def finalizar_processamento_e_gerar_arquivos(df_processed, input_file):
    """Gera o arquivo Excel de resumo e chama a função de plotagem."""
    
    # Gerar carimbo de data/hora
    timestamp = datetime.now().strftime('%y%m%d_%H%M%S')
    
    # 1. Definir o caminho de saída para o Excel
    diretorio_entrada = os.path.dirname(input_file)
    nome_excel_saida = f'TCGA_Resumo_Estatisticas_{timestamp}.xlsx'
    caminho_excel_saida = os.path.join(diretorio_entrada, nome_excel_saida)
    
    total_samples = len(df_processed)
    
    # 2. Criar o DataFrame de Resumo
    summary_df = df_processed.groupby(['Sample_Code', 'Sample_Type_Desc']).size().reset_index(name='Contagem Absoluta')
    summary_df['Porcentagem (%)'] = (summary_df['Contagem Absoluta'] / total_samples * 100).round(2)
    summary_df.columns = ['Código do Tipo', 'Descrição da Amostra', 'Contagem Absoluta', 'Porcentagem (%)']
    summary_df = summary_df.sort_values(by='Código do Tipo').reset_index(drop=True)

    print(f"\nEscrevendo a planilha de resumo estatístico no arquivo Excel: {caminho_excel_saida}...")
    
    # 3. Processo de Escrita no Excel
    try:
        with pd.ExcelWriter(caminho_excel_saida, engine='openpyxl') as writer:
            sheet_name_resumo = "Resumo Estatístico"
            summary_df.to_excel(writer, sheet_name=sheet_name_resumo, index=False)
            print(f" -> Planilha '{sheet_name_resumo}' criada com sucesso.")

    except Exception as e:
        print(f"ERRO: Ocorreu um erro ao escrever no Excel. Verifique se o arquivo está aberto ou se há permissões de escrita. Erro: {e}")
        return 0
    
    # 4. Geração do Gráfico
    if not summary_df.empty:
        gerar_grafico_barras_estatistico(summary_df, total_samples, caminho_excel_saida)

    print(f"\nProcessamento concluído! O arquivo Excel está em: {caminho_excel_saida}")
    return 1

def processar_dados_tcga(input_file):
    """Lê, extrai e sumariza os tipos de amostra."""
    
    # 1. Tentar ler o arquivo
    df, separador = tentar_ler_csv_ou_tsv(input_file)

    if df is None or df.empty:
        print("\nERRO CRÍTICO: Não foi possível carregar os dados. Verifique se o arquivo está no formato CSV/TSV e não está vazio.")
        return None

    # 2. Identificar a coluna de barcode (a primeira coluna)
    barcode_column = df.columns[0]
    df.rename(columns={barcode_column: 'TCGA_Barcode'}, inplace=True)
    barcode_column = 'TCGA_Barcode'
    
    print(f"Coluna de Barcode identificada como: '{barcode_column}' (Separador usado: '{separador}')")

    # 3. Extrair os códigos de tipo de amostra
    df['Sample_Code'] = df[barcode_column].apply(extract_sample_type_code)
    
    # Mapear o código para uma descrição legível
    df['Sample_Type_Desc'] = df['Sample_Code'].map(SAMPLE_TYPE_CODES).fillna('Outros/Desconhecido')

    print("\n--- Contagem de Amostras por Código TCGA ---")
    total_samples = len(df)
    summary = df.groupby(['Sample_Code', 'Sample_Type_Desc']).size().reset_index(name='Contagem')
    summary['Porcentagem'] = (summary['Contagem'] / total_samples * 100).round(2).astype(str) + '%'
    summary.columns = ['Código', 'Descrição da Amostra', 'Contagem', 'Porcentagem']
    summary = summary[['Código', 'Descrição da Amostra', 'Contagem', 'Porcentagem']]
    print(summary.to_string(index=False))
    print(f"\nTOTAL DE AMOSTRAS PROCESSADAS: {total_samples}")
    
    return df

def selecionar_arquivo_bruto():
    """Abre uma caixa de diálogo para o usuário selecionar o arquivo de dados."""
    # Inicializa o Tkinter e esconde a janela principal
    root = tk.Tk()
    root.withdraw()

    print("--- Seletor de Arquivo ---")
    
    nome_arquivo_entrada = filedialog.askopenfilename(
        title="Selecione o arquivo de dados brutos combinados (CSV ou TSV)",
        filetypes=(("Arquivos de Dados", "*.tsv *.csv"), ("Todos os Arquivos", "*.*"))
    )
    
    if not nome_arquivo_entrada:
        print("\nNenhum arquivo selecionado. Encerrando o programa.")
        return None
    else:
        print(f"\nArquivo selecionado com sucesso: {nome_arquivo_entrada}")
        return nome_arquivo_entrada

if __name__ == "__main__":
    caminho_arquivo = selecionar_arquivo_bruto()
    
    if caminho_arquivo:
        df_processado = processar_dados_tcga(caminho_arquivo)

        if df_processado is not None and not df_processado.empty:
            finalizar_processamento_e_gerar_arquivos(df_processado, caminho_arquivo)
        else:
            print("\nO processamento falhou ou o DataFrame está vazio. Não foi possível gerar o Excel ou o Gráfico.")