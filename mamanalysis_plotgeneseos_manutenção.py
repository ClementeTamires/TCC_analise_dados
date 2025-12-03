import pandas as pd
import matplotlib.pyplot as plt
import io
import tkinter as tk
from tkinter import filedialog
import sys

def load_and_analyze_data():
    """
    Carrega o arquivo EXCEL de dados de pacientes (PAM50 e expressão gênica),
    processa os dados para calcular a frequência e a contagem absoluta de expressão 
    de genes individuais e plota o gráfico de barras empilhadas por subtipo PAM50.
    """
    
    # ----------------------------------------------------------------------
    # GENES DE INTERESSE - ADAPTADO PARA EOSINÓFILOS: IL5, IL33, IL25, TSLP
    GENES_OF_INTEREST = ['IL5', 'IL33', 'IL25', 'TSLP']
    # Nome da aba do Excel que contém os dados de expressão e PAM50
    SHEET_NAME = 'Dados_PAM50' 
    # Ordem dos subtipos PAM50 para o gráfico
    SUBTYPE_ORDER = ['LumA', 'LumB', 'Her2', 'Basal']
    
    # CORES PERSONALIZADAS PARA OS GENES (4 CORES) - NOVA PALETA DE ALTO CONTRASTE
    # 1. IL5: Vermelho ('#e41a1c')
    # 2. IL33: Azul ('#377eb8')
    # 3. IL25: Verde ('#4daf4a') 
    # 4. TSLP: Roxo/Violeta ('#984ea3')
    GENE_COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
    # Cores que requerem texto branco para contraste:
    DARK_COLORS = ['#377eb8', '#984ea3']
    # ----------------------------------------------------------------------
    
    # 1. Seleção do Arquivo com Tkinter filedialog
    
    # Inicializa o Tkinter e esconde a janela principal
    root = tk.Tk()
    root.withdraw() 
    
    # MENSAGEM ATUALIZADA PARA EXCEL
    print("Aguardando a seleção do arquivo Excel (ex: 'Analise_Genes_20251123_2036.xlsx')...")
    
    # Abre o diálogo para seleção de arquivo EXCEL
    file_path = filedialog.askopenfilename(
        title="Selecione o arquivo Excel de Dados de Pacientes",
        # TIPOS DE ARQUIVO ATUALIZADOS PARA EXCEL
        filetypes=[("Arquivos Excel", "*.xlsx"), ("Todos os arquivos", "*.*")]
    )
    
    # Verifica se o usuário selecionou um arquivo
    if not file_path:
        print("Nenhum arquivo selecionado. Encerrando a análise.")
        # Se um root foi criado, ele deve ser destruído.
        if 'root' in locals():
            root.destroy()
        return

    # Certifica-se de destruir a instância do Tk após o uso do filedialog
    if 'root' in locals():
        root.destroy()

    # 2. Carregamento e Pré-processamento dos Dados
    try:
        # Carrega o arquivo Excel, especificando a aba. O cabeçalho é na linha 1 (header=0).
        df = pd.read_excel(file_path, sheet_name=SHEET_NAME, header=0)
        
        # Renomeia a coluna PAM50.
        if len(df.columns) > 1:
            # Assumimos que a coluna 1 (índice 1) é o subtipo PAM50
            df = df.rename(columns={df.columns[1]: 'PAM50Subtype'})
        else:
            raise ValueError("O arquivo Excel/aba tem dados insuficientes. Verifique a estrutura.")
        
        # Garante que as colunas essenciais estão presentes
        required_cols = ['PAM50Subtype'] + GENES_OF_INTEREST
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            # Mensagem de erro atualizada
            raise ValueError(f"A aba '{SHEET_NAME}' está faltando as colunas essenciais: {missing_cols}. Colunas presentes: {list(df.columns)}")
            
        # Filtra e converte colunas de expressão para numérico
        df = df[required_cols].copy()
        df[GENES_OF_INTEREST] = df[GENES_OF_INTEREST].apply(pd.to_numeric, errors='coerce')
        
        # Remove pacientes sem subtipo PAM50 (NaN na coluna Subtype)
        df = df.dropna(subset=['PAM50Subtype'])
        
        # FILTRAGEM: Exclui o subtipo 'Normal' da análise
        df = df[df['PAM50Subtype'] != 'Normal'].copy()
        print("Subtipo 'Normal' excluído da análise.")
        
        # REORDENAMENTO: Define a ordem desejada dos subtipos
        df['PAM50Subtype'] = pd.Categorical(
            df['PAM50Subtype'], 
            categories=SUBTYPE_ORDER, 
            ordered=True
        )
        # Remove quaisquer linhas que possam ter se tornado NaN devido à categorização (embora a linha anterior de filtragem já devesse ter lidado com isso)
        df = df.dropna(subset=['PAM50Subtype'])
        print(f"Subtipos reordenados no gráfico: {SUBTYPE_ORDER}")
        
        # 3. Binarização da Expressão (Gene Expresso se valor > 0)
        # Cria um novo DataFrame binário (0 ou 1)
        df_expressed = df[GENES_OF_INTEREST].apply(lambda x: (x > 0).astype(int))
        df_expressed['PAM50Subtype'] = df['PAM50Subtype']

    except Exception as e:
        # Mensagens de erro atualizadas
        print(f"Erro ao carregar ou processar o arquivo Excel: {e}")
        print(f"Verifique se o arquivo é um .xlsx válido e se a aba '{SHEET_NAME}' existe e tem a estrutura esperada (PAM50 na segunda coluna, genes nas colunas seguintes).")
        return

    # Se o DataFrame estiver vazio após o carregamento, saímos.
    if df_expressed.empty:
        print("Não foi possível processar os dados de expressão. Encerrando.")
        return
        
    # 4. Cálculo de Contagens e Frequências
    
    # A. Contagem absoluta de pacientes em cada subtipo (para normalização)
    # O groupby agora usará a ordem categórica definida
    total_patients_per_subtype = df_expressed.groupby('PAM50Subtype').size()
    
    # B. Contagem absoluta de pacientes que expressam cada gene, por subtipo
    df_absolute_counts = df_expressed.groupby('PAM50Subtype')[GENES_OF_INTEREST].sum()

    # C. Cálculo da Frequência Média (Porcentagem) - Base do novo eixo Y
    df_percentage = df_absolute_counts.div(total_patients_per_subtype, axis=0) * 100
    
    # 5. Plotagem do Gráfico de Barras Empilhadas
    fig, ax = plt.subplots(figsize=(12, 8))

    # Cores e ponto de início para o empilhamento
    
    # Reverte para acumular a Porcentagem para o empilhamento (bottoms)
    bottoms = pd.Series([0.0] * len(df_percentage.index), index=df_percentage.index)

    # Plotar cada gene como uma "camada"
    for i, gene in enumerate(GENES_OF_INTEREST):
        # Valores de altura: Porcentagem (Eixo Y)
        heights = df_percentage[gene] 
        # Contagem Absoluta (para rótulo)
        absolute_count = df_absolute_counts[gene]
        
        # Cria a barra empilhada (usando a Porcentagem como altura)
        bars = ax.bar(
            df_percentage.index, # Eixo X: Subtipos PAM50
            heights, # Altura: Frequência Média (%)
            bottom=bottoms, # Ponto de início (Porcentagem)
            label=gene, # Legenda (o gene individual)
            color=GENE_COLORS[i], # Usa a cor definida na lista personalizada
            edgecolor='black'
        )
        
        # Adicionar os rótulos de contagem (absoluta) e porcentagem na barra
        for bar, absolute_count_val, percentage_val in zip(bars, absolute_count, heights):
            # Encontra o subtipo correspondente. 
            x_labels = df_percentage.index
            
            # Calcula o índice da barra
            bar_index = list(ax.get_xticks()).index(bar.get_x() + bar.get_width()/2)
            subtype = x_labels[bar_index]

            # A posição Y para o texto é o ponto médio da barra
            yval = bottoms.loc[subtype] + bar.get_height() / 2
            
            # Formatar o texto: Contagem Absoluta (Porcentagem%)
            label_text = f'{int(absolute_count_val)}\n({percentage_val:.1f}%)'
            
            # Apenas adicionar o rótulo se a barra for alta o suficiente (ex: > 3% para clareza)
            if bar.get_height() > 3: 
                # Determina a cor do texto com base na cor da barra
                current_color = GENE_COLORS[i]
                # Usa texto branco nas cores mais escuras (Roxo e Azul)
                text_color = 'white' if current_color in DARK_COLORS else 'black'
                
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    yval,
                    label_text,
                    ha='center',
                    va='center',
                    color=text_color, 
                    # AUMENTA A FONTE DOS RÓTULOS INTERNOS
                    fontsize=11, 
                    fontweight='bold'
                )
                
        # Atualiza o ponto de início (bottom) para a próxima barra empilhada
        bottoms += heights

    # 5. Finalização do Gráfico
    genes_str = ', '.join(GENES_OF_INTEREST)
    
    # O limite Y total é baseado na soma máxima das frequências
    y_max = bottoms.max() * 1.05
    
    ax.set_title(f'Frequência de Expressão Individual dos Genes de Eosinófilos ({genes_str}) por Subtipo PAM50 (Excluído: Normal)', fontsize=14, pad=20)
    # AUMENTA A FONTE DO TÍTULO DO EIXO X
    ax.set_xlabel('Subtipo PAM50', fontsize=14) 
    
    # AUMENTA A FONTE DOS NOMES DOS SUBTIPOS (TICK LABELS)
    ax.tick_params(axis='x', labelsize=12)
    
    # Reverte para o rótulo de Frequência
    ax.set_ylabel('Frequência de Expressão (%)', fontsize=12) 
    
    # Define o intervalo dos ticks com base no valor máximo (incrementos de 10%)
    ax.set_yticks(range(0, int(y_max) + 10, 10))
    ax.set_ylim(0, y_max)

    # Configuração da legenda
    ax.legend(title='Gene Individual Expresso', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()

    print("\nAnálise concluída. O gráfico de barras empilhadas para os genes de eosinófilos foi gerado e exibido.")
    print("O subtipo 'Normal' foi excluído da análise e do gráfico.")
    print(f"A ordem dos subtipos no eixo X é: {SUBTYPE_ORDER}")
    print("Eixo Y: Frequência de Expressão (Porcentagem de pacientes do subtipo que expressam o gene).")
    print("Rótulos nas barras: Contagem Absoluta (Porcentagem de pacientes do subtipo que expressam o gene).")
    print("O total da barra é a soma das frequências individuais, indicando a média de genes expressos por paciente no subtipo.")

# Executa a função principal
if __name__ == '__main__':
    # Adiciona um try/except para capturar exceções do tkinter se a interface não estiver disponível
    try:
        load_and_analyze_data()
    except Exception as e:
        print(f"\nOcorreu um erro na execução da função principal: {e}")