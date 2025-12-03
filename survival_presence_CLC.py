import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog
import sys

# Função para configurar a interface gráfica (oculta) e solicitar o arquivo
def load_data_file():
    """Abre uma janela de diálogo para selecionar o arquivo de dados (CSV ou Excel)."""
    # Configura a janela principal do Tkinter, mas a mantém oculta
    root = Tk()
    root.withdraw() 
    
    file_path = filedialog.askopenfilename(
        title="Selecione o arquivo de Dados de Sobrevida (CSV ou Excel)",
        filetypes=[
            ("Arquivos Excel", "*.xlsx"),
            ("Arquivos CSV", "*.csv"), 
            ("Todos os arquivos", "*.*")
        ]
    )
    
    # Verifica se o usuário cancelou a seleção
    if not file_path:
        print("Nenhum arquivo selecionado. Encerrando o programa.")
        sys.exit(0)
    
    # Tenta carregar os dados
    try:
        if file_path.lower().endswith('.csv'):
            # Leitura de arquivo CSV
            data = pd.read_csv(file_path, index_col=0)
        elif file_path.lower().endswith(('.xlsx', '.xls')):
            # Leitura de arquivo Excel
            # Tenta ler a planilha chamada 'Dados_Sobrevida', se existir
            sheet_name_to_load = 'Dados_Sobrevida'
            try:
                data = pd.read_excel(file_path, sheet_name=sheet_name_to_load, index_col=0)
            except ValueError:
                print(f"Aviso: Planilha '{sheet_name_to_load}' não encontrada. Tentando carregar a primeira planilha (índice 0).")
                data = pd.read_excel(file_path, sheet_name=0, index_col=0)
        else:
            print("Formato de arquivo não suportado. Por favor, selecione .csv ou .xlsx.")
            sys.exit(1)
            
        # Limpar nomes de colunas (remover espaços em branco) para evitar KeyErrors comuns
        data.columns = data.columns.str.strip()
            
        print(f"Arquivo carregado com sucesso: {file_path}")
        return data
    except Exception as e:
        print(f"Erro ao carregar o arquivo: {e}")
        # Uma saída limpa sem usar messagebox
        sys.exit(1)

def create_survival_groups(df: pd.DataFrame, gene_list: list, group_title: str) -> pd.DataFrame:
    """
    Calcula a SOMA da expressão dos genes na lista e divide as amostras
    em grupos de Alta e Baixa Expressão com base na mediana.
    
    As regras de corte (cutoff) são:
    - Alta Expressão: Expressão Agregada >= Mediana
    - Baixa Expressão: Expressão Agregada < Mediana
    """
    # 1. Calcular a expressão AGREGADA (SOMA) dos genes no grupo
    # A coluna de tempo e evento não devem ser incluídas no cálculo
    expression_cols = [gene for gene in gene_list if gene in df.columns]
    
    if not expression_cols:
        print(f"Erro: Nenhum dos genes {gene_list} foi encontrado no DataFrame.")
        return None

    # Usando a SOMA (sum) para agregar a expressão, conforme solicitado pelo usuário.
    df['sum_expression'] = df[expression_cols].sum(axis=1)
    
    # 2. Definir o ponto de corte (mediana)
    # O ponto de corte é a mediana da soma das expressões
    median_expression = df['sum_expression'].median()
    
    # 3. Criar a coluna de grupo de sobrevida
    # Corte: >= Mediana para Alta Expressão, < Mediana para Baixa Expressão.
    df[group_title] = df['sum_expression'].apply(
        lambda x: 'Alta Expressão' if x >= median_expression else 'Baixa Expressão'
    )
    
    return df

def plot_survival(df: pd.DataFrame, group_column: str, gene_list_str: str, time_col: str, event_col: str, fig_title: str):
    """
    Plota as curvas de Kaplan-Meier e realiza o teste log-rank.
    Cria uma nova figura para ser exibida separadamente.
    """
    # Cria uma nova figura e eixo para cada plotagem
    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    
    # Colunas de Sobrevida
    T = df[time_col]
    E = df[event_col]
    
    # Separar os grupos
    high_exp = (df[group_column] == 'Alta Expressão')
    low_exp = (df[group_column] == 'Baixa Expressão')

    # Verificar se há dados suficientes nos grupos
    if high_exp.sum() < 2 or low_exp.sum() < 2:
        ax.set_title(f"Sobrevida Agregada: {gene_list_str} - Dados Insuficientes", fontsize=10)
        ax.text(0.5, 0.5, "Dados insuficientes para um ou ambos os grupos.", 
                transform=ax.transAxes, ha='center', va='center')
        ax.grid(False)
        ax.legend([], [], frameon=False) # Remove a legenda se não houver dados
        return # Não retorna figura, permite plt.show() no main
        

    # 1. Ajustar o modelo Kaplan-Meier (KM)
    kmf = KaplanMeierFitter()
    
    # Grupo de Alta Expressão
    kmf.fit(T[high_exp], E[high_exp], label=f'Alta Expressão (n={high_exp.sum()})')
    kmf.plot_survival_function(ax=ax, ci_show=True, color='red', linewidth=2)
    
    # Grupo de Baixa Expressão
    kmf.fit(T[low_exp], E[low_exp], label=f'Baixa Expressão (n={low_exp.sum()})')
    kmf.plot_survival_function(ax=ax, ci_show=True, color='blue', linewidth=2, linestyle='--')
    
    # 2. Teste Log-Rank
    results = logrank_test(T[high_exp], T[low_exp], E[high_exp], E[low_exp], alpha=.99)
    p_value = results.p_value
    
    # 3. Configuração do Gráfico
    ax.set_title(fig_title, fontsize=12, fontweight='bold')
    ax.set_xlabel("Tempo (Dias)")
    ax.set_ylabel("Probabilidade de Sobrevida Global")
    
    # Adicionar o P-valor (Fontsize ajustado de 10 para 12)
    ax.text(0.05, 0.15, f"Log-rank P-value: {p_value:.4f}", 
            transform=ax.transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    ax.grid(True, linestyle='--', alpha=0.6)
    # Garante a legenda de cores (Fontsize ajustado de 10 para 12)
    ax.legend(loc='lower left', frameon=True, fontsize=12) 
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    plt.tight_layout()

def main():
    # 1. Carregar os dados
    print("Iniciando a análise...")
    df = load_data_file()
    
    # 2. DEFINIR AS COLUNAS DE TEMPO E EVENTO (Ajuste aqui se o nome for diferente no seu arquivo)
    # Estes nomes de coluna foram confirmados pelo usuário.
    TIME_COL = 'OS_Time_nature2012'
    EVENT_COL = 'OS_event_nature2012'
    
    # 3. VERIFICAÇÃO DE ERRO: Garante que as colunas de sobrevida existem
    if TIME_COL not in df.columns or EVENT_COL not in df.columns:
        print("\nERRO CRÍTICO: As colunas de Sobrevida não foram encontradas.")
        print(f"O script está procurando por: Tempo='{TIME_COL}' e Evento='{EVENT_COL}'")
        print(f"Colunas encontradas no seu arquivo: {list(df.columns)}")
        print("\n** Ação necessária: ** Por favor, atualize as variáveis TIME_COL e EVENT_COL no código acima (na função main) com os nomes exatos das colunas do seu arquivo.")
        sys.exit(1)
        
    # Normalizando nomes de colunas (caso o índice esteja sem nome)
    if df.index.name is None:
        df.index.name = 'Amostra'
    df = df.reset_index()

    # Garantir que o nome da primeira coluna (índice original) não seja usado por engano
    if df.columns[0] == 'Unnamed: 0':
        df = df.rename(columns={'Unnamed: 0': 'Amostra'})

    # 4. Definir os grupos de genes
    GENES_G1 = ['CLC', 'EPX', 'IL5RA', 'PRG2']
    GENES_G2 = ['EPX', 'IL5RA', 'PRG2']
    GROUP_COL_1 = 'Grupo_G1'
    GROUP_COL_2 = 'Grupo_G2'
    
    plt.style.use('ggplot')
    
    # --- GRÁFICO 1: Análise do Grupo 1: [CLC, EPX, IL5RA, PRG2] ---
    df_g1 = create_survival_groups(df.copy(), GENES_G1, GROUP_COL_1)
    if df_g1 is not None:
        title_g1 = f"Curvas de Sobrevida de Kaplan-Meier\nExpressão Agregada (Soma) do Grupo: {', '.join(GENES_G1)}"
        plot_survival(df_g1, GROUP_COL_1, ", ".join(GENES_G1), TIME_COL, EVENT_COL, title_g1)
        # Mostrar o primeiro gráfico imediatamente
        plt.show() 
    
    # --- GRÁFICO 2: Análise do Grupo 2: [EPX, IL5RA, PRG2] ---
    df_g2 = create_survival_groups(df.copy(), GENES_G2, GROUP_COL_2)
    if df_g2 is not None:
        title_g2 = f"Curvas de Sobrevida de Kaplan-Meier\nExpressão Agregada (Soma) do Grupo: {', '.join(GENES_G2)}"
        plot_survival(df_g2, GROUP_COL_2, ", ".join(GENES_G2), TIME_COL, EVENT_COL, title_g2)
        # Mostrar o segundo gráfico imediatamente
        plt.show()

if __name__ == "__main__":
    # Garantir que o Tkinter não inicie a janela principal, apenas o filedialog
    try:
        main()
    except Exception as e:
        print(f"Ocorreu um erro no programa principal: {e}")