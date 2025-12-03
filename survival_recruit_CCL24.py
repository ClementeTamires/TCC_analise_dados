import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog
import sys
import os # Importar o módulo os para manipulação de caminhos/arquivos

# Configuração global de estilo para os gráficos
plt.style.use('ggplot')

# Função para configurar a interface gráfica (oculta) e solicitar o arquivo
def load_data_file():
    """Abre uma janela de diálogo para selecionar o arquivo de dados (CSV ou Excel)."""
    # Configura a janela principal do Tkinter, mas a mantém oculta
    root = Tk()
    # Impede que a janela principal apareça
    root.withdraw() 
    
    # Abre a caixa de diálogo para seleção do arquivo
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
        # Uma saída limpa
        sys.exit(1)

def create_survival_groups(df: pd.DataFrame, gene_list: list, group_title: str) -> pd.DataFrame:
    """
    Calcula a SOMA da expressão dos genes na lista e divide as amostras
    em grupos de Alta e Baixa Expressão com base na mediana.
    
    Esta função foi ajustada para:
    1. Usar a SOMA da expressão (em vez da média) conforme solicitado.
    2. Usar a Mediana como ponto de corte.
    3. 'Alta Expressão' = SOMA >= Mediana (Incluindo a mediana para consistência).
    4. 'Baixa Expressão' = SOMA < Mediana.
    """
    # 1. Calcular a expressão AGREGADA (SOMA) dos genes no grupo
    # Filtra apenas os genes que realmente existem no DataFrame
    expression_cols = [gene for gene in gene_list if gene in df.columns]
    
    if not expression_cols:
        print(f"Erro: Nenhum dos genes {gene_list} foi encontrado no DataFrame.")
        return None

    # NOVO/CONFIRMADO: Usando a SOMA (sum) da expressão dos genes
    df['sum_expression'] = df[expression_cols].sum(axis=1)
    
    # 2. Definir o ponto de corte (mediana)
    # O ponto de corte é a mediana da soma das expressões (Corte Mediano)
    median_expression = df['sum_expression'].median()
    
    # 3. Criar a coluna de grupo de sobrevida
    # Alta Expressão: Amostras com soma maior ou igual à mediana (inclui as medianas)
    # Baixa Expressão: Amostras com soma menor que a mediana
    df[group_title] = df['sum_expression'].apply(
        lambda x: 'Alta Expressão' if x >= median_expression else 'Baixa Expressão'
    )
    
    # Remove a coluna temporária de soma
    df = df.drop(columns=['sum_expression'])
    
    return df

def plot_survival(df: pd.DataFrame, group_column: str, gene_list_str: str, time_col: str, event_col: str, fig_title: str, save_path: str):
    """
    Plota as curvas de Kaplan-Meier, realiza o teste log-rank e salva o gráfico.
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
        ax.set_title(f"Sobrevida Agregada: {gene_list_str} - Dados Insuficientes", fontsize=12, fontweight='bold')
        ax.text(0.5, 0.5, "Dados insuficientes (n < 2 em um ou ambos os grupos).", 
                  transform=ax.transAxes, ha='center', va='center', color='red')
        ax.grid(False)
        ax.axis('off') # Oculta os eixos
        plt.tight_layout()
        
        # Tenta salvar mesmo o aviso
        try:
            plt.savefig(save_path)
            print(f"Aviso: O gráfico (com alerta de dados insuficientes) foi salvo em: {os.path.abspath(save_path)}")
        except Exception as e:
            print(f"Erro ao salvar o gráfico: {e}")
            
        plt.show() # Tenta exibir
        return
        
    # 1. Ajustar o modelo Kaplan-Meier (KM)
    kmf = KaplanMeierFitter()
    
    # Grupo de Alta Expressão (Vermelho sólido)
    kmf.fit(T[high_exp], E[high_exp], label=f'Alta Expressão (n={high_exp.sum()})')
    kmf.plot_survival_function(ax=ax, ci_show=True, color='red', linewidth=2)
    
    # Grupo de Baixa Expressão (Azul tracejado)
    kmf.fit(T[low_exp], E[low_exp], label=f'Baixa Expressão (n={low_exp.sum()})')
    kmf.plot_survival_function(ax=ax, ci_show=True, color='blue', linewidth=2, linestyle='--')
    
    # 2. Teste Log-Rank
    results = logrank_test(T[high_exp], T[low_exp], E[high_exp], E[low_exp], alpha=.99)
    p_value = results.p_value
    
    # 3. Configuração do Gráfico
    ax.set_title(fig_title, fontsize=14, fontweight='bold')
    ax.set_xlabel("Tempo (Dias)")
    ax.set_ylabel("Probabilidade de Sobrevida Global (S(t))")
    
    # Adicionar o P-valor (POSIÇÃO AJUSTADA para 0.20 para subir mais a legenda)
    ax.text(0.05, 0.20, f"Teste Log-rank P-value: {p_value:.4f}", 
              transform=ax.transAxes, fontsize=12, 
              bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round,pad=0.5'))
    
    ax.grid(True, linestyle=':', alpha=0.7)
    # Garante a legenda de cores
    ax.legend(loc='lower left', frameon=True, fontsize=12, title='Grupos de Expressão') 
    
    # Configuração estética (Limites do eixo Y, linhas de eixo)
    ax.set_ylim(0.0, 1.05)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    plt.tight_layout()

    # NOVO: SALVAR O GRÁFICO
    try:
        plt.savefig(save_path)
        print(f"Gráfico salvo com sucesso em: {os.path.abspath(save_path)}")
    except Exception as e:
        print(f"Erro ao salvar o gráfico: {e}")
            
    # Tenta exibir o gráfico (pode falhar em ambientes headless)
    plt.show() 

def main():
    # 1. Carregar os dados
    print("Iniciando a análise de sobrevida...")
    df = load_data_file()
    
    # 2. DEFINIR AS COLUNAS DE TEMPO E EVENTO (Ajuste aqui se o nome for diferente no seu arquivo)
    # Estas são as colunas de sobrevida esperadas
    TIME_COL = 'OS_Time_nature2012'
    EVENT_COL = 'OS_event_nature2012'
    
    # 3. VERIFICAÇÃO DE ERRO: Garante que as colunas de sobrevida existem
    if TIME_COL not in df.columns or EVENT_COL not in df.columns:
        print("\nERRO CRÍTICO: As colunas de Sobrevida não foram encontradas.")
        print(f"O script está procurando por: Tempo='{TIME_COL}' e Evento='{EVENT_COL}'")
        print(f"Colunas encontradas no seu arquivo: {list(df.columns)}")
        print("\n** Ação necessária: ** Por favor, atualize as variáveis TIME_COL e EVENT_COL no código com os nomes exatos das colunas do seu arquivo.")
        sys.exit(1)
        
    # Normalizando nomes de colunas (caso o índice esteja sem nome)
    if df.index.name is None:
        df.index.name = 'Amostra'
    df = df.reset_index(drop=True)

    # 4. Definir os NOVOS grupos de genes conforme solicitado
    # Grupo 1: CCL11, CCL24, CCL26
    GENES_G1 = ['CCL11', 'CCL24', 'CCL26']
    GROUP_COL_1 = 'Grupo_CCL11_CCL24_CCL26'
    FILE_G1 = 'Sobrevida_Grupo1_CCL11_CCL24_CCL26.png'
    
    # Grupo 2: CCL11, CCL26
    GENES_G2 = ['CCL11', 'CCL26']
    GROUP_COL_2 = 'Grupo_CCL11_CCL26'
    FILE_G2 = 'Sobrevida_Grupo2_CCL11_CCL26.png'
    
    
    # --- GRÁFICO 1: Análise do Grupo 1: [CCL11, CCL24, CCL26] ---
    print(f"\n--- Processando Grupo 1: {GENES_G1} ---")
    df_g1 = create_survival_groups(df.copy(), GENES_G1, GROUP_COL_1)
    if df_g1 is not None:
        title_g1 = f"Curvas de Sobrevida de Kaplan-Meier\nExpressão Agregada do Grupo: {', '.join(GENES_G1)}"
        # Passa o caminho do arquivo para salvar
        plot_survival(df_g1, GROUP_COL_1, ", ".join(GENES_G1), TIME_COL, EVENT_COL, title_g1, FILE_G1)
    
    # --- GRÁFICO 2: Análise do Grupo 2: [CCL11, CCL26] ---
    print(f"\n--- Processando Grupo 2: {GENES_G2} ---")
    df_g2 = create_survival_groups(df.copy(), GENES_G2, GROUP_COL_2)
    if df_g2 is not None:
        title_g2 = f"Curvas de Sobrevida de Kaplan-Meier\nExpressão Agregada do Grupo: {', '.join(GENES_G2)}"
        # Passa o caminho do arquivo para salvar
        plot_survival(df_g2, GROUP_COL_2, ", ".join(GENES_G2), TIME_COL, EVENT_COL, title_g2, FILE_G2)

if __name__ == "__main__":
    # Garante que o Tkinter não inicie a janela principal, apenas o filedialog
    try:
        main()
    except Exception as e:
        print(f"Ocorreu um erro no programa principal: {e}")