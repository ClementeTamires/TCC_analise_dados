# Script 2: Análise Detalhada de Genes (PAM50, Sobrevida, Combinações)

import pandas as pd
import os
import itertools
from tkinter import Tk, filedialog
from scipy.stats import kruskal, chi2_contingency
from lifelines import CoxPHFitter
import warnings

# --- Configurações Iniciais ---
print("Iniciando a análise detalhada...")

# Suprimir avisos comuns do lifelines e pandas
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# --- 0. Carregamento dos Dados (com Tkinter) ---

# Esconder a janela raiz do Tkinter
Tk().withdraw() 

print("Selecione o arquivo CSV de dados brutos")

# Pedir o arquivo de dados brutos gerado pelo Script 1
arquivo_csv_path = filedialog.askopenfilename(
    title="Selecione o arquivo CSV de dados brutos (ex: ..._combined_raw_data.csv)"
)

if not arquivo_csv_path:
    print("Seleção de arquivo cancelada. Encerrando o script.")
    exit()

print(f"Arquivo selecionado: {arquivo_csv_path}")

# Obter o nome da pasta e do arquivo para o Log
nome_arquivo = os.path.basename(arquivo_csv_path)
caminho_pasta = os.path.dirname(arquivo_csv_path)
nome_saida_excel = f"Analise_Genes_{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.xlsx"
caminho_saida_excel = os.path.join(caminho_pasta, nome_saida_excel)

print(f"Carregando {nome_arquivo}...")
try:
    df_merged = pd.read_csv(arquivo_csv_path, index_col=0)
except FileNotFoundError:
    print(f"ERRO: O arquivo '{arquivo_csv_path}' não foi encontrado.")
    exit()
except Exception as e:
    print(f"ERRO ao ler o arquivo CSV: {e}")
    exit()

print(f"Dados carregados com sucesso. Temos {df_merged.shape[0]} amostras e {df_merged.shape[1]} colunas.")

# --- 1. Preparação dos Dados e Definição de Grupos ---

# Estes são os nomes oficiais (HUGO) dos genes de interesse
# PRG2  -> MBP (Major Basic Protein)
# EPX   -> EPO (Eosinophil Peroxidase)
# CLC -> Gal10 (Galectin-10)
# IL5RA -> IL5R1 (IL-5 Receptor alpha)
genes_of_interest = ['PRG2', 'EPX', 'CLC', 'IL5RA']
genes_presentes = [gene for gene in genes_of_interest if gene in df_merged.columns]
genes_ausentes = [gene for gene in genes_of_interest if gene not in df_merged.columns]

if genes_ausentes:
    print(f"Atenção: Os seguintes genes não foram encontrados: {genes_ausentes}")
if not genes_presentes:
    print("ERRO: Nenhum dos genes de assinatura foi encontrado. Encerrando script.")
    exit()

print(f"Analisando os genes: {genes_presentes}")

# Definição das colunas de análise
col_pam50 = 'PAM50Call_RNAseq'
col_tempo = 'OS_Time_nature2012'
col_evento = 'OS_event_nature2012'

# Verificar se as colunas de análise existem
colunas_necessarias = {col_pam50: "PAM50", col_tempo: "Sobrevida (Tempo)", col_evento: "Sobrevida (Evento)"}
colunas_faltando = [col for col in colunas_necessarias if col not in df_merged.columns]

if colunas_faltando:
    print("\n--- ATENÇÃO ---")
    for col in colunas_faltando:
        print(f"A coluna '{col}' (necessária para {colunas_necessarias[col]}) não foi encontrada.")
    print("O script continuará, mas as análises dependentes dessas colunas falharão ou serão puladas.")
    print("-----------------\n")


# --- DEFINIÇÃO DE GRUPOS DE EXPRESSÃO (Binário: Expresso > 0) ---
# Esta é a etapa central para responder suas perguntas
# Assumimos que "expresso" significa valor de expressão > 0

print("Definindo grupos de expressão (expresso se > 0)...")
binary_cols = []
for gene in genes_presentes:
    col_name = f"{gene}_expresso"
    df_merged[col_name] = (df_merged[gene] > 0)
    binary_cols.append(col_name)

# Função para criar o nome do grupo de combinação
def get_expression_group(row):
    expressed_genes = []
    for gene in genes_presentes:
        if row[f"{gene}_expresso"]:
            expressed_genes.append(gene)
    
    count = len(expressed_genes)
    
    if count == 0:
        group_name = "Nenhum"
    else:
        # Ordena para garantir que CLC_IL5RA e IL5RA_CLC sejam o mesmo grupo
        group_name = "_".join(sorted(expressed_genes))
        
    return pd.Series([count, group_name], index=['contagem_genes_expressos', 'grupo_expressao'])

# Aplica a função para criar as colunas de grupo
df_merged[['contagem_genes_expressos', 'grupo_expressao']] = df_merged.apply(get_expression_group, axis=1)

print("Grupos de combinação definidos.")


# --- 2. Preparação do Arquivo Excel ---
print(f"Preparando o arquivo Excel de saída: '{caminho_saida_excel}'...")
writer_excel = pd.ExcelWriter(caminho_saida_excel, engine='openpyxl')
workbook = writer_excel.book

# Função auxiliar para escrever no Excel
def write_to_excel(df, sheet_name, startrow=0, startcol=0, header=True, index=True):
    # Função para ajustar a largura das colunas
    df.to_excel(writer_excel, sheet_name=sheet_name, startrow=startrow, startcol=startcol, header=header, index=index)
    worksheet = writer_excel.sheets[sheet_name]
    for i, col in enumerate(df.columns if header else []):
        col_width = max(len(str(col)), df[col].astype(str).str.len().max()) + 2
        worksheet.column_dimensions[chr(65 + i + startcol + (1 if index else 0))].width = col_width
    if index:
         idx_width = max(len(str(df.index.name)), df.index.astype(str).str.len().max()) + 2
         worksheet.column_dimensions[chr(65 + startcol)].width = idx_width

# --- Aba: Dados_PAM50 (Raw) ---
if col_pam50 in df_merged.columns:
    print("Escrevendo Aba: Dados_PAM50...")
    df_pam50_raw = df_merged.dropna(subset=[col_pam50] + genes_presentes)
    df_pam50_raw = df_pam50_raw[[col_pam50] + genes_presentes]
    write_to_excel(df_pam50_raw, 'Dados_PAM50', index=True)

# --- Aba: Dados_Sobrevida (Raw) ---
if col_tempo in df_merged.columns and col_evento in df_merged.columns:
    print("Escrevendo Aba: Dados_Sobrevida...")
    df_sobrevida_raw = df_merged.dropna(subset=[col_tempo, col_evento] + genes_presentes)
    df_sobrevida_raw = df_sobrevida_raw[[col_tempo, col_evento] + genes_presentes]
    write_to_excel(df_sobrevida_raw, 'Dados_Sobrevida', index=True)

# --- Aba: Analise_Estatistica (Contagens) ---
print("Escrevendo Aba: Analise_Estatistica (Contagens)...")
# Esta aba é um relatório, vamos construir as tabelas uma por uma
current_row = 0

# 1. Totais da Amostra
total_pacientes = len(df_merged)
df_total = pd.DataFrame({'Total de Pacientes': [total_pacientes]})
write_to_excel(df_total, 'Analise_Estatistica', startrow=current_row, index=False)
current_row += 3

# 2. Contagem por PAM50
if col_pam50 in df_merged.columns:
    df_pam50_counts = df_merged[col_pam50].value_counts().to_frame('Absoluto')
    df_pam50_counts['Porcentagem'] = (df_pam50_counts['Absoluto'] / df_pam50_counts['Absoluto'].sum() * 100).round(2)
    df_pam50_counts.index.name = "PAM50 Subtipo"
    write_to_excel(df_pam50_counts, 'Analise_Estatistica', startrow=current_row)
    current_row += df_pam50_counts.shape[0] + 3

# 3. Contagem por Sobrevida (Evento)
if col_evento in df_merged.columns:
    df_surv_counts = df_merged[col_evento].map({0:'Vivo', 1:'Morto'}).value_counts().to_frame('Absoluto')
    df_surv_counts['Porcentagem'] = (df_surv_counts['Absoluto'] / df_surv_counts['Absoluto'].sum() * 100).round(2)
    df_surv_counts.index.name = "Status Sobrevida (Evento)"
    write_to_excel(df_surv_counts, 'Analise_Estatistica', startrow=current_row)
    current_row += df_surv_counts.shape[0] + 3

# 4. Contagem por Grupos de Expressão (Geral)
df_group_counts = df_merged['grupo_expressao'].value_counts().to_frame('Absoluto')
df_group_counts['Porcentagem'] = (df_group_counts['Absoluto'] / total_pacientes * 100).round(2)
df_group_counts.index.name = "Grupo de Expressão (Geral)"
write_to_excel(df_group_counts, 'Analise_Estatistica', startrow=current_row)
current_row += df_group_counts.shape[0] + 3

# 5. Específico do PAM50 (Crosstab)
if col_pam50 in df_merged.columns:
    df_pam50_crosstab = pd.crosstab(df_merged['grupo_expressao'], df_merged[col_pam50])
    df_pam50_crosstab_pct = df_pam50_crosstab.apply(lambda r: (r/r.sum() * 100).round(1), axis=0) # Porcentagem por subtipo
    
    # Escreve o cabeçalho para esta seção
    pd.DataFrame(["Específico PAM50: Contagem Absoluta por Grupo de Expressão"]).to_excel(
        writer_excel, sheet_name='Analise_Estatistica', startrow=current_row, index=False, header=False)
    current_row += 1
    write_to_excel(df_pam50_crosstab, 'Analise_Estatistica', startrow=current_row)
    current_row += df_pam50_crosstab.shape[0] + 2

    pd.DataFrame(["Específico PAM50: Porcentagem por Subtipo (%)"]).to_excel(
        writer_excel, sheet_name='Analise_Estatistica', startrow=current_row, index=False, header=False)
    current_row += 1
    write_to_excel(df_pam50_crosstab_pct, 'Analise_Estatistica', startrow=current_row)
    current_row += df_pam50_crosstab_pct.shape[0] + 3

# 6. Específico da Sobrevida (Crosstab)
if col_evento in df_merged.columns:
    df_surv_crosstab = pd.crosstab(df_merged['grupo_expressao'], df_merged[col_evento].map({0:'Vivo', 1:'Morto'}))
    
    # Teste Chi-Square
    try:
        chi2, p_value, dof, expected = chi2_contingency(df_surv_crosstab)
        p_value_text = f"P-Valor Chi-Square (Associação entre Grupo e Status): {p_value:.4e}"
    except ValueError:
        p_value_text = "Teste Chi-Square falhou (provavelmente poucos dados)."
        
    pd.DataFrame(["Específico Sobrevida: Contagem Absoluta por Grupo de Expressão"]).to_excel(
        writer_excel, sheet_name='Analise_Estatistica', startrow=current_row, index=False, header=False)
    current_row += 1
    write_to_excel(df_surv_crosstab, 'Analise_Estatistica', startrow=current_row)
    current_row += df_surv_crosstab.shape[0] + 2
    
    pd.DataFrame([p_value_text]).to_excel(
        writer_excel, sheet_name='Analise_Estatistica', startrow=current_row, index=False, header=False)
    current_row += 3

# --- Aba: Correlacao_Continua (Nível de Expressão) ---
print("Escrevendo Aba: Correlacao_Continua (Nível de Expressão)...")
# Responde: "dentro dos que expressam CLC e IL5ra, quanto maior a expressão..."

results_continua = []

# Iterar por cada grupo de expressão (exceto 'Nenhum')
for group_name, group_df in df_merged[df_merged['contagem_genes_expressos'] > 0].groupby('grupo_expressao'):
    
    genes_no_grupo = group_name.split('_')
    
    for gene in genes_no_grupo:
        
        # 1. Relação com PAM50 (Kruskal-Wallis) - REMOVIDO
        
        # 2. Relação com Sobrevida (Cox Proportional Hazards)
        if col_tempo in group_df.columns and col_evento in group_df.columns:
            df_cox = group_df[[col_tempo, col_evento, gene]].dropna()
            if df_cox.shape[0] > 10 and df_cox[col_evento].sum() > 1: # Mínimo de dados
                try:
                    cph = CoxPHFitter()
                    cph.fit(df_cox, duration_col=col_tempo, event_col=col_evento, formula=gene)
                    summary = cph.summary
                    p_val_cox = summary.loc[gene, 'p']
                    hr = summary.loc[gene, 'exp(coef)']
                    ci = f"[{summary.loc[gene, 'exp(coef) lower 95%']:.2f}-{summary.loc[gene, 'exp(coef) upper 95%']:.2f}]"
                except Exception as e:
                    p_val_cox = pd.NA
                    hr = pd.NA
                    ci = str(e) # Registrar o erro
            else:
                p_val_cox = pd.NA
                hr = pd.NA
                ci = "Dados insuficientes"

            results_continua.append({
                'Grupo Expressão': group_name,
                'Gene Analisado': gene,
                'Teste': 'Nível Expressão vs Sobrevida',
                'Estatística': 'Cox PH P-Value',
                'Valor': p_val_cox,
                'Hazard Ratio (HR)': hr,
                'HR (IC 95%)': ci
            })

if results_continua:
    df_results_continua = pd.DataFrame(results_continua)
    write_to_excel(df_results_continua, 'Correlacao_Continua', index=False)
else:
    pd.DataFrame(["Nenhuma análise contínua foi executada."]).to_excel(
        writer_excel, sheet_name='Correlacao_Continua', index=False, header=False)

# --- Aba: Log ---
print("Escrevendo Aba: Log...")
log_data = [
    ("Data da Análise", pd.Timestamp.now().isoformat()),
    ("Arquivo de Entrada", nome_arquivo),
    ("Pasta de Entrada", caminho_pasta),
    ("Arquivo de Saída", caminho_saida_excel),
    ("Genes de Interesse", ", ".join(genes_of_interest)),
    ("Genes Encontrados", ", ".join(genes_presentes)),
    ("Genes Ausentes", ", ".join(genes_ausentes) if genes_ausentes else "Nenhum"),
    ("", ""),
    ("Metodologia de Análise", ""),
    ("Definição de 'Gene Expresso'", "Valor de expressão > 0 (usado nas abas 'Analise_Estatistica')"),
    ("Grupos de Expressão", "Baseado na combinação de quais genes tinham expressão > 0 para cada amostra."),
    ("Aba 'Analise_Estatistica'", "Contagens (absolutas e %) e Teste Chi-Square para associação entre 'grupos de expressão' e 'status de sobrevida'."),
    ("Aba 'Correlacao_Continua'", "Análise do NÍVEL de expressão (valor numérico) DENTRO de cada 'grupo de expressão'."),
    ("... vs Sobrevida", "Modelo de Riscos Proporcionais de Cox (CoxPHFitter)."),
    ("... (Interpretação HR)", "Hazard Ratio (HR) > 1 sugere maior risco (pior sobrevida) com o aumento da expressão; HR < 1 sugere menor risco (melhor sobrevida).")
]
df_log = pd.DataFrame(log_data, columns=['Item', 'Descrição'])
write_to_excel(df_log, 'Log', index=False)


# --- Finalização ---
try:
    writer_excel.close()
    print("\n--- Análise Concluída ---")
    print(f"Arquivo Excel '{caminho_saida_excel}' salvo com sucesso.")
except Exception as e:
    print(f"\nERRO ao salvar o Excel: {e}")
    print("Verifique se o arquivo já está aberto ou se você tem permissão para escrever no diretório.")