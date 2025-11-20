# Script: Processamento de Dados Brutos (Versão com Tkinter)

import pandas as pd
import os
from tkinter import Tk, filedialog # Importar o necessário do Tkinter

# --- 0. Configurar o Tkinter e Obter informações do usuário ---

# Esconder a janela raiz do Tkinter
Tk().withdraw() 

print("Abrindo caixas de diálogo para seleção de arquivos...")

# Pedir os nomes dos arquivos usando o seletor de arquivos gráfico
print("Selecione o arquivo de RNAseq (.gz)")
arquivo_rnaseq = filedialog.askopenfilename(title="Selecione o arquivo de RNAseq (.gz)")
print("Selecione o arquivo Clinical")
arquivo_clinical = filedialog.askopenfilename(title="Selecione o arquivo Clinical (.txt, .tsv, etc.)")
print("Selecione o arquivo Phenotypes/Survival")
arquivo_phenotypes = filedialog.askopenfilename(title="Selecione o arquivo Phenotypes/Survival (.txt, .tsv, etc.)")

# Verificar se o usuário cancelou alguma seleção
if not all([arquivo_rnaseq, arquivo_clinical, arquivo_phenotypes]):
    print("\nSeleção de arquivo cancelada. Encerrando o script.")
    exit()

print(f"Arquivo RNAseq selecionado: {arquivo_rnaseq}")
print(f"Arquivo Clinical selecionado: {arquivo_clinical}")
print(f"Arquivo Phenotypes selecionado: {arquivo_phenotypes}")

# Pedir o nome do banco de dados para criar a pasta (ainda via console)
nome_banco_dados = input("\nDigite o nome do banco de dados (ex: TCGA-BRCA): ")

if not nome_banco_dados:
    print("Nome do banco de dados não fornecido. Encerrando.")
    exit()

# Criar a pasta para o banco de dados
try:
    os.makedirs(nome_banco_dados, exist_ok=True)
    print(f"Pasta '{nome_banco_dados}' criada/verificada com sucesso.")
except OSError as error:
    print(f"Erro ao criar a pasta '{nome_banco_dados}': {error}")
    # Encerra o script se não for possível criar a pasta
    exit()


# --- 1. Carregar os Três Arquivos ---

print("\nCarregando arquivos...")
try:
    # Carregar dados de expressão (RNAseq)
    df_expr = pd.read_csv(arquivo_rnaseq, 
                            sep='\t', 
                            index_col=0)
    print(f"Arquivo RNAseq carregado.")

    # Carregar dados clínicos
    df_pheno = pd.read_csv(arquivo_clinical, 
                             sep='\t', 
                             index_col=0)
    print(f"Arquivo Clinical carregado.")

    # Carregar dados de fenótipo/sobrevida
    df_surv = pd.read_csv(arquivo_phenotypes, 
                          sep='\t', 
                          index_col=0)
    print(f"Arquivo Phenotypes carregado.")

except FileNotFoundError as e:
    print(f"\nErro: Arquivo não encontrado.")
    print(f"Detalhe: {e}")
    print("Por favor, verifique se os nomes dos arquivos estão corretos e no mesmo diretório.")
    exit()
except Exception as e:
    print(f"\nOcorreu um erro ao ler os arquivos: {e}")
    exit()


# --- 2. Preparar e Juntar os Dados ---

print("\nProcessando e juntando os dados...")

# Transpor a matriz de expressão para que as amostras fiquem nas linhas
df_expr_T = df_expr.T

# Antes de juntar, removemos a coluna '_PATIENT' duplicada do df_surv, se ela existir.
# Vamos manter a coluna '_PATIENT' que já existe no df_pheno (que é o índice).
df_surv_cleaned = df_surv.copy()
if '_PATIENT' in df_surv_cleaned.columns:
    df_surv_cleaned = df_surv_cleaned.drop(columns=['_PATIENT'])
    print("Coluna '_PATIENT' removida do arquivo de phenotypes para evitar duplicata.")

# Agora, junte o df_pheno (clinical) com o df_surv_cleaned (phenotypes)
df_clinical_full = df_pheno.join(df_surv_cleaned, how='inner')

# Juntar os dados clínicos completos com os dados de expressão
df_merged = df_clinical_full.join(df_expr_T, how='inner')

if df_merged.empty:
    print("\nAVISO: O DataFrame final está vazio.")
    print("Isso geralmente acontece se os IDs das amostras (índice) não correspondem entre os arquivos.")
    print("Verifique se os arquivos são compatíveis e se 'index_col=0' está correto para seus dados.")
else:
    print("Dados juntados com sucesso.")


# --- 3. Salvar o Arquivo ---

# Definir o nome do arquivo de saída e o caminho completo
nome_arquivo_saida = f"{nome_banco_dados}_combined_raw_data.csv"
caminho_saida = os.path.join(nome_banco_dados, nome_arquivo_saida)

print(f"\nFormato do DataFrame final: {df_merged.shape}")

# Salvar o DataFrame combinado no caminho especificado
try:
    print(f"Salvando dados combinados em '{caminho_saida}'...")
    df_merged.to_csv(caminho_saida)
    print("Arquivo salvo com sucesso!")
except Exception as e:
    print(f"Ocorreu um erro ao salvar o arquivo: {e}")