Este repositório é uma coleção de scripts gerados para o estudo de monografia "Exploração de dados públcios para investigação in silico de genes relacionados a eosinófilos no câncer de mama." de Tamires Avila de Souza Clemente, defendido em 08 de Dezembro de 2025, pelo Instituto de Microbiologia Paulo de Góes, na Universidade Federal do Rio de Janeiro (UFRJ), campus Fundão. A descrição dos scripts usados se encontra a seguir, utilizados para gerar as figuras de resultados na seguinte ordem:

s1_mamanalysis.py
Usa as bibliotecas pandas, tkinter e os para importar e juntar os dados dos arquivos de dados de “gene expression RNAseq - IlluminalHiseq”, “phenotype - curated survival data” e “phenotype - Phenotypes” em um único arquivo .csv. 
O gene expression RNA - IlluminaHiseq contém dados de expressão gênica de câncer de mama invasivo, sendo obtidos de amostras de tumor sólido primário (código 01), tumor metastático (código 06) e tecido sólido normal (código 011) adjacentes ao tumor. Por sua vez, o dataset “phenotype - Curated Survival Data” incluí uma tabela curada para análise estatística de dados de sobrevivência, enquanto que “phenotypes - Phenotypes” contém informações sobre o estado clínico do paciente, como informações demográficas, de diagnóstico, patologia, sobrevida e tratamento. 

sample_type_plotandexcel.py
Usa as bibliotecas pandas, tkinter, os, matplotlib.pyplot, sci.py, numpy, datetime e math para fazer um arquivo excel e plotar um gráfico que contenha a contagem absoluta e porcentagem de pacientes correlacionado aos tipos de amostra. O código foi gerado com auxílio de IA (Gemini).
O código cria uma planilha com o nome “Resumo Estatístico” e o gráfico com a representação gráfica da distribuição tipo de amostra-paciente.

mamanalysis_PAM50_sample.py
Usa as bibliotecas matplotlib, numpy e pandas para plotar um gráfico que analisa a distribuição de amostras de acordo com os subtipos PAM50, selecionando o arquivo gerado em s1_mamanalysis.py.

s4_mamanalysis_plotgenesandeos_mama.py
Utiliza as bibliotecas pandas e matplotlib.pyplot para consruir gráficos de heatmap, colunas empilhadas com combinações de genes e colunas empilhadas com a expressão gênica simultânea do grupo de genes escolhido. Nesse caso, o grupo de genes selecionados foi CLC, EPX, IL5RA, PRG2.

survival_presence_CL.py
Utiliza as biliotecas pandas, lifelines e matplotlib para construção do gráfico de análise de sobrevida por Kaplan-Meier, comparando a alta expressão contra baixa expressão. Os grupos de alta e baixa expressão são definidos a partir da soma dos dados de genes selecionadospor paciente, que gera a expressão média agregada. Em seguida se faz um cálculo da mediana dessa expressão média agregada, e por fim se classifica as amostras em alta expressão para aquelas que forem maior que a mediana e baixa expressão para aquelas qu forem menor ou igual a mediana. Nesse script se analisou os genes CLC, EPX, IL5RA, PRG2, em que dois gráficos são gerados: um com CLC e outro sem CLC.

s3_mamanalysis.py
Utiliza as bibliotecas pandas, matplotlb, scipy, lifelines e seaborn para construir dois tipos de gráficos: o primeiro tipo se refere aos gráficos de sobrevida Kaplan-Meier comparando os grupos de expressão simultanea dos genes contra o grupo "demais", que é a combinação dos genes (3 genes, 2 genes, 1 gene, nenhum gene), e o segundo tipo se refere aos gráficos de frequência de expressão gênica por subtipo molecular de câncer de mama, em um gráfico de barras empilhado, seguido de gráficos da expressão indiidual de cada gene do grupo. Nesse script foram usados os genes CLC, EPX, IL5RA, PRG2.

mamanalysis_plotgeneseos_recuit.py
Segue a mesma lógica do script s4_mamanalysis_plotgenesandeos_mama.py, porém voltado para análise dos genes de recrutamento de eosinófilos: CCL11, CCL24, CCL26.

survival_recruit_CCL24.py
Segue a mesma logica do script survival_presence_CL.py, porém analisando os genes de recrutamento CCL11, CCL24, CCL26, em que um dos gráficos foi gerado sem CCL24.

s3_mamanalysis_survivalrecruit.py
Segue a mesma lógica do script s3_mamanalysis.py, porém analisando os genes de recrutamento CCL11, CCL24, CCL26.

mamanalysis_plotgeneseos_manutenção.py
Segue a mesma lógica do script s4_mamanalysis_plotgenesandeos_mama.py, porém analisa os genes de manutenção de sobrevida dos eosinófilos: IL5, IL33, IL25, TSLP.

survival_manut_IL5.py
Segue a mesma lógica do script survival_presence_CL.py, porém analisa os genes referentes à manutenção dos eosinófilos: IL5, IL33, IL25, TSLP, em que um dos gráficos é gerado sem IL5.

s3_mamanalysis_manutenção
Segue a mesma lógica do script s3_mamanalysis.py, porém analisando os genes de manutenção dos eosinófilos: IL5, IL33, IL25, TSLP.
