# Análise de Choques de Política Monetária sobre Preços de Ativos no Brasil

Este repositório implementa a metodologia de **Alessi & Kerssenfischer (2019)** para analisar os efeitos da política monetária sobre preços de ativos no Brasil, utilizando um **Modelo de Fatores Dinâmicos Estrutural (SDFM)**.

## Objetivo

Investigar como choques de política monetária afetam diferentes classes de ativos financeiros brasileiros (ações, títulos, câmbio, spreads) através de uma abordagem de alta dimensionalidade que supera as limitações dos modelos VAR tradicionais.

## Metodologia

### Modelo de Fatores Dinâmicos Estrutural (SDFM)

Seguindo **Alessi & Kerssenfischer (2019)**, implementamos um SDFM baseado na padronização de **Barigozzi, Lippi & Luciani (2016)** que:

1. **Extrai fatores estáticos** de um painel de 70+ variáveis econômicas brasileiras
2. **Modela a dinâmica** dos fatores através de um VAR com correção de viés de Kilian
3. **Identifica choques estruturais** de política monetária 
4. **Calcula funções de impulso-resposta** para todas as variáveis do sistema

### Vantagens sobre Modelos VAR Tradicionais

- **Alta dimensionalidade**: Incorpora informação de 70+ variáveis vs 3-8 no VAR
- **Robustez**: Menos sensível à seleção específica de variáveis
- **Riqueza informacional**: Captura dinâmicas de múltiplas classes de ativos simultaneamente
- **Não-fundamentalidade**: Fatores capturam informação não observável aos agentes

## Estrutura do Projeto

```
├── R/                       # Funções de estimação e análise
│   ├── data_download/       # Coleta de dados (BCB, B3, ANBIMA)
│   ├── modeling/           # Estimação SDFM e cálculo de IRFs
│   └── preprocessing/      # Análise de dados
├── script/                 # Scripts principais de execução
│   ├── model_alessi.R     # Estimação completa do modelo
│   ├── clean.R            # Preparação dos dados
│   └── download.R         # Coleta automática de dados
├── data/                   # Dados brutos e processados
│   ├── raw/               # Dados originais (BCB, B3, ANBIMA)
│   └── processed/         # Dados limpos para análise
└── img/                   # Gráficos e resultados
```

## Base de Dados

### Variáveis Utilizadas (70+ séries)

**Setor Real**
- Indicadores de consumo, vendas e produção industrial
- Utilização da capacidade instalada
- Mercado de trabalho

**Preços**
- IPCA e componentes setoriais
- Índice de Preços ao Produtor (IPP)
- Preços de commodities (agro e energia)

**Setor Externo**
- Taxas de câmbio (USD, EUR, GBP, ARS, JPY)
- Preços internacionais de commodities

**Mercado Financeiro**
- **Juros**: Selic, DI, yields (3M, 1A, 2A, 3A, 5A)
- **Spreads**: Corporativo e bancário
- **Ações**: IBRx-100 e índices setoriais
- **Crédito**: Volume e taxas
- **Expectativas**: Índice de Decisão de Aplicação (IDA)

## Implementação da Metodologia

### Padronização BLL (Barigozzi, Lippi & Luciani, 2016)

O modelo segue a padronização específica proposta por Barigozzi, Lippi & Luciani (2016) e aplicada por Alessi & Kerssenfischer (2019):

1. **Cálculo do desvio padrão** das primeiras diferenças
2. **Remoção de tendência linear** dos dados em nível
3. **Normalização** pelos desvios padrão das diferenças
4. **Extração de fatores** via decomposição de componentes principais

### Estimação do Modelo

**Etapa 1: Seleção do Número de Fatores**
- Fatores estáticos: Critérios de Bai & Ng (2002)
- Fatores dinâmicos: Método de Amengual & Watson (2007)

**Etapa 2: Estimação Sequencial**
- Extração de fatores estáticos via PCA
- Modelagem VAR dos fatores com correção de Kilian (1998)
- Identificação estrutural via decomposição espectral
- Cálculo de funções de impulso-resposta

**Etapa 3: Análise de Robustez**
- Verificação de estabilidade do sistema
- Análise de sensibilidade ao número de fatores
- Diagnósticos de especificação

## Principais Funcionalidades

### Scripts de Análise

**`model_alessi.R`** - Script principal
- Implementa a metodologia completa de Alessi & Kerssenfischer (2019)
- Estima o modelo SDFM com diagnósticos automáticos
- Gera análises de robustez e visualizações

**`factor_estimation.R`** - Funções centrais
- Padronização BLL conforme Barigozzi, Lippi & Luciani (2016)
- Estimação de fatores estáticos e dinâmicos
- Correção de viés de Kilian (1998)
- Cálculo de funções de impulso-resposta

**`impulse_responde.R`** - Análise de respostas
- Funções de impulso-resposta estruturais
- Bandas de confiança via bootstrap
- Visualizações customizáveis

### Scripts de Dados

**`download.R`** - Coleta automática
- Download via APIs oficiais (BCB, Yahoo Finance)
- Dados de expectativas (ANBIMA)
- Índices setoriais (B3)

**`clean.R`** - Preparação
- Limpeza e organização dos dados
- Aplicação da padronização BLL
- Controle de qualidade

## Tratamento dos Dados

### Abordagem Metodológica

Seguindo **Alessi & Kerssenfischer (2019)**, a preparação dos dados prioriza:

- **Padronização BLL**: Método específico que combina detrending linear com normalização pelas diferenças
- **Sem ajustes adicionais**: Não são aplicados testes de raiz unitária individuais ou ajustes sazonais explícitos
- **Confiança na metodologia**: A padronização BLL é robusta para dados não-estacionários
- **Transformações pós-estimação**: Aplicadas apenas na interpretação das IRFs

### Características dos Dados

**Período**: 2010-2024 (frequência mensal)
**Dimensão**: 70+ variáveis econômicas brasileiras
**Fontes**: Banco Central do Brasil, B3, ANBIMA
**Cobertura**: Setor real, preços, mercado financeiro, setor externo

## Resultados Principais

### Impactos de um Choque Contracionista (50 bp)

**Mercado Acionário**
- Queda imediata de ~3% no IBRx-100
- Efeitos mais pronunciados em setores cíclicos
- Recuperação gradual ao longo de 12-18 meses

**Mercado de Câmbio**
- Apreciação de ~8% do Real frente ao Dólar
- Transmissão rápida (1-2 meses)
- Efeitos persistentes

**Mercado de Renda Fixa**
- Aumentos em yields de todas as maturidades
- Curva de juros se inclina (maior impacto no longo prazo)
- Ampliação de spreads corporativos

**Setor Imobiliário**
- Respostas mais fortes e persistentes
- Reflexo da sensibilidade aos juros
- Impactos duradouros nos preços

### Propriedades do Modelo Estimado

- **Estabilidade**: Sistema VAR estável (autovalores < 1)
- **Capacidade explicativa**: 56% da variância pelos fatores estáticos
- **Ortogonalidade**: Fatores dinâmicos apropriadamente identificados
- **Robustez**: Resultados consistentes em diferentes especificações

## Instalação e Uso

### Instalação
```r
# Instale os pacotes necessários
install.packages(c("MASS", "readr", "dplyr", "ggplot2"))
```

### Execução
```r
# Carregar funções principais
source("R/modeling/factor_estimation.R")
source("script/model_alessi.R")

# Executar modelo
resultado <- main_sdfm()
```

## Estrutura dos Dados

### Variáveis Incluídas (70+ séries)

**Setor Real**
- Produção industrial, vendas, consumo
- Utilização da capacidade instalada
- Indicadores do mercado de trabalho

**Preços**
- IPCA e componentes, IPP
- Preços de commodities

**Setor Externo**
- Taxas de câmbio principais
- Commodities internacionais

**Mercado Financeiro**
- Taxas de juros (Selic, DI, yields)
- Spreads corporativo e financeiro
- Índices acionários setoriais
- Indicadores de expectativas

## Referências

**Alessi, L. & Kerssenfischer, M.** (2019). The response of asset prices to monetary policy shocks: Stronger than thought. *Journal of Applied Econometrics*, 34(5), 661–672.

**Barigozzi, M., Lippi, M. & Luciani, M.** (2016). Non-stationary dynamic factor models for large datasets. *Federal Reserve Bank of New York Staff Reports*, no. 741.

**Bai, J. & Ng, S.** (2002). Determining the Number of Factors in Approximate Factor Models. *Econometrica*, 70(1), 191-221.

**Kilian, L.** (1998). Small-Sample Confidence Intervals for Impulse Response Functions. *Review of Economics and Statistics*, 80(2), 218-230.

---

**Implementação**: Gabriel Arruda  
**Disciplina**: Macroeconomia II  
**Professor**: João Caldeira  
**Universidade Federal de Santa Catarina (UFSC)**  
**2024**