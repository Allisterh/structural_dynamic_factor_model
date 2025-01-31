# Análise de Choques de Política Monetária usando SDFM

Este repositório contém o código e os dados utilizados no artigo das disciplinas de Macroeconomia II e Econometria II, ministradas pelos professores João Caldeira e Guilherme Moura na UFSC.

## Visão Geral

O projeto investiga os efeitos da política monetária sobre diversos preços de ativos no Brasil através de um Modelo de Fatores Dinâmicos Estrutural. A análise captura a resposta de diferentes ativos financeiros (ações, títulos, câmbio) a choques de política monetária, abordando o problema da não-fundamentalidade comum em modelos VAR tradicionais.

## Estrutura do Projeto

```
├── R/
│   ├── data_download/        # Scripts para coleta de dados
│   │   ├── bcb.R            # Download de dados do Banco Central
│   │   └── exchange.R       # Download de taxas de câmbio
│   ├── preprocessing/        # Limpeza e transformação de dados
│   │   ├── seasonality.R    # Testes e ajustes de sazonalidade
│   │   └── stationarity.R   # Testes de raiz unitária
│   ├── modeling/            # Scripts principais de modelagem
│   │   ├── factor_estimation.R     # Funções de estimação SDFM
│   │   ├── impulse_response.R      # Cálculo e plotagem de IRF
│   │   └── svensson_model.R        # Estimação da curva de juros
│   └── clean.R              # Script principal de preparação
├── data/
│   ├── raw/                 # Dados originais baixados
│   └── processed/           # Dados limpos e transformados
└── img/                     # Gráficos e figuras gerados
```

## Descrição Detalhada dos Scripts

### Scripts de Download de Dados (data_download/)

#### bcb.R
- Função `download_bcb_data()`: 
  - Download de séries temporais do BCB usando a API
  - Permite download paralelo para múltiplas séries
  - Organiza dados em formato wide
  - Parâmetros principais:
    - `id`: Códigos das séries do BCB
    - `start_date`: Data inicial
    - `end_date`: Data final
    - `parallel`: Opção de processamento paralelo

#### exchange.R
- Função `download_cambio()`:
  - Download de taxas de câmbio via Yahoo Finance
  - Calcula médias mensais
  - Formata dados para análise
  - Parâmetros:
    - `currency_tickers`: Códigos das moedas
    - `base_currency`: Moeda base (default: BRL)

### Scripts de Pré-processamento (preprocessing/)

#### seasonality.R
- Função `check_seasonality()`:
  - Testa sazonalidade em séries temporais
  - Usa testes QS, Friedman e Kruskal-Wallis
  - Retorna diagnóstico detalhado
  - Identifica variáveis que precisam de ajuste sazonal

#### stationarity.R
- Funções principais:
  - `adf_test()`: Implementa teste ADF com diferentes especificações
  - `remove_unit_root()`: Aplica diferenciação sequencial
  - Implementa testes de raiz unitária em painel
  - Retorna controle de transformações aplicadas

### Scripts de Modelagem (modeling/)

#### factor_estimation.R
- Implementa funções essenciais para SDFM:
  - `bai_ng_criteria()`: Determina número de fatores estáticos
  - `amengual_watson()`: Estima fatores dinâmicos
  - `estimate_dfm()`: Estima o modelo completo
  - Inclui correção de Kilian para VAR

#### impulse_response.R
- Funções para análise de impulso-resposta:
  - `compute_irf_dfm()`: Calcula IRFs do modelo
  - `plot_irf()`: Gera gráficos com bandas de confiança
  - Implementa bootstrap para intervalos de confiança

#### svensson_model.R
- Implementa modelo de Svensson para curva de juros:
  - `svensson_rate()`: Calcula taxas para maturidades específicas
  - `fit_svensson()`: Estima parâmetros do modelo
  - `generate_fixed_maturity_series()`: Gera séries de taxas fixas

### Script Principal (clean.R)
- Orquestra todo o processo de limpeza:
  - Carrega e organiza dados brutos
  - Aplica transformações necessárias
  - Trata sazonalidade e estacionariedade
  - Prepara dados para estimação do modelo

## Metodologia Detalhada

### 1. Preparação dos Dados

#### 1.1 Teste e Ajuste de Sazonalidade
- Aplicação de testes combinados:
  - Teste QS para periodicidade sazonal
  - Teste de Friedman para variação sazonal não-paramétrica
  - Teste de Kruskal-Wallis para diferenças sazonais
- Ajuste via procedimento X-11 ARIMA quando necessário

#### 1.2 Tratamento de Estacionariedade
- Testes de raiz unitária:
  - Teste ADF individual para cada série
  - Testes em painel (Maddala-Wu, Choi, Levin-Lin-Chu)
- Diferenciação sequencial quando necessário
- Controle de transformações aplicadas

#### 1.3 Construção da Curva de Juros
- Implementação do modelo de Svensson:
  - Estimação de parâmetros via otimização L-BFGS-B
  - Interpolação para maturidades fixas
  - Construção de séries temporais consistentes

### 2. Estimação do Modelo

#### 2.1 Determinação do Número de Fatores
- Fatores estáticos:
  - Critérios de informação de Bai-Ng (IC1, IC2, IC3)
  - Análise de scree plot
  - Decomposição da variância explicada

- Fatores dinâmicos:
  - Procedimento de Amengual-Watson
  - Análise de robustez com diferentes especificações

#### 2.2 Estimação do SDFM
1. Extração de fatores estáticos via PCA
2. Modelagem VAR dos fatores com correção de Kilian
3. Identificação de choques estruturais:
   - Decomposição espectral
   - Normalização de efeito unitário
   - Decomposição de Cholesky

#### 2.3 Análise de Impulso-Resposta
- Cálculo de IRFs estruturais
- Bootstrap wild com 800 replicações
- Construção de bandas de confiança
- Horizonte de análise de 50 períodos

### 3. Análise e Diagnóstico

#### 3.1 Testes de Robustez
- Análise de sensibilidade ao número de fatores
- Testes de estacionariedade em painel
- Diagnóstico de especificação do modelo

#### 3.2 Visualização de Resultados
- Gráficos de impulso-resposta com bandas de confiança
- Decomposição da variância
- Análise das cargas fatoriais

## Resultados Principais

- Queda imediata de ~3% no mercado acionário após choque contracionista
- Apreciação de 8% na taxa USD/BRL
- Aumentos significativos nos yields de diferentes maturidades
- Efeitos mais fortes em ativos mais arriscados e setor imobiliário
- Respostas assimétricas entre diferentes classes de ativos

## Instalação

1. Clone o repositório:
```bash
git clone https://github.com/seuperfil/choques-monetarios.git
```

2. Instale os pacotes R necessários:
```R
install.packages(c(
    "tidyverse", "GetBCBData", "tidyquant", "seasonal",
    "tseries", "zoo", "moments", "ggplot2", "patchwork"
))
```

## Como Usar

1. Coleta de Dados:
```R
source("R/data_download/bcb.R")
source("R/data_download/exchange.R")
```

2. Pré-processamento:
```R
source("R/clean.R")
```

3. Estimação do Modelo:
```R
source("R/modeling/factor_estimation.R")
```


## Referências 

- ALESSI, L.; KERSSENFISCHER, M. The response of asset prices to monetary policy shocks: Stronger than thought. Journal of Applied Econometrics, v. 34, n. 5, p. 661–672, 2019. Disponível em: <https://onlinelibrary.wiley.com/doi/abs/10.1002/jae.2706>.

- Stock, J. H., & Watson, M. W. (2016). Dynamic factor models, factor-augmented vector autoregressions, and structural vector autoregressions in macroeconomics. In Handbook of macroeconomics (Vol. 2, pp. 415-525). Elsevier.