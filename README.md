# Modelo de Fator Dinâmico Estrutural em R

Este projeto implementa um modelo de fator dinâmico estrutural (SDFM - Structural Dynamic Factor Model) em R para analisar choques da politica monetária no preço dos ativos, seguindo a metodologia de Stock & Watson (2016). 

## Dados

- Dados brutos são baixados de múltiplas fontes em `download.R`:
  - Taxas de câmbio do Yahoo Finance
  - Várias séries do Banco Central do Brasil 
  - Índices de commodities
  - Dados de emprego do SIDRA
- Os dados são pré-processados em `clean.R`:
  - Sazonalidade é removida usando X-13ARIMA-SEATS
  - Variáveis nominais são transformadas em log
  - Testes de raiz unitária são realizados, diferenciação aplicada para alcançar estacionariedade
  - Dados processados finais salvos em `final_data.csv`

## Estimação do Modelo de Fatores

- `factor_estimation.R` contém funções para:
  - Determinar o número de fatores via critérios de Bai & Ng (2002)
  - Estimar o número de fatores dinâmicos via Amengual & Watson (2007)
  - Estimar o modelo de fator com restrições na matriz de cargas
- O modelo é estimado em `model.R`:
  - 6 fatores, 2 defasagens no VAR
  - Restrições de identificação impõem uma estrutura de blocos nas cargas
  - Separa atividade real, preços, variáveis monetárias
  - Estimado via mínimos quadrados iterativos
- Usa as normalizações de "unit effect" e "named factor" de Stock & Watson para identificar os choques estruturais e escala dos fatores:
  - "Unit effect" (Seção 4.1.3): define que um choque estrutural unitário causa um aumento contemporâneo de uma unidade em uma variável observada específica. Equação (32) formaliza isso como H_jj = 1.  
  - "Named factor" (Seção 2.1.3.2): associa cada fator a uma variável específica, normalizando as cargas dos fatores. Equação (12) mostra a forma da matriz Lambda sob essa normalização.

## Análise Estrutural

- Choques estruturais identificados via restrições zero na matriz de impacto
- Normalização de "unit effect" usada para definir a escala
- Quando o número de fatores dinâmicos (q) < fatores estáticos (r), estima a matriz 'G'
- IRFs estruturais computadas por `compute_SIRF()` em `impulse_responde.R`, seguindo a equação (58) do artigo: 
  - X_t = Lambda * Phi(L)^(-1) * G * H * eps_t + e_t
  - Lambda é a matriz de cargas, Phi(L) são os coeficientes do VAR dos fatores, G é a matriz que relaciona os choques eta aos fatores, e H é a matriz de impacto
- Alguns bugs permanecem no cálculo da SIRF: 
  - Todos os choques têm o mesmo impacto nos fatores independente das restrições
  - O problema provavelmente está no cálculo da representação de médias móveis (coeficientes C_h na função `compute_SIRF()`)
  - Tentei calcular as irf com base no método utilizado pelo pacote `vars`

## Próximos Passos

- Debugar o cálculo da SIRF, focar na recursão para os coeficientes MA
- Considerar esquemas de identificação alternativos
- Reportar decomposições de variância
- Analisar transmissões de choques específicos de interesse

## Referências 

- Stock, J. H., & Watson, M. W. (2016). Dynamic factor models, factor-augmented vector autoregressions, and structural vector autoregressions in macroeconomics. In Handbook of macroeconomics (Vol. 2, pp. 415-525). Elsevier.
- Códigos em MATLAB dos autores: http://www.princeton.edu/~mwatson/ddisk/Stock_Watson_DFM_HOM_replication_files_20160312.zip