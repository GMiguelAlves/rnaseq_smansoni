# 025-parse — Processamento de Metadados com metaQC

Este modulo utiliza o **metaQC** para padronizar, validar e normalizar metadados de RNA‑seq provenientes do ENA, de tabelas fornecidas por autores ou de uma combinação de ambas as fontes.
O resultado é uma tabela de metadados canônica única, pronta para análises posteriores (DESeq2, organização de FASTQs, QC, integração multiestudo).

## Schema canônico

Após o processamento, todos os estudos compartilham exatamente as seguintes colunas:
sample_id, stage, sex, condition, replicate, batch, dataset, lane, run_accession

Todos os demais campos utilizados durante o processamento são mapeados para esse schema fixo, pronto para integração.

## Visão geral do pipeline
01_raw_metadata/
│
▼
metaqc validate
│
▼
015_intermediate/
│
(opcional) metaqc enrich
│
▼
metaqc parse
│
▼
02_parsed/
│
▼
metaqc merge
│
▼
03_final/


## Detalhamento dos componentes

### 1. `validate`
Validação inicial dos metadados brutos (normalmente vindos do ENA).

- Mapeamento de aliases de colunas (ex.: `Sample` → `sample_id`)
- Detecção de identificadores de amostra duplicados
- Verificação do schema mínimo
- Geração de relatório

**Saídas:**
- `metadata_clean.csv`
- `validation_report.txt`

### 2. `enrich` (opcional)
Integra os metadados do ENA com metadados fornecidos pelo autor para resolver problemas comuns:

- `sample_id` ausente ou inconsistente na tabela do ENA
- Anotações biológicas ausentes (condição, estágio, sexo)
- Múltiplas lanes de sequenciamento por amostra biológica
- Inconsistências entre corridas (runs)

**Saídas:**
- `*_enriched.csv`

### 3. `parse`
Etapa central do pipeline. Responsável por:

- Padronização semântica (ex.: estágio de vida, sexo, condição)
- Extração de características biológicas a partir de texto livre ou campos codificados
- Normalização de categorias
- **Construção do `sample_id` canônico**

#### Geração do `sample_id`
O `sample_id` **nunca** é herdado diretamente do ENA ou do autor.
Ele é sempre reconstruído deterministicamente utilizando o padrão:
SM_<ESTÁGIO><SEXO><CONDIÇÃO><REPLICATA><BATCH>


Exemplos:

| Entrada (nome bruto)                                           | `sample_id` parseado         |
|----------------------------------------------------------------|------------------------------|
| `X_Eggs_R1` (ENA, PRJEB32839)                                  | `SM_EGGS_MIX_UNKN_R1_B1`     |
| `M_Cercariae_R2` (ENA, PRJEB32839)                             | `SM_CERC_MAL_UNKN_R2_B1`     |
| `Schistosoma_mansoni__testes_from_single_sex_infection-sc-1878522` (metadados do autor) | `SM_ADUL_MAL_SING_RNA_B2` |

**Regras:**
- `stage`, `sex` e `condition` são normalizados a partir dos campos fonte configurados.
- `replicate` recebe o prefixo `R` (ex.: `R1`, `R2`).
- `batch` é preservado conforme informado ou atribuído via configuração YAML.
- Valores ausentes recebem um valor de fallback controlado (ex.: `UNKN` para desconhecido).

#### Tratamento de múltiplas lanes
Quando uma amostra biológica foi sequenciada em várias lanes, o parser gera uma tabela no formato longo (long format):

| sample_id | lane | run_accession |
|-----------|------|---------------|
| …         | L1   | ERRxxxx       |
| …         | L2   | ERRyyyy       |

Isso preserva o mapeamento entre cada corrida e sua lane, mantendo um único `sample_id` canônico.

### 4. `merge`
Combina vários estudos já parseados em uma única tabela final.

- Concatenação vertical de todos os estudos
- Imposição do conjunto canônico de colunas
- Validação de consistência entre os datasets
- Remoção opcional de corridas não mapeadas (`--drop-unmatched`)

**Saída:** `03_final/merged_metadata.csv`

## Configuração (YAML)

Cada estudo pode ter seu próprio arquivo YAML de configuração. Exemplo:

```yaml
columns:
  run_field: run_accession           # coluna que contém o identificador da corrida
  sample_field: sample_alias         # coluna usada para parsing

defaults:
  dataset: PRJEB32839                # identificador curto do estudo
  batch: B1                          # batch padrão, caso não seja encontrado

parsing:
  stage:
    source: Sample Name              # campo de onde extrair o estágio
    regex_map:
      eggs: eggs
      cercariae: cercariae
  sex:
    source: Sample Name
    regex_map:
      "^X_": mixed
      "^M_": male
      "^F_": female

```
Para enriquecimento (metadados do autor):

```yaml
enrichment:
  join_key: run_accession
  lane_expansion: true
  column_mapping:
    author_sample_id: sample_id
    author_condition: condition
  cleanup:
    drop_unmatched: false
```

Qualidade e validação
Após o parsing, cada linha é verificada:

Todas as colunas obrigatórias devem existir ('sample_id, stage, sex, condition, replicate, batch, dataset, lane, run_accession')

sample_id nunca pode estar vazio

O campo dataset é obrigatório para rastreabilidade

Cenários encontrados
Cenário Etapas  Observações
Apenas ENA      validate → parse        sample_id construído integralmente a partir dos campos do ENA
Apenas autor    enrich → parse  Tabela do autor é obrigatória; enriquecimento é mandatório
Híbrido (ENA + autor + lanes)   validate → enrich → melt → parse        Combina todas as funcionalidades para máxima completude

Uso
```` bash
# 1. Validar um arquivo bruto de metadados
metaqc validate PRJEB32839_raw.csv --schema schemas/meu_estudo.yaml

# 2. Opcional: enriquecer com metadados do autor
metaqc enrich PRJEB32839_clean.csv --author author_meta.csv --config configs/PRJEB32839_enrich.yaml

# 3. Parsear para o schema canônico
metaqc parse PRJEB32839_enriched.csv --config configs/PRJEB32839_parse.yaml

# 4. Mesclar vários estudos parseados
metaqc merge 02_parsed/ -o 03_final/merged.csv
```

Estrutura de diretórios
025-parse/
├── 01_raw_metadata/          # Arquivos originais de metadados
├── 015_intermediate/         # Arquivos limpos/enriquecidos
├── 02_parsed/                # Tabelas canônicas por estudo
├── 03_final/                 # Metadados finais mesclados
├── configs/                  # Configurações YAML por estudo
│   ├── PRJEB32839.yaml
│   └── PRJEB14695.yaml
├── logs/                     # Relatórios de validação e parsing
└── README.md

Dependências
metaQc v0.2 ou superior
