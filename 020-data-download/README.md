# 020-data-download

Módulo responsável pelo download, organização e preparação de arquivos FASTQ brutos para análises transcriptômicas.
Suporta múltiplos datasets de forma modular, com configuração individual por projeto.

---

## Objetivo

Esta etapa realiza:

- download automatizado de FASTQ a partir de ENA/SRA;
- organização de arquivos por projeto;
- preparação de dados para QC e alinhamento.

---

## Estrutura

```bash
020-data-download/
├── download_fastq.sh
├── datasets/
│   ├── PRJEB14695/
│   │   ├── config.yaml
│   │   ├── ena-file-download-read_run-PRJEB14695-fastq_ftp.sh
│   │   └── metadata.csv
│   └── PRJXXXX/
│       ├── config.yaml
│       ├── ftp_links.txt
│       └── metadata.csv
├── logs/
```

---

## Filosofia

Cada dataset deve ser autocontido dentro de `datasets/<PROJECT_ID>/`.

Isso permite o isolamento entre projetos e fácil manutenção com escalabilidade para múltiplos estudos.

Cada dataset contém:

- `config.yaml`
- arquivo de links de download
- metadata associada (optional)

---

## Configuração

Exemplo de `config.yaml`:

```yaml
project:
  id: PRJEB14695
  organism: Schistosoma mansoni
  source: ENA
  accession: PRJEB14695

download:
  link_file: ena-file-download-read_run-PRJEB14695-fastq_ftp.sh
  protocol: ftp
  threads: 8
  retry: true
  continue_download: true

library:
  layout: paired
  compression: gz
  read_suffix:
    r1: "_1.fastq.gz"
    r2: "_2.fastq.gz"

validation:
  check_md5: false
  expected_files: auto

metadata:
  metadata_file: metadata.csv
  parser: ena
```

---

## Arquivos de links suportados

### Script original do ENA

Exemplo:

```bash
wget -nc ftp://ftp.sra.ebi.ac.uk/.../ERR1117826_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/.../ERR1117826_2.fastq.gz
```

### Lista simples

```text
ftp://ftp.sra.ebi.ac.uk/.../ERR1117826_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/.../ERR1117826_2.fastq.gz
```

O script extrai URLs automaticamente.

---

## Execução

### Execução local

```bash
bash download_fastq.sh PRJEB14695
```

---

## Output

Arquivos baixados:

```bash
/OUTPUTDIR/<PROJECT_ID>/fastq_ftp/
```

Exemplo:

```bash
/OUTPUTDIR/PRJEB14695/fastq_ftp/
├── ERR1117826_1.fastq.gz
├── ERR1117826_2.fastq.gz
└── ...
```

## Manifesto de renomeacao

Depois que `025-parse/030-metadata_final/AllProjects_metadata.csv` for gerado,
crie um manifesto run-level antes de renomear FASTQs:

```bash
python generate_rename_manifest.py \
  --metadata ../025-parse/030-metadata_final/AllProjects_metadata.csv \
  --project PRJEB14695 \
  --output PRJEB14695_rename_manifest.csv
```

O manifesto preserva o `run_accession` no nome novo, evitando colisoes quando
uma mesma amostra biologica possui varios runs tecnicos:

```text
ERR506074_1.fastq.gz -> SM_ADUL_OVAR_FEM_MIXS_SC1878518_B2_ERR506074_R1.fastq.gz
```

Logs:

```bash
020-data-download/logs/
├── fastq_download_<jobid>.log
└── fastq_download_<jobid>.err
```

---

## Adicionando novo dataset

1. Criar diretório:

```bash
mkdir datasets/PRJXXXX
```

2. Adicionar:

- `config.yaml`
- arquivo de links (`ftp_links.txt` ou script ENA)

3. Executar:

```bash
sbatch download_fastq.sh PRJXXXX
```

---

## Debug

Caso haja erro de configuração ou paths:

```bash
cat logs/fastq_download_<jobid>.log
cat logs/fastq_download_<jobid>.err
```

Verificar especialmente:

- `BASE_DIR`
- `DATASET_DIR`
- `CONFIG`
- `OUTPUT`

---
