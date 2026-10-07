# Data

## Источник

**GSE173634** — Breast Cancer Single Cell Atlas  
Ссылка: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173634)

## Описание

Single-cell RNA-seq профилирование **35,276 клеток** из **32 клеточных линий** 
рака груди. Данные охватывают все основные подтипы:
- Luminal A
- Luminal B
- HER2-enriched
- Basal-like

Платформа: Illumina NovaSeq 6000  
Протокол: Drop-seq

## Файлы

| Файл | Размер | Описание |
|------|--------|----------|
| `GSE173634_RAW_UMI_counts.txt.gz` | 327 MB | Исходная матрица UMI-счетов |
| `GSE173634_Human_BreastCancer_32CellLines.h5ad` | ~1 GB | AnnData-объект (после конвертации) |
| `MT.csv` | < 1 KB | Список митохондриальных генов |

## Как скачать

1. Перейдите на страницу [GSE173634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173634).
2. В разделе **Supplementary file** скачайте `GSE173634_RAW_UMI_counts.txt.gz`.
3. Поместите файл в эту папку (`data/`).

## Конвертация в AnnData

Исходный TXT-файл (327 MB) не открывается стандартными текстовыми редакторами 
из-за размера. Используйте следующий скрипт:

```python
import scanpy as sc
import pandas as pd

# Загрузка UMI-счетов
counts = pd.read_csv("GSE173634_RAW_UMI_counts.txt.gz", sep="\t", index_col=0)

# Создание AnnData (транспонирование: клетки × гены)
adata = sc.AnnData(counts.T)
adata.var_names_make_unique()

# Сохранение
adata.write("GSE173634_Human_BreastCancer_32CellLines.h5ad")