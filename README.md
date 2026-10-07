markdown
# scRNA-seq Analysis Pipeline for Breast Cancer

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Scanpy](https://img.shields.io/badge/scanpy-1.9+-green.svg)](https://scanpy.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Пайплайн для анализа single-cell RNA-seq данных, полученных из образцов 
рака молочной железы. Включает контроль качества, удаление двойных клеток, 
нормализацию, кластеризацию, аннотацию типов клеток и дифференциальную 
экспрессию.

## Biological Context

**Рак груди** — гетерогенное заболевание, которое остаётся самой 
распространённой формой рака среди женщин (25% всех случаев). 
Стандартная стратификация пациентов основана на статусе рецепторов (ER, PR, HER2) 
и гистологическом подтипе. Однако около 20% пациентов не имеют 
действенных биомаркеров, что затрудняет выбор терапии.

**Single-cell RNA-seq (scRNA-seq)** позволяет исследовать транскриптом 
отдельных клеток, выявляя гетерогенность опухоли. Это критично для понимания:
- Какие субпопуляции клеток присутствуют в опухоли.
- Как они различаются по экспрессии клинически значимых генов.
- Какие клетки могут быть чувствительны или устойчивы к терапии.

Данный пайплайн работает с **Breast Cancer Single Cell Atlas** 
(GSE173634) — набором из **35,276 клеток** из **32 клеточных линий**, 
охватывающих все основные подтипы рака груди (Luminal A/B, HER2-enriched, 
Basal-like). Атлас был создан для автоматической 
диагностики и предсказания лекарственной чувствительности.

## Pipeline Steps

```
┌─────────────────┐
│  Загрузка данных │  ← GSE173634 (h5ad)
└────────┬────────┘
         ↓
┌─────────────────┐
│  QC метрики     │  ← n_genes, n_counts, percent_mito
└────────┬────────┘
         ↓
┌─────────────────┐
│  Фильтрация     │  ← min/max genes, counts, mito%
└────────┬────────┘
         ↓
┌─────────────────┐
│  Scrublet       │  ← удаление двойных клеток
└────────┬────────┘
         ↓
┌─────────────────┐
│  Нормализация   │  ← log-normalization, HVG
└────────┬────────┘
         ↓
┌─────────────────┐
│  PCA + Leiden   │  ← кластеризация
└────────┬────────┘
         ↓
┌─────────────────┐
│  Wilcoxon DE    │  ← маркерные гены по кластерам
└────────┬────────┘
         ↓
┌─────────────────┐
│  HTML отчёт     │  ← визуализация результатов
└─────────────────┘
```

## Data

**Источник:** [GSE173634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173634) (GEO, NCBI)

**Описание:** Single-cell RNA-seq профилирование 35,276 клеток из 
32 клеточных линий рака груди. Данные получены с платформы Illumina 
NovaSeq 6000, протокол Drop-seq.

**Исходный формат:** `GSE173634_RAW_UMI_counts.txt.gz` (327 MB) — матрица UMI-счетов.

### Preprocessing

Исходный TXT-файл (327 MB) не открывается стандартными текстовыми редакторами 
из-за размера. Для конвертации в AnnData (`.h5ad`) использовался следующий подход:

```python
import scanpy as sc
import pandas as pd

# Загрузка UMI-счетов (разреженная матрица)
counts = pd.read_csv("GSE173634_RAW_UMI_counts.txt.gz", sep="\t", index_col=0)

# Создание AnnData
adata = sc.AnnData(counts.T)  # транспонирование: клетки × гены
adata.var_names_make_unique()

# Сохранение в h5ad
adata.write("GSE173634_Human_BreastCancer_32CellLines.h5ad")
Результат: AnnData-объект с 35,276 клетками и ~20,000 генов.
```

**Список митохондриальных генов: data/MT.csv — используется для удаления
MT-генов перед анализом.**

## Tools & Technologies
Python 3.10+

**Scanpy** — анализ scRNA-seq данных

**Scrublet** — удаление двойных клеток

**Pandas, NumPy** — работа с данными

**Matplotlib, Seaborn** — визуализация

**PyYAML** — конфигурация

### Installation
```bash
git clone https://github.com/diatrimma/scRNA-seq-Breast-Cancer-Pipeline.git
cd scRNA-seq-Breast-Cancer-Pipeline
pip install -r requirements.txt
```
### Usage
Скачайте данные GSE173634 и поместите .h5ad-файл в data/.

Настройте параметры в config.yaml.

Запустите пайплайн:

```bash
python pipeline.py --config config.yaml
```
Сгенерируйте HTML-отчёт:

```bash
python generate_report.py --results results/
Results
```
## Пайплайн генерирует:

**QC-визуализации:** violin plots, scatter plots до и после фильтрации.

**Анализ двойников:** гистограмма doublet scores.

**PCA-визуализация:** проекция клеток на главные компоненты.

**Кластеры Leiden:** таблица с топ-маркерными генами.

**HTML-отчёт:** интерактивный отчёт со всеми результатами.

## Skills Demonstrated
Работа с scRNA-seq данными (AnnData, Scanpy).

Контроль качества и фильтрация клеток.

Удаление двойных клеток (Scrublet).

Кластеризация (Leiden) и дифференциальная экспрессия (Wilcoxon).

Визуализация данных и генерация HTML-отчётов.

Воспроизводимость анализа (config-файлы, логирование).

**Author:** Анастасия Попова

**GitHub:** @diatrimma

**Email:** popova.ai@phystech.edu

**License:** MIT License
