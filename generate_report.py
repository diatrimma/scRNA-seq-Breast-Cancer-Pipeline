"""
HTML Report Generator for scRNA-seq Analysis
=============================================

Скрипт принимает результаты выполнения pipeline.py и генерирует 
интерактивный HTML-отчёт со всеми визуализациями и таблицами.

Usage:
    python generate_report.py --adata results/adata_processed.h5ad \
                              --markers results/top_markers_leiden.csv \
                              --output results/scRNA_seq_analysis_report.html

Author: Anastasia Popova
Date: 2025
"""

import argparse
import base64
import logging
from datetime import datetime
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def fig_to_base64(fig) -> str:
    """
    Преобразует matplotlib-фигуру в base64-строку для встраивания в HTML.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Фигура для конвертации.

    Returns
    -------
    str
        Base64-кодированное изображение в формате PNG.
    """
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=100)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def generate_qc_plots(adata: sc.AnnData) -> dict:
    """
    Генерирует QC-визуализации: violin plots и scatter plots.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData с QC-метриками в .obs.

    Returns
    -------
    dict
        Словарь с base64-изображениями.
    """
    images = {}

    # Violin plots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    metrics = ["n_genes", "n_counts", "percent_mito"]
    titles = [
        "Количество генов",
        "Количество counts",
        "Процент митохондриальных генов",
    ]
    for i, (metric, title) in enumerate(zip(metrics, titles)):
        if metric in adata.obs.columns:
            sns.violinplot(y=adata.obs[metric], ax=axes[i])
            axes[i].set_title(title)
            axes[i].set_ylabel("")
    plt.tight_layout()
    images["qc_violin"] = fig_to_base64(fig)
    plt.close(fig)

    # Scatter: n_counts vs percent_mito
    if "n_counts" in adata.obs.columns and "percent_mito" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(adata.obs["n_counts"], adata.obs["percent_mito"], alpha=0.5)
        ax.set_xlabel("n_counts")
        ax.set_ylabel("percent_mito")
        ax.set_title("n_counts vs percent_mito")
        images["qc_scatter1"] = fig_to_base64(fig)
        plt.close(fig)

    # Scatter: n_counts vs n_genes
    if "n_counts" in adata.obs.columns and "n_genes" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(adata.obs["n_counts"], adata.obs["n_genes"], alpha=0.5)
        ax.set_xlabel("n_counts")
        ax.set_ylabel("n_genes")
        ax.set_title("n_counts vs n_genes")
        images["qc_scatter2"] = fig_to_base64(fig)
        plt.close(fig)

    return images


def generate_pca_plot(adata: sc.AnnData) -> str:
    """Генерирует PCA-визуализацию с окраской по кластерам."""
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.pca(adata, color=["leiden"], ax=ax, show=False)
    img = fig_to_base64(fig)
    plt.close(fig)
    return img


def generate_doublet_plot(adata: sc.AnnData) -> str:
    """Генерирует гистограмму doublet scores."""
    fig, ax = plt.subplots(figsize=(8, 6))
    if "doublet_score" in adata.obs.columns:
        ax.hist(adata.obs["doublet_score"], bins=50, alpha=0.7, color="blue")
        ax.set_xlabel("Doublet score")
        ax.set_ylabel("Частота")
        ax.set_title("Распределение doublet scores")
        ax.grid(True, alpha=0.3)
    img = fig_to_base64(fig)
    plt.close(fig)
    return img


def create_html_report(
    adata: sc.AnnData,
    markers_df: pd.DataFrame,
    images: dict,
    qc_stats: dict = None,
    doublet_stats: dict = None,
) -> str:
    """
    Создаёт HTML-отчёт со всеми визуализациями и таблицами.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData с кластерами.
    markers_df : pd.DataFrame
        Таблица маркерных генов.
    images : dict
        Словарь с base64-изображениями.
    qc_stats : dict, optional
        Статистика QC.
    doublet_stats : dict, optional
        Статистика двойников.

    Returns
    -------
    str
        HTML-содержимое отчёта.
    """
    # Статистика по кластерам
    cluster_stats = (
        adata.obs.groupby("leiden")
        .agg(
            {
                "n_genes": ["mean", "std"],
                "n_counts": ["mean", "std"],
                "percent_mito": ["mean", "std"],
            }
        )
        .round(2)
    )
    cluster_stats_html = cluster_stats.to_html(classes="table table-striped")

    # Топ-5 маркеров для каждого кластера
    top_markers_per_cluster = {}
    for cluster in sorted(adata.obs["leiden"].unique(), key=int):
        cluster_markers = markers_df[markers_df["group"] == cluster]
        top_5 = cluster_markers.head(5)[["names", "scores", "pvals_adj"]]
        top_markers_per_cluster[cluster] = top_5.to_html(
            classes="table table-sm", index=False
        )

    # Опциональные секции
    qc_section = ""
    if qc_stats:
        qc_section = f"""
        <div class="row">
            <div class="col-md-3 mb-3">
                <div class="metric-card">
                    <div class="metric-value">{qc_stats.get('initial_cells', '—')}</div>
                    <div class="metric-label">Клеток до фильтрации</div>
                </div>
            </div>
            <div class="col-md-3 mb-3">
                <div class="metric-card">
                    <div class="metric-value">{qc_stats.get('after_qc', '—')}</div>
                    <div class="metric-label">После QC</div>
                </div>
            </div>
            <div class="col-md-3 mb-3">
                <div class="metric-card">
                    <div class="metric-value">{qc_stats.get('after_doublets', '—')}</div>
                    <div class="metric-label">После удаления двойников</div>
                </div>
            </div>
            <div class="col-md-3 mb-3">
                <div class="metric-card">
                    <div class="metric-value">{adata.n_vars}</div>
                    <div class="metric-label">Генов</div>
                </div>
            </div>
        </div>
        """

    html = f"""<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Анализ одноклеточных данных RNA-Seq</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background-color: #f8f9fa;
            color: #333;
            line-height: 1.6;
        }}
        .header {{
            background: linear-gradient(135deg, #6a11cb 0%, #2575fc 100%);
            color: white;
            padding: 2rem 0;
            margin-bottom: 2rem;
            text-align: center;
        }}
        .card {{
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            margin-bottom: 1.5rem;
            border: none;
        }}
        .card-header {{
            background-color: #f1f8ff;
            border-bottom: 1px solid #e3f2fd;
            font-weight: 600;
            padding: 1rem 1.5rem;
        }}
        .card-body {{ padding: 1.5rem; }}
        .section-title {{
            color: #2c3e50;
            border-bottom: 2px solid #3498db;
            padding-bottom: 0.5rem;
            margin-bottom: 1.5rem;
        }}
        .plot-img {{
            width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
            margin-bottom: 1rem;
        }}
        .metric-card {{
            text-align: center;
            padding: 1rem;
            background-color: #fff;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
        }}
        .metric-value {{
            font-size: 2rem;
            font-weight: 700;
            color: #2980b9;
        }}
        .metric-label {{
            font-size: 0.9rem;
            color: #7f8c8d;
        }}
        .footer {{
            background-color: #2c3e50;
            color: white;
            padding: 1.5rem 0;
            margin-top: 2rem;
            text-align: center;
        }}
    </style>
</head>
<body>
    <div class="header">
        <div class="container">
            <h1>Анализ одноклеточных данных RNA-Seq</h1>
            <p class="lead">Отчёт по анализу данных рака молочной железы (GSE173634)</p>
        </div>
    </div>

    <div class="container">
        <div class="card">
            <div class="card-header">Общая информация</div>
            <div class="card-body">
                <h3 class="section-title">Статистика данных</h3>
                {qc_section}
            </div>
        </div>

        <div class="card">
            <div class="card-header">Контроль качества (QC)</div>
            <div class="card-body">
                <h3 class="section-title">Распределение QC метрик</h3>
                <img src="data:image/png;base64,{images.get('qc_violin', '')}" class="plot-img" alt="QC violin">
                <div class="row">
                    <div class="col-md-6">
                        <h4>n_counts vs percent_mito</h4>
                        <img src="data:image/png;base64,{images.get('qc_scatter1', '')}" class="plot-img" alt="QC scatter 1">
                    </div>
                    <div class="col-md-6">
                        <h4>n_counts vs n_genes</h4>
                        <img src="data:image/png;base64,{images.get('qc_scatter2', '')}" class="plot-img" alt="QC scatter 2">
                    </div>
                </div>
            </div>
        </div>

        <div class="card">
            <div class="card-header">Анализ двойных клеток</div>
            <div class="card-body">
                <h3 class="section-title">Результаты Scrublet</h3>
                <img src="data:image/png;base64,{images.get('doublet_plot', '')}" class="plot-img" alt="Doublet plot">
            </div>
        </div>

        <div class="card">
            <div class="card-header">Кластеризация и визуализация</div>
            <div class="card-body">
                <h3 class="section-title">PCA анализ</h3>
                <img src="data:image/png;base64,{images.get('pca_plot', '')}" class="plot-img" alt="PCA plot">

                <h3 class="section-title">Статистика по кластерам</h3>
                {cluster_stats_html}

                <h3 class="section-title">Топ маркерные гены по кластерам</h3>
                <div class="accordion" id="markersAccordion">
"""

    for i, (cluster, table_html) in enumerate(top_markers_per_cluster.items()):
        collapsed = "collapsed" if i > 0 else ""
        show = "show" if i == 0 else ""
        expanded = "true" if i == 0 else "false"
        html += f"""
                    <div class="accordion-item">
                        <h2 class="accordion-header" id="heading{cluster}">
                            <button class="accordion-button {collapsed}" type="button"
                                    data-bs-toggle="collapse" data-bs-target="#collapse{cluster}"
                                    aria-expanded="{expanded}" aria-controls="collapse{cluster}">
                                Кластер {cluster}
                            </button>
                        </h2>
                        <div id="collapse{cluster}" class="accordion-collapse collapse {show}"
                             aria-labelledby="heading{cluster}" data-bs-parent="#markersAccordion">
                            <div class="accordion-body">
                                <div class="table-responsive">{table_html}</div>
                            </div>
                        </div>
                    </div>
"""

    html += f"""
                </div>
            </div>
        </div>
    </div>

    <div class="footer">
        <div class="container">
            <p>Анализ одноклеточных RNA-Seq данных | Сгенерировано: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
            <p>Использованы библиотеки: Scanpy, Scrublet, Matplotlib, Pandas, NumPy</p>
        </div>
    </div>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
"""
    return html


def main(args):
    """Основная логика генерации отчёта."""
    logger.info("Загрузка AnnData из %s", args.adata)
    adata = sc.read_h5ad(args.adata)

    logger.info("Загрузка маркеров из %s", args.markers)
    markers_df = pd.read_csv(args.markers)

    logger.info("Генерация графиков")
    images = generate_qc_plots(adata)
    images["pca_plot"] = generate_pca_plot(adata)
    images["doublet_plot"] = generate_doublet_plot(adata)

    logger.info("Создание HTML-отчёта")
    html = create_html_report(adata, markers_df, images)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    logger.info("Отчёт сохранён: %s", output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate HTML report for scRNA-seq analysis")
    parser.add_argument("--adata", type=str, required=True, help="Path to processed .h5ad")
    parser.add_argument("--markers", type=str, required=True, help="Path to markers CSV")
    parser.add_argument("--output", type=str, default="results/report.html", help="Output HTML path")
    args = parser.parse_args()
    main(args)