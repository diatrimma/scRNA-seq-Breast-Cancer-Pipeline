"""
scRNA-seq Analysis Pipeline for Breast Cancer
==============================================

Пайплайн для анализа single-cell RNA-seq данных (10x Genomics) 
из образцов рака молочной железы.

Шаги:
    1. Загрузка данных (h5ad)
    2. Контроль качества (QC)
    3. Удаление двойных клеток (Scrublet)
    4. Нормализация и отбор HVG
    5. PCA и кластеризация (Leiden)
    6. Дифференциальная экспрессия (Wilcoxon)

Author: Anastasia Popova
Date: 2025
"""

import argparse
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import yaml

warnings.filterwarnings("ignore")

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def load_config(config_path: str) -> dict:
    """Загружает YAML-конфиг с параметрами анализа."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def load_data(h5ad_path: str, mt_genes_path: str) -> sc.AnnData:
    """
    Загружает AnnData и удаляет митохондриальные гены.

    Parameters
    ----------
    h5ad_path : str
        Путь к файлу .h5ad с данными.
    mt_genes_path : str
        Путь к CSV со списком митохондриальных генов.

    Returns
    -------
    sc.AnnData
        Объект AnnData без митохондриальных генов.
    """
    logger.info("Загрузка данных из %s", h5ad_path)
    adata = sc.read_h5ad(h5ad_path)
    logger.info("Исходное число клеток: %d, генов: %d", adata.n_obs, adata.n_vars)

    # Удаление митохондриальных генов
    mt_genes = pd.read_csv(mt_genes_path, index_col=0, header=0)
    mask = ~adata.var_names.isin(mt_genes.squeeze())
    adata = adata[:, mask].copy()
    adata.var_names_make_unique()
    logger.info("После удаления MT-генов: %d генов", adata.n_vars)

    return adata


def compute_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """
    Вычисляет QC-метрики: n_genes, n_counts, percent_mito.

    Parameters
    ----------
    adata : sc.AnnData
        Объект AnnData.

    Returns
    -------
    sc.AnnData
        AnnData с добавленными QC-метриками в .obs.
    """
    logger.info("Вычисление QC-метрик")
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, "A1") else (adata.X > 0).sum(axis=1)
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)

    mito_genes = adata.var_names.str.startswith("MT-")
    if mito_genes.sum() > 0:
        mito_sum = adata.X[:, mito_genes].sum(axis=1)
        mito_sum_1d = mito_sum.A1 if hasattr(mito_sum, "A1") else mito_sum
        adata.obs["percent_mito"] = mito_sum_1d / adata.obs["n_counts"] * 100
    else:
        adata.obs["percent_mito"] = 0

    return adata


def filter_cells(adata: sc.AnnData, qc_params: dict) -> sc.AnnData:
    """
    Фильтрует клетки по QC-порогам.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData с QC-метриками.
    qc_params : dict
        Словарь с порогами: min_genes, max_genes, min_counts,
        max_counts, max_percent_mito.

    Returns
    -------
    sc.AnnData
        Отфильтрованный AnnData.
    """
    logger.info("Фильтрация клеток по QC-порогам")
    mask = (
        (adata.obs["n_genes"] >= qc_params["min_genes"])
        & (adata.obs["n_genes"] <= qc_params["max_genes"])
        & (adata.obs["n_counts"] >= qc_params["min_counts"])
        & (adata.obs["n_counts"] <= qc_params["max_counts"])
        & (adata.obs["percent_mito"] <= qc_params["max_percent_mito"])
    )
    n_before = adata.n_obs
    adata = adata[mask].copy()
    logger.info("Удалено клеток: %d, осталось: %d", n_before - adata.n_obs, adata.n_obs)
    return adata


def remove_doublets(adata: sc.AnnData) -> tuple:
    """
    Ищет и удаляет двойные клетки с помощью Scrublet.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData после QC-фильтрации.

    Returns
    -------
    tuple
        (adata_clean, doublet_stats, scrub) — очищенный AnnData,
        статистика и объект Scrublet.
    """
    logger.info("Поиск двойных клеток (Scrublet)")
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    doublets_found = int(predicted_doublets.sum())
    doublet_percentage = doublets_found / len(predicted_doublets) * 100

    adata = adata[~adata.obs["predicted_doublet"]].copy()
    logger.info("Найдено двойников: %d (%.2f%%)", doublets_found, doublet_percentage)

    doublet_stats = {
        "doublets_found": doublets_found,
        "doublet_percentage": doublet_percentage,
    }
    return adata, doublet_stats, scrub


def normalize_and_cluster(adata: sc.AnnData, n_top_genes: int = 2000) -> sc.AnnData:
    """
    Нормализация, отбор HVG, PCA, neighbors, Leiden-кластеризация.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData после удаления двойников.
    n_top_genes : int
        Число высоковариабельных генов.

    Returns
    -------
    sc.AnnData
        AnnData с кластерами в .obs['leiden'].
    """
    logger.info("Нормализация и отбор HVG")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)

    logger.info("PCA")
    sc.tl.pca(adata, use_highly_variable=True)

    logger.info("Кластеризация (Leiden)")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.leiden(adata)

    return adata


def find_markers(adata: sc.AnnData, groupby: str = "leiden") -> pd.DataFrame:
    """
    Дифференциальная экспрессия (Wilcoxon).

    Parameters
    ----------
    adata : sc.AnnData
        AnnData с кластерами.
    groupby : str
        Колонка в .obs для группировки.

    Returns
    -------
    pd.DataFrame
        Таблица с маркерными генами.
    """
    logger.info("Поиск маркерных генов (Wilcoxon)")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", key_added="wilcoxon")
    markers_df = sc.get.rank_genes_groups_df(adata, group=None, key="wilcoxon")
    return markers_df


def main(args):
    """Основной пайплайн."""
    config = load_config(args.config)

    # 1. Загрузка
    adata = load_data(config["data"]["h5ad_path"], config["data"]["mt_genes_path"])
    initial_cells = adata.n_obs

    # 2. QC
    adata = compute_qc_metrics(adata)
    adata = filter_cells(adata, config["qc"])
    after_qc = adata.n_obs

    # 3. Двойники
    adata, doublet_stats, scrub = remove_doublets(adata)
    after_doublets = adata.n_obs

    # 4. Нормализация и кластеризация
    adata = normalize_and_cluster(adata, config["clustering"]["n_top_genes"])

    # 5. Маркеры
    markers_df = find_markers(adata)

    # 6. Сохранение
    output_dir = Path(config["output"]["dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    markers_df.to_csv(output_dir / "top_markers_leiden.csv", index=False)
    adata.write(output_dir / "adata_processed.h5ad")
    logger.info("Результаты сохранены в %s", output_dir)

    # 7. Статистика
    qc_stats = {
        "initial_cells": initial_cells,
        "after_qc": after_qc,
        "after_doublets": after_doublets,
        **config["qc"],
    }

    logger.info("Пайплайн завершён успешно")
    return adata, qc_stats, doublet_stats, markers_df, scrub


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="scRNA-seq analysis pipeline")
    parser.add_argument("--config", type=str, default="config.yaml", help="Path to config file")
    args = parser.parse_args()
    main(args)