from __future__ import annotations
from typing import TYPE_CHECKING

import pandas as pd

from .basic_cleaners import (
    apply_mutations_to_sequences,
)
from ..core.pipeline import multiout_step, pipeline_step

if TYPE_CHECKING:
    from typing import Dict, List, Optional, Tuple

__all__ = [
    "prepare_rbd_records",
    "add_reference_sequences_by_target",
    "apply_mutations_preserving_wild_type",
]


def __dir__() -> List[str]:
    """Return exported names."""

    return __all__


def _normalize_target_name(
    raw_target_name: str,
    target_name_aliases: Dict[str, str],
) -> str:
    """Resolve one raw target label to a canonical target name."""

    raw_target_name = str(raw_target_name).strip()
    return target_name_aliases.get(raw_target_name, raw_target_name)


def _prepare_rbd_source_table(
    dataset: pd.DataFrame,
    *,
    label_column: str,
) -> pd.DataFrame:
    """Run the shared raw-table preparation used by RBD source cleaners."""

    result = dataset.copy()
    result["target"] = result["target"].fillna("Wuhan-Hu-1").astype(str)
    result["variant_class"] = result["variant_class"].fillna("")

    for column in ("pass_pre_count_filter", "pass_ACE2bind_expr_filter"):
        if column in result.columns:
            result = result.loc[result[column].fillna(False)].copy()
    result["aa_substitutions"] = (
        result["aa_substitutions"].fillna("").astype(str).str.strip()
    )
    result[label_column] = pd.to_numeric(result[label_column], errors="coerce")
    return result.loc[
        result[label_column].notna()
        & ~result["aa_substitutions"].str.contains(r"\*", regex=True, na=False)
    ].copy()

@multiout_step(main="main", summary="summary")
def prepare_rbd_records(
    dataset: pd.DataFrame,
    *,
    mode: str,
    reference_sequences: Dict[str, str],
    target_name_aliases: Dict[str, str],
    column_mapping: Dict[str, str],
    label_column: str = "label",
    fallback_reference_sequence: Optional[str] = None,
    require_known_reference_sequence: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare RBD source tables for ACE2 or antibody cleaners."""

    result = dataset.rename(columns=column_mapping)
    raw_reference_id = (
        str(result["target"].fillna("Wuhan-Hu-1").iloc[0]).strip()
        if not result.empty
        else "Wuhan-Hu-1"
    )
    result = _prepare_rbd_source_table(
        result,
        label_column=label_column,
    )
    variant_class = result["variant_class"].astype(str).str.strip().str.lower()
    wt_mask = variant_class.eq("wildtype")
    result.loc[wt_mask, "aa_substitutions"] = ""
    keep_mask = wt_mask | ~(
        variant_class.isin(["synonymous", "stop"])
        | result["aa_substitutions"].str.contains("-", regex=False)
        | result["aa_substitutions"].eq("")
    )
    result = result.loc[keep_mask].copy()
    result["label"] = result[label_column]
    result["mut_info"] = result["aa_substitutions"]
    result.loc[wt_mask, "mut_info"] = "WT"
    summary = pd.DataFrame()

    if mode == "ace2":
        result["name"] = result["target"].map(
            lambda value: _normalize_target_name(
                raw_target_name=str(value),
                target_name_aliases=target_name_aliases,
            )
        )
        standardized_df = result[["name", "mut_info", "variant_class", "label"]]
    elif mode == "antibody":
        reference_id = _normalize_target_name(
            raw_target_name=raw_reference_id,
            target_name_aliases=target_name_aliases,
        )
        reference_sequence = str(
            reference_sequences.get(reference_id)
            or (
                fallback_reference_sequence
                if not require_known_reference_sequence
                else reference_sequences[reference_id]
            )
        ).strip()
        result["reference_id"] = reference_id
        result["antibody_name"] = result["name"].astype(str)
        standardized_df = result[["reference_id", "antibody_name", "mut_info", "label"]]
        summary = pd.DataFrame(
            [{"reference_id": reference_id, "reference_sequence": reference_sequence}]
        )
    return standardized_df.reset_index(drop=True), summary.reset_index(drop=True)


@pipeline_step
def add_reference_sequences_by_target(
    dataset: pd.DataFrame,
    reference_sequences: Dict[str, str],
    name_column: str = "name",
    sequence_column: str = "sequence",
    fallback_reference_sequence: Optional[str] = None,
) -> pd.DataFrame:
    """Attach reference sequences to standardized RBD rows."""

    result = dataset.copy()
    result[sequence_column] = result[name_column].map(reference_sequences)
    if fallback_reference_sequence is not None:
        result[sequence_column] = result[sequence_column].fillna(
            str(fallback_reference_sequence).strip()
        )
    return result


@multiout_step(main="success", failed="failed")
def apply_mutations_preserving_wild_type(
    dataset: pd.DataFrame,
    sequence_column: str = "sequence",
    name_column: str = "name",
    mutation_column: str = "mut_info",
    wt_identifier: str = "WT",
    mutation_sep: str = ",",
    is_zero_based: bool = True,
    sequence_type: str = "protein",
    num_workers: int = 4,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Materialize RBD mutant sequences while keeping WT rows intact."""

    mutation_text = dataset[mutation_column].fillna("").astype(str).str.upper()
    wt_mask = mutation_text.eq(str(wt_identifier).upper())

    wt_rows = dataset.loc[wt_mask].copy()
    wt_rows["mut_seq"] = wt_rows[sequence_column].astype(str)

    non_wt_rows = dataset.loc[~wt_mask].copy()
    if non_wt_rows.empty:
        return wt_rows.reset_index(drop=True), pd.DataFrame(
            columns=[*dataset.columns, "error_message"]
        )

    mutation_result = apply_mutations_to_sequences(
        non_wt_rows,
        sequence_column=sequence_column,
        name_column=name_column,
        mutation_column=mutation_column,
        mutation_sep=mutation_sep,
        is_zero_based=is_zero_based,
        sequence_type=sequence_type,
        num_workers=num_workers,
    )
    successful_rows = pd.concat([wt_rows, mutation_result.main], axis=0).sort_index(
        kind="stable"
    )
    return (
        successful_rows.reset_index(drop=True),
        mutation_result.side.get("failed", pd.DataFrame()).reset_index(drop=True),
    )
