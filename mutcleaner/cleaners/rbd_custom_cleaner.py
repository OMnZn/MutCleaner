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
    "standardize_rbd_ace2_records",
    "standardize_rbd_antibody_records",
    "add_reference_sequences_by_target",
    "apply_mutations_preserving_wild_type",
]


def __dir__() -> List[str]:
    """Return exported names."""

    return __all__


@pipeline_step
def standardize_rbd_ace2_records(
    dataset: pd.DataFrame,
    target_name_aliases: Dict[str, str],
) -> pd.DataFrame:
    """Standardize RBD-ACE2 rows after generic extraction and filtering."""

    result = dataset.copy()
    result["target"] = result["target"].astype(str).str.strip()
    result["variant_class"] = result["variant_class"].fillna("")
    result["aa_substitutions"] = (
        result["aa_substitutions"].fillna("").astype(str).str.strip()
    )
    variant_class = result["variant_class"].astype(str).str.strip().str.lower()
    wt_mask = variant_class.eq("wildtype")
    result["name"] = result["target"].map(
        lambda value: target_name_aliases.get(value, value)
    )
    result["mut_info"] = result["aa_substitutions"]
    result.loc[wt_mask, "mut_info"] = "WT"
    return result[["name", "mut_info", "variant_class", "label"]].reset_index(
        drop=True
    )


@pipeline_step
def standardize_rbd_antibody_records(
    dataset: pd.DataFrame,
    target_name_aliases: Dict[str, str],
) -> pd.DataFrame:
    """Standardize RBD-antibody rows after generic extraction and filtering."""

    result = dataset.copy()
    result["target"] = result["target"].astype(str).str.strip()
    result["variant_class"] = result["variant_class"].fillna("")
    result["aa_substitutions"] = (
        result["aa_substitutions"].fillna("").astype(str).str.strip()
    )

    reference_ids = result["target"].map(
        lambda value: target_name_aliases.get(value, value)
    )
    unique_reference_ids = reference_ids.dropna().unique()
    if len(unique_reference_ids) != 1:
        raise ValueError(
            "RBD antibody source tables must contain exactly one target "
            f"reference, got {sorted(map(str, unique_reference_ids))}"
        )

    variant_class = result["variant_class"].astype(str).str.strip().str.lower()
    wt_mask = variant_class.eq("wildtype")
    result["reference_id"] = unique_reference_ids[0]
    result["antibody_name"] = result["name"].astype(str)
    result["mut_info"] = result["aa_substitutions"]
    result.loc[wt_mask, "mut_info"] = "WT"
    return result[["reference_id", "antibody_name", "mut_info", "label"]].reset_index(
        drop=True
    )


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


@multiout_step(main="successful", failed="failed")
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
