from __future__ import annotations
from typing import TYPE_CHECKING

import pandas as pd

from ..core.pipeline import pipeline_step

if TYPE_CHECKING:
    from typing import Dict, List, Optional

__all__ = [
    "standardize_rbd_target_names",
    "standardize_rbd_ace2_records",
    "standardize_rbd_antibody_records",
    "add_reference_sequences_by_target",
]


def __dir__() -> List[str]:
    """Return exported names."""

    return __all__


@pipeline_step
def standardize_rbd_target_names(
    dataset: pd.DataFrame,
    target_name_aliases: Dict[str, str],
    name_column: str = "name",
) -> pd.DataFrame:
    """Canonicalize RBD target/reference names."""

    result = dataset.copy()

    def canonicalize(value: object) -> object:
        if pd.isna(value):
            return value
        stripped_value = str(value).strip()
        return target_name_aliases.get(stripped_value, stripped_value)

    result[name_column] = result[name_column].map(canonicalize)
    return result.reset_index(drop=True)


@pipeline_step
def standardize_rbd_ace2_records(
    dataset: pd.DataFrame,
    mutation_column: str = "mut_info",
    variant_class_column: str = "variant_class",
) -> pd.DataFrame:
    """Mark wild-type RBD-ACE2 records in the mutation column."""

    def _mark_wild_type_mutations(dataset: pd.DataFrame) -> pd.DataFrame:
        result = dataset.copy()
        variant_class = (
            result[variant_class_column].astype("string").str.strip().str.lower()
        )
        wt_mask = variant_class.eq("wildtype").fillna(False)
        result.loc[wt_mask, mutation_column] = "WT"
        return result.reset_index(drop=True)

    return _mark_wild_type_mutations(dataset)


@pipeline_step
def standardize_rbd_antibody_records(
    dataset: pd.DataFrame,
    mutation_column: str = "mut_info",
    variant_class_column: str = "variant_class",
) -> pd.DataFrame:
    """Mark wild-type RBD-antibody records in the mutation column."""

    def _mark_wild_type_mutations(dataset: pd.DataFrame) -> pd.DataFrame:
        result = dataset.copy()
        variant_class = (
            result[variant_class_column].astype("string").str.strip().str.lower()
        )
        wt_mask = variant_class.eq("wildtype").fillna(False)
        result.loc[wt_mask, mutation_column] = "WT"
        return result.reset_index(drop=True)

    return _mark_wild_type_mutations(dataset)


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
