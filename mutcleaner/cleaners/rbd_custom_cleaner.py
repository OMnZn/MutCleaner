from __future__ import annotations

import re
from typing import TYPE_CHECKING

import pandas as pd

from .basic_cleaners import apply_mutations_to_sequences
from ..core.pipeline import multiout_step, pipeline_step

if TYPE_CHECKING:
    from typing import Any, Dict, List, Optional, Tuple

__all__ = [
    "prepare_ace2_binding_records",
    "prepare_rbd_antibody_records",
    "add_reference_sequences_by_target",
    "apply_mutations_preserving_wild_type",
    "capture_rbd_antibody_variants",
]
RAW_ANTIBODY_REQUIRED_COLUMNS = ("name", "target", "label", "aa_substitutions")
MUTATION_PATTERN = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")


def __dir__() -> List[str]:
    """Return exported names."""

    return __all__


def _normalize_bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)

    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
        "yes": True,
        "no": False,
        "t": True,
        "f": False,
    }
    return series.astype(str).str.strip().str.lower().map(mapping).fillna(False)


def _normalize_target_name(
    raw_target_name: str,
    target_name_aliases: Dict[str, str],
    reference_sequences: Optional[Dict[str, str]] = None,
    allow_unknown: bool = False,
) -> str:
    """Resolve one raw target label to a canonical target name."""

    def normalize_key(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())

    normalized_target_name = str(raw_target_name).strip()
    if reference_sequences and normalized_target_name in reference_sequences:
        return normalized_target_name

    normalized_key = normalize_key(normalized_target_name)
    normalized_aliases = {
        normalize_key(alias): canonical_name
        for alias, canonical_name in target_name_aliases.items()
    }
    if normalized_key in normalized_aliases:
        return normalized_aliases[normalized_key]

    if reference_sequences:
        direct_matches = [
            canonical_name
            for canonical_name in reference_sequences
            if normalize_key(canonical_name) == normalized_key
        ]
        if len(direct_matches) == 1:
            return direct_matches[0]

    if allow_unknown or not reference_sequences:
        return normalized_target_name

    raise ValueError(
        f"Unknown RBD target '{raw_target_name}'. Add it to target_name_aliases "
        "or reference_sequences in the cleaner config."
    )


def _available_column_mapping(
    dataset: pd.DataFrame,
    column_mapping: Dict[str, str],
) -> Dict[str, str]:
    return {
        raw_name: standardized_name
        for raw_name, standardized_name in column_mapping.items()
        if raw_name in dataset.columns and raw_name != standardized_name
    }


def _select_reference_sequence(
    *,
    reference_id: str,
    raw_target_name: str,
    reference_sequences: Dict[str, str],
    fallback_reference_sequence: Optional[str],
    require_known_reference_sequence: bool,
    expected_reference_length: Optional[int],
) -> str:
    if reference_id in reference_sequences:
        reference_sequence = str(reference_sequences[reference_id]).strip()
    elif raw_target_name in reference_sequences:
        reference_sequence = str(reference_sequences[raw_target_name]).strip()
    elif not require_known_reference_sequence and fallback_reference_sequence:
        reference_sequence = str(fallback_reference_sequence).strip()
    else:
        raise ValueError(
            f"No reference sequence configured for target '{reference_id}'. "
            "Add it to reference_sequences or disable "
            "require_known_reference_sequence."
        )

    if (
        expected_reference_length is not None
        and int(expected_reference_length) > 0
        and len(reference_sequence) != int(expected_reference_length)
    ):
        raise ValueError(
            f"Reference sequence for target '{reference_id}' must have length "
            f"{expected_reference_length}, got {len(reference_sequence)}"
        )
    return reference_sequence


def _normalize_rbd_mutations(
    mutation_text: str,
    reference_sequence: str,
) -> List[str]:
    normalized_mutations: List[str] = []
    for token in str(mutation_text).split():
        match = MUTATION_PATTERN.fullmatch(token.strip())
        if match is None:
            raise ValueError(f"Unable to parse mutation token: {token}")

        _, site_text, mutant_aa = match.groups()
        site_zero_based = int(site_text) - 1
        if site_zero_based < 0 or site_zero_based >= len(reference_sequence):
            raise ValueError(
                f"Mutation position {site_text} is outside the reference RBD "
                f"length ({len(reference_sequence)}): {token}"
            )

        reference_wildtype_aa = reference_sequence[site_zero_based]
        if mutant_aa == reference_wildtype_aa:
            continue
        normalized_mutations.append(
            f"{reference_wildtype_aa}{site_zero_based}{mutant_aa}"
        )
    return normalized_mutations


@pipeline_step
def prepare_ace2_binding_records(
    dataset: pd.DataFrame,
    reference_sequences: Dict[str, str],
    target_name_aliases: Dict[str, str],
    column_mapping: Dict[str, str],
    target_column: str = "target",
    mutation_column: str = "aa_substitutions",
    label_column: str = "log10Ka",
    variant_class_column: str = "variant_class",
    default_target_name: str = "Wuhan-Hu-1",
    output_name_column: str = "name",
    output_mutation_column: str = "mut_info",
    output_label_column: str = "label",
) -> pd.DataFrame:
    """Standardize one ACE2 table into cleaner-ready columns."""

    result = dataset.copy()
    available_column_mapping = _available_column_mapping(result, column_mapping)
    if available_column_mapping:
        result = result.rename(columns=available_column_mapping)

    if mutation_column not in result.columns:
        raise ValueError(
            f"{mutation_column} column is missing in ACE2 raw input data"
        )
    if label_column not in result.columns:
        raise ValueError(f"{label_column} column is missing in ACE2 raw input data")

    if target_column not in result.columns:
        result[target_column] = default_target_name
    else:
        result[target_column] = (
            result[target_column].fillna(default_target_name).astype(str)
        )
    if variant_class_column not in result.columns:
        result[variant_class_column] = ""
    if "pass_pre_count_filter" in result.columns:
        result = result.loc[
            _normalize_bool_series(result["pass_pre_count_filter"])
        ].copy()
    if "pass_ACE2bind_expr_filter" in result.columns:
        result = result.loc[
            _normalize_bool_series(result["pass_ACE2bind_expr_filter"])
        ].copy()

    result[mutation_column] = (
        result[mutation_column].fillna("").astype(str).str.strip()
    )
    result = result.loc[
        ~result[mutation_column].str.contains(r"\*", regex=True, na=False)
    ].copy()

    result[output_name_column] = result[target_column].map(
        lambda value: _normalize_target_name(
            raw_target_name=str(value),
            target_name_aliases=target_name_aliases,
            reference_sequences=reference_sequences,
        )
    )
    result[output_mutation_column] = result[mutation_column]
    result[output_label_column] = pd.to_numeric(result[label_column], errors="coerce")

    standardized = result[
        [
            output_name_column,
            output_mutation_column,
            variant_class_column,
            output_label_column,
        ]
    ].rename(columns={variant_class_column: "variant_class"})
    standardized = standardized.loc[standardized[output_label_column].notna()].copy()

    mutation_text = (
        standardized[output_mutation_column].fillna("").astype(str).str.strip()
    )
    variant_class = (
        standardized["variant_class"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )
    synonymous_mask = variant_class.eq("synonymous")
    stop_mask = variant_class.eq("stop")
    deletion_mask = mutation_text.str.contains("-", regex=False)
    empty_mutation_mask = mutation_text.eq("")
    wt_mask = variant_class.eq("wildtype")
    standardized.loc[wt_mask, output_mutation_column] = "WT"
    keep_mask = wt_mask | ~(
        synonymous_mask | stop_mask | deletion_mask | empty_mutation_mask
    )
    standardized = standardized.loc[keep_mask].copy()
    return standardized.reset_index(drop=True)


@multiout_step(main="main", variants="variants", summary="summary")
def prepare_rbd_antibody_records(
    dataset: pd.DataFrame,
    reference_sequences: Dict[str, str],
    target_name_aliases: Dict[str, str],
    column_mapping: Dict[str, str],
    fallback_reference_sequence: Optional[str] = None,
    require_known_reference_sequence: bool = True,
    expected_reference_length: Optional[int] = 201,
    name_column: str = "name",
    target_column: str = "target",
    label_column: str = "label",
    mutation_column: str = "aa_substitutions",
    variant_class_column: str = "variant_class",
    default_target_name: str = "Wuhan-Hu-1",
    output_reference_column: str = "reference_id",
    output_group_column: str = "group_name",
    output_mutation_column: str = "mut_info",
    output_label_column: str = "label",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Standardize one antibody raw table into cleaner-ready records."""

    result = dataset.copy()
    available_column_mapping = _available_column_mapping(result, column_mapping)
    if available_column_mapping:
        result = result.rename(columns=available_column_mapping)

    missing_columns = sorted(
        set(RAW_ANTIBODY_REQUIRED_COLUMNS) - set(result.columns)
    )
    if missing_columns:
        raise ValueError(
            f"Raw antibody table is missing required columns {missing_columns}"
        )

    if target_column not in result.columns:
        result[target_column] = default_target_name
    else:
        result[target_column] = (
            result[target_column].fillna(default_target_name).astype(str)
        )
    if variant_class_column not in result.columns:
        result[variant_class_column] = ""

    unique_targets = sorted(
        result[target_column].dropna().astype(str).str.strip().unique().tolist()
    )
    if len(unique_targets) != 1:
        raise ValueError(
            f"Expected exactly one target per raw score table, got {unique_targets}"
        )

    raw_reference_id = unique_targets[0]
    reference_id = _normalize_target_name(
        raw_target_name=raw_reference_id,
        target_name_aliases=target_name_aliases,
        reference_sequences=reference_sequences,
        allow_unknown=not require_known_reference_sequence,
    )
    reference_sequence = _select_reference_sequence(
        reference_id=reference_id,
        raw_target_name=raw_reference_id,
        reference_sequences=reference_sequences,
        fallback_reference_sequence=fallback_reference_sequence,
        require_known_reference_sequence=require_known_reference_sequence,
        expected_reference_length=expected_reference_length,
    )

    result[mutation_column] = result[mutation_column].fillna("").astype(str).str.strip()
    result[label_column] = pd.to_numeric(result[label_column], errors="coerce")
    result = result.loc[result[label_column].notna()].copy()

    if "pass_pre_count_filter" in result.columns:
        result = result.loc[
            _normalize_bool_series(result["pass_pre_count_filter"])
        ].copy()
    if "pass_ACE2bind_expr_filter" in result.columns:
        result = result.loc[
            _normalize_bool_series(result["pass_ACE2bind_expr_filter"])
        ].copy()

    variant_class = (
        result[variant_class_column].fillna("").astype(str).str.strip().str.lower()
    )
    result["is_wildtype_variant"] = variant_class.eq("wildtype") | result[
        mutation_column
    ].eq("")
    result = result.loc[
        ~result[mutation_column].str.contains(r"\*", regex=True, na=False)
    ].copy()
    aggregation_spec: Dict[str, Tuple[str, str]] = {
        output_label_column: (label_column, "mean"),
        "replicate_count": (label_column, "count"),
        "is_wildtype_variant": ("is_wildtype_variant", "any"),
    }

    aggregated = (
        result.groupby(
            [name_column, target_column, mutation_column],
            as_index=False,
            sort=False,
        )
        .agg(**aggregation_spec)
        .rename(
            columns={
                name_column: "antibody_name",
                mutation_column: "aa_substitutions",
                output_label_column: output_label_column,
            }
        )
    )

    records: List[Dict[str, Any]] = []
    for row in aggregated.itertuples(index=False):
        aa_substitutions = str(row.aa_substitutions).strip()
        label = float(getattr(row, output_label_column))
        if pd.isna(label):
            continue

        if bool(row.is_wildtype_variant) or not aa_substitutions:
            normalized_mutation_text = "WT"
        else:
            normalized_mutations = _normalize_rbd_mutations(
                mutation_text=aa_substitutions,
                reference_sequence=reference_sequence,
            )
            if not normalized_mutations:
                continue
            normalized_mutation_text = ",".join(normalized_mutations)

        records.append(
            {
                output_reference_column: reference_id,
                output_group_column: f"{reference_id}__{row.antibody_name}",
                "sequence": reference_sequence,
                output_mutation_column: normalized_mutation_text,
                output_label_column: label,
                "antibody_name": str(row.antibody_name),
                "replicate_count": int(row.replicate_count),
            }
        )

    standardized_df = pd.DataFrame.from_records(
        records,
        columns=[
            output_reference_column,
            output_group_column,
            "sequence",
            output_mutation_column,
            output_label_column,
            "antibody_name",
            "replicate_count",
        ],
    )
    variant_df = standardized_df.loc[
        standardized_df[output_mutation_column].ne("WT")
    ].reset_index(drop=True)
    summary = pd.DataFrame.from_records(
        [
            {
                "reference_id": reference_id,
                "source_target_name": raw_reference_id,
                "reference_sequence": reference_sequence,
                "n_rows": int(len(variant_df)),
                "n_antibodies": (
                    int(standardized_df["antibody_name"].nunique())
                    if not standardized_df.empty
                    else 0
                ),
                "label_definition": (
                    "label is copied directly from the standardized antibody label "
                    "column and then WT-subtracted within each antibody group"
                ),
            }
        ]
    )
    return standardized_df.reset_index(drop=True), variant_df, summary


@pipeline_step
def add_reference_sequences_by_target(
    dataset: pd.DataFrame,
    reference_sequences: Dict[str, str],
    name_column: str = "name",
    sequence_column: str = "sequence",
) -> pd.DataFrame:
    """Attach reference sequences to standardized RBD rows."""

    result = dataset.copy()
    result[sequence_column] = result[name_column].map(reference_sequences)
    missing_targets = sorted(
        {
            str(target)
            for target in result.loc[
                result[sequence_column].isna(),
                name_column,
            ].dropna()
        }
    )
    if missing_targets:
        raise ValueError(
            "Missing reference sequences for targets: "
            + ", ".join(missing_targets)
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
    if not wt_rows.empty:
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


@multiout_step(main="main", variants="variants")
def capture_rbd_antibody_variants(
    dataset: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Keep WT-subtracted antibody rows in the pipeline and as a side artifact."""

    result = dataset.reset_index(drop=True)
    return result, result.copy()
