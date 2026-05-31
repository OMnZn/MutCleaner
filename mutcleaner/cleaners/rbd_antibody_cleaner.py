from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    average_labels_by_name,
    convert_to_mutation_dataset_format,
    merge_columns,
    read_dataset,
    subtract_labels_by_wt,
    validate_mutations,
)
from .rbd_custom_cleaner import (
    add_reference_sequences_by_target,
    apply_mutations_preserving_wild_type,
    prepare_rbd_records,
)
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline, multiout_step
from ..core.sequence import ProteinSequence

if TYPE_CHECKING:
    from typing import Any, Dict, List, Optional, Tuple, Union

__all__ = [
    "RBDAntibodyConfig",
    "create_rbd_antibody_cleaner",
    "clean_rbd_antibody_dataset",
]

logger = logging.getLogger(__name__)

DEFAULT_RBD_REFERENCE_SEQUENCES = {
    "Wuhan-Hu-1": "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST",
}

DEFAULT_RBD_TARGET_NAME_ALIASES = {
    "SARS-CoV-2": "Wuhan-Hu-1",
    "SARS_CoV_2": "Wuhan-Hu-1",
    "Wuhan_Hu_1": "Wuhan-Hu-1",
}


@multiout_step(main="main", variants="variants")
def capture_rbd_antibody_variants(
    dataset: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Keep WT-subtracted antibody rows in the pipeline and as a side artifact."""

    result = dataset.reset_index(drop=True)
    return result, result.copy()


def __dir__() -> List[str]:
    """Return exported names."""

    return __all__


@dataclass
class RBDAntibodyConfig(BaseCleanerConfig):
    """Configuration for the RBD antibody cleaner."""

    reference_sequences: Dict[str, str] = field(
        default_factory=lambda: deepcopy(DEFAULT_RBD_REFERENCE_SEQUENCES)
    )
    target_name_aliases: Dict[str, str] = field(
        default_factory=lambda: deepcopy(DEFAULT_RBD_TARGET_NAME_ALIASES)
    )
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "name": "name",
            "target": "target",
            "score": "label",
            "aa_substitutions": "aa_substitutions",
            "variant_class": "variant_class",
            "n_aa_substitutions": "n_aa_substitutions",
            "pass_pre_count_filter": "pass_pre_count_filter",
            "pass_ACE2bind_expr_filter": "pass_ACE2bind_expr_filter",
        }
    )
    fallback_reference_sequence: Optional[str] = None
    require_known_reference_sequence: bool = True
    validate_mut_workers: int = 16
    process_workers: int = 16
    label_columns: List[str] = field(default_factory=lambda: ["label"])
    primary_label_column: str = "label"
    pipeline_name: str = "RBDAntibody"

    @classmethod
    def from_dict(
        cls,
        config_dict: Dict[str, Any],
    ) -> "RBDAntibodyConfig":
        config_fields = {one_field.name for one_field in fields(cls)}
        payload = {
            key: value for key, value in config_dict.items() if key in config_fields
        }
        config = cls(**{**payload, "validate_config": False})
        config.validate()
        return config

    def validate(self) -> None:
        super().validate()

        if self.require_known_reference_sequence and not self.reference_sequences:
            raise ValueError(
                "reference_sequences cannot be empty when "
                "require_known_reference_sequence=True"
            )
        if (
            not self.require_known_reference_sequence
            and not self.reference_sequences
            and self.fallback_reference_sequence is None
        ):
            raise ValueError(
                "Provide reference_sequences or fallback_reference_sequence"
            )
        if not self.label_columns:
            raise ValueError("label_columns cannot be empty")
        if self.primary_label_column not in self.label_columns:
            raise ValueError(
                f"primary_label_column '{self.primary_label_column}' "
                f"must be in label_columns {self.label_columns}"
            )

        required_standard_columns = {
            "name",
            "target",
            "label",
            "aa_substitutions",
        }
        missing_standard_columns = required_standard_columns - set(
            self.column_mapping.values()
        )
        if missing_standard_columns:
            raise ValueError(
                "column_mapping must provide standardized columns "
                f"{sorted(required_standard_columns)}, missing "
                f"{sorted(missing_standard_columns)}"
            )

        for target_name, sequence in self.reference_sequences.items():
            if len(str(sequence).strip()) <= 0:
                raise ValueError(
                    f"Reference sequence for target '{target_name}' cannot be empty"
                )

        if (
            self.fallback_reference_sequence is not None
            and len(str(self.fallback_reference_sequence).strip()) <= 0
        ):
            raise ValueError(
                "fallback_reference_sequence cannot be empty when provided"
            )


def _resolve_rbd_antibody_config(
    config: Optional[Union[RBDAntibodyConfig, Dict[str, Any], str, Path]] = None,
) -> RBDAntibodyConfig:
    default_config = RBDAntibodyConfig()
    if config is None:
        final_config = default_config
    elif isinstance(config, RBDAntibodyConfig):
        final_config = config
    elif isinstance(config, dict):
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        final_config = RBDAntibodyConfig.from_json(config)
    else:
        raise TypeError(
            "config must be RBDAntibodyConfig, dict, str, Path or None, "
            f"got {type(config)}"
        )
    final_config.validate()
    return final_config


def create_rbd_antibody_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[RBDAntibodyConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create the RBD antibody cleaning pipeline."""

    final_config = _resolve_rbd_antibody_config(config)

    logger.info(
        "RBD antibody dataset will be cleaned with pipeline: %s",
        final_config.pipeline_name,
    )
    logger.debug("Configuration:\n%s", final_config.get_summary())

    pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)
    pipeline = (
        pipeline.delayed_then(
            prepare_rbd_records,
            mode="antibody",
            reference_sequences=final_config.reference_sequences,
            target_name_aliases=final_config.target_name_aliases,
            column_mapping=final_config.column_mapping,
            fallback_reference_sequence=final_config.fallback_reference_sequence,
            require_known_reference_sequence=(
                final_config.require_known_reference_sequence
            ),
        )
        .delayed_then(
            validate_mutations,
            mutation_column="mut_info",
            format_mutations=True,
            mutation_sep=" ",
            is_zero_based=False,
            exclude_patterns="WT",
            cache_results=False,
            num_workers=final_config.validate_mut_workers,
        )
        .delayed_then(
            average_labels_by_name,
            name_columns=["reference_id", "antibody_name", "mut_info"],
            label_columns=final_config.label_columns,
        )
        .delayed_then(
            add_reference_sequences_by_target,
            reference_sequences=final_config.reference_sequences,
            name_column="reference_id",
            sequence_column="sequence",
            fallback_reference_sequence=final_config.fallback_reference_sequence,
        )
        # `group_name` is only an internal WT-subtraction key. For the current
        # antibody tables one file is usually one reference, so `antibody_name`
        # would often be enough on its own, but we keep the combined key to
        # preserve the existing downstream grouping behavior.
        .delayed_then(
            merge_columns,
            columns_to_merge=["reference_id", "antibody_name"],
            new_column_name="group_name",
            separator="_",
        )
        .delayed_then(
            apply_mutations_preserving_wild_type,
            sequence_column="sequence",
            name_column="group_name",
            mutation_column="mut_info",
            wt_identifier="WT",
            mutation_sep=",",
            is_zero_based=True,
            sequence_type="protein",
            num_workers=final_config.process_workers,
        )
        .delayed_then(
            subtract_labels_by_wt,
            name_column="group_name",
            label_columns=final_config.label_columns,
            mutation_column="mut_info",
            wt_identifier="WT",
            in_place=True,
            drop_wt_row=True,
        )
        .delayed_then(capture_rbd_antibody_variants)
        .delayed_then(
            convert_to_mutation_dataset_format,
            name_column="antibody_name",
            mutation_column="mut_info",
            sequence_column="sequence",
            mutated_sequence_column="mut_seq",
            label_column=final_config.primary_label_column,
            include_wild_type=False,
            is_zero_based=True,
        )
    )

    if isinstance(dataset_or_path, (str, Path)):
        pipeline.add_delayed_step(read_dataset, 0, file_format="csv")
    elif dataset_or_path is not None and not isinstance(dataset_or_path, pd.DataFrame):
        raise TypeError(
            "dataset_or_path must be pd.DataFrame, str, Path, or None, "
            f"got {type(dataset_or_path)}"
        )

    return pipeline


def clean_rbd_antibody_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Execute the RBD antibody cleaning pipeline."""

    pipeline.execute()
    dataset_df, reference_sequences = pipeline.data
    summary_df = pipeline.get_artifact("prepare_rbd_records.summary")
    summary = summary_df.iloc[0].to_dict()
    source_reference_id = str(summary["reference_id"])

    if isinstance(dataset_df, pd.DataFrame) and dataset_df.empty:
        dataset = MutationDataset(name=pipeline.name)
        reference_sequence = str(summary["reference_sequence"])
        dataset.add_reference_sequence(
            source_reference_id,
            ProteinSequence(reference_sequence, name=source_reference_id),
        )
        logger.info(
            "Successfully cleaned RBD antibody dataset: 0 mutation sets from 1 reference"
        )
        return pipeline, dataset

    dataset = MutationDataset.from_dataframe(dataset_df, reference_sequences)
    for one_reference_sequence in dataset.reference_sequences.values():
        one_reference_sequence.name = source_reference_id
    logger.info(
        "Successfully cleaned RBD antibody dataset: %s mutation sets from %s references",
        len(dataset.mutation_sets),
        len(dataset.reference_sequences),
    )
    return pipeline, dataset
