#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwschema
~~~~~~~~~~~~~~

This module provides schema definitions for different sections of the
``mwTab`` Metabolomics Workbench format.
"""

import sys

from schema import Schema, Optional, Or

if sys.version_info.major == 2:
    str = unicode


metabolomics_workbench_schema = Schema(
    {
        "VERSION": str,
        "CREATED_ON": str,
        Optional("STUDY_ID"): str,
        Optional("ANALYSIS_ID"): str,
        Optional("PROJECT_ID"): str,
        Optional("HEADER"): str,
        Optional("DATATRACK_ID"): str
    }
)

project_schema = Schema(
    {
        "PROJECT_TITLE": str,
        Optional("PROJECT_TYPE"): str,
        "PROJECT_SUMMARY": str,
        "INSTITUTE": str,
        Optional("DEPARTMENT"): str,
        Optional("LABORATORY"): str,
        "LAST_NAME": str,
        "FIRST_NAME": str,
        "ADDRESS": str,
        "EMAIL": str,
        "PHONE": str,
        Optional("FUNDING_SOURCE"): str,
        Optional("PROJECT_COMMENTS"): str,
        Optional("PUBLICATIONS"): str,
        Optional("CONTRIBUTORS"): str,
        Optional("DOI"): str
    }
)

study_schema = Schema(
    {
        "STUDY_TITLE": str,
        Optional("STUDY_TYPE"): str,
        "STUDY_SUMMARY": str,
        "INSTITUTE": str,
        Optional("DEPARTMENT"): str,
        Optional("LABORATORY"): str,
        "LAST_NAME": str,
        "FIRST_NAME": str,
        "ADDRESS": str,
        "EMAIL": str,
        "PHONE": str,
        Optional("NUM_GROUPS"): str,
        Optional("TOTAL_SUBJECTS"): str,
        Optional("NUM_MALES"): str,
        Optional("NUM_FEMALES"): str,
        Optional("STUDY_COMMENTS"): str,
        Optional("PUBLICATIONS"): str,  # assumed
        Optional("SUBMIT_DATE"): str  # assumed
    }
)

subject_schema = Schema(
    {
        "SUBJECT_TYPE": str,
        "SUBJECT_SPECIES": str,
        Optional("TAXONOMY_ID"): str,
        Optional("GENOTYPE_STRAIN"): str,
        Optional("AGE_OR_AGE_RANGE"): str,
        Optional("WEIGHT_OR_WEIGHT_RANGE"): str,
        Optional("HEIGHT_OR_HEIGHT_RANGE"): str,
        Optional("GENDER"): str,
        Optional("HUMAN_RACE"): str,
        Optional("HUMAN_ETHNICITY"): str,
        Optional("HUMAN_TRIAL_TYPE"): str,
        Optional("HUMAN_LIFESTYLE_FACTORS"): str,
        Optional("HUMAN_MEDICATIONS"): str,
        Optional("HUMAN_PRESCRIPTION_OTC"): str,
        Optional("HUMAN_SMOKING_STATUS"): str,
        Optional("HUMAN_ALCOHOL_DRUG_USE"): str,
        Optional("HUMAN_NUTRITION"): str,
        Optional("HUMAN_INCLUSION_CRITERIA"): str,
        Optional("HUMAN_EXCLUSION_CRITERIA"): str,
        Optional("ANIMAL_ANIMAL_SUPPLIER"): str,
        Optional("ANIMAL_HOUSING"): str,
        Optional("ANIMAL_LIGHT_CYCLE"): str,
        Optional("ANIMAL_FEED"): str,
        Optional("ANIMAL_WATER"): str,
        Optional("ANIMAL_INCLUSION_CRITERIA"): str,
        Optional("CELL_BIOSOURCE_OR_SUPPLIER"): str,
        Optional("CELL_STRAIN_DETAILS"): str,
        Optional("SUBJECT_COMMENTS"): str,
        Optional("CELL_PRIMARY_IMMORTALIZED"): str,
        Optional("CELL_PASSAGE_NUMBER"): str,
        Optional("CELL_COUNTS"): str,
        Optional("SPECIES_GROUP"): str
    }
)

subject_sample_factors_schema = Schema(
    [
        {
            "Subject ID": str,
            "Sample ID": str,
            "Factors": dict,
            Optional("Additional sample data"): {
                Optional("RAW_FILE_NAME"): str,
                Optional(str): str
            }
        }
    ]
)

collection_schema = Schema(
    {
        "COLLECTION_SUMMARY": str,
        Optional("COLLECTION_PROTOCOL_ID"): str,
        Optional("COLLECTION_PROTOCOL_FILENAME"): str,
        Optional("COLLECTION_PROTOCOL_COMMENTS"): str,
        Optional("SAMPLE_TYPE"): str,  # assumed optional due to large number of files without
        Optional("COLLECTION_METHOD"): str,
        Optional("COLLECTION_LOCATION"): str,
        Optional("COLLECTION_FREQUENCY"): str,
        Optional("COLLECTION_DURATION"): str,
        Optional("COLLECTION_TIME"): str,
        Optional("VOLUMEORAMOUNT_COLLECTED"): str,
        Optional("STORAGE_CONDITIONS"): str,
        Optional("COLLECTION_VIALS"): str,
        Optional("STORAGE_VIALS"): str,
        Optional("COLLECTION_TUBE_TEMP"): str,
        Optional("ADDITIVES"): str,
        Optional("BLOOD_SERUM_OR_PLASMA"): str,
        Optional("TISSUE_CELL_IDENTIFICATION"): str,
        Optional("TISSUE_CELL_QUANTITY_TAKEN"): str
    }
)

treatment_schema = Schema(
    {
        "TREATMENT_SUMMARY": str,
        Optional("TREATMENT_PROTOCOL_ID"): str,
        Optional("TREATMENT_PROTOCOL_FILENAME"): str,
        Optional("TREATMENT_PROTOCOL_COMMENTS"): str,
        Optional("TREATMENT"): str,
        Optional("TREATMENT_COMPOUND"): str,
        Optional("TREATMENT_ROUTE"): str,
        Optional("TREATMENT_DOSE"): str,
        Optional("TREATMENT_DOSEVOLUME"): str,
        Optional("TREATMENT_DOSEDURATION"): str,
        Optional("TREATMENT_VEHICLE"): str,
        Optional("ANIMAL_VET_TREATMENTS"): str,
        Optional("ANIMAL_ANESTHESIA"): str,
        Optional("ANIMAL_ACCLIMATION_DURATION"): str,
        Optional("ANIMAL_FASTING"): str,
        Optional("ANIMAL_ENDP_EUTHANASIA"): str,
        Optional("ANIMAL_ENDP_TISSUE_COLL_LIST"): str,
        Optional("ANIMAL_ENDP_TISSUE_PROC_METHOD"): str,
        Optional("ANIMAL_ENDP_CLINICAL_SIGNS"): str,
        Optional("HUMAN_FASTING"): str,
        Optional("HUMAN_ENDP_CLINICAL_SIGNS"): str,
        Optional("CELL_STORAGE"): str,
        Optional("CELL_GROWTH_CONTAINER"): str,
        Optional("CELL_GROWTH_CONFIG"): str,
        Optional("CELL_GROWTH_RATE"): str,
        Optional("CELL_INOC_PROC"): str,
        Optional("CELL_MEDIA"): str,
        Optional("CELL_ENVIR_COND"): str,
        Optional("CELL_HARVESTING"): str,
        Optional("PLANT_GROWTH_SUPPORT"): str,
        Optional("PLANT_GROWTH_LOCATION"): str,
        Optional("PLANT_PLOT_DESIGN"): str,
        Optional("PLANT_LIGHT_PERIOD"): str,
        Optional("PLANT_HUMIDITY"): str,
        Optional("PLANT_TEMP"): str,
        Optional("PLANT_WATERING_REGIME"): str,
        Optional("PLANT_NUTRITIONAL_REGIME"): str,
        Optional("PLANT_ESTAB_DATE"): str,
        Optional("PLANT_HARVEST_DATE"): str,
        Optional("PLANT_GROWTH_STAGE"): str,
        Optional("PLANT_METAB_QUENCH_METHOD"): str,
        Optional("PLANT_HARVEST_METHOD"): str,
        Optional("PLANT_STORAGE"): str,
        Optional("CELL_PCT_CONFLUENCE"): str,
        Optional("CELL_MEDIA_LASTCHANGED"): str
    }
)

sampleprep_schema = Schema(
    {
        "SAMPLEPREP_SUMMARY": str,
        Optional("SAMPLEPREP_PROTOCOL_ID"): str,
        Optional("SAMPLEPREP_PROTOCOL_FILENAME"): str,
        Optional("SAMPLEPREP_PROTOCOL_COMMENTS"): str,
        Optional("PROCESSING_METHOD"): str,
        Optional("PROCESSING_STORAGE_CONDITIONS"): str,
        Optional("EXTRACTION_METHOD"): str,
        Optional("EXTRACT_CONCENTRATION_DILUTION"): str,
        Optional("EXTRACT_ENRICHMENT"): str,
        Optional("EXTRACT_CLEANUP"): str,
        Optional("EXTRACT_STORAGE"): str,
        Optional("SAMPLE_RESUSPENSION"): str,
        Optional("SAMPLE_DERIVATIZATION"): str,
        Optional("SAMPLE_SPIKING"): str,
        Optional("ORGAN"): str,
        Optional("ORGAN_SPECIFICATION"): str,
        Optional("CELL_TYPE"): str,
        Optional("SUBCELLULAR_LOCATION"): str
    }
)

chromatography_schema = Schema(
    {
        Optional("CHROMATOGRAPHY_SUMMARY"): str,
        "CHROMATOGRAPHY_TYPE": str,
        "INSTRUMENT_NAME": str,
        "COLUMN_NAME": str,
        Optional("FLOW_GRADIENT"): str,
        Optional("FLOW_RATE"): str,
        Optional("COLUMN_TEMPERATURE"): str,
        Optional("METHODS_FILENAME"): str,
        Optional("SOLVENT_A"): str,
        Optional("SOLVENT_B"): str,
        Optional("METHODS_ID"): str,
        Optional("COLUMN_PRESSURE"): str,
        Optional("INJECTION_TEMPERATURE"): str,
        Optional("INTERNAL_STANDARD"): str,
        Optional("INTERNAL_STANDARD_MT"): str,
        Optional("RETENTION_INDEX"): str,
        Optional("RETENTION_TIME"): str,
        Optional("SAMPLE_INJECTION"): str,
        Optional("SAMPLING_CONE"): str,
        Optional("ANALYTICAL_TIME"): str,
        Optional("CAPILLARY_VOLTAGE"): str,
        Optional("MIGRATION_TIME"): str,
        Optional("OVEN_TEMPERATURE"): str,
        Optional("PRECONDITIONING"): str,
        Optional("RUNNING_BUFFER"): str,
        Optional("RUNNING_VOLTAGE"): str,
        Optional("SHEATH_LIQUID"): str,
        Optional("TIME_PROGRAM"): str,
        Optional("TRANSFERLINE_TEMPERATURE"): str,
        Optional("WASHING_BUFFER"): str,
        Optional("WEAK_WASH_SOLVENT_NAME"): str,
        Optional("WEAK_WASH_VOLUME"): str,
        Optional("STRONG_WASH_SOLVENT_NAME"): str,
        Optional("STRONG_WASH_VOLUME"): str,
        Optional("TARGET_SAMPLE_TEMPERATURE"): str,
        Optional("SAMPLE_LOOP_SIZE"): str,
        Optional("SAMPLE_SYRINGE_SIZE"): str,
        Optional("RANDOMIZATION_ORDER"): str,
        Optional("CHROMATOGRAPHY_COMMENTS"): str
    }
)

analysis_schema = Schema(
    {
        "ANALYSIS_TYPE": str,
        Optional("LABORATORY_NAME"): str,
        Optional("OPERATOR_NAME"): str,
        Optional("DETECTOR_TYPE"): str,
        Optional("SOFTWARE_VERSION"): str,
        Optional("ACQUISITION_DATE"): str,
        Optional("ANALYSIS_PROTOCOL_FILE"): str,
        Optional("ACQUISITION_PARAMETERS_FILE"): str,
        Optional("PROCESSING_PARAMETERS_FILE"): str,
        Optional("DATA_FORMAT"): str,

        # not specified in mwTab specification (assumed)
        Optional("ACQUISITION_ID"): str,
        Optional("ACQUISITION_TIME"): str,
        Optional("ANALYSIS_COMMENTS"): str,
        Optional("ANALYSIS_DISPLAY"): str,
        Optional("INSTRUMENT_NAME"): str,
        Optional("INSTRUMENT_PARAMETERS_FILE"): str,
        Optional("NUM_FACTORS"): str,
        Optional("NUM_METABOLITES"): str,
        Optional("PROCESSED_FILE"): str,
        Optional("RANDOMIZATION_ORDER"): str,
        Optional("RAW_FILE"): str,
    }
)

ms_schema = Schema(
    {
        "INSTRUMENT_NAME": str,
        "INSTRUMENT_TYPE": str,
        "MS_TYPE": str,
        "ION_MODE": str,
        "MS_COMMENTS": str,  # changed to required
        Optional("CAPILLARY_TEMPERATURE"): str,
        Optional("CAPILLARY_VOLTAGE"): str,
        Optional("COLLISION_ENERGY"): str,
        Optional("COLLISION_GAS"): str,
        Optional("DRY_GAS_FLOW"): str,
        Optional("DRY_GAS_TEMP"): str,
        Optional("FRAGMENT_VOLTAGE"): str,
        Optional("FRAGMENTATION_METHOD"): str,
        Optional("GAS_PRESSURE"): str,
        Optional("HELIUM_FLOW"): str,
        Optional("ION_SOURCE_TEMPERATURE"): str,
        Optional("ION_SPRAY_VOLTAGE"): str,
        Optional("IONIZATION"): str,
        Optional("IONIZATION_ENERGY"): str,
        Optional("IONIZATION_POTENTIAL"): str,
        Optional("MASS_ACCURACY"): str,
        Optional("PRECURSOR_TYPE"): str,
        Optional("REAGENT_GAS"): str,
        Optional("SOURCE_TEMPERATURE"): str,
        Optional("SPRAY_VOLTAGE"): str,
        Optional("ACTIVATION_PARAMETER"): str,
        Optional("ACTIVATION_TIME"): str,
        Optional("ATOM_GUN_CURRENT"): str,
        Optional("AUTOMATIC_GAIN_CONTROL"): str,
        Optional("BOMBARDMENT"): str,
        Optional("CDL_SIDE_OCTOPOLES_BIAS_VOLTAGE"): str,
        Optional("CDL_TEMPERATURE"): str,
        Optional("DATAFORMAT"): str,
        Optional("DESOLVATION_GAS_FLOW"): str,
        Optional("DESOLVATION_TEMPERATURE"): str,
        Optional("INTERFACE_VOLTAGE"): str,
        Optional("IT_SIDE_OCTOPOLES_BIAS_VOLTAGE"): str,
        Optional("LASER"): str,
        Optional("MATRIX"): str,
        Optional("NEBULIZER"): str,
        Optional("OCTPOLE_VOLTAGE"): str,
        Optional("PROBE_TIP"): str,
        Optional("RESOLUTION_SETTING"): str,
        Optional("SAMPLE_DRIPPING"): str,
        Optional("SCAN_RANGE_MOVERZ"): str,
        Optional("SCANNING"): str,
        Optional("SCANNING_CYCLE"): str,
        Optional("SCANNING_RANGE"): str,
        Optional("SKIMMER_VOLTAGE"): str,
        Optional("TUBE_LENS_VOLTAGE"): str,
        Optional("MS_RESULTS_FILE"): Or(str, dict)
    }
)

nmr_schema = Schema(
    {
        "INSTRUMENT_NAME": str,
        "INSTRUMENT_TYPE": str,
        "NMR_EXPERIMENT_TYPE": str,
        Optional("NMR_COMMENTS"): str,
        Optional("FIELD_FREQUENCY_LOCK"): str,
        Optional("STANDARD_CONCENTRATION"): str,
        "SPECTROMETER_FREQUENCY": str,
        Optional("NMR_PROBE"): str,
        Optional("NMR_SOLVENT"): str,
        Optional("NMR_TUBE_SIZE"): str,
        Optional("SHIMMING_METHOD"): str,
        Optional("PULSE_SEQUENCE"): str,
        Optional("WATER_SUPPRESSION"): str,
        Optional("PULSE_WIDTH"): str,
        Optional("POWER_LEVEL"): str,
        Optional("RECEIVER_GAIN"): str,
        Optional("OFFSET_FREQUENCY"): str,
        Optional("PRESATURATION_POWER_LEVEL"): str,
        Optional("CHEMICAL_SHIFT_REF_CPD"): str,
        Optional("TEMPERATURE"): str,
        Optional("NUMBER_OF_SCANS"): str,
        Optional("DUMMY_SCANS"): str,
        Optional("ACQUISITION_TIME"): str,
        Optional("RELAXATION_DELAY"): str,
        Optional("SPECTRAL_WIDTH"): str,
        Optional("NUM_DATA_POINTS_ACQUIRED"): str,
        Optional("REAL_DATA_POINTS"): str,
        Optional("LINE_BROADENING"): str,
        Optional("ZERO_FILLING"): str,
        Optional("APODIZATION"): str,
        Optional("BASELINE_CORRECTION_METHOD"): str,
        Optional("CHEMICAL_SHIFT_REF_STD"): str,
        Optional("BINNED_INCREMENT"): str,
        Optional("BINNED_DATA_NORMALIZATION_METHOD"): str,
        Optional("BINNED_DATA_PROTOCOL_FILE"): str,
        Optional("BINNED_DATA_CHEMICAL_SHIFT_RANGE"): str,
        Optional("BINNED_DATA_EXCLUDED_RANGE"): str
    }
)

data_schema = Schema(
    [
        {
            Or("Metabolite", "Bin range(ppm)", only_one=True): str,
            Optional(str): str,
        },
    ]
)

extended_schema = Schema(
    [
        {
            "Metabolite": str,
            Optional(str): str,
            "sample_id": str
        },
    ]
)

ms_metabolite_data_schema = Schema(
    {
        "Units": str,
        "Data": data_schema,
        "Metabolites": data_schema,
        Optional("Extended"): extended_schema
    }
)

nmr_binned_data_schema = Schema(
    {
        "Units": str,
        "Data": data_schema
    }
)

section_schema_mapping = {
    "METABOLOMICS WORKBENCH": metabolomics_workbench_schema,
    "PROJECT": project_schema,
    "STUDY": study_schema,
    "ANALYSIS": analysis_schema,
    "SUBJECT": subject_schema,
    "SUBJECT_SAMPLE_FACTORS": subject_sample_factors_schema,
    "COLLECTION": collection_schema,
    "TREATMENT": treatment_schema,
    "SAMPLEPREP": sampleprep_schema,
    "CHROMATOGRAPHY": chromatography_schema,
    "MS": ms_schema,
    "NM": nmr_schema,
    "MS_METABOLITE_DATA": ms_metabolite_data_schema,
    "NMR_METABOLITE_DATA": ms_metabolite_data_schema,
    "NMR_BINNED_DATA": nmr_binned_data_schema,
}
