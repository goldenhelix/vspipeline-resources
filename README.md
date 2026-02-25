# VSPipeline Resources

This repository provides tasks for creating and managing VarSeq projects. It serves two purposes:

1. **Standalone user tasks** — tasks at the root of this repository that users can run directly from the UI.
2. **Pipeline integration tasks** — tasks in the `pipeline/` subdirectory used as building blocks in secondary analysis workflows (e.g., PacBio WDL Somatic, Sentieon Secondary Analysis). These are not intended for independent use.

---

## Standalone User Tasks

### VSPipeline Create Project (`create_project_from_vcfs.task.yaml`)

Creates a VarSeq project from a folder of VCF.gz files using a project template.

**Parameters:**
- **Input Folder** — directory containing `*.vcf.gz` files to import
- **Project Directory** — where the project will be created (default: `AppData/Projects`)
- **Project Name** — name for the project; defaults to the input folder name if blank
- **VarSeq Project Template** — template file to use for project creation
- **Overwrite Existing VarSeq Projects** — whether to overwrite if the project already exists

---

### Update Sample Catalog from Manifest (`update_catalog_from_manifest.task.yaml`)

Populates or updates the SampleCatalog from a TSV or CSV manifest file. Uses fuzzy header matching to map manifest columns to catalog fields.

**Parameters:**
- **Manifest File** — a `.tsv`, `.csv`, or `.txt` file with sample metadata

---

### Update Sample Catalog with PanelApp AUS Genes (`update_catalog_with_panelapp_genes.task.yaml`)

Fetches a gene panel from the PanelApp Australia API and updates the SampleCatalog with the panel name and gene list for a given sample.

**Parameters:**
- **Sample Name** — the sample to update in the catalog
- **Panel Name** — the PanelApp panel name (e.g., `Achromatopsia`)

---

## Pipeline Integration Tasks (`pipeline/`)

These tasks are designed to be referenced as steps in larger workflow definitions. They are used as a git submodule in secondary analysis workflow repositories.

- **`pipeline/generate_vsbatch.task.yaml`** — generates VSBatch files from per-sample JSON metadata, querying the SampleCatalog for sample relationships
- **`pipeline/vspipeline.task.yaml`** — runs VSPipeline against a pre-generated `.vsbatch` file
- **`pipeline/run_varseq_reports.task.yaml`** — runs VSPGx and/or STR reports for a VarSeq project
- **`pipeline/vspipeline_utilities.py`** — Python utility library used by the pipeline tasks
