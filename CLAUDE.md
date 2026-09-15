# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`plant-food-research-open/assemblyqc` is an nf-core-style Nextflow (DSL2) pipeline that runs a large set of
independent QC tools over one or more genome assemblies (optionally with GFF3 annotations, HiC/long/short reads,
and parental data) and renders every tool's results into a single unified `report.html`. Almost every tool is
optional and gated behind a `<tool>_skip` parameter.

## Common commands

```bash
# Run the built-in minimal test profile (needs Nextflow + Docker/Singularity)
nextflow run main.nf -profile docker,test --outdir results

# Run pipeline-level nf-test suite (tests/ dir), matching CI
nf-test test --tag test --profile +docker --verbose

# Run a single pipeline-level test
nf-test test tests/orthofinder/main.nf.test --profile +docker --verbose

# Re-record snapshots after an intentional output change
nf-test test <path> --profile +docker --verbose --update-snapshots

# Lint everything (prettier, mypy/ruff on bin/*.py, version consistency, Nextflow lint) before committing
pre-commit run --all-files

# nf-core pipeline lint (schema/module/template conformance)
nf-core pipelines lint .

# Regenerate nextflow_schema.json after adding/changing a param in nextflow.config
nf-core pipelines schema build
```

Module-level nf-test suites under `modules/nf-core/**/tests`, `modules/gallvp/**/tests`,
`subworkflows/nf-core/**/tests` and `subworkflows/gallvp/**/tests` are excluded from the pipeline-level run
(see `ignore` in `nf-test.config`) — they belong to the upstream module repos, not this pipeline.

`nextflow config -o json .` currently fails locally with `No such variable: meta` (reproduces on a clean
checkout, unrelated to any in-progress change). This breaks anything that shells out to it, e.g.
`nf-core pipelines bump-version`'s regeneration of `ro-crate-metadata.json` — that file has to be regenerated
in an environment where the command succeeds.

## Architecture

### Entry point and per-run flow

`main.nf` wraps everything: `PIPELINE_INITIALISATION` (parses/validates `--input` against
`assets/schema_input.json`, resolves reads/xref-assembly channels) → the `ASSEMBLYQC` workflow in
`workflows/assemblyqc.nf` (the actual pipeline logic, ~1100 lines, one section per tool) →
`PIPELINE_COMPLETION` (email/summary). `workflows/assemblyqc.nf` is the file to read to understand how any two
tools' channels connect.

### Module/subworkflow provenance — three distinct trees

- `modules/nf-core/`, `subworkflows/nf-core/`: vendored verbatim from nf-core/modules via `nf-core modules install`.
- `modules/gallvp/`, `subworkflows/gallvp/`: vendored the same way but from the community repo
  `GallVp/nxf-components` (a second repo tracked in `modules.json`). Used for tools that don't have an
  nf-core module (e.g. `plotsr`, `ltrretriever`, `syri`). Install/update with:
  `nf-core modules --git-remote https://github.com/GallVp/nxf-components.git install <module>`.
- `modules/local/`, `subworkflows/local/`: bespoke to this pipeline, usually a thin process wrapping a script in
  `bin/` (Perl for the legacy `assemblathon_stats`/Circos-plotting tools, Python for everything newer). This is
  also where the report is built (`modules/local/createreport`).

Do not hand-edit files under the `nf-core`/`gallvp` trees — patch via `nf-core modules patch` if a diff is
unavoidable (see `*.diff` files in a few module dirs), otherwise changes will be silently lost on the next
`nf-core modules update`.

Per-process CLI args, `publishDir`, and resource overrides live in `conf/modules.config`, keyed by
`withName: '.*:ASSEMBLYQC:<PROCESS>'` (or a subworkflow-qualified pattern for a process reused in several
places, e.g. `.*:FASTA_SYNTENY:MINIMAP2_ALIGN`).

### Adding a pipeline step / parameter

Follow `docs/CONTRIBUTING.md`'s checklist. The two easy-to-miss bits:

- Channel naming: `ch_output_from_<process>` for a process's own output, `ch_<previousprocess>_for_<nextprocess>`
  for something threaded between two steps.
- New/changed params must be added to **both** `nextflow.config`'s `params { }` block (with a default) and
  `nextflow_schema.json` (regenerate with `nf-core pipelines schema build`, not by hand) — keep them in sync:
  `nf-core pipelines lint` flags a param missing from the schema, and nf-schema rejects unrecognised params
  at runtime.
- Version strings: most modules emit a Nextflow "topic" channel (`emit: versions_x, topic: versions`) collected
  automatically near the end of `assemblyqc.nf` via `channel.topic("versions")`; a few older/local modules still
  use the classic `path("versions.yml")` + explicit `ch_versions = ch_versions.mix(X.out.versions)` pattern —
  match whichever style the module you're touching already uses.

### The report

`modules/local/createreport` stages every enabled tool's published outputs into its own subdirectory (e.g.
`busco_outputs/`, `tidk_outputs/`) and runs `bin/assemblyqc.py`, which:

1. Calls one `parse_<tool>_folder()` from `bin/report_modules/parsers/<tool>_parser.py` per tool — each reads
   the raw files in that tool's staged subdirectory and returns a `{TOOL_KEY: ...}` dict (a flat table, or a
   list of per-genome dicts for tools that report per-assembly).
2. Merges all of those dicts into one `all_stats_dicts` and hands it to `ReportPrinter`
   (`bin/report_modules/report_printer.py`), which renders `bin/report_modules/templates/base.html` with Jinja2.
3. `base.html` conditionally includes a `{% include 'tool/tool.html' %}` block per tool, guarded by
   `{% if 'TOOL_KEY' in all_stats_dicts %}` — both the top nav button and the tab body have to be added there
   for a new tool to show up.

Tools that report per-genome detail (BUSCO, TIDK, ...) follow the same dropdown pattern:
`templates/<tool>/dropdown.html` (a `<select>` switching visible `.tabcontent-<TOOL>` divs via the shared
`showContent()` JS in `templates/js.html`), plus `summary_contents.html` (cross-genome table) and
`report_contents.html` (per-genome detail, one div per assembly). Copy an existing tool's three files as the
template for a new one rather than inventing a new layout.

### Test layout

`tests/<name>/main.nf.test` are full-pipeline nf-tests (nf-test's `nextflow_pipeline` blocks), each asserting
against a checked-in `.snap` snapshot of `report.json` plus stable output paths. Several of these snapshots
embed the pipeline's own version string (`Workflow.workflow.manifest.version`) — see the "Pipeline version
bumping" note in `docs/CONTRIBUTING.md` for why a version bump must update them in the same commit.
