# Changelog

All notable changes to `protools4py` are recorded here. Format:
[Keep a Changelog](https://keepachangelog.com/) style, sections per
`git tag v*` release (newest first). Unreleased changes stay under
`## Development` until the next tag is cut.

## Development

Unreleased changes (will move under the next `v*` tag when cut).

### Added
- `skills/` — agent skills per module group (sequence, structure, cluster,
  dataset, database, dock), installable via `npx skills add GCS-ZHN/protools`.
- `docs/architecture/overview.md` — single source of truth for project
  layout, module boundaries, and key entrypoints.
- `AGENTS.md` — agent-facing notes: workflow, hard rules, build /
  install / test commands, conventions to read & update the
  architecture doc instead of restating it.

### Changed
- `README.md` — slimmed to install + module index; defers layout
  details to `docs/architecture/overview.md`.
- `README.md` — added skills installation instructions
  (`npx skills add GCS-ZHN/protools`).
