# Repository Tour

This page explains the role of the main publication-facing files and folders in the repository.

## Public Entry Points

- `README.md`: the front door for GitHub visitors and the shortest path into installation, quick start, and repo navigation
- `docs/`: the authoritative public documentation site
- `examples/`: small, copyable files for first runs and command examples
- `tutorials/`: interactive notebook tutorials for users who want a hands-on path

## Support And Maintenance Material

- `scripts/`: optional setup and data-preparation helpers
- `LICENSE`: repository license
- `THIRD_PARTY_SOFTWARE.md`: attribution for adapted or bundled third-party components

## Unsupported Development Material

- `src/cofolder/ui/`, `.streamlit/`, and `run_ui.sh` retain an obsolete Streamlit
  prototype for development reference only. They are not installed, tested, or
  supported in the initial release, and their configuration and command construction
  may not match the current recipe contracts.
- Use the supported `bias`, `validate`, `screen`, and `oracle` commands or the
  [Python API execution references](../reference/python-api/index.md) instead.

## Archived Material

- `legacy/`: historical files retained for traceability and migration support, not the recommended path for new users

## Recommended Reading Order

1. `README.md`
2. `getting-started/installation.md`
3. `getting-started/quickstart.md`
4. `user-guide/overview.md`
5. `tutorials/index.md`

## Recommended Hands-On Order

1. copy from `examples/`
2. run a first command from the quick start
3. continue with the written tutorials in `docs/tutorials/`
4. switch to the notebooks in `tutorials/` when you want interactive exploration
