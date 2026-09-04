# Screen Recipe

`Screen` accepts the screen-only analysis arguments `ifp_filter_threshold`,
`ifp_ligand_chain`, `cluster_ifps=False`, and
`ifp_cluster_similarity_threshold=0.5`. A single ligand chain is selected
automatically; pass `ifp_ligand_chain` for a multi-ligand system.

`Screen.run()` writes `screen_results.csv` and `screen_results_with_scores.csv`,
returns the latter as a `pandas.DataFrame`, and, when clustering is enabled, writes
`ifp_cluster_summary.csv`. Both consolidated outputs always include
`ifp_cluster_id` and `ifp_cluster_status`; disabled clustering uses `not_applied`.

See the [Screen guide](../../user-guide/screen.md) for the full stable score schema,
filter audit contract, clustering summary fields, and complete commands.

::: cofolder.recipes.screen
    options:
      show_root_heading: true
      show_source: true
      heading_level: 2
