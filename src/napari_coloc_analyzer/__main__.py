"""Convenience launcher: ``python -m napari_coloc_analyzer``."""

import napari

viewer = napari.Viewer()
viewer.window.add_plugin_dock_widget(
    "napari-coloc-analyzer", "Colocalization Analyzer"
)
napari.run()
