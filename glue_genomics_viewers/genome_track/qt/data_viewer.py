import logging
import os

from echo import delay_callback
from matplotlib.axes._base import _TransformedBoundsLocator
import numpy as np
import pandas as pd
from glue.core import Subset
from glue.utils import nonpartial, defer_draw, decorate_all_methods
from glue.viewers.matplotlib.qt.data_viewer import MatplotlibDataViewer
import coolbox.api as cb
from coolbox.utilities.bed import ReadBed

from qtpy import QtWidgets
from ...data import BedgraphData, BedGraph, GenomeRange
from ...subsets import GenomicRangeSubsetState

from .layer_style_editor import GenomeTrackLayerStyleEditor
from .options_widget import GenomeTrackOptionsWidget
from ..layer_artist import GenomeProfileLayerArtist, GenomeLoopLayerArtist, GenomeTrackLayerArtist
from ..state import GenomeTrackState
from ..utils import PanTrackerMixin


__all__ = ['GenomeTrackViewer']


@decorate_all_methods(defer_draw)
class GenomeTrackViewer(MatplotlibDataViewer, PanTrackerMixin):

    LABEL = 'Genome Track Viewer'

    _layer_style_widget_cls = GenomeTrackLayerStyleEditor
    _options_cls = GenomeTrackOptionsWidget
    _state_cls = GenomeTrackState

    large_data_size = 2e7

    tools = ['select:xrange']

    def __init__(self, session, parent=None, state=None):
        super().__init__(session, parent=parent, state=state)
        self._layer_artist_container.on_changed(self.reflow_tracks)

        self.init_pan_tracking(self.axes)
        self._setup_annotation_track()

        self.state.add_callback('chr', self._redraw_annotation)
        self.state.add_callback('start', self._redraw_annotation)
        self.state.add_callback('end', self._redraw_annotation)
        self.state.add_callback('show_annotations', self.reflow_tracks)
        self._setup_zoom_to_layer_action()

    def _setup_zoom_to_layer_action(self):
        layer_artist_view = self._view.layer_list
        act = QtWidgets.QAction('Zoom to Layer', layer_artist_view)
        act.triggered.connect(nonpartial(self._zoom_to_layer))
        layer_artist_view.addAction(act)

    def _zoom_to_layer(self):
        layer = self._view.layer_list.current_artist().layer

        if not isinstance(layer, Subset):
            return
        try:
            chr, start, end = layer.subset_state.extent()
        except AttributeError:
            return

        with delay_callback(self.state, 'chr', 'start', 'end'):
            self.state.chr = chr.lstrip('chr')
            self.state.start = start
            self.state.end = end

    def get_data_layer_artist(self, layer=None, layer_state=None):
        cls = GenomeProfileLayerArtist if isinstance(layer, BedgraphData) else GenomeLoopLayerArtist
        result = self.get_layer_artist(cls, layer=layer, layer_state=layer_state)
        result.state.add_callback('zorder', self.reflow_tracks)
        result.state.add_callback('visible', self.reflow_tracks)
        return result

    def get_subset_layer_artist(self, layer=None, layer_state=None):
        data = layer.data if layer else None
        cls = GenomeProfileLayerArtist if isinstance(data, BedgraphData) else GenomeLoopLayerArtist
        result = self.get_layer_artist(cls, layer=layer, layer_state=layer_state)
        result.state.add_callback('zorder', self.reflow_tracks)
        result.state.add_callback('visible', self.reflow_tracks)
        return result

    def on_pan_end(self):
        self.reflow_tracks()

    def apply_roi(self, roi, override_mode=None):
        if len(self.layers) == 0:
            return

        x = roi.to_polygon()[0]
        chr = self.state.chr
        if not chr.startswith('chr'):
            chr = 'chr' + chr
        start, end = min(x), max(x)
        state = GenomicRangeSubsetState(chr, start, end)
        self.apply_subset_state(state, override_mode=override_mode)

    def _setup_annotation_track(self):
        self.annotations_ax = GenomeTrackLayerArtist._setup_track_axes(
            self.axes, 'annotations', self.state)
        annotation_path = os.environ.get('GLUEGENES_GENE_FILE', 'gencodeVM23_bed12.bed.bgz')
        if os.path.exists(annotation_path) and annotation_path is not None:
            pass
        else:
            annotation_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),'..','data','gencodeVM23_bed12.bed.bgz')
        print(f'ANNOTATION PATH: {annotation_path}')
        self.annotation = cb.BED(
            annotation_path,
            num_rows=None,
            gene_style='simple',
            bed_type='bed12',
            fontsize=8,
        )
        self.annotations_ax.get_yaxis().set_visible(False)
        self.annotations_ax.spines['left'].set_visible(False)

    def _redraw_annotation(self, *args):
        if self.panning:
            return
        self.annotations_ax.set_visible(self.state.show_annotations)

        if not self.state.show_annotations:
            self.axes.figure.canvas.draw_idle()
            return

        self.annotations_ax.clear()

        start = max(int(self.state.start), 1)
        end = max(int(self.state.end), 1)
        start, end = min(start, end), max(start, end)
        gr = GenomeRange('chr' + self.state.chr, start, end)
        intervals = BedGraph._tabix_query(self.annotation.properties['file'], gr)

        if intervals:
            parsed = list(ReadBed('\t'.join(l) for l in intervals))
            bed_type = 'bed12'
            df = self.annotation.intervals2dataframe(parsed, bed_type)
            df.bed_type = bed_type
        else:
            df = pd.DataFrame()

        self.annotation.is_draw_labels = True
        self.annotation.plot_genes(
            self.annotations_ax,
            cb.GenomeRange(f"{self.state.chr}:{start}-{end}"),
            df
        )
        self.axes.figure.canvas.draw_idle()

    def reflow_tracks(self, *args):
        """
        Reorder each track top->bottom sorted by zorder, removing any gaps of non-visible tracks.
        """
        axes = {}
        for artist in self._layer_artist_container.artists:
            if not artist.state.visible:
                continue
            axes[artist.track_axes] = min(axes.get(artist.track_axes, np.inf), artist.zorder)

        axes = list(sorted(axes, key=axes.get))
        if self.state.show_annotations:
            axes = [self.annotations_ax] + axes

        track_cnt = len(axes)
        height = 0.9 / track_cnt
        for i, ax in enumerate(axes):
            bounds = _TransformedBoundsLocator([0, i / track_cnt, 1, height], ax.axes.transAxes)
            ax.set_axes_locator(bounds)

        self._redraw_annotation()
