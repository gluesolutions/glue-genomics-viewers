import os
import operator
import math
#import numpy as np

from glue.core.subset import roi_to_subset_state, InequalitySubsetState, MultiOrState, combine_multiple#, ElementSubsetState
from glue.core.coordinates import Coordinates, LegacyCoordinates
from glue.core.coordinate_helpers import dependent_axes
from glue.core.util import update_ticks

from glue.viewers.scatter.layer_artist import ScatterLayerArtist
from glue.viewers.image.compat import update_image_viewer_state #This will not work

from glue.viewers.image.frb_artist import imshow
from glue.viewers.image.composite_array import CompositeArray

from .layer_artist import HeatmapLayerArtist, HeatampSubsetLayerArtist


__all__ = ['MatplotlibHeatmapMixin']


EXTRA_FOOTER = """
# Set tick label size - for now tick_params (called lower down) doesn't work
# properly, but these lines won't be needed in future.
ax.coords[{x_att_axis}].set_ticklabel(size={x_ticklabel_size})
ax.coords[{y_att_axis}].set_ticklabel(size={y_ticklabel_size})
""".strip()


class MatplotlibHeatmapMixin(object):

    def setup_callbacks(self):
        self._wcs_set = False
        self._changing_slice_requires_wcs_update = None
        self.axes.set_adjustable('datalim')
        self.state.add_callback('x_att', self._update_axes)
        self.state.add_callback('y_att', self._update_axes)
        self.state.add_callback('reference_data', self._set_wcs)
        self.axes._composite = CompositeArray()
        self.axes._composite_image = imshow(self.axes, self.axes._composite, aspect='auto',
                                            origin='lower', interpolation='nearest')
        #self._set_wcs()
    def _update_axes(self, *args):
        if self.state.x_att_world is not None:
            self.state.x_axislabel = self.state.x_att_world.label
            x_ticks = self.state.reference_data.coords.get_tick_labels(self.state.x_axislabel) 
            update_ticks(self.axes, 'x', ['categorical'], False, x_ticks)
            #self.axes.xaxis.set_tick_params(rotation=0)
        if self.state.y_att_world is not None:
            self.state.y_axislabel = self.state.y_att_world.label
            y_ticks =self.state.reference_data.coords.get_tick_labels(self.state.y_axislabel)
            update_ticks(self.axes, 'y', ['categorical'], False, y_ticks)

        self.state.reset_limits()
        self.axes.figure.canvas.draw_idle()

    def add_data(self, data):
        result = super(MatplotlibHeatmapMixin, self).add_data(data)
        return result

    def _set_wcs(self, event=None, relim=True):
    
        if self.state.x_att is None or self.state.y_att is None or self.state.reference_data is None:
            return
    
        ref_coords = getattr(self.state.reference_data, 'coords', None)
        
        #if ref_coords is None or isinstance(ref_coords, LegacyCoordinates):
        #    self.axes.reset_wcs(slices=self.state.wcsaxes_slice,
        #                        wcs=get_identity_wcs(self.state.reference_data.ndim))
        #else:
        #    self.axes.reset_wcs(slices=self.state.wcsaxes_slice, wcs=ref_coords)

        
        # Reset the axis labels to match the fact that the new axes have no labels
        self.state.x_axislabel = ''
        self.state.y_axislabel = ''
    
        self._update_appearance_from_settings()
        self._update_axes()

        #self.update_x_ticklabel()
        #self.update_y_ticklabel()
    
        if relim:
            self.state.reset_limits()
        

        # Determine whether changing slices requires changing the WCS
        if ref_coords is None or type(ref_coords) == Coordinates:
            self._changing_slice_requires_wcs_update = False
        else:
            ix = self.state.x_att.axis
            iy = self.state.y_att.axis
            x_dep = list(dependent_axes(ref_coords, ix))
            y_dep = list(dependent_axes(ref_coords, iy))
            if ix in x_dep:
                x_dep.remove(ix)
            if iy in x_dep:
                x_dep.remove(iy)
            if ix in y_dep:
                y_dep.remove(ix)
            if iy in y_dep:
                y_dep.remove(iy)
            self._changing_slice_requires_wcs_update = bool(x_dep or y_dep)
    
        self._wcs_set = True
    
    def apply_roi(self, roi, override_mode=None):

        # Force redraw to get rid of ROI. We do this because applying the
        # subset state below might end up not having an effect on the viewer,
        # for example there may not be any layers, or the active subset may not
        # be one of the layers. So we just explicitly redraw here to make sure
        # a redraw will happen after this method is called.
        self.redraw()

        if len(self.layers) == 0:
            return

        if self.state.x_att is None or self.state.y_att is None or self.state.reference_data is None:
            return

        # We do NOT want to define a subset on the pixel coordinates, since they are not
        # meaningful in a data matrix/heatmap. Instead, we will do the following:
        #   1) Get extent of ROI
        #   2) Use coords to get a reverse mapping from pixel x back to 
        #   1) Get a list of all ids in a roi
        #   2) Create an OrState subset concatenating them
        
        # FIXME: A subset defined this way works well in the HeatmapViewer, and propagates correctly
        # through join_on_key, but I think join_on_key gets broken by clustering...
        
        cmin = round(roi.min) # this is not exactly the same logic as roi_to_subset, but it works okay
        cmax = math.ceil(roi.max)

        if roi.ori == 'x':
            selection_component_id = self.state.x_linked_attribute or self.state.x_real_attribute #select on linked_attribute if available, otherwise on real one
            ticks = self.state.reference_data.coords.x_axis_ticks[cmin:cmax]
            #selection_component_id = self.state.x_linked_attribute#.id['orsam_id']

        elif roi.ori == 'y':
            selection_component_id = self.state.y_linked_attribute or self.state.y_real_attribute
            ticks = self.state.reference_data.coords.y_axis_ticks[cmin:cmax] 

        # By far the most stable solution is to define the subset ON the metadata data set, which means
        # that we need a reference to it...
        

        states = [InequalitySubsetState(selection_component_id, int(t), operator.eq) for t in ticks]
        subset_state = MultiOrState(states)
        #subset_state = ElementSubsetState(indices = np.array([0,30000]), data = self.state.reference_data)

        #from glue.core.command import ApplySubsetState

        #cmd = ApplySubsetState(data_collection=self._data,
        #                        subset_state=subset_state,
        #                        override_mode=override_mode)
        #self._session.command_stack.do(cmd)

        self.apply_subset_state(subset_state, override_mode=override_mode)

    def _scatter_artist(self, axes, state, layer=None, layer_state=None):
        if len(self._layer_artist_container) == 0:
            raise Exception("Can only add a scatter plot overlay once an image is present")
        return ScatterLayerArtist(axes, state, layer=layer, layer_state=None)

    def get_data_layer_artist(self, layer=None, layer_state=None):
        if layer.ndim == 1:
            cls = self._scatter_artist
        else:
            cls = HeatmapLayerArtist
        return self.get_layer_artist(cls, layer=layer, layer_state=layer_state)

    def get_subset_layer_artist(self, layer=None, layer_state=None):
        if layer.ndim == 1:
            cls = self._scatter_artist
        else:
            cls = HeatampSubsetLayerArtist
        return self.get_layer_artist(cls, layer=layer, layer_state=layer_state)

    @staticmethod
    def update_viewer_state(rec, context):
        return update_image_viewer_state(rec, context)

    def show_crosshairs(self, x, y):

        if getattr(self, '_crosshairs', None) is not None:
            self._crosshairs.remove()

        self._crosshairs, = self.axes.plot([x], [y], '+', ms=12,
                                           mfc='none', mec='#d32d26',
                                           mew=1, zorder=100)

        self.axes.figure.canvas.draw_idle()

    def hide_crosshairs(self):
        if getattr(self, '_crosshairs', None) is not None:
            self._crosshairs.remove()
            self._crosshairs = None
            self.axes.figure.canvas.draw_idle()

    def _script_header(self):

        imports = []
        imports.append('import matplotlib.pyplot as plt')
        imports.append('from glue.viewers.matplotlib.mpl_axes import init_mpl')
        imports.append('from glue.viewers.image.composite_array import CompositeArray')
        imports.append('from glue.viewers.image.frb_artist import imshow')

        script = ""
        script += "fig, ax = init_mpl(wcs=True)\n"
        script += "ax.set_aspect('{0}')\n".format(self.state.aspect)

        script += '\ncomposite = CompositeArray()\n'
        script += "image = imshow(ax, composite, origin='lower', interpolation='nearest', aspect='{0}')\n\n".format(self.state.aspect)

        dindex = self.session.data_collection.index(self.state.reference_data)

        script += "ref_data = data_collection[{0}]\n".format(dindex)

        ref_coords = self.state.reference_data.coords

        if hasattr(ref_coords, 'wcs'):
            script += "ax.reset_wcs(slices={0}, wcs=ref_data.coords.wcs)\n".format(self.state.wcsaxes_slice)
        elif hasattr(ref_coords, 'wcsaxes_dict'):
            raise NotImplementedError()
        else:
            imports.append('from glue.viewers.image.viewer import get_identity_wcs')
            script += "ax.reset_wcs(slices={0}, wcs=get_identity_wcs(ref_data.ndim))\n".format(self.state.wcsaxes_slice)

        script += "# for the legend\n"
        script += "legend_handles = []\n"
        script += "legend_labels = []\n"
        script += "legend_handler_dict = dict()\n\n"

        return imports, script

    def _script_footer(self):
        imports, script = super(MatplotlibImageMixin, self)._script_footer()
        options = dict(x_att_axis=0 if self.state.x_att is None else self.state.reference_data.ndim - self.state.x_att.axis - 1,
                       y_att_axis=1 if self.state.y_att is None else self.state.reference_data.ndim - self.state.y_att.axis - 1,
                       x_ticklabel_size=self.state.x_ticklabel_size,
                       y_ticklabel_size=self.state.y_ticklabel_size)
        return [], EXTRA_FOOTER.format(**options) + os.linesep * 2 + script
