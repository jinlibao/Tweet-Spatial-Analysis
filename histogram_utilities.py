import pandas as pd
import numpy as np

from bokeh.models import ColumnDataSource, CustomJS, HoverTool, TapTool
from bokeh.palettes import plasma
from bokeh.plotting import figure


class HistogramData:

    def __init__(self, data_df, bins_text_list=None, bins_list=20):
        # If bins is an int, it defines the number of equal-width bins in the given range (10, by default).
        # If bins is a sequence, it defines the bin edges, including the rightmost edge, allowing for non-uniform bin widths.
        hist, edges = np.histogram(data_df, bins=bins_list)

        bin_list_len = len(hist)
        bottom_list = [0] * bin_list_len
        #fill_color_list = ['red'] * bin_list_len
        fill_color_list = plasma(bin_list_len)
        fill_color_list.reverse()
        line_color_list = ['black'] * bin_list_len
        id_list = [-1] * bin_list_len

        if bins_text_list is None:
            bins_text_list = [""] * bin_list_len

        data_histogram_cds = {
            'bottom': bottom_list,
            'top': hist,
            'left': edges[:-1],
            'right': edges[1:],
            'fill_color': fill_color_list,
            'line_color': line_color_list,
            'bins_text': bins_text_list,
            'id': id_list
        }

        self.df_histogram_cds = pd.DataFrame(data_histogram_cds)


class HistogramController:

    def __init__(self, df_cds):
        self.chr = None

        self.cds = ColumnDataSource(data=dict(bottom=[], top=[], left=[], right=[], fill_color=[], line_color=[]))
        self.cds.data = self.cds.from_df(df_cds)

        self.hover_idx = ColumnDataSource(data=dict(idx=[]))
        data_hover_idx = {
            'idx': [-1]
        }
        df_hover_idx = pd.DataFrame(data_hover_idx)
        self.hover_idx.data = self.hover_idx.from_df(df_hover_idx)

        self.selected = ColumnDataSource(data=dict(x=[], y=[], idx=[]))
        self.selected.on_change('data', self.selected_changed)

    def selected_changed(self, attrname, old, new):
        if len(new['idx']) > 0:
            print("Selected Count: idx: " + str(new['idx']))


class HistogramPlot:

    def __init__(   self,
                    histogram_controller, callback_hover, callback_tap,
                    title, x_axis_label, y_axis_label
                 ):

        self.histogram_controller = histogram_controller

        self.p = figure(    plot_height=200, plot_width=400,
                            title=title,
                            x_axis_label=x_axis_label,
                            y_axis_label=y_axis_label,
                            tools=["pan, wheel_zoom"], )

        self.r = self.p.quad(   bottom='bottom',
                                top='top',
                                left='left',
                                right='right',
                                source=self.histogram_controller.cds,
                                fill_color='fill_color', line_color='line_color')

        self.dummy_r = self.p.circle(   x='x', y='y',
                                        source=self.histogram_controller.selected,
                                        fill_color='white', line_color=None, fill_alpha=0.0, size=1)

        hover_tool = HoverTool(tooltips=[('count', "@top"), ('bin', "@bins_text")],
                                callback=CustomJS.from_py_func(callback_hover), renderers=[self.r])

        tap_tool = TapTool(callback=CustomJS.from_py_func(callback_tap), renderers=[self.r])
        self.p.add_tools(tap_tool, hover_tool)
