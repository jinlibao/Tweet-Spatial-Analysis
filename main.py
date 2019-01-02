
from bokeh import events
from bokeh.io import curdoc
from bokeh.layouts import Column, Row
from bokeh.models import Circle, CustomJS, HoverTool, Range1d, TapTool, Panel, Tabs
from bokeh.plotting import figure
from bokeh.tile_providers import CARTODBPOSITRON

from configuration_utilities import *
from data_preprocessing import TweetDataPreProcessing
from file_utilities import *
from tweet_data import *
from user_profile_data import UserProfileDetails
from widget_utilities import *

tweet_spatial_analysis_config = TweetSpatialAnalysisConfig("Tweet-Spatial-Analysis/tweet_spatial_analysis.ini")
logger.info(tweet_spatial_analysis_config)

pre_processor = TweetDataPreProcessing(None)
pre_processor.read_from_json(   "Tweet-Spatial-Analysis/data/tweet_mean_all.json",
                                "Tweet-Spatial-Analysis/data/tweets_median_working.json",
                                "Tweet-Spatial-Analysis/data/tweets_median_non_working.json",
                                tweet_spatial_analysis_config)

file_open = FileOpen("Tweet-Spatial-Analysis/data", "user-info.csv")
user_info = UserProfileDetails(file_open)
user_info.process()

tweet_data_controller = TweetDataController(pre_processor, tweet_spatial_analysis_config, user_info)
map_widgets = MapWidgets(tweet_data_controller, tweet_spatial_analysis_config)
tweet_data_controller.selection_details = map_widgets.text_selection_details

east_min, north_min = lon_lat_to_east_north(tweet_spatial_analysis_config.longitude[0], tweet_spatial_analysis_config.latitude[0])
east_max, north_max = lon_lat_to_east_north(tweet_spatial_analysis_config.longitude[1], tweet_spatial_analysis_config.latitude[1])

# Configuring Plot Tools: https://bokeh.pydata.org/en/latest/docs/user_guide/tools.html
p = figure(plot_width=800, plot_height=800,
           tools=["pan, wheel_zoom"],
           x_range=Range1d(east_min, east_max), y_range=Range1d(north_min, north_max),
           title="Tweet Spatial Analysis", toolbar_location="above",
           output_backend="webgl",
           lod_factor=10, lod_threshold=100,
           x_axis_type="mercator", y_axis_type="mercator")
p.add_tile(CARTODBPOSITRON)
p.xaxis[0].axis_label = 'West - East: WSG 84 Web Mercator'
p.yaxis[0].axis_label = 'South - North: WSG 84 Web Mercator'

circles_renderer = p.circle(x='x', y='y', source=tweet_data_controller.circles, fill_color='color', line_color=None, fill_alpha=0.9, size=5)

# There is either something I'm missing, or potentially a problem with bokeh:
# Related to task 065: BBug with circle_renderer.selection_glyph and point 2429397887 - cheated by setting alpha to 0.0 to work:
# To replicate: Find 2429397887, Switch to working, Switch back to all - erronous point appears with id 2429397887 ?
# Cheated by making alpha 0.0.
circles_renderer.selection_glyph = Circle(fill_alpha=0.0, fill_color="olivedrab", line_color=None)
circles_renderer.nonselection_glyph = Circle(fill_alpha=0.2, fill_color="blue", line_color=None)

# This is a dummy renderer needed to enable changes on the tweet_data_controller.lod_dummy source to be emitted.
# and its related on_change callback to be called.
lod_dummy_renderer = p.circle(x='x', y='y', source=tweet_data_controller.lod_dummy, fill_color='white', line_color=None, fill_alpha=0.0, size=1)

patch_dissolve_renderer = p.patch(x='x', y='y', source=tweet_data_controller.patch_dissolve, fill_color="wheat", line_color="wheat", line_alpha=0.4, fill_alpha=0.4)
sde_ellipse_renderer = p.ellipse(x='x', y='y', width='width', height='height', angle='angle', source=tweet_data_controller.sde_ellipse, fill_color="#cab2d6", line_alpha=0.5, fill_alpha=0.5)
sibling_ellipses_renderer = p.ellipse(x='x', y='y', width='width', height='height', angle='angle', source=tweet_data_controller.sibling_ellipses, line_color="darkmagenta", fill_alpha=0.0)
siblings_renderer = p.circle(x='x', y='y', source=tweet_data_controller.siblings, fill_color="orange", line_color="orange", fill_alpha=0.8, size=3)
selected_circle_renderer = p.circle(x='x', y='y', source=tweet_data_controller.selected_circle, fill_color="olivedrab", line_color=None, fill_alpha=1, size=5)

patch_dissolve_blend_renderer = p.patch(x='x', y='y', source=tweet_data_controller.patch_dissolve_blend, fill_color="wheat", line_color="wheat", line_alpha=0.4, fill_alpha=0.4)
sde_ellipse_blend_renderer = p.ellipse(x='x', y='y', width='width', height='height', angle='angle', source=tweet_data_controller.sde_ellipse_blend, fill_color="#cab2d6", line_alpha=0.5, fill_alpha=0.5)
sibling_ellipses_blend_renderer = p.ellipse(x='x', y='y', width='width', height='height', angle='angle', source=tweet_data_controller.sibling_ellipses_blend, line_color="darkmagenta", fill_alpha=0.0)
siblings_blend_renderer = p.circle(x='x', y='y', source=tweet_data_controller.siblings_blend, fill_color="orange", line_color="orange", fill_alpha=0.8, size=3)
selected_circle_blend_renderer = p.circle(x='x', y='y', source=tweet_data_controller.selected_circle_blend, fill_color="indigo", line_color=None, fill_alpha=1, size=5)

find_circle_renderer = p.circle(x='x', y='y', source=tweet_data_controller.find_circle, line_color="#410967", fill_color="orange", fill_alpha=0.0, size=15)

tweet_data_controller.cr = circles_renderer
tweet_data_controller.csr = selected_circle_renderer
tweet_data_controller.csbr = selected_circle_blend_renderer
tweet_data_controller.er = sde_ellipse_renderer
tweet_data_controller.ebr = sde_ellipse_blend_renderer
tweet_data_controller.ser = sibling_ellipses_renderer
tweet_data_controller.sebr = sibling_ellipses_blend_renderer
tweet_data_controller.pdr = patch_dissolve_renderer
tweet_data_controller.pdbr = patch_dissolve_blend_renderer
tweet_data_controller.sr = siblings_renderer
tweet_data_controller.sbr = siblings_blend_renderer

# Fake a Legend
for idx in range(0, len(tweet_data_controller.histogram_controller_count.cds.data['fill_color'])):
    legend_text = tweet_data_controller.histogram_controller_count.cds.data['bins_text'][idx]
    legend_color = tweet_data_controller.histogram_controller_count.cds.data['fill_color'][idx]
    p.circle(0, 0, legend=legend_text, fill_color=legend_color, alpha = 0.75, size = 1, line_color=None)
p.legend.background_fill_alpha = 0.25

# Causes:
# Uncaught Error: reference {"id":"1005","type":"ColumnDataSource"} isn't known (not in Document?)
# Thrown when running serve, but not when running directly!!
def callback_hover(hover_idx = tweet_data_controller.hover_idx):
    indices = cb_data.index["1d"].indices
    if len(indices) > 0:
        hover_idx.data['idx'] = [indices[0]]

hover_tool = HoverTool(tooltips=[('id', "@id"), ('area', '@area'), ('count', '@count'), ('distance', '@distance'), ('ratio', '@ratio')], callback=CustomJS.from_py_func(callback_hover), renderers=[circles_renderer])

# This is only called when a glyph is tapped on.
# Tapping on the 'bare' map isn't registered.
def callback_tap(   hover_idx = tweet_data_controller.hover_idx,
                    circles = tweet_data_controller.circles,
                    selected_circle = tweet_data_controller.selected_circle,
                    toggle_sde_e = map_widgets.toggle_sde_ellipse,
                    toggle_se = map_widgets.toggle_sibling_ellipses,
                    toggle_d = map_widgets.toggle_dissolve,
                    toggle_b = map_widgets.toggle_blend,
                    ba = tweet_data_controller.blend_active
                    ):
    idx = hover_idx.data['idx']

    new_data = dict()
    new_data['x'] = [circles.data['x'][idx]]
    new_data['y'] = [circles.data['y'][idx]]
    new_data['id'] = [circles.data['id'][idx]]
    selected_circle.data = new_data

    toggle_sde_e.active = False
    toggle_se.active = False
    toggle_d.active = False
    toggle_b.active = False
    ba.active = False

tap_tool = TapTool(callback = CustomJS.from_py_func(callback_tap), renderers=[circles_renderer])
p.add_tools(hover_tool, tap_tool)


def callback_hover_count(
        hhci = tweet_data_controller.histogram_controller_count.hover_idx
    ):
    indices = cb_data.index["1d"].indices
    if len(indices) > 0:
        print([indices[0]])
        hhci.data['idx'] = [indices[0]]

def callback_tap_count(
        hhci = tweet_data_controller.histogram_controller_count.hover_idx,
        sc = tweet_data_controller.histogram_controller_count.selected
    ):

    idx = hhci.data['idx']

    new_data = dict()
    new_data['x'] = [0]
    new_data['y'] = [0]
    new_data['idx'] = [idx[0]]
    sc.data = new_data

histogram_plot_count = HistogramPlot(tweet_data_controller.histogram_controller_count, callback_hover_count, callback_tap_count,
                                     "Number of Siblings within a Tweet Data Point", "Number of Siblings", "Tweets with Count")
tweet_data_controller.hr_count = histogram_plot_count.r



def callback_hover_area(
        hhci = tweet_data_controller.histogram_controller_area.hover_idx
    ):
    indices = cb_data.index["1d"].indices
    if len(indices) > 0:
        print([indices[0]])
        hhci.data['idx'] = [indices[0]]

def callback_tap_area(
        hhci = tweet_data_controller.histogram_controller_area.hover_idx,
        sc = tweet_data_controller.histogram_controller_area.selected
    ):

    idx = hhci.data['idx']

    new_data = dict()
    new_data['x'] = [0]
    new_data['y'] = [0]
    new_data['idx'] = [idx[0]]
    sc.data = new_data

histogram_plot_area = HistogramPlot(tweet_data_controller.histogram_controller_area, callback_hover_area, callback_tap_area, "Area of SDE", "Area", "Count")
tweet_data_controller.hr_area = histogram_plot_area.r


def callback_hover_distance(
        hhci = tweet_data_controller.histogram_controller_distance.hover_idx
    ):
    indices = cb_data.index["1d"].indices
    if len(indices) > 0:
        print([indices[0]])
        hhci.data['idx'] = [indices[0]]

def callback_tap_distance(
        hhci = tweet_data_controller.histogram_controller_distance.hover_idx,
        sc = tweet_data_controller.histogram_controller_distance.selected
    ):

    idx = hhci.data['idx']

    new_data = dict()
    new_data['x'] = [0]
    new_data['y'] = [0]
    new_data['idx'] = [idx[0]]
    sc.data = new_data

histogram_plot_distance = HistogramPlot(tweet_data_controller.histogram_controller_distance, callback_hover_distance, callback_tap_distance, "Distance between W and Non-W", "Distance", "Count")
tweet_data_controller.hr_distance = histogram_plot_distance.r


def callback_hover_ratio(
        hhci = tweet_data_controller.histogram_controller_ratio.hover_idx
    ):
    indices = cb_data.index["1d"].indices
    if len(indices) > 0:
        print([indices[0]])
        hhci.data['idx'] = [indices[0]]

def callback_tap_ratio(
        hhci = tweet_data_controller.histogram_controller_ratio.hover_idx,
        sc = tweet_data_controller.histogram_controller_ratio.selected
    ):

    idx = hhci.data['idx']

    new_data = dict()
    new_data['x'] = [0]
    new_data['y'] = [0]
    new_data['idx'] = [idx[0]]
    sc.data = new_data

histogram_plot_ratio = HistogramPlot(tweet_data_controller.histogram_controller_ratio, callback_hover_ratio, callback_tap_ratio, "X/Y Ratio", "Ratio", "Count")
tweet_data_controller.hr_ratio = histogram_plot_ratio.r


def lod_start(ld = tweet_data_controller.lod_dummy):
    print("LOD Start:")
    new_data = dict()
    new_data['x'] = [0]
    new_data['y'] = [0]
    new_data['lod'] = [0]
    ld.data = new_data
    ld.change.emit()

p.js_on_event(events.LODStart, CustomJS.from_py_func(lod_start))

def lod_end(ld = tweet_data_controller.lod_dummy):
    print("LOD End:")
    new_data = dict()
    new_data['x'] = [0]
    new_data['y'] = [0]
    new_data['lod'] = [1]
    ld.data = new_data
    ld.change.emit()

p.js_on_event(events.LODEnd, CustomJS.from_py_func(lod_end))

def pan_start(ld = tweet_data_controller.lod_dummy):
    print("Pan Start:")
    new_data = dict()
    new_data['x'] = [0]
    new_data['y'] = [0]
    new_data['lod'] = [0]
    ld.data = new_data


p.js_on_event(events.PanStart, CustomJS.from_py_func(pan_start))

def pan_end(ld = tweet_data_controller.lod_dummy):
    print("Pan End:")
    new_data = dict()
    new_data['x'] = [0]
    new_data['y'] = [0]
    new_data['lod'] = [1]
    ld.data = new_data

p.js_on_event(events.PanEnd, CustomJS.from_py_func(pan_end))

#def mouse_move(ld = tweet_data_controller.lod_dummy):
#    print("Mouse Move:")

#p.js_on_event(events.MouseMove, CustomJS.from_py_func(mouse_move))

#def mouse_wheel(ld = tweet_data_controller.lod_dummy):
#    print("Mouse Wheel:")

#p.js_on_event(events.MouseWheel, CustomJS.from_py_func(mouse_wheel))

lhs = Column(   map_widgets.radio_button_data_type, map_widgets.toggle_data, map_widgets.text_selection_details,
                map_widgets.toggle_sde_ellipse, map_widgets.toggle_sibling_ellipses, map_widgets.toggle_dissolve,
                map_widgets.toggle_blend, map_widgets.slider_blend,
                map_widgets.text_input, map_widgets.button_find)

filter_sliders \
    = Column(   Row(    Column(map_widgets.button_count_start_minus, map_widgets.button_count_start_plus, width=50),
                        map_widgets.range_slider_count,
                        Column(map_widgets.button_count_end_minus, map_widgets.button_count_end_plus)),
                Row(    Column(map_widgets.button_area_start_minus, map_widgets.button_area_start_plus, width=50),
                        map_widgets.range_slider_area,
                        Column(map_widgets.button_area_end_minus, map_widgets.button_area_end_plus)),
                Row(    Column(map_widgets.button_distance_start_minus, map_widgets.button_distance_start_plus, width=50),
                        map_widgets.range_slider_distance,
                        Column(map_widgets.button_distance_end_minus, map_widgets.button_distance_end_plus)),
                Row(    Column(map_widgets.button_ratio_start_minus, map_widgets.button_ratio_start_plus, width=50),
                        map_widgets.range_slider_ratio,
                        Column(map_widgets.button_ratio_end_minus, map_widgets.button_ratio_end_plus)),
                map_widgets.text_count,
                map_widgets.filters_active)


filter_histograms = Column(histogram_plot_count.p, histogram_plot_area.p, histogram_plot_distance.p, histogram_plot_ratio.p)

tab_filter_sliders = Panel(child=filter_sliders, title="Filter by Sliders")
tab_filter_histograms = Panel(child=filter_histograms, title="Filter by Histograms")

tabs_rhs = Tabs(tabs = [tab_filter_sliders, tab_filter_histograms, ])

l = Row(lhs, p, tabs_rhs)

curdoc().add_root(l)
curdoc().title = "Tweet Spatial Analysis"


