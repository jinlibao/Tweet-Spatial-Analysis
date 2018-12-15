from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import Circle, CustomJS, HoverTool, Range1d, TapTool
from bokeh.plotting import figure
from bokeh.tile_providers import CARTODBPOSITRON

from data_preprocessing import TweetDataPreProcessing
from file_utilities import *
from tweet_data import *
from user_profile_data import UserProfileDetails
from widget_utilities import *

pre_processor = TweetDataPreProcessing(None)
pre_processor.read_from_json(   "Tweet-Spatial-Analysis/data/tweet_mean_all.json",
                                "Tweet-Spatial-Analysis/data/tweets_median_working.json",
                                "Tweet-Spatial-Analysis/data/tweets_median_non_working.json")

file_open = FileOpen("Tweet-Spatial-Analysis/data", "user-info.csv")
user_info = UserProfileDetails(file_open)
user_info.process()

tweet_data_controller = TweetDataController(pre_processor, user_info)
map_widgets = MapWidgets(tweet_data_controller, user_info)
tweet_data_controller.selection_details = map_widgets.text_selection_details

latitude_range = [37, 43]
longitude_range = [-78.0, -72.0]

east_min, north_min = lon_lat_to_east_north(longitude_range[0], latitude_range[0])
east_max, north_max = lon_lat_to_east_north(longitude_range[1], latitude_range[1])
# print(str(east_min) + ", " + str(north_min) + " to " + str(east_max) + ", " + str(north_max))
# -8682920.281875338, 4439106.787250583 to -8015003.337115697, 5311971.846945471

# Configuring Plot Tools: https://bokeh.pydata.org/en/latest/docs/user_guide/tools.html
p = figure(plot_width=800, plot_height=800,
           tools=["pan,wheel_zoom,save"],
           x_range=Range1d(east_min, east_max), y_range=Range1d(north_min, north_max),
           title="Tweet Spatial Analysis", toolbar_location="above",
           output_backend="webgl",
           lod_factor=10, lod_threshold=100)
p.add_tile(CARTODBPOSITRON)

circles_renderer = p.circle(x='x', y='y', source=tweet_data_controller.circles, fill_color='color', line_color=None, fill_alpha=0.9, size=5)
# There is either something I'm missing, or potentially a problem with bokeh:
# Related to task 065: BBug with circle_renderer.selection_glyph and point 2429397887 - cheated by setting alpha to 0.0 to work:
# To replicate: Find 2429397887, Switch to working, Switch back to all - erronous point appears with id 2429397887 ?
# Cheated by making alpha 0.0.
circles_renderer.selection_glyph = Circle(fill_alpha=0.0, fill_color="olivedrab", line_color=None)
circles_renderer.nonselection_glyph = Circle(fill_alpha=0.2, fill_color="blue", line_color=None)

sde_ellipse_renderer = p.ellipse(x='x', y='y', width='width', height='height', angle='angle', source=tweet_data_controller.sde_ellipse, fill_color="#cab2d6", fill_alpha=0.5)
siblings_renderer = p.circle(x='x', y='y', source=tweet_data_controller.siblings, fill_color="orange", line_color="orange", fill_alpha=0.8, size=3)
sibling_ellipses_renderer = p.ellipse(x='x', y='y', width='width', height='height', angle='angle', source=tweet_data_controller.sibling_ellipses, line_color="darkmagenta", fill_alpha=0.0)
patch_dissolve_renderer = p.patch(x='x', y='y', source=tweet_data_controller.patch_dissolve, fill_color="wheat", line_color="wheat", fill_alpha=0.4)
selected_circle_renderer = p.circle(x='x', y='y', source=tweet_data_controller.selected_circle, fill_color="olivedrab", line_color=None, fill_alpha=1, size=5)
find_circle_renderer = p.circle(x='x', y='y', source=tweet_data_controller.find_circle, line_color="#410967", fill_color="orange", fill_alpha=0.0, size=15)


# Causes:
# Uncaught Error: reference {"id":"1005","type":"ColumnDataSource"} isn't known (not in Document?)
# Thrown when running serve, but not when running directly!!
def callback_hover(hover_idx = tweet_data_controller.hover_idx):
    indices = cb_data.index["1d"].indices
    if len(indices) > 0:
        hover_idx.data['idx'] = [indices[0]]

hover_tool = HoverTool(tooltips=[('id', "@id"), ('area', '@area'), ('count', '@count')], callback=CustomJS.from_py_func(callback_hover), renderers=[circles_renderer])

# This is only called when a glyph is tapped on.
# Tapping on the 'bare' map isn't registered.
def callback_tap(   hover_idx = tweet_data_controller.hover_idx,
                    circles = tweet_data_controller.circles,
                    selected_circle=tweet_data_controller.selected_circle,
                    toggle_sde_e = map_widgets.toggle_sde_ellipse,
                    toggle_se = map_widgets.toggle_sibling_ellipses,
                    toggle_d = map_widgets.toggle_dissolve
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

tap_tool = TapTool(callback = CustomJS.from_py_func(callback_tap), renderers=[circles_renderer])
p.add_tools(hover_tool, tap_tool)

lhs = column(   map_widgets.radio_button_data_type, map_widgets.text_selection_details,
                map_widgets.toggle_sde_ellipse, map_widgets.toggle_sibling_ellipses, map_widgets.toggle_dissolve,
                map_widgets.text_input, map_widgets.button_find,
                map_widgets.toggle_blend, map_widgets.slider_blend)

rhs = column(   row(    column(map_widgets.button_count_start_minus, map_widgets.button_count_start_plus, width=50),
                        map_widgets.range_slider_count,
                        column(map_widgets.button_count_end_minus, map_widgets.button_count_end_plus)),
                row(    column(map_widgets.button_area_start_minus, map_widgets.button_area_start_plus, width=50),
                        map_widgets.range_slider_area,
                        column(map_widgets.button_area_end_minus, map_widgets.button_area_end_plus)),
                map_widgets.text_count,
                map_widgets.filters_active)

l = row(lhs, p, rhs)

curdoc().add_root(l)
curdoc().title = "Tweet Spatial Analysis"
