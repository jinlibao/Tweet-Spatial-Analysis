import json
import logging.config

import pandas as pd

from bokeh.models import ColumnDataSource

from analysis_utilities import *
from file_utilities import *
from shapely_utilities import *
from histogram_utilities import *

logger = logging.getLogger()

def print_progress_to_console(idx):
    if idx == 0:
        print("", end="")
    elif idx % 1000 == 0:
        print("*", end="")
    elif idx % 100 == 0:
        print(".", end="")

class TweetData:

    def __init__(self, data_name):
        self.name = data_name
        self.df = None
        self.max_count = -1.0
        self.max_area = -1.0
        #self.max_xy_ratio = -1.0
        self.max_a = -1.0
        self.max_b = -1.0

    def read_from_json(self, filename):
        self.df = pd.read_json(filename, orient='records', lines=True)

    # Cleaning:
    # An area of 0 is allowed - potentially implying they don't move - always tweeting from the same location.
    # An area of NaN? Mark as '-1' to identify undefined.
    # For both cases mark s and b as -1
    # Angle ranges from -1.57 to 1.57- set nan values to 0
    def create_dataframe(self,
                         data_df,
                         id_column_name,
                         latitude_column_name,
                         longitude_column_name,
                         area_column_name,
                         x_y_ratio,
                         angle_column_name,
                         distance_name):

        self.df = pd.DataFrame(columns=['id', 'lat', 'lon', 'x', 'y', 'area', 'a', 'b', 'angle', 'distance', 'ratio'])

        # Keep an eye open for: A value is trying to be set on a copy of a slice from a DataFrame
        # Caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
        # Recommends using df.loc[row_idx, col_idx]
        self.df['id'] = data_df[id_column_name].copy(True)
        self.df['lat'] = data_df[latitude_column_name].copy(True)
        self.df['lon'] = data_df[longitude_column_name].copy(True)
        self.df['area'] = data_df[area_column_name].copy(True)
        self.df['angle'] = data_df[angle_column_name].copy(True)
        self.df['distance'] = data_df[distance_name].copy(True)

        for row_idx in range(0, data_df.shape[0]):
            # If the area is nan, then set it to -1 to indicate it is undefined.
            if math.isnan(self.df['area'][row_idx]):
                logger.warning("idx: %s : id %s has area NaN.", row_idx, self.df['id'][row_idx])
                self.df.loc[row_idx, 'area'] = -1.0

            if math.isnan(self.df['distance'][row_idx]):
                logger.warning("idx: %s : id %s has distance NaN.", row_idx, self.df['id'][row_idx])
                self.df.loc[row_idx, 'distance'] = -1.0

            # If any angles are nan, then clean to 0.0
            #if (math.isnan(self.df['angle'][row_idx])):
            #    self.df.loc[row_idx, 'angle'] = 0.0

            self.df.loc[row_idx, 'x'], self.df.loc[row_idx, 'y'] = lon_lat_to_east_north(
                self.df['lon'][row_idx],
                self.df['lat'][row_idx])

            self.df.loc[row_idx, 'a'], self.df.loc[row_idx, 'b'], self.df.loc[row_idx, 'ratio'] = \
                self.calculate_ellipse_properties(
                    self.df['id'][row_idx],
                    row_idx,
                    self.df['area'][row_idx],
                    data_df[x_y_ratio][row_idx])

            print_progress_to_console(row_idx)

        print("")
        self.max_a = self.df['a'].max()
        self.max_b = self.df['b'].max()
        self.max_area = self.df['area'].max()

        self.find_datapoints_enclosed_with_ellipse()
        self.calculate_colour()

        logger.info("\n%s", self.df.head())

    def find_datapoints_enclosed_with_ellipse(self):
        logger.info("Find Datapoints Enclosed within Ellipse and Dissolve.")
        print("Find Datapoints Enclosed within Ellipse and Dissolve.")
        row_idx = 0
        counts = []
        dissolves = []

        for row in self.df.itertuples():
            siblings_data = dict()
            selected_point_id = row.id
            xo = row.x
            yo = row.y
            a = row.a
            b = row.b
            theta = row.angle

            grouped_ids = []
            grouped_idxs = []
            row2_idx = 0
            for row2 in self.df.itertuples():
                if row2.id != selected_point_id:
                    if is_point_inside_ellipse(xo, yo, a, b, theta, row2.x, row2.y):
                        grouped_ids.append(row2.id)
                        grouped_idxs.append(row2_idx)
                row2_idx += 1

            count = len(grouped_ids)
            counts.append(count)

            dissolve_data = self.create_dissolve(selected_point_id, row_idx, grouped_idxs)
            dissolves.append(dissolve_data['area'])

            siblings_data['id'] = selected_point_id
            siblings_data['idx'] = row_idx
            siblings_data['count'] = count
            siblings_data['count'] = count
            siblings_data['ids'] = grouped_ids
            siblings_data['idxs'] = grouped_idxs
            siblings_data['dissolve'] = dissolve_data
            #print(siblings_data)

            #filename = self.name + "_" + str(selected_point_id) + "_" + str(row_idx) + ".json"
            filename = self.name + "_" + str(selected_point_id) + ".json"
            with open("sibling_data/" + filename, 'w') as fp:
                json.dump(siblings_data, fp)

            if count > self.max_count:
                self.max_count = count

            print_progress_to_console(row_idx)
            row_idx += 1

        print("")
        self.df['count'] = counts
        self.df['dissolve_area'] = dissolves

    def create_dissolve(self, id, idx, grouped_idxs):
        #print("Create Dissolve: " + str(circle_idx))
        #print(self.tweet_data_df['idxs'][circle_idx])
        ellipses = dict()
        ellipses['x'] = []
        ellipses['y'] = []
        ellipses['a'] = []
        ellipses['b'] = []
        ellipses['angle'] = []
        for row_idx in grouped_idxs:
            if self.df['a'][row_idx] > 0.0:
                ellipses['x'].append(self.df['x'][row_idx])
                ellipses['y'].append(self.df['y'][row_idx])
                ellipses['a'].append(self.df['a'][row_idx])
                ellipses['b'].append(self.df['b'][row_idx])
                ellipses['angle'].append(self.df['angle'][row_idx])

        # Need to remember to include the datapoints own ellipse details. Otherwise the following error occurs:
        # AttributeError("'MultiPolygon' object has no attribute 'exterior'")
        # Probably caused by the sibling ellipses NOT overlapping i.e. separate.
        if self.df['a'][idx] > 0.0:
            ellipses['x'].append(self.df['x'][idx])
            ellipses['y'].append(self.df['y'][idx])
            ellipses['a'].append(self.df['a'][idx])
            ellipses['b'].append(self.df['b'][idx])
            ellipses['angle'].append(self.df['angle'][idx])

        if any(ellipses['x']):
            new_data = dissolve_ellipses(ellipses)
        else:
            new_data = dict()
            new_data['x'] = []
            new_data['y'] = []
            new_data['area'] = 0.0
            logger.warning("Empty Ellipses: Point ID: %s : a:%s : b:%s : angle:%s : IDXs", id, self.df['a'][idx], self.df['b'][idx], self.df['angle'][idx], str(grouped_idxs))

        return new_data

    def calculate_colour(self):
        logger.info("Calculating Color.")
        print("Calculating Color.")

        max_count = self.max_count
        #count_max = 600.0

        #start_r = 178.0
        #start_g = 34.0
        #start_b = 34.0
        #end_r = 255.0
        #end_g = 215.0
        #end_b = 0.0
        end_r = 178.0
        end_g = 34.0
        end_b = 34.0
        start_r = 255.0
        start_g = 215.0
        start_b = 0.0

        #interval_r = (end_r - start_r) / count_max
        #interval_g = (end_g - start_g) / count_max
        #interval_b = (end_b - start_b) / count_max

        interval_r = (end_r - start_r)
        interval_g = (end_g - start_g)
        interval_b = (end_b - start_b)

        row_idx = 0
        colours = []
        for row in self.df.itertuples():
            count = row.count
            scale_factor = 1.0 - 1.0/(1.0 + count)

            point_r = hex(int(start_r + scale_factor * interval_r)).split('x')[-1]
            point_g = hex(int(start_g + scale_factor * interval_g)).split('x')[-1]
            point_b = hex(int(start_b + scale_factor * interval_b)).split('x')[-1]
            colour_hex_str = "#" + str(point_r) + str(point_g) + str(point_b)
            colours.append(colour_hex_str)

            print_progress_to_console(row_idx)
            row_idx += 1

        print("")
        self.df['color'] = colours

    # Using the x/y ratio and area calculate the major and minor axes for the ellipse.
    # If the area is 0.0 or -1, then we can not calculate these values
    def calculate_ellipse_properties(self, row_id, idx, area, xy_ratio):
        # Remember: Viewing the data in excel or notepad++, the first row is the headers, and the row count starts from 1
        # So, where the pandas index starts from 0, the starting row in excel or notepad++ is actually 2.
        if area == 0.0 or area == -1:
            logger.warning("id: %s (%s): The ellipse area is %s - can not calculate 'a' or 'b'.", row_id, idx, area)
            return -1, -1, -1
        else:
            if xy_ratio == "Infinity":
                logger.warning("id %s (%s): The ellipse x/y ratio is 'Infinity' - can not calculate 'a' or 'b'.", row_id, idx)
                return -1, -1, -1

            xy_ratio_f = float(xy_ratio)
            if math.isnan(xy_ratio_f):
                logger.warning("id %s (%s): The ellipse x/y ratio is 'nan' - can not calculate 'a' or 'b'.", row_id, idx)
                return -1, -1, -1

            if math.isinf(xy_ratio_f):
                logger.warning("id %s (idx): The ellipse x/y ratio is 'inf' - can not calculate 'a' or 'b'.", row_id, idx)
                return -1, -1, -1

            return calculate_ellipse_a_and_b_from_area_and_ratio(area, xy_ratio_f)

    def __str__(self):
        df_str = "TweetData:\nName: " + str(self.name)
        if self.df is not None:
            df_str += "\nShape: "+ str(self.df.shape)
            df_str += "\nColumn Names:\n" + str(self.df.columns)
            df_str += "\nColumn dtypes:\n" + str(self.df.dtypes)
        else:
            df_str += "\nDataframe NOT defined."

        df_str += "\nMax Count: " + str(self.max_count)
        df_str += "\nMax Area: " + str(self.max_area)
        df_str += "\nMax a: " + str(self.max_a)
        df_str += "\nMax b: " + str(self.max_b)

        return df_str

class MapTweetData:

    def __init__(self, tweet_data, config):
        self.name = tweet_data.name
        print(self.name)
        self.tweet_data_df = tweet_data.df

        self.histogram_count = HistogramData(self.tweet_data_df['count'], bins_text_list=config.bins_count_text, bins_list=config.bins_count)
        self.histogram_area = HistogramData(self.tweet_data_df['area'], bins_list=100)
        self.histogram_distance = HistogramData(self.tweet_data_df['distance'], bins_list=50)
        self.histogram_ratio = HistogramData(self.tweet_data_df['ratio'], bins_text_list=config.bins_ratio_text, bins_list=config.bins_ratio)
        self.histogram_dissolve = HistogramData(self.tweet_data_df['dissolve_area'], bins_list=100)

        self.sibling_data_manager = SiblingDataManager(config.sibling_data_folder, self.name)

    def clear_selected_circle(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['id'] = []
        return new_data

    def update_selected_circle(self, circle_id):
        id_df = self.tweet_data_df.loc[self.tweet_data_df['id'] == circle_id]
        idx = id_df.index
        circle_idx = idx[0]

        new_data = dict()
        new_data['x'] = [self.tweet_data_df['x'][circle_idx]]
        new_data['y'] = [self.tweet_data_df['y'][circle_idx]]
        new_data['id'] = [circle_id]
        return new_data, circle_idx

    def clear_sde_ellipse(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['width'] = []
        new_data['height'] = []
        new_data['angle'] = []
        return new_data

    def update_sde_ellipse(self, circle_idx, alpha = 0.5):
        #print("Update SDE Ellipse: " + str(circle_idx) + " : alpha: " + str(alpha))
        new_data = dict()
        new_data['x'] = [self.tweet_data_df['x'][circle_idx]]
        new_data['y'] = [self.tweet_data_df['y'][circle_idx]]
        new_data['width'] = [self.tweet_data_df['a'][circle_idx] * 2.0]
        new_data['height'] = [self.tweet_data_df['b'][circle_idx] * 2.0]
        new_data['angle'] = [self.tweet_data_df['angle'][circle_idx]]
        return new_data

    def clear_siblings(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        return new_data

    def update_siblings(self, id):
        #print("Update Siblings: " + str(id))
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []

        idxs_list = self.sibling_data_manager.get_sibling_idxs(id)
        for row_idx in idxs_list:
            new_data['x'].append(self.tweet_data_df['x'][row_idx])
            new_data['y'].append(self.tweet_data_df['y'][row_idx])
        return new_data

    def clear_sibling_ellipses(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['width'] = []
        new_data['height'] = []
        new_data['angle'] = []
        return new_data

    def update_sibling_ellipses(self, id):
        #print("Update Sibling Ellipses: " + str(id))
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['width'] = []
        new_data['height'] = []
        new_data['angle'] = []

        idxs_list = self.sibling_data_manager.get_sibling_idxs(id)
        for row_idx in idxs_list:
            if self.tweet_data_df['a'][row_idx] > 0.0:
                new_data['x'].append(self.tweet_data_df['x'][row_idx])
                new_data['y'].append(self.tweet_data_df['y'][row_idx])
                new_data['width'].append(self.tweet_data_df['a'][row_idx] * 2.0)
                new_data['height'].append(self.tweet_data_df['b'][row_idx] * 2.0)
                new_data['angle'].append(self.tweet_data_df['angle'][row_idx])
        return new_data

    def clear_dissolve(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        return new_data

    def update_dissolve(self, id):
        #print("Update Dissolve: " + str(id))
        dissolve_data = self.sibling_data_manager.get_dissolve(id)
        new_data = dict()
        new_data['x'] = dissolve_data['x']
        new_data['y'] = dissolve_data['y']
        return new_data

    def update_find_circle(self, find_circle_idx):
        new_data = dict()
        new_data['x'] = [self.tweet_data_df['x'][find_circle_idx]]
        new_data['y'] = [self.tweet_data_df['y'][find_circle_idx]]
        return new_data


class FilterSettings:

    def __init__(self, config):
        self.count_start = config.count[0]
        self.count_end = config.count[1]
        self.count_active = True
        self.area_start = config.area[0]
        self.area_end = config.area[1]
        self.area_active = True
        self.distance_start = config.distance[0]
        self.distance_end = config.distance[1]
        self.distance_active = True
        self.ratio_start = config.ratio[0]
        self.ratio_end = config.ratio[1]
        self.ratio_max_unbounded = config.ratio[1] - 1.0
        self.ratio_active = True
        self.dissolve_start = config.dissolve[0]
        self.dissolve_end = config.dissolve[1]
        self.dissolve_active = True

class TweetDataController:

    def __init__(self, pre_processor, config, user_info):
        self.all = MapTweetData(pre_processor.tweet_data_all, config)
        self.working = MapTweetData(pre_processor.tweet_data_working, config)
        self.non_working = MapTweetData(pre_processor.tweet_data_non_working, config)
        self.active_dataset = self.working
        self.blend_dataset = self.non_working

        self.config = config
        self.user_info = user_info
        self.selection_details = None

        self.hover_idx = ColumnDataSource(data=dict(id=[], idx=[]))
        data_hover_idx = {
            'id': [-1],
            'idx': [-1]
        }
        df_hover_idx = pd.DataFrame(data_hover_idx)
        self.hover_idx.data = self.hover_idx.from_df(df_hover_idx)

        self.circles = ColumnDataSource(data=dict(x=[], y=[], id=[], distance=[]))
        self.circles.data = self.circles.from_df(self.active_dataset.tweet_data_df)

        self.lod_dummy = ColumnDataSource(data=dict(x=[], y=[], lod=[]))
        self.lod_dummy.on_change('data', self.lod_dummy_changed)

        self.selected_circle = ColumnDataSource(data=dict(x=[], y=[], id=[]))
        self.selected_circle.on_change('data', self.selected_circle_changed)

        self.sde_ellipse = ColumnDataSource(data=dict(x=[], y=[], width=[], height=[], angle=[]))
        self.siblings = ColumnDataSource(data=dict(x=[], y=[]))
        self.sibling_ellipses = ColumnDataSource(data=dict(x=[], y=[], width=[], height=[], angle=[]))
        self.patch_dissolve = ColumnDataSource(data=dict(x=[], y=[]))

        self.selected_circle_blend = ColumnDataSource(data=dict(x=[], y=[], id=[]))
        self.sde_ellipse_blend = ColumnDataSource(data=dict(x=[], y=[], width=[], height=[], angle=[], alpha=[]))
        self.siblings_blend = ColumnDataSource(data=dict(x=[], y=[]))
        self.sibling_ellipses_blend = ColumnDataSource(data=dict(x=[], y=[], width=[], height=[], angle=[]))
        self.patch_dissolve_blend = ColumnDataSource(data=dict(x=[], y=[]))

        self.find_circle = ColumnDataSource(data=dict(id=[], x=[], y=[]))

        self.filter_settings = FilterSettings(self.config)

        self.circle_id = -1
        self.circle_idx = -1
        self.find_circle_idx = -1

        self.blend_active = False
        self.circles_active = True

        self.csr = None
        self.csbr = None
        self.er = None
        self.ebr = None
        self.ser = None
        self.sebr = None
        self.pdr = None
        self.pdbr = None
        self.sr = None
        self.sbr = None
        self.hr_count = None
        self.hr_area = None
        self.hr_distance = None
        self.hr_ratio = None
        self.hr_dissolve = None

        self.histogram_controller_count = HistogramController(self.active_dataset.histogram_count.df_histogram_cds)
        self.histogram_controller_area = HistogramController(self.active_dataset.histogram_area.df_histogram_cds)
        self.histogram_controller_distance = HistogramController(self.active_dataset.histogram_distance.df_histogram_cds)
        self.histogram_controller_ratio = HistogramController(self.active_dataset.histogram_ratio.df_histogram_cds)
        self.histogram_controller_dissolve = HistogramController(self.active_dataset.histogram_dissolve.df_histogram_cds)

    def find_id(self, id_value):
        print("Find: ID: " + str(id_value))
        # This is a dataframe of the selected (single) row.
        id_df = self.active_dataset.tweet_data_df.loc[self.active_dataset.tweet_data_df['id'] == id_value]
        if id_df.empty:
            print("No row with ID: " + str(id_value))
            self.clear_find_circle()
        else:
            idx = id_df.index
            # We can assume that there will be only one index value in this list, as IDs are unique.
            self.find_circle_idx = idx[0]
            self.update_find_circle()

    def reset_renderer_properties(self):
        self.csr.glyph.fill_alpha = 1.0
        self.er.glyph.line_alpha = 0.5
        self.er.glyph.fill_alpha = 0.5
        self.sr.glyph.fill_alpha = 0.8
        self.ser.glyph.line_alpha = 1.0
        self.pdr.glyph.line_alpha = 0.4
        self.pdr.glyph.fill_alpha = 0.4

    def clear_circles(self):
        print("Clear Circles:")
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['id'] = []
        new_data['distance'] = []
        new_data['color'] = []
        self.circles.data = new_data

    def toggle_data(self, toggle_on):
        if toggle_on:
            print("Toggle Data: On")
            self.clear_circles()
            self.circles_active = False
        else:
            print("Toggle Data: Off")
            self.circles_active = True
            self.apply_filters()

    def lod_dummy_changed(self, attrname, old, new):
        #print("LOD Dummy Change: " + str(new))
        if len(new['lod']) > 0:
            lod_value = new['lod'][0]

            if lod_value is 0:
                print("LOD Dummy Callback: Start")
                self.clear_circles()

            if lod_value is 1:
                print("LOD Dummy Callback: End")
                self.apply_filters()

    def selected_circle_changed(self, attrname, old, new):
        #print("Selected Circle Change: " + str(new))
        if len(new['id']) > 0:
            self.circle_id = new['id'][0]
            print("Selected Circle Change: id: " + str(self.circle_id))

            self.working.sibling_data_manager.get_sibling_data(self.circle_id)
            self.non_working.sibling_data_manager.get_sibling_data(self.circle_id)

            self.update_selection_details()

            if self.blend_dataset is not None:
                self.selected_circle_blend.data = self.blend_dataset.clear_selected_circle()

            self.clear_sde_ellipse()
            self.clear_siblings()
            self.clear_sibling_ellipses()
            self.clear_dissolve()
            self.reset_renderer_properties()
            self.blend_active = False

    def selected_count_changed(self, attrname, old, new):
        if len(new['idx']) > 0:
            print("Selected Histogram Count: idx: " + str(new['idx']))

    def selected_area_changed(self, attrname, old, new):
        if len(new['idx']) > 0:
            print("Selected Histogram Area: idx: " + str(new['idx']))

    def update_selection_details(self):
        if self.circle_id > -1:
            id_df = self.active_dataset.tweet_data_df.loc[self.active_dataset.tweet_data_df['id'] == self.circle_id]
            self.circle_idx = id_df.index[0]
            print("Selected Circle Change: " + str(self.circle_id) + " : " + str(self.circle_idx))

            area_working = self.working.tweet_data_df['area'][self.circle_idx]
            count_working = self.working.tweet_data_df['count'][self.circle_idx]
            ratio_working = self.working.tweet_data_df['ratio'][self.circle_idx]
            area_non_working = self.non_working.tweet_data_df['area'][self.circle_idx]
            count_non_working = self.non_working.tweet_data_df['count'][self.circle_idx]
            ratio_non_working = self.non_working.tweet_data_df['ratio'][self.circle_idx]
            distance = self.working.tweet_data_df['distance'][self.circle_idx]

            username, profile_text = self.user_info.find_user_profile(int(self.circle_id))
            details_str = "<b>Selected ID</b>: " + str(self.circle_id) + "<br/>"
            details_str += "<b>Username</b>: " + str(username) + "<br/>"
            details_str += "<b>Profile</b>: " + str(profile_text) + "<br/>"
            details_str += "<b>Working</b>:<br/>"
            details_str += "<b>&nbsp;&nbsp;Area</b>: " + str(area_working) + "<br/>"
            details_str += "<b>&nbsp;&nbsp;Count</b>: " + str(count_working) + "<br/>"
            details_str += "<b>&nbsp;&nbsp;Ratio</b>: " + str(ratio_working) + "<br/>"
            details_str += "<b>Non Working</b>:<br/>"
            details_str += "<b>&nbsp;&nbsp;Area</b>: " + str(area_non_working) + "<br/>"
            details_str += "<b>&nbsp;&nbsp;Count</b>: " + str(count_non_working) + "<br/>"
            details_str += "<b>&nbsp;&nbsp;Ratio</b>: " + str(ratio_non_working) + "<br/>"
            details_str += "<b>Distance</b>: " + str(distance)

            self.selection_details.text = details_str

    def clear_selected_circle(self):
        self.selected_circle.data = self.active_dataset.clear_selected_circle()
        if self.blend:
            if self.blend_dataset is not None:
                self.selected_circle_blend.data = self.blend_dataset.clear_selected_circle()

    def update_selected_circle(self):
        if self.circle_id > -1:
            print("Updating Selected Circle: " + str(self.circle_id))
            self.selected_circle.data, self.circle_idx = self.active_dataset.update_selected_circle(self.circle_id)
            if self.blend:
                if self.blend_dataset is not None:
                    self.selected_circle_blend.data, dummy_circle_idx = self.blend_dataset.update_selected_circle(self.circle_id)

    def clear_sde_ellipse(self):
        self.sde_ellipse.data = self.active_dataset.clear_sde_ellipse()
        if self.blend:
            if self.blend_dataset is not None:
                self.sde_ellipse_blend.data = self.blend_dataset.clear_sde_ellipse()

    def update_sde_ellipse(self):
        print("Update SDE Ellipse: " + str(self.circle_idx))
        if self.circle_idx > -1:
            self.sde_ellipse.data = self.active_dataset.update_sde_ellipse(self.circle_idx)
            if self.blend_active:
                if self.blend_dataset is not None:
                    self.sde_ellipse_blend.data = self.blend_dataset.update_sde_ellipse(self.circle_idx)

    def clear_siblings(self):
        self.siblings.data = self.active_dataset.clear_siblings()
        if self.blend_dataset is not None:
            self.siblings_blend.data = self.blend_dataset.clear_siblings()

    def update_siblings(self):
        print("Update Siblings: " + str(self.circle_id) + " : " + str(self.circle_idx))
        if self.circle_idx > -1:
            self.siblings.data = self.active_dataset.update_siblings(self.circle_id)
            if self.blend_active:
                if self.blend_dataset is not None:
                    self.siblings_blend.data = self.blend_dataset.update_siblings(self.circle_id)

    def clear_sibling_ellipses(self):
        self.sibling_ellipses.data = self.active_dataset.clear_sibling_ellipses()
        if self.blend_dataset is not None:
            self.sibling_ellipses_blend.data = self.blend_dataset.clear_sibling_ellipses()

    def update_sibling_ellipses(self):
        print("Update Sibling Ellipses: " + str(self.circle_id) + " : " + str(self.circle_idx))
        if self.circle_idx > -1:
            self.sibling_ellipses.data = self.active_dataset.update_sibling_ellipses(self.circle_id)
            if self.blend_active:
                if self.blend_dataset is not None:
                    self.sibling_ellipses_blend.data = self.blend_dataset.update_sibling_ellipses(self.circle_id)

    def clear_dissolve(self):
        self.patch_dissolve.data = self.active_dataset.clear_dissolve()
        if self.blend_dataset is not None:
            self.patch_dissolve_blend.data = self.blend_dataset.clear_dissolve()

    def update_dissolve(self):
        print("Update Dissolve: " + str(self.circle_id) + " : " + str(self.circle_idx))
        if self.circle_idx > -1:
            self.patch_dissolve.data = self.active_dataset.update_dissolve(self.circle_id)
            if self.blend_active:
                if self.blend_dataset is not None:
                    self.patch_dissolve_blend.data = self.blend_dataset.update_dissolve(self.circle_id)

    def clear_find_circle(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        self.find_circle.data = new_data

    def update_find_circle(self):
        if self.find_circle_idx > -1:
            self.find_circle.data = self.active_dataset.update_find_circle(self.find_circle_idx)

    def clear_all(self):
        self.clear_selected_circle()
        self.clear_sde_ellipse()
        self.clear_siblings()
        self.clear_sibling_ellipses()
        self.clear_dissolve()
        self.clear_find_circle()

    def update_histograms(self):
        print(self.active_dataset.name)
        self.hr_count.data_source.data = self.histogram_controller_count.cds.from_df(self.active_dataset.histogram_count.df_histogram_cds)
        self.hr_area.data_source.data = self.histogram_controller_area.cds.from_df(self.active_dataset.histogram_area.df_histogram_cds)
        self.hr_distance.data_source.data = self.histogram_controller_distance.cds.from_df(self.active_dataset.histogram_distance.df_histogram_cds)
        self.hr_ratio.data_source.data = self.histogram_controller_ratio.cds.from_df(self.active_dataset.histogram_ratio.df_histogram_cds)

    def switch_tweet_dataset(self, value):
        if value == 0:
            print("Switching to 'mean all'.")
            self.active_dataset = self.all
            self.blend_dataset = None
        elif value == 1:
            print("Switching to 'median working'.")
            self.active_dataset = self.working
            self.blend_dataset = self.non_working
        else:
            print("Switching to 'median non_working'.")
            self.active_dataset = self.non_working
            self.blend_dataset = self.working

        self.apply_filters()
        self.update_selected_circle()
        self.update_selection_details()
        self.update_histograms()

    def num_of_points_active(self):
        active = len(self.circles.data['id'])
        total = len(self.active_dataset.tweet_data_df['id'])
        return active, total

    def filter_circles_by_count(self, count_start, count_end):
        self.filter_settings.count_start = count_start
        self.filter_settings.count_end = count_end
        self.apply_filters()

    def filter_circles_by_area(self, area_start, area_end):
        self.filter_settings.area_start = area_start
        self.filter_settings.area_end = area_end
        self.apply_filters()

    def filter_circles_by_distance(self, distance_start, distance_end):
        self.filter_settings.distance_start = distance_start
        self.filter_settings.distance_end = distance_end
        self.apply_filters()

    def filter_circles_by_ratio(self, ratio_start, ratio_end):
        self.filter_settings.ratio_start = ratio_start
        self.filter_settings.ratio_end = ratio_end
        self.apply_filters()

    def filter_circles_by_dissolve(self, dissolve_start, dissolve_end):
        self.filter_settings.dissolve_start = dissolve_start
        self.filter_settings.dissolve_end = dissolve_end
        self.apply_filters()

    def filters_active(self, active_list):
        self.filter_settings.count_active = False
        self.filter_settings.area_active = False
        self.filter_settings.distance_active = False
        self.filter_settings.ratio_active = False
        self.filter_settings.dissolve_active = False

        for filter_check_value in active_list:
            if filter_check_value is 0:
                self.filter_settings.count_active = True

            if filter_check_value is 1:
                self.filter_settings.area_active = True

            if filter_check_value is 2:
                self.filter_settings.distance_active = True

            if filter_check_value is 3:
                self.filter_settings.ratio_active = True

            if filter_check_value is 4:
                self.filter_settings.dissolve_active = True

        self.apply_filters()

    def apply_filters(self):
        if self.circles_active:
            subset_df = self.active_dataset.tweet_data_df
            if self.filter_settings.count_active:
                subset_df = subset_df.loc[
                            (subset_df['count'] >= self.filter_settings.count_start)
                        &   (subset_df['count'] <= self.filter_settings.count_end)]

            if self.filter_settings.area_active:
                subset_df = subset_df.loc[
                            (subset_df['area'] >= self.filter_settings.area_start)
                        &   (subset_df['area'] <= self.filter_settings.area_end)]

            if self.filter_settings.distance_active:
                subset_df = subset_df.loc[
                            (subset_df['distance'] >= self.filter_settings.distance_start)
                        &   (subset_df['distance'] <= self.filter_settings.distance_end)]

            if self.filter_settings.ratio_active:
                if self.filter_settings.ratio_end > self.filter_settings.ratio_max_unbounded:
                    subset_df = subset_df.loc[(subset_df['ratio'] >= self.filter_settings.ratio_start)]
                else:
                    subset_df = subset_df.loc[
                                (subset_df['ratio'] >= self.filter_settings.ratio_start)
                            &   (subset_df['ratio'] <= self.filter_settings.ratio_end)]

            if self.filter_settings.dissolve_active:
                subset_df = subset_df.loc[
                            (subset_df['dissolve_area'] >= self.filter_settings.dissolve_start)
                        &   (subset_df['dissolve_area'] <= self.filter_settings.dissolve_end)]

            self.circles.data = self.circles.from_df(subset_df)

    def turn_blend_on(self, ellipse_active, siblings_active, dissolve_active):
        print("Turn Blend On:")
        self.blend_active = True
        if self.circle_idx > -1:
            if self.blend_dataset is not None:
                self.selected_circle_blend.data, dummy_circle_idx = self.blend_dataset.update_selected_circle(self.circle_id)
                if ellipse_active:
                    self.sde_ellipse_blend.data = self.blend_dataset.update_sde_ellipse(self.circle_idx)
                    self.siblings_blend.data = self.blend_dataset.update_siblings(self.circle_id)

                if siblings_active:
                    self.sibling_ellipses_blend.data = self.blend_dataset.update_sibling_ellipses(self.circle_id)

                if dissolve_active:
                    self.patch_dissolve_blend.data = self.blend_dataset.update_dissolve(self.circle_id)

    def turn_blend_off(self):
        print("Turn Blend Off:")
        self.blend_active = False
        if self.blend_dataset is not None:
            self.selected_circle_blend.data = self.blend_dataset.clear_selected_circle()
            self.sde_ellipse_blend.data = self.blend_dataset.clear_sde_ellipse()
            self.siblings_blend.data = self.blend_dataset.clear_siblings()
            self.sibling_ellipses_blend.data = self.blend_dataset.clear_sibling_ellipses()
            self.patch_dissolve_blend.data = self.blend_dataset.clear_dissolve()
            self.reset_renderer_properties()

    def blend(self, blend_ratio, ellipse_active, siblings_active, dissolve_active):
        #print("Blend: "+ str(blend_ratio))
        self.csr.glyph.fill_alpha = blend_ratio
        self.er.glyph.line_alpha = blend_ratio
        self.er.glyph.fill_alpha = blend_ratio
        self.sr.glyph.line_alpha = blend_ratio
        self.sr.glyph.fill_alpha = blend_ratio
        self.ser.glyph.line_alpha = blend_ratio
        self.pdr.glyph.line_alpha = blend_ratio
        self.pdr.glyph.fill_alpha = blend_ratio

        if self.blend_dataset is not None:
            self.csbr.glyph.fill_alpha = 1.0 - blend_ratio

            if ellipse_active:
                self.ebr.glyph.line_alpha = 1.0 - blend_ratio
                self.ebr.glyph.fill_alpha = 1.0 - blend_ratio
                self.sbr.glyph.line_alpha = 1.0 - blend_ratio
                self.sbr.glyph.fill_alpha = 1.0 - blend_ratio

            if siblings_active:
                self.sebr.glyph.line_alpha = 1.0 - blend_ratio

            if dissolve_active:
                self.pdbr.glyph.line_alpha = 1.0 - blend_ratio
                self.pdbr.glyph.fill_alpha = 1.0 - blend_ratio

class SiblingDataManager:

    def __init__(self, sibling_data_folder, dataset_name):
        self.sibling_data_folder = sibling_data_folder
        self.dataset_name = dataset_name

        self.sibling_data = dict()

    def get_sibling_data(self, id):
        filename = str(self.dataset_name) + "_" + str(id) + ".json"
        logger.info("Get Sibling Data: " + filename)
        if id not in self.sibling_data:
            self.sibling_data[id] = self.load_sibling_data(filename)

        return self.sibling_data[id]

    def load_sibling_data(self, filename):
        logger.info("Load Sibling Data: " + filename)
        file_open = FileOpen(self.sibling_data_folder, filename)
        logger.info(file_open)

        data = dict()
        with open(file_open.absolute, 'r') as sibling_file:
            sibling_json = json.load(sibling_file)

        logger.debug(data)
        return sibling_json

    def get_sibling_idxs(self, id):
        sibling_data = self.get_sibling_data(id)
        return sibling_data['idxs']

    def get_dissolve(self, id):
        sibling_data = self.get_sibling_data(id)
        return sibling_data['dissolve']

    def __str__(self):
        manager_str = "Sibling Data Manager: "
        manager_str += "\nFolder: " + str(self.sibling_data_folder)
        manager_str += "\nDataset Name: " + str(self.dataset_name)
        manager_str += "\nNumber of Sibling Sets: " + str(len(self.sibling_data))
        return manager_str


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
import os

def main():
    logging.basicConfig(
        filename='logs/tweet_data.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s')

    logger.info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    logger.info("Tweet Data:")

    cwd = os.getcwd()
    sibling_data_folder = cwd + "/sibling_data/"
    sibling_data_manager = SiblingDataManager(sibling_data_folder, "working")
    logger.info(sibling_data_manager)

    sibling_data = sibling_data_manager.get_sibling_data(379005817)
    logger.info(sibling_data)
    logger.info(sibling_data_manager)

    sibling_data = sibling_data_manager.get_sibling_data(464331848)
    logger.info(sibling_data)
    logger.info(sibling_data_manager)

    sibling_data = sibling_data_manager.get_sibling_data(379005817)
    logger.info(sibling_data)
    logger.info(sibling_data_manager)

if __name__ == '__main__':
    main()