import logging.config
import pandas as pd

from bokeh.models import ColumnDataSource

from analysis_utilities import *
from shapely_utilities import *

logger = logging.getLogger()

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
                         angle_column_name):

        self.df = pd.DataFrame(columns=['id', 'lat', 'lon', 'x', 'y', 'area', 'a', 'b', 'angle'])

        # Keep an eye open for: A value is trying to be set on a copy of a slice from a DataFrame
        # Caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
        # Recommends using df.loc[row_idx, col_idx]
        self.df['id'] = data_df[id_column_name].copy(True)
        self.df['lat'] = data_df[latitude_column_name].copy(True)
        self.df['lon'] = data_df[longitude_column_name].copy(True)
        self.df['area'] = data_df[area_column_name].copy(True)
        self.df['angle'] = data_df[angle_column_name].copy(True)

        for row_idx in range(0, data_df.shape[0]):
            # If the area is nan, then set it to -1 to indicate it is undefined.
            if math.isnan(self.df['area'][row_idx]):
                self.df.loc[row_idx, 'area'] = -1.0

            # If any angles are nan, then clean to 0.0
            #if (math.isnan(self.df['angle'][row_idx])):
            #    self.df.loc[row_idx, 'angle'] = 0.0

            self.df.loc[row_idx, 'x'], self.df.loc[row_idx, 'y'] = lon_lat_to_east_north(
                self.df['lon'][row_idx],
                self.df['lat'][row_idx])

            self.df.loc[row_idx, 'a'], self.df.loc[row_idx, 'b'] = self.calculate_ellipse_properties(
                self.df['id'][row_idx],
                row_idx,
                self.df['area'][row_idx],
                data_df[x_y_ratio][row_idx])

        self.max_a = self.df['a'].max()
        self.max_b = self.df['b'].max()
        self.max_area = self.df['area'].max()

        self.find_datapoints_enclosed_with_ellipse()
        self.calculate_colour()

        logger.info("\n%s", self.df.head())

    def find_datapoints_enclosed_with_ellipse(self):
        row_idx = 0
        counts = []
        grouped_ids_list = []
        grouped_idxs_list = []
        for row in self.df.itertuples():
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

            row_idx += 1

            count = len(grouped_ids)
            counts.append(count)
            grouped_ids_list.append(grouped_ids)
            grouped_idxs_list.append(grouped_idxs)

            if count > self.max_count:
                self.max_count = count

        self.df['count'] = counts
        self.df['ids'] = grouped_ids_list
        self.df['idxs'] = grouped_idxs_list

        logger.info("\n%s", self.df.head())

    def calculate_colour(self):
        #max_count = self.max_count
        count_max = 600.0

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

        colours = []
        for row in self.df.itertuples():
            count = row.count
            scale_factor = 1.0 - 1.0/(1.0 + count)

            point_r = hex(int(start_r + scale_factor * interval_r)).split('x')[-1]
            point_g = hex(int(start_g + scale_factor * interval_g)).split('x')[-1]
            point_b = hex(int(start_b + scale_factor * interval_b)).split('x')[-1]
            colour_hex_str = "#" + str(point_r) + str(point_g) + str(point_b)
            colours.append(colour_hex_str)

        self.df['color'] = colours

    # Using the x/y ratio and area calculate the major and minor axes for the ellipse.
    # If the area is 0.0 or -1, then we can not cacluate these values
    def calculate_ellipse_properties(self, row_id, idx, area, xy_ratio):
        # Remember: Viewing the data in excel or notepad++, the first row is the headers, and the row count starts from 1
        # So, where the pandas index starts from 0, the starting row in excel or notepad++ is actually 2.
        if area == 0.0 or area == -1:
            logger.warning("id: %s (%s): The ellipse area is %s - can not calculate 'a' or 'b'.", row_id, idx, area)
            return -1, -1
        else:
            if xy_ratio == "Infinity":
                logger.warning("id %s (%s): The ellipse x/y ratio is 'Infinity' - can not calculate 'a' or 'b'.", row_id, idx)
                return -1, -1

            xy_ratio_f = float(xy_ratio)
            if math.isnan(xy_ratio_f):
                logger.warning("id %s (%s): The ellipse x/y ratio is 'nan' - can not calculate 'a' or 'b'.", row_id, idx)
                return -1, -1

            if math.isinf(xy_ratio_f):
                logger.warning("id %s (idx): The ellipse x/y ratio is 'inf' - can not calculate 'a' or 'b'.", row_id, idx)
                return -1, -1

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

    def __init__(self, tweet_data_df):
        self.tweet_data_df = tweet_data_df

    def update_selected_circle(self, circle_id):
        id_df = self.tweet_data_df.loc[self.tweet_data_df['id'] == circle_id]
        idx = id_df.index
        circle_idx = idx[0]

        new_data = dict()
        new_data['x'] = [self.tweet_data_df['x'][idx[0]]]
        new_data['y'] = [self.tweet_data_df['y'][idx[0]]]
        new_data['id'] = [circle_id]

        return new_data, circle_idx

    def update_sde_ellipse(self, circle_idx):
        print("Update SDE Ellipse: " + str(circle_idx))

        new_data = dict()
        new_data['x'] = [self.tweet_data_df['x'][circle_idx]]
        new_data['y'] = [self.tweet_data_df['y'][circle_idx]]
        new_data['width'] = [self.tweet_data_df['a'][circle_idx] * 2.0]
        new_data['height'] = [self.tweet_data_df['b'][circle_idx] * 2.0]
        new_data['angle'] = [self.tweet_data_df['angle'][circle_idx]]

        return new_data

    def update_siblings(self, circle_idx):
        print("Update Siblings: " + str(circle_idx))
        print(str(circle_idx) + " : " + str(self.tweet_data_df['idxs'][circle_idx]))
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        for row_idx in self.tweet_data_df['idxs'][circle_idx]:
            new_data['x'].append(self.tweet_data_df['x'][row_idx])
            new_data['y'].append(self.tweet_data_df['y'][row_idx])

        return new_data

    def update_sibling_ellipses(self, circle_idx):
        print("Update Sibling Ellipses: " + str(circle_idx))

        print(self.tweet_data_df['idxs'][circle_idx])
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['width'] = []
        new_data['height'] = []
        new_data['angle'] = []
        for row_idx in self.tweet_data_df['idxs'][circle_idx]:
            new_data['x'].append(self.tweet_data_df['x'][row_idx])
            new_data['y'].append(self.tweet_data_df['y'][row_idx])
            new_data['width'].append(self.tweet_data_df['a'][row_idx] * 2.0)
            new_data['height'].append(self.tweet_data_df['b'][row_idx] * 2.0)
            new_data['angle'].append(self.tweet_data_df['angle'][row_idx])

        return new_data

    def update_dissolve(self, circle_idx):
        print("Update Dissolve: " + str(circle_idx))

        print(self.tweet_data_df['idxs'][circle_idx])
        ellipses = dict()
        ellipses['x'] = []
        ellipses['y'] = []
        ellipses['a'] = []
        ellipses['b'] = []
        ellipses['angle'] = []
        for row_idx in self.tweet_data_df['idxs'][circle_idx]:
            ellipses['x'].append(self.tweet_data_df['x'][row_idx])
            ellipses['y'].append(self.tweet_data_df['y'][row_idx])
            ellipses['a'].append(self.tweet_data_df['a'][row_idx])
            ellipses['b'].append(self.tweet_data_df['b'][row_idx])
            ellipses['angle'].append(self.tweet_data_df['angle'][row_idx])

        # Need to remember to include the datapoints own ellipse details. Otherwise the following error occurs:
        # AttributeError("'MultiPolygon' object has no attribute 'exterior'")
        # Probably caused by the sibling ellipses NOT overlapping i.e. separate.
        ellipses['x'].append(self.tweet_data_df['x'][circle_idx])
        ellipses['y'].append(self.tweet_data_df['y'][circle_idx])
        ellipses['a'].append(self.tweet_data_df['a'][circle_idx])
        ellipses['b'].append(self.tweet_data_df['b'][circle_idx])
        ellipses['angle'].append(self.tweet_data_df['angle'][circle_idx])

        new_data = dict()
        new_data['x'], new_data['y'] = dissolve_ellipses(ellipses)

        return new_data

    def update_find_circle(self, find_circle_idx):
        new_data = dict()
        new_data['x'] = [self.tweet_data_df['x'][find_circle_idx]]
        new_data['y'] = [self.tweet_data_df['y'][find_circle_idx]]

        return new_data

class FilterSettings:

    def __init__(self, count_start, count_end, area_start, area_end):
        self.count_start = count_start
        self.count_end = count_end
        self.count_active = True
        self.area_start = area_start
        self.area_end = area_end
        self.area_active = True


class TweetDataController:

    def __init__(self, pre_processor, user_info):
        self.all = MapTweetData(pre_processor.tweet_data_all.df)
        self.working = MapTweetData(pre_processor.tweet_data_working.df)
        self.non_working = MapTweetData(pre_processor.tweet_data_non_working.df)
        self.active_dataset = self.all

        self.user_info = user_info
        self.selection_details = None
        #self.text_username = None
        #self.text_profile = None

        self.hover_idx = ColumnDataSource(data=dict(id=[], idx=[]))
        data_hover_idx = {
            'id': [-1],
            'idx': [-1]
        }
        df_hover_idx = pd.DataFrame(data_hover_idx)
        self.hover_idx.data = self.hover_idx.from_df(df_hover_idx)

        self.circles = ColumnDataSource(data=dict(x=[], y=[], id=[]))
        self.circles.data = self.circles.from_df(self.active_dataset.tweet_data_df)

        self.selected_circle = ColumnDataSource(data=dict(x=[], y=[], id=[]))
        sc_data = {
            'x': [],
            'y': [],
            'id': []
        }
        df_selected_circle = pd.DataFrame(sc_data)
        self.selected_circle.data = self.selected_circle.from_df(df_selected_circle)
        self.selected_circle.on_change('data', self.selected_circle_changed)

        self.sde_ellipse = ColumnDataSource(data=dict(x=[], y=[], width=[], height=[], angle=[]))
        sde_e_data = {
            'x': [],
            'y': [],
            'width': [],
            'height': [],
            'angle': [],
        }
        df_sde_ellipse = pd.DataFrame(sde_e_data)
        self.sde_ellipse.data = self.sde_ellipse.from_df(df_sde_ellipse)

        self.siblings = ColumnDataSource(data=dict(x=[], y=[]))
        s_data = {
            'x': [],
            'y': []
        }
        df_siblings = pd.DataFrame(s_data)
        self.siblings.data = self.siblings.from_df(df_siblings)

        self.sibling_ellipses = ColumnDataSource(data=dict(x=[], y=[], width=[], height=[], angle=[]))
        se_data = {
            'x': [],
            'y': [],
            'width': [],
            'height': [],
            'angle': [],
        }
        df_sibling_ellipses = pd.DataFrame(se_data)
        self.sibling_ellipses.data = self.sibling_ellipses.from_df(df_sibling_ellipses)

        self.patch_dissolve = ColumnDataSource(data=dict(x=[], y=[]))
        patch_data = {
            'x': [],
            'y': []
        }
        df_patch = pd.DataFrame(patch_data)
        self.patch_dissolve.data = self.patch_dissolve.from_df(df_patch)

        self.find_circle = ColumnDataSource(data=dict(id=[], x=[], y=[]))

        self.filter_settings = FilterSettings(0, 600, 0, 131000)

        self.circle_id = -1
        self.circle_idx = -1
        self.find_circle_idx = -1

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

    def selected_circle_changed(self, attrname, old, new):
        #print("Selected Circle Change: " + str(new))
        if len(new['id']) > 0:
            self.circle_id = new['id'][0]
            print("Selected Circle Change: id: " + str(self.circle_id))
            self.update_selection_details()

            self.clear_sde_ellipse()
            self.clear_siblings()
            self.clear_sibling_ellipses()
            self.clear_dissolve()

    def update_selection_details(self):
        id_df = self.active_dataset.tweet_data_df.loc[self.active_dataset.tweet_data_df['id'] == self.circle_id]
        print(id_df)
        self.circle_idx = id_df.index[0]

        area = self.active_dataset.tweet_data_df['area'][self.circle_idx]
        count = self.active_dataset.tweet_data_df['count'][self.circle_idx]
        print("Selected Circle Change: idx: " + str(self.circle_idx))

        username, profile_text = self.user_info.find_user_profile(int(self.circle_id))
        details_str = "<b>Selected ID</b>: " + str(self.circle_id) + "<br/>"
        details_str += "<b>Username</b>: " + str(username) + "<br/>"
        details_str += "<b>Profile</b>: " + str(profile_text) + "<br/>"
        details_str += "<b>Area</b>: " + str(area) + "<br/>"
        details_str += "<b>Count</b>: " + str(count)

        self.selection_details.text = details_str

    def clear_selected_circle(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['id'] = []
        self.selected_circle.data = new_data

    def update_selected_circle(self):
        if self.circle_id > -1:
            print("Updating Selected Circle: " + str(self.circle_id))
            self.selected_circle.data, self.circle_idx = self.active_dataset.update_selected_circle(self.circle_id)

    def clear_sde_ellipse(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['width'] = []
        new_data['height'] = []
        new_data['angle'] = []
        self.sde_ellipse.data = new_data

    def update_sde_ellipse(self):
        print("Update SDE Ellipse: " + str(self.circle_idx))
        if self.circle_idx > -1:
            self.sde_ellipse.data = self.active_dataset.update_sde_ellipse(self.circle_idx)

    def clear_siblings(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        self.siblings.data = new_data

    def update_siblings(self):
        print("Update Siblings: " + str(self.circle_idx))
        if self.circle_idx > -1:
            self.siblings.data = self.active_dataset.update_siblings(self.circle_idx)

    def clear_sibling_ellipses(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        new_data['width'] = []
        new_data['height'] = []
        new_data['angle'] = []
        self.sibling_ellipses.data = new_data

    def update_sibling_ellipses(self):
        print("Update Sibling Ellipses: " + str(self.circle_idx))
        if self.circle_idx > -1:
            self.sibling_ellipses.data = self.active_dataset.update_sibling_ellipses(self.circle_idx)

    def clear_dissolve(self):
        new_data = dict()
        new_data['x'] = []
        new_data['y'] = []
        self.patch_dissolve.data = new_data

    def update_dissolve(self):
        print("Update Dissolve: " + str(self.circle_idx))
        if self.circle_idx > -1:
            self.patch_dissolve.data = self.active_dataset.update_dissolve(self.circle_idx)

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

    def switch_tweet_dataset(self, value):
        if value == 0:
            self.active_dataset = self.all
            print("Switching to 'mean all'.")
        elif value == 1:
            self.active_dataset = self.working
            print("Switching to 'median working'.")
        else:
            self.active_dataset = self.non_working
            print("Switching to 'median non-working'.")

        self.circles.data = self.circles.from_df(self.active_dataset.tweet_data_df)

        self.update_selected_circle()
        self.update_selection_details()

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

    def filters_active(self, active_list):
        self.filter_settings.count_active = False
        self.filter_settings.area_active = False
        for filter_check_value in active_list:
            if filter_check_value is 0:
                self.filter_settings.count_active = True

            if filter_check_value is 1:
                self.filter_settings.area_active = True

        self.apply_filters()

    def apply_filters(self):
        subset_df = self.active_dataset.tweet_data_df
        if self.filter_settings.count_active:
            subset_df = subset_df.loc[
                        (subset_df['count'] >= self.filter_settings.count_start)
                    &   (subset_df['count'] <= self.filter_settings.count_end)]

        if self.filter_settings.area_active:
            subset_df = subset_df.loc[
                        (subset_df['area'] >= self.filter_settings.area_start)
                    &   (subset_df['area'] <= self.filter_settings.area_end)]

        self.circles.data = self.circles.from_df(subset_df)

    def turn_blend_on(self):
        print("Turn Blend On:")

    def turn_blend_off(self):
        print("Turn Blend Off:")

    def blend(self, blend_ratio):
        print("Blend: "+ str(blend_ratio))
