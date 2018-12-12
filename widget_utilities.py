
from bokeh.models import RangeSlider
from bokeh.models.widgets import Button, CheckboxGroup, Paragraph, RadioButtonGroup, Slider, TextInput, Toggle

class MapWidgets:

    def __init__(self, map_tweet_data, user_info):
        self.map_tweet_data = map_tweet_data
        self.user_info = user_info

        self.toggle_sde_ellipse = Toggle(label="Toggle Ellipse", active=False, width=150)
        self.toggle_sde_ellipse.on_click(self.toggle_sde_ellipse_callback)

        self.toggle_sibling_ellipses = Toggle(label="Toggle Siblings", active=False, width=150)
        self.toggle_sibling_ellipses.on_click(self.toggle_sibling_ellipses_callback)

        self.toggle_dissolve = Toggle(label="Toggle Dissolve", active=False, width=150)
        self.toggle_dissolve.on_click(self.toggle_dissolve_callback)

        self.toggle_user_info = Toggle(label="Toggle User Info", active=False, width=150)
        self.toggle_user_info.on_click(self.toggle_user_info_callback)

        self.text_id = Paragraph(text='Selected ID: -1 (-1)')
        self.text_username = Paragraph(text='Username: ')
        self.text_profile = Paragraph(text='Profile: ')

        self.text_input = TextInput(value="-1", title="Enter ID:", width=150)

        self.button_find = Button(label="Find ID", width=150)
        self.button_find.on_click(self.button_find_callback)

        self.radio_button_data_type = RadioButtonGroup(labels=['all', 'working', 'non-working'], active=0, width=300)
        self.radio_button_data_type.on_change('active', self.radio_button_data_type_change)

        self.range_slider_count = RangeSlider(start=0, end=600, value=(0, 600), step=1, title="Sibling Count")
        self.range_slider_count.on_change('value', self.range_slider_count_change)
        self.button_count_start_minus = Button(label="-", width=15)
        self.button_count_start_minus.on_click(self.button_count_start_minus_callback)
        self.button_count_start_plus = Button(label="+", width=15)
        self.button_count_start_plus.on_click(self.button_count_start_plus_callback)
        self.button_count_end_minus = Button(label="-", width=15)
        self.button_count_end_minus.on_click(self.button_count_end_minus_callback)
        self.button_count_end_plus = Button(label="+", width=15)
        self.button_count_end_plus.on_click(self.button_count_end_plus_callback)

        self.range_slider_area = RangeSlider(start=-100, end=131000, value=(0, 131000), step=100, title="Area Size")
        self.range_slider_area.on_change('value', self.range_slider_area_change)
        self.button_area_start_minus = Button(label="-", width=15)
        self.button_area_start_minus.on_click(self.button_area_start_minus_callback)
        self.button_area_start_plus = Button(label="+", width=15)
        self.button_area_start_plus.on_click(self.button_area_start_plus_callback)
        self.button_area_end_minus = Button(label="-", width=15)
        self.button_area_end_minus.on_click(self.button_area_end_minus_callback)
        self.button_area_end_plus = Button(label="+", width=15)
        self.button_area_end_plus.on_click(self.button_area_end_plus_callback)

        self.filters_active = CheckboxGroup(labels=["Filter by Sibling Count", "Filter by Area"], active=[0, 1])
        self.filters_active.on_change('active', self.filters_active_change)

        self.toggle_blend = Toggle(label="Toggle Blend", active=False, width=150)
        self.toggle_blend.on_click(self.toggle_blend_callback)

        self.slider_blend = Slider(start=0.0, end=1.0, value=0.5, step=0.025, title="Blend Ratio")
        self.slider_blend.on_change('value', self.slider_blend_change)

        self.text_count = Paragraph(text="")
        self.update_text_count()

    def parse_text_to_find_idx(self):
        print(self.text_id.text)
        start_bracket = self.text_id.text.find("(")
        end_bracket = self.text_id.text.find(")")
        idx = int(self.text_id.text[start_bracket + 1:end_bracket])
        return idx

    def parse_text_to_find_id(self):
        print(self.text_id.text)
        start_bracket = self.text_id.text.find(":")
        end_bracket = self.text_id.text.find("(")
        select_id = int(self.text_id.text[start_bracket + 1:end_bracket])
        return select_id

    def toggle_sde_ellipse_callback(self, arg):
        print("Toggle Ellipse: Callback: " + str(self.map_tweet_data.circle_id) + " : " + str(self.map_tweet_data.circle_idx))
        if arg:
            self.map_tweet_data.update_sde_ellipse()
            self.map_tweet_data.update_siblings()
        else:
            self.map_tweet_data.clear_sde_ellipse()
            self.map_tweet_data.clear_siblings()

    def toggle_sibling_ellipses_callback(self, arg):
        if arg:
            self.map_tweet_data.update_sibling_ellipses()
        else:
            self.map_tweet_data.clear_sibling_ellipses()

    def toggle_dissolve_callback(self, arg):
        print("Toggle Dissolve: Callback: " + str(self.map_tweet_data.circle_id) + " : " + str(self.map_tweet_data.circle_idx))
        if arg:
            self.map_tweet_data.update_dissolve()
        else:
            self.map_tweet_data.clear_dissolve()

    def toggle_user_info_callback(self, arg):
        print("Toggle User Info: Callback: " + str(self.map_tweet_data.circle_id))

        if arg:
            username, profile_text = self.user_info.find_user_profile(int(self.map_tweet_data.circle_id))
            self.text_username.text = "Username: " + str(username)
            self.text_profile.text = "Profile: " + str(profile_text)
        else:
            self.text_username.text = "Username:"
            self.text_profile.text = "Profile:"

    def radio_button_data_type_change(self, attrname, old, new):
        self.map_tweet_data.switch_tweet_data(new)

        self.map_tweet_data.apply_filters()

        if self.toggle_sde_ellipse.active:
            self.map_tweet_data.update_sde_ellipse()
            self.map_tweet_data.update_siblings()
        else:
            self.map_tweet_data.clear_sde_ellipse()
            self.map_tweet_data.clear_siblings()

        if self.toggle_sibling_ellipses.active:
            self.map_tweet_data.update_sibling_ellipses()
        else:
            self.map_tweet_data.clear_sibling_ellipses()

        if self.toggle_dissolve.active:
            self.map_tweet_data.update_dissolve()
        else:
            self.map_tweet_data.clear_dissolve()

        self.map_tweet_data.clear_find_circle()

    def button_find_callback(self):
        id_value = int(self.text_input.value)
        print("Button Find: Callback: " + str(id_value))

        # This is a dataframe of the selected (single) row.
        id_df = self.map_tweet_data.tweet_data_df.loc[self.map_tweet_data.tweet_data_df['id'] == id_value]
        if id_df.empty:
            print("No row with ID: " + str(id_value))
            self.map_tweet_data.clear_find_circle()
        else:
            idx = id_df.index
            # We can assume that there will be only one index value in this list, as IDs are unique.
            self.map_tweet_data.update_find_circle_index(idx[0])
            self.map_tweet_data.update_find_circle()

    def range_slider_count_change(self, attrname, old, new):
        start_count = int(new[0])
        end_count = int(new[1])
        self.map_tweet_data.filter_circles_by_count(start_count, end_count)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def button_count_start_minus_callback(self):
        value_start = int(self.range_slider_count.value[0])
        value_end = int(self.range_slider_count.value[1])

        if value_start > self.range_slider_count.start:
            value_start -= 1

        new_values = [value_start, value_end]
        self.range_slider_count.value = new_values

        self.map_tweet_data.filter_circles_by_count(value_start, value_end)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def button_count_start_plus_callback(self):
        value_start = int(self.range_slider_count.value[0])
        value_end = int(self.range_slider_count.value[1])

        if value_start < value_end:
            value_start += 1

        new_values = [value_start, value_end]
        self.range_slider_count.value = new_values

        self.map_tweet_data.filter_circles_by_count(value_start, value_end)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def button_count_end_minus_callback(self):
        value_start = int(self.range_slider_count.value[0])
        value_end = int(self.range_slider_count.value[1])

        if value_end > value_start:
            value_end -= 1

        new_values = [value_start, value_end]
        self.range_slider_count.value = new_values

        self.map_tweet_data.filter_circles_by_count(value_start, value_end)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def button_count_end_plus_callback(self):
        value_start = int(self.range_slider_count.value[0])
        value_end = int(self.range_slider_count.value[1])

        if value_end < self.range_slider_count.end:
            value_end += 1

        new_values = [value_start, value_end]
        self.range_slider_count.value = new_values

        self.map_tweet_data.filter_circles_by_count(value_start, value_end)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def range_slider_area_change(self, attrname, old, new):
        start_area = int(new[0])
        end_area = int(new[1])
        self.map_tweet_data.filter_circles_by_area(start_area, end_area)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def button_area_start_minus_callback(self):
        value_start = int(self.range_slider_area.value[0])
        value_end = int(self.range_slider_area.value[1])

        if value_start > self.range_slider_area.start:
            value_start -= 100

        new_values = [value_start, value_end]
        self.range_slider_area.value = new_values

        self.map_tweet_data.filter_circles_by_area(value_start, value_end)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def button_area_start_plus_callback(self):
        value_start = int(self.range_slider_area.value[0])
        value_end = int(self.range_slider_area.value[1])

        if value_start < value_end:
            value_start += 100

        new_values = [value_start, value_end]
        self.range_slider_area.value = new_values

        self.map_tweet_data.filter_circles_by_area(value_start, value_end)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def button_area_end_minus_callback(self):
        value_start = int(self.range_slider_area.value[0])
        value_end = int(self.range_slider_area.value[1])

        if value_end > value_start:
            value_end -= 100

        new_values = [value_start, value_end]
        self.range_slider_area.value = new_values

        self.map_tweet_data.filter_circles_by_area(value_start, value_end)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def button_area_end_plus_callback(self):
        value_start = int(self.range_slider_area.value[0])
        value_end = int(self.range_slider_area.value[1])

        if value_end < self.range_slider_area.end:
            value_end += 100

        new_values = [value_start, value_end]
        self.range_slider_area.value = new_values

        self.map_tweet_data.filter_circles_by_area(value_start, value_end)
        self.map_tweet_data.clear_all()
        self.update_text_count()

    def update_text_count(self):
        count_active, count_total = self.map_tweet_data.num_of_points_active()
        percentage = (count_active/count_total) * 100.0
        text_count_str = "Count: Total: " + str(count_active) + " : " + str(percentage) + "%"
        self.text_count.text = text_count_str

    def filters_active_change(self, attrname, old, new):
        self.map_tweet_data.filters_active(new)
        self.update_text_count()

    def toggle_blend_callback(self, arg):
        print("Toggle Blend: Callback: " + str(self.map_tweet_data.circle_id) + " : " + str(self.map_tweet_data.circle_idx))
        if arg:
            self.map_tweet_data.turn_blend_on()
            self.slider_blend.disabled = False
        else:
            self.map_tweet_data.turn_blend_off()
            self.slider_blend.disabled = True

    def slider_blend_change(self, attrname, old, new):
        print(new)
        self.map_tweet_data.blend(new)
