
from bokeh.models import RangeSlider
from bokeh.models.widgets import Button, CheckboxGroup, Div, Paragraph, RadioButtonGroup, Slider, TextInput, Toggle


class MapWidgets:

    def __init__(self, tweet_data_controller, user_info):
        self.tweet_data_controller = tweet_data_controller

        self.toggle_sde_ellipse = Toggle(label="Toggle Ellipse", active=False, width=150)
        self.toggle_sde_ellipse.on_click(self.toggle_sde_ellipse_callback)

        self.toggle_sibling_ellipses = Toggle(label="Toggle Siblings", active=False, width=150)
        self.toggle_sibling_ellipses.on_click(self.toggle_sibling_ellipses_callback)

        self.toggle_dissolve = Toggle(label="Toggle Dissolve", active=False, width=150)
        self.toggle_dissolve.on_click(self.toggle_dissolve_callback)

        self.text_selection_details = Div(text='Selected ID:')

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
        self.toggle_blend.disabled = True

        self.slider_blend = Slider(start=0.0, end=1.0, value=0.5, step=0.025, title="Blend Ratio")
        self.slider_blend.on_change('value', self.slider_blend_change)
        self.slider_blend.disabled = True

        self.text_count = Paragraph(text="")
        self.update_text_count()

    def toggle_sde_ellipse_callback(self, arg):
        print("Toggle Ellipse: Callback: " + str(self.tweet_data_controller.circle_id) + " : " + str(self.tweet_data_controller.circle_idx))
        if arg:
            self.tweet_data_controller.update_sde_ellipse()
            self.tweet_data_controller.update_siblings()
        else:
            self.tweet_data_controller.clear_sde_ellipse()
            self.tweet_data_controller.clear_siblings()

    def toggle_sibling_ellipses_callback(self, arg):
        if arg:
            self.tweet_data_controller.update_sibling_ellipses()
        else:
            self.tweet_data_controller.clear_sibling_ellipses()

    def toggle_dissolve_callback(self, arg):
        print("Toggle Dissolve: Callback: " + str(self.tweet_data_controller.circle_id) + " : " + str(self.tweet_data_controller.circle_idx))
        if arg:
            self.tweet_data_controller.update_dissolve()
        else:
            self.tweet_data_controller.clear_dissolve()

    def radio_button_data_type_change(self, attrname, old, new):
        if new == 0:
            self.toggle_blend.disabled = True
            self.toggle_blend.active = False
            self.tweet_data_controller.blend_active = False
        else:
            self.toggle_blend.disabled = False
            self.toggle_blend.active = False
            self.tweet_data_controller.blend_active = False

        self.tweet_data_controller.switch_tweet_dataset(new)
        self.tweet_data_controller.apply_filters()

        if self.toggle_sde_ellipse.active:
            self.tweet_data_controller.update_sde_ellipse()
            self.tweet_data_controller.update_siblings()
        else:
            self.tweet_data_controller.clear_sde_ellipse()
            self.tweet_data_controller.clear_siblings()

        if self.toggle_sibling_ellipses.active:
            self.tweet_data_controller.update_sibling_ellipses()
        else:
            self.tweet_data_controller.clear_sibling_ellipses()

        if self.toggle_dissolve.active:
            self.tweet_data_controller.update_dissolve()
        else:
            self.tweet_data_controller.clear_dissolve()

        self.tweet_data_controller.clear_find_circle()

    def button_find_callback(self):
        id_value = int(self.text_input.value)
        print("Button Find: Callback: " + str(id_value))
        self.tweet_data_controller.find_id(id_value)

    def range_slider_count_change(self, attrname, old, new):
        start_count = int(new[0])
        end_count = int(new[1])
        self.tweet_data_controller.filter_circles_by_count(start_count, end_count)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def button_count_start_minus_callback(self):
        value_start = int(self.range_slider_count.value[0])
        value_end = int(self.range_slider_count.value[1])

        if value_start > self.range_slider_count.start:
            value_start -= 1

        new_values = [value_start, value_end]
        self.range_slider_count.value = new_values

        self.tweet_data_controller.filter_circles_by_count(value_start, value_end)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def button_count_start_plus_callback(self):
        value_start = int(self.range_slider_count.value[0])
        value_end = int(self.range_slider_count.value[1])

        if value_start < value_end:
            value_start += 1

        new_values = [value_start, value_end]
        self.range_slider_count.value = new_values

        self.tweet_data_controller.filter_circles_by_count(value_start, value_end)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def button_count_end_minus_callback(self):
        value_start = int(self.range_slider_count.value[0])
        value_end = int(self.range_slider_count.value[1])

        if value_end > value_start:
            value_end -= 1

        new_values = [value_start, value_end]
        self.range_slider_count.value = new_values

        self.tweet_data_controller.filter_circles_by_count(value_start, value_end)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def button_count_end_plus_callback(self):
        value_start = int(self.range_slider_count.value[0])
        value_end = int(self.range_slider_count.value[1])

        if value_end < self.range_slider_count.end:
            value_end += 1

        new_values = [value_start, value_end]
        self.range_slider_count.value = new_values

        self.tweet_data_controller.filter_circles_by_count(value_start, value_end)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def range_slider_area_change(self, attrname, old, new):
        start_area = int(new[0])
        end_area = int(new[1])
        self.tweet_data_controller.filter_circles_by_area(start_area, end_area)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def button_area_start_minus_callback(self):
        value_start = int(self.range_slider_area.value[0])
        value_end = int(self.range_slider_area.value[1])

        if value_start > self.range_slider_area.start:
            value_start -= 100

        new_values = [value_start, value_end]
        self.range_slider_area.value = new_values

        self.tweet_data_controller.filter_circles_by_area(value_start, value_end)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def button_area_start_plus_callback(self):
        value_start = int(self.range_slider_area.value[0])
        value_end = int(self.range_slider_area.value[1])

        if value_start < value_end:
            value_start += 100

        new_values = [value_start, value_end]
        self.range_slider_area.value = new_values

        self.tweet_data_controller.filter_circles_by_area(value_start, value_end)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def button_area_end_minus_callback(self):
        value_start = int(self.range_slider_area.value[0])
        value_end = int(self.range_slider_area.value[1])

        if value_end > value_start:
            value_end -= 100

        new_values = [value_start, value_end]
        self.range_slider_area.value = new_values

        self.tweet_data_controller.filter_circles_by_area(value_start, value_end)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def button_area_end_plus_callback(self):
        value_start = int(self.range_slider_area.value[0])
        value_end = int(self.range_slider_area.value[1])

        if value_end < self.range_slider_area.end:
            value_end += 100

        new_values = [value_start, value_end]
        self.range_slider_area.value = new_values

        self.tweet_data_controller.filter_circles_by_area(value_start, value_end)
        self.tweet_data_controller.clear_all()
        self.update_text_count()

    def update_text_count(self):
        count_active, count_total = self.tweet_data_controller.num_of_points_active()
        percentage = (count_active/count_total) * 100.0
        text_count_str = "Count: Total: " + str(count_active) + " : " + str(percentage) + "%"
        self.text_count.text = text_count_str

    def filters_active_change(self, attrname, old, new):
        self.tweet_data_controller.filters_active(new)
        self.update_text_count()

    def toggle_blend_callback(self, arg):
        print("Toggle Blend: Callback: " + str(self.tweet_data_controller.circle_id) + " : " + str(self.tweet_data_controller.circle_idx))
        if arg:
            self.tweet_data_controller.turn_blend_on(self.toggle_sde_ellipse.active, self.toggle_sibling_ellipses.active, self.toggle_dissolve.active)
            self.slider_blend.disabled = False
        else:
            self.tweet_data_controller.turn_blend_off()
            self.slider_blend.disabled = True

    def slider_blend_change(self, attrname, old, new):
        self.tweet_data_controller.blend(new, self.toggle_sde_ellipse.active, self.toggle_sibling_ellipses.active, self.toggle_dissolve.active)
