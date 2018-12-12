def initialise_cds():
    global tweet_data

    tweet_data = pre_processor.df_tweet_data_all.df

    source_circle.data = source_circle.from_df(tweet_data)

    data_hover_idx = {
        'id' : [-1],
        'idx' : [-1]
    }

    df_hover_idx = pd.DataFrame(data_hover_idx)
    source_hover_idx.data = source_hover_idx.from_df(df_hover_idx)

    e_data = {
        'x': [],
        'y': [],
        'width': [],
        'height': [],
        'angle': [],
    }

    df_ellipse = pd.DataFrame(e_data)
    source_ellipse.data = source_ellipse.from_df(df_ellipse)

    es_data = {
        'x': [],
        'y': [],
        'width': [],
        'height': [],
        'angle': []
    }

    df_ellipse_siblings = pd.DataFrame(es_data)
    source_ellipse_siblings.data = source_ellipse_siblings.from_df(df_ellipse_siblings)

    cie_data = {
        'x': [],
        'y': []
    }

    df_circle_inside_ellipse = pd.DataFrame(cie_data)
    source_circle_inside_ellipse.data = source_circle_inside_ellipse.from_df(df_circle_inside_ellipse)

    patch_data = {
        'x': [],
        'y': []
    }

    df_patch = pd.DataFrame(patch_data)
    source_patch_dissolve.data = source_patch_dissolve.from_df(df_patch)

    def toggle_ellipse_callback(arg):
        global ellipse_idx
        global interface_properties
        global map_tweet_data

        # idx = parse_text_to_find_idx()
        id = parse_text_to_find_id()
        df = tweet_data.loc[tweet_data['id'] == id]
        idx = df.index[0]
        print("Toggle Ellipse: Callback: " + str(id) + " : " + str(idx))

        if arg:
            map_tweet_data.update_ellipse(idx)
            map_tweet_data.update_circle_inside_ellipse(idx)
        else:
            map_tweet_data.clear_ellipse()
            map_tweet_data.clear_circle_inside_ellipse()

        interface_properties.toggle_ellipse_active = arg
        ellipse_idx = idx



