
import logging
import logging.config
import numpy as np

from bokeh.palettes import plasma

from file_utilities import *
from tweet_data import *

logger = logging.getLogger()

pd.set_option('expand_frame_repr', False)
pd.options.display.max_rows = 999
# Options and Settings: https://pandas.pydata.org/pandas-docs/stable/options.html

num_of_rows_to_process = 10000
# Define as None to process all rows.


class TweetDataPreProcessing:

    def __init__(self, file_details, encoding="utf-8", separator=","):
        self.file_details = file_details
        self.encoding = encoding
        self.separator = separator
        self.df = None
        self.tweet_data_all = None
        self.tweet_data_working = None
        self.tweet_data_non_working = None

    def trim_columns(self, df):
        # Trimming / stripping the white space can create an empty cell. This is then NOT picked up by isspace()
        # and not converted to nan.
        # So, we convert cells made purely of white space to nan first, and then trim any other cells.
        df = df.applymap(lambda x: np.nan if isinstance(x, str) and x.isspace() else x)
        # Replacing blank values (white space) with NaN in pandas
        # https://stackoverflow.com/questions/13445241/replacing-blank-values-white-space-with-nan-in-pandas
        df = df.applymap(lambda x: x.strip() if type(x) is str else x)
        # Strip / trim all strings of a dataframe
        # https://stackoverflow.com/questions/40950310/strip-trim-all-strings-of-a-dataframe
        return df

    def read_data(self):
        df = None
        alternative_encodings = ['latin-1']

        date_columns = []
        # Lines below were commented out to switch off using the metadata relating to dates.
        # By applying the 'parse_dates' functionality in the 'read_csv' the underlying data format is lost.
        # i.e. YYYY/MM/DD is converted to the underlying python date-time format.
        # if self.metadata is not None:
        #     date_columns = self.metadata.date

        try:
            df = pd.read_csv(self.file_details.absolute, encoding=self.encoding, sep=self.separator, nrows=num_of_rows_to_process)
            # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html

            if df.shape[0] == 0:
                logger.warning("No rows in file:%s", self.file_details.filename)
                return None

            if df is not None:
                df = self.trim_columns(df)
                return df

        except UnicodeDecodeError as ude:
            logger.warning("UnicodeDecodeError:%s: Reading:%s", ude, self.file_details.absolute)
        except pd.errors.EmptyDataError as ede:
            logger.error("pandas.errors.EmptyDataError:%s: Reading:%s", ede, self.file_details.absolute)
            # This file is completely empty. No column headers, no data.
            return None

        for encoding in alternative_encodings:
            try:
                self.encoding = encoding
                df = pd.read_csv(self.file_details.absolute, encoding=self.encoding,
                                 sep=self.separator, parse_dates=date_columns)

                if df.shape[0] == 0:
                    logger.warning("No rows in file:%s", self.file_details.filename)
                    return None

                if df is not None:
                    self.trim_columns(df)
                    return df

            except UnicodeDecodeError as ude:
                logger.error("UnicodeDecodeError:%s: Reading:%s", ude, self.file_details.absolute)

        return df

    def process(self):
        self.df = self.read_data()
        self.df_details()
        #self.identify_rows_with_infinity_string()
        #self.identify_rows_with_nan()
        #self.df_details()

        logger.info("Processing: All")
        print("Processing: All")
        self.tweet_data_all = TweetData("all")
        self.tweet_data_all.create_dataframe(self.df, 'User-ID', 'latitude-mean-all-tweets',
                                                'longitude-mean-all-tweets', 'area-all-tweets', 'x/y-all-tweets',
                                                'theta-all-tweets', 'medians-distance')
        logger.info(self.tweet_data_all)
        self.write_to_json(self.tweet_data_all.df, "tweet_mean_all.json")
        logger.info("Processed.")
        print("Processed.")

        logger.info("Processing: Working")
        print("Processing: Working")
        self.tweet_data_working = TweetData("working")
        self.tweet_data_working.create_dataframe(self.df, 'User-ID', 'latitude-median-working-tweets',
                                            'longitude-median-working-tweets', 'area-working-tweets',
                                            'x/y-working-tweets', 'theta-working-tweets', 'medians-distance')
        logger.info(self.tweet_data_working)
        self.write_to_json(self.tweet_data_working.df, "tweets_median_working.json")
        logger.info("Processed.")
        print("Processed.")

        logger.info("Processing: Non-Working")
        print("Processing: Non-Working")
        self.tweet_data_non_working = TweetData("non-working")
        self.tweet_data_non_working.create_dataframe(self.df, 'User-ID', 'latitude-median-nonworking-tweets',
                                                'longitude-median-nonworking-tweets', 'area-nonworking-tweets',
                                                'x/y-nonworking-tweets', 'theta-nonworking-tweets', 'medians-distance')
        logger.info(self.tweet_data_non_working)
        self.write_to_json(self.tweet_data_non_working.df, "tweets_median_non_working.json")
        logger.info("Processed.")
        print("Processed.")

    def process_color(self, tweet_data_df, config):

        fill_color_list = plasma(len(config.bins_count))
        fill_color_list.reverse()

        for idx in range(0, tweet_data_df.shape[0]):
            count = tweet_data_df['count'][idx]
            for idx2 in range(0, len(config.bins_count)-1):
                if count < config.bins_count[idx2]:
                    tweet_data_df.loc[idx, 'color'] = fill_color_list[idx2 - 1]
                    break
                else:
                    tweet_data_df.loc[idx, 'color'] = fill_color_list[idx2]

    def read_from_json(self, mean_all, working, non_working, config):
        self.tweet_data_all = TweetData("all")
        self.tweet_data_all.read_from_json(mean_all)
        self.process_color(self.tweet_data_all.df, config)

        self.tweet_data_working = TweetData("working")
        self.tweet_data_working.read_from_json(working)
        self.process_color(self.tweet_data_working.df, config)

        self.tweet_data_non_working = TweetData("non-working")
        self.tweet_data_non_working.read_from_json(non_working)
        self.process_color(self.tweet_data_non_working.df, config)

    def identify_rows_with_nan(self):
        logger.info("Identify rows with NaN values:")
        nan_rows_df = self.df.loc[self.df.isnull().any(axis=1)]
        logger.info("nan_rows: Shape: %s", nan_rows_df.shape)
        if nan_rows_df.shape[0] > 0:
            logger.info("There are rows with null values: ")
            idx_str = "Num of Rows: " + str(nan_rows_df.shape[0]) + " : (index, id) :"
            for index, row in nan_rows_df.iterrows():
                idx_str += " (" + str(index) + ", " + str(row['User-ID']) + ")"
            logger.info(idx_str)

            self.df = self.df.dropna()
            logger.info("Drop nan: Shape: %s", self.df.shape)
        else:
            logger.info("There are NO rows with null values.")

    def identify_rows_with_infinity_string(self):
        logger.info("Identify rows with 'Infinity' string values:")
        #has_inf_rows = self.df.isin(["Infinity"]).any(1)
        has_inf_rows = self.df[self.df['x/y-nonworking-tweets'] == "Infinity"]
        logger.info(has_inf_rows.shape)

    def df_details(self, log_column_details=True):
        logger.info("TweetDataPreProcessing:")
        logger.info("Filename: %s", self.file_details.absolute)
        if self.df is not None:
            logger.info("Shape: %s", self.df.shape)
            if log_column_details:
                logger.info("Column Names:\n%s", self.df.columns)
                logger.info("Column Types:\n%s", self.df.dtypes)
        else:
            logger.info("\nDataframe NOT defined.")

    def write_to_json(self, df, filename):
        logger.info("Writing to JSON: " + filename)
        print("Writing to JSON: " + filename)
        df.to_json(filename, orient='records', lines=True)
        logger.info("Written.")
        print("Written.")

    def __str__(self):
        df_str = "TweetDataPreProcessing:\nFilename: " + str(self.file_details.absolute)
        if self.df is not None:
            df_str += "\nShape: "+ str(self.df.shape)
            df_str += "\nColumn Names:\n" + str(self.df.columns)
            df_str += "\nColumn dtypes:\n" + str(self.df.dtypes)
        else:
            df_str += "\nDataframe NOT defined."

        return df_str


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Confirming that the median distance within the dataset follows the haversine between lat/lon points
# and follows the Earth's (curved) surface.
def calculate_distance_non_projected_havershine_test(df, num_of_rows):
    for row_idx in range(0, num_of_rows):
        dist = calculate_distance_non_projected_havershine(
            df['latitude-median-working-tweets'][row_idx],
            df['longitude-median-working-tweets'][row_idx],
            df['latitude-median-nonworking-tweets'][row_idx],
            df['longitude-median-nonworking-tweets'][row_idx])
        print(str(row_idx) + " : " + str(dist))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def main():
    logging.basicConfig(
        filename='logs/data_preprocessing.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s')

    logger.info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    logger.info("Tweet Data: Pre-Processing:")

    file_open = FileOpen("data", "tweet-data.csv")
    logger.info(file_open)

    pre_processor = TweetDataPreProcessing(file_open)
    logger.info(pre_processor)
    pre_processor.process()


if __name__ == '__main__':
    main()

