import pandas as pd
import numpy as np
import logging.config

from utils import file_utilities as fu

logger = logging.getLogger()

pd.set_option('expand_frame_repr', False)
pd.options.display.max_rows = 999
# See: Options and Settings: https://pandas.pydata.org/pandas-docs/stable/options.html


class UserProfileDetails:

    def __init__(self, file_details, encoding="utf-8", separator=","):
        self.file_details = file_details
        self.encoding = encoding
        self.separator = separator
        self.df = None

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
            df = pd.read_csv(self.file_details.absolute, encoding=self.encoding, sep=self.separator)
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

    def find_user_profile(self, point_id):
        user_profile_df = self.df.loc[self.df['User-ID'] == point_id]

        if user_profile_df.empty:
            return None, None

        idx = user_profile_df.index[0]
        #userid = user_profile_df['User-ID'][idx]
        username = user_profile_df['UserName'][idx]
        profile_text = user_profile_df['Profile-Text'][idx]
        return username, profile_text

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

    def __str__(self):
        df_str = "UserProfileDetails:\nFilename: " + str(self.file_details.absolute)
        if self.df is not None:
            df_str += "\nShape: "+ str(self.df.shape)
            df_str += "\nColumn Names:\n" + str(self.df.columns)
            df_str += "\nColumn dtypes:\n" + str(self.df.dtypes)
        else:
            df_str += "\nDataframe NOT defined."

        return df_str

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def main():
    logging.basicConfig(
        filename='logs/user_profile_details.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s')

    logger.info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    logger.info("User Profile Details::")

    file_open = fu.FileOpen("data", "user-info.csv")
    logger.info(file_open)

    user_profile_details = UserProfileDetails(file_open)
    user_profile_details.process()
    logger.info(user_profile_details)

    user_id = 681473
    username, profile_text = user_profile_details.find_user_profile(user_id)
    print("UserID: " + str(user_id) + " : " + str(username) + " : " + str(profile_text))

    user_id = -1
    username, profile_text = user_profile_details.find_user_profile(user_id)
    print("UserID: " + str(user_id) + " : " + str(username) + " : " + str(profile_text))


if __name__ == '__main__':
    main()
