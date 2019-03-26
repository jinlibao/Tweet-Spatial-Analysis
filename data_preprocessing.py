#!/usr/bin/env python3
import math
import timeit
import logging
import logging.config
import numpy as np
import pandas as pd
import datetime

from argparse import ArgumentParser
from utils import analysis_utilities as au
from utils import configuration_utilities as cu
from utils import file_utilities as fu
from utils import tweet_data_utilities as tdu
from mpi4py import MPI

logger = logging.getLogger()

pd.set_option('expand_frame_repr', False)
pd.options.display.max_rows = 999
# Options and Settings: https://pandas.pydata.org/pandas-docs/stable/options.html

num_of_rows_to_process = None
# Define as None to process all rows.


class TweetDataPreProcessing:

    def __init__(self, file_details, config=None, encoding="utf-8", separator=","):
        self.file_details = file_details
        self.config = config
        self.encoding = encoding
        self.separator = separator
        self.df = None
        self.tweet_data_all = None
        self.tweet_data_working = None
        self.tweet_data_non_working = None
        self.tweet_data_working_overlap = None

    def build_overlap_matrix(self, df, rows, filename='./overlap_matrix.csv'):
        start_time = timeit.default_timer()
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        if rank == 0:
            # rows, cols = df.shape
            rows, cols = 400, 20
            tweet_data_working_overlap = np.zeros((rows, rows))

            for i in range(rows):
                for j in range(rows):
                    if i == j:
                        tweet_data_working_overlap[i, j] = 1
                    elif j < i:
                        if df['a'][i] == 0 or df['b'][i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                            tweet_data_working_overlap[i, j] = -1
                        else:
                            tweet_data_working_overlap[i, j] = au.are_two_ellipses_overlapping(
                                df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i],
                                df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j],
                            )
                if (i + 1) % 10 == 0:
                    print('CPU {:04d}: {:6.2f}% accomplished'.format(rank,  (i + 1) / rows * 100))

            for i in range(rows):
                tweet_data_working_overlap[0:i, i] = tweet_data_working_overlap[i, 0:i]
            print('CPU {:04d}: Assigned transposed tril to triu (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))

            tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap.astype(int), columns=df['id'][0:rows], index=df['id'][0:rows])
            print('CPU {:04d}: Converted numpy array to pandas DataFrame (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
            # print(tweet_data_working_overlap)
            print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(rank, filename, timeit.default_timer() - start_time))
            tweet_data_working_overlap.to_csv(filename, sep=',', header=True, index=True)
            self.tweet_data_working_overlap = tweet_data_working_overlap

            elapsed_time = timeit.default_timer() - start_time
            hour = math.floor(elapsed_time / 3600)
            minute = math.floor((elapsed_time - hour * 3600) / 60)
            second = elapsed_time - 3600 * hour - 60 * minute
            print('Time elapsed: {:d} hours {:d} minutes {:.2f} seconds'.format(hour, minute, second))

    def build_overlap_matrix_parallel(self, df, rows, filename='./overlap_matrix_parallel.csv'):
        start_time = timeit.default_timer()
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        status = MPI.Status()
        # rows, cols = df.shape
        rows, cols = 400, 10

        row_start, rows_local = self.find_range_1D(rows, size, rank)
        # print('Total CPUs: {:4d}, Rank {:4d}, rows: {:5d}, rows_per_cpu: {:5d}, row_start: {:5d}'.format(size, rank, rows, rows_local, row_start))

        tweet_data_working_overlap_local = np.zeros((rows_local, rows))
        for i in range(row_start, row_start + rows_local):
            for j in range(rows):
                if i == j:
                    tweet_data_working_overlap_local[i - row_start, j] = 1
                elif j < i:
                    if df['a'][i] == 0 or df['b'][i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                        tweet_data_working_overlap_local[i - row_start, j] = -1
                    else:
                        tweet_data_working_overlap_local[i - row_start, j] = au.are_two_ellipses_overlapping(
                            df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i],
                            df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j],
                        )

            if (i - row_start + 1) % 10 == 0:
                print('CPU {:04d}: {:6.2f}% accomplished'.format(rank,  (i - row_start + 1) / rows_local * 100))
        # print(tweet_data_working_overlap, rank)

        if rank != 0:
            comm.send(tweet_data_working_overlap_local, dest=0)

        if rank == 0:
            tweet_data_working_overlap = np.zeros((rows, rows))
            tweet_data_working_overlap[0:rows_local] = tweet_data_working_overlap_local
            for i in range(row_start, row_start + rows_local):
                tweet_data_working_overlap[0:i, i] = tweet_data_working_overlap[i, 0:i]
            number_of_sources = 1
            while number_of_sources < size:
                data = comm.recv(source=MPI.ANY_SOURCE, status=status)
                k = status.Get_source()
                print('CPU {:04d}: Received data from CPU {:04d} (Time: {:.2f} seconds)'.format(rank,  k, timeit.default_timer() - start_time))
                row_start, rows_local = self.find_range_1D(rows, size, k)
                tweet_data_working_overlap[row_start:row_start + rows_local] = data + tweet_data_working_overlap[row_start:row_start + rows_local]
                for i in range(row_start, row_start + rows_local):
                    tweet_data_working_overlap[0:i, i] = tweet_data_working_overlap[i, 0:i]
                number_of_sources += 1
            print('CPU {:04d}: Merged received data to global numpy array (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))

            tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap.astype(int), columns=df['id'][0:rows], index=df['id'][0:rows])
            print('CPU {:04d}: Converted numpy array to pandas DataFrame (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
            # print(tweet_data_working_overlap)
            print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(rank, filename, timeit.default_timer() - start_time))
            tweet_data_working_overlap.to_csv(filename, sep=',', header=True, index=True)
            self.tweet_data_working_overlap = tweet_data_working_overlap

            elapsed_time = timeit.default_timer() - start_time
            hour = math.floor(elapsed_time / 3600)
            minute = math.floor((elapsed_time - hour * 3600) / 60)
            second = elapsed_time - 3600 * hour - 60 * minute
            print('Time elapsed: {:d} hours {:d} minutes {:.2f} seconds'.format(hour, minute, second))

    def build_overlap_matrix_parallel_block(self, df, rows, filename='overlap_matrix_parallel_block.csv'):
        start_time = timeit.default_timer()
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        status = MPI.Status()
        # rows, cols = df.shape
        rows, cols = 400, 10
        cpu_required, row_start, col_start, rows_local, cols_local = self.find_range_2D(rows, size, rank)
        # print('Total CPUs: {:4d}, Rank {:4d}, rows: {:5d}, rows_per_cpu: {:5d}, row_start: {:5d}, cols_per_cpu: {:5d}, col_start: {:5d}'.format(size, rank, rows, rows_local, row_start, cols_local, col_start))

        if rank < cpu_required:
            tweet_data_working_overlap_local = np.zeros((rows_local, cols_local))
            for i in range(row_start, row_start + rows_local):
                for j in range(col_start, col_start + cols_local):
                    if i == j:
                        tweet_data_working_overlap_local[i - row_start, j - col_start] = 1
                    elif j < i:
                        if df['a'][i] == 0 or df['b'][i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                            tweet_data_working_overlap_local[i - row_start, j - col_start] = -1
                        else:
                            tweet_data_working_overlap_local[i - row_start, j - col_start] = au.are_two_ellipses_overlapping(
                                df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i],
                                df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j]
                            )
                if row_start == col_start:
                    tweet_data_working_overlap_local[0:i - row_start + 1, i - row_start] = tweet_data_working_overlap_local[i - row_start, 0:i - row_start + 1]

                if (i - row_start + 1) % 10 == 0:
                    print('CPU {:04d}: {:6.2f}% accomplished'.format(rank,  (i - row_start + 1) / rows_local * 100))

            # print(tweet_data_working_overlap, rank)

            if rank != 0:
                comm.send(tweet_data_working_overlap_local, dest=0)
                # print('CPU {:04d}: Sent data to CPU {:04d}'.format(rank,  0))

            if rank == 0:
                tweet_data_working_overlap = np.zeros((rows, rows))
                tweet_data_working_overlap[0:rows_local, 0:cols_local] = tweet_data_working_overlap_local

                number_of_sources = 1
                while number_of_sources < cpu_required:
                    data = comm.recv(source=MPI.ANY_SOURCE, status=status)
                    k = status.Get_source()
                    print('CPU {:04d}: Received data from CPU {:04d} (Time: {:.2f} seconds)'.format(rank,  k, timeit.default_timer() - start_time))
                    cpu_required, row_start, col_start, rows_local, cols_local = self.find_range_2D(rows, size, k)
                    tweet_data_working_overlap[row_start:row_start + rows_local, col_start:col_start + cols_local] = data
                    if row_start != col_start:
                        tweet_data_working_overlap[col_start:col_start + cols_local, row_start:row_start + rows_local] = data.T
                    number_of_sources += 1
                print('CPU {:04d}: Merged received data to global numpy array (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))


                tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap.astype(int), columns=df['id'][0:rows], index=df['id'][0:rows])
                print('CPU {:04d}: Converted numpy array to pandas DataFrame (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
                # print(tweet_data_working_overlap)
                print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(rank, filename, timeit.default_timer() - start_time))
                tweet_data_working_overlap.to_csv(filename, sep=',', header=True, index=True)
                self.tweet_data_working_overlap = tweet_data_working_overlap

                elapsed_time = timeit.default_timer() - start_time
                hour = math.floor(elapsed_time / 3600)
                minute = math.floor((elapsed_time - hour * 3600) / 60)
                second = elapsed_time - 3600 * hour - 60 * minute
                print('Time elapsed: {:d} hours {:d} minutes {:.2f} seconds'.format(hour, minute, second))


    def find_range_1D(self, number, number_of_cpu, rank):
        number_per_cpu = number / number_of_cpu
        number_per_cpu_ceil = math.ceil(number_per_cpu)
        number_per_cpu_floor = math.floor(number_per_cpu)
        threshold = number - number_per_cpu_floor * number_of_cpu

        if rank < threshold:
            start = rank * number_per_cpu_ceil
            number_local = number_per_cpu_ceil
        else:
            start = threshold * number_per_cpu_ceil + (rank - threshold) * number_per_cpu_floor
            number_local = number_per_cpu_floor

        return (start, number_local)

    def find_range_2D(self, rows, cpu_available, rank):
        cols = rows
        m = math.floor(math.sqrt(2 * cpu_available + 1 / 4) - 1 / 2)
        cpu_required = int(m * (m + 1) / 2)
        i = math.floor(math.sqrt(2 * rank + 1 / 4) - 1 / 2)
        j = int(rank - i * (i + 1) / 2)

        row_start, rows_local = self.find_range_1D(rows, m, i)
        col_start, cols_local = self.find_range_1D(cols, m, j)

        return (cpu_required, row_start, col_start, rows_local, cols_local)


    def test_overlap(self):

        dt = datetime.datetime.now()
        year = dt.year
        month = dt.month
        day= dt.day
        hour = dt.hour
        minute = dt.minute
        second = dt.second

        # algorithm = 'serial'
        # algorithm = 'parallel'
        algorithm = 'parallel_block'

        parser = ArgumentParser()
        parser.add_argument('-l', '--loc', dest='loc', help='Location to save the output file')
        parser.add_argument('-n', '--name', dest='name', help='Name for the output file')
        parser.add_argument('-o', '--outfile', dest='filename', help="Path for saving the output file")
        parser.add_argument('-r', '--rows', dest='rows', help="Path for saving the output file")
        args = parser.parse_args()

        if args.loc:
            loc = args.loc
        else:
            loc = '.'

        if args.name:
            name = args.name
        else:
            if platform.system() == 'Linux':
                node = 'teton'
                # node = 'moran'
                # partition = 'hugemem'
                partition = 'regular'
                # partition = 'knl'
                # partition = 'gpu'
            elif platform.system() == 'Darwin':
                node = 'libaoutrage'
                partition = 'macOS'
            else:
                node = 'intel'
                partition = 'partition'

            name = 'overlap_matrix_{:s}_{:s}_{:s}_{:d}_{:02d}_{:02d}_{:02d}_{:02d}_{:02d}.csv'.format(node, partition, algorithm, year, month, day, hour, minute, second)

        if args.filename:
            filename = args.filename
        else:
            filename = '{:s}/{:s}'.format(loc, name)

        # rows, _ = tdp.tweet_data_working.df.shape
        if args.rows:
            rows = int(args.rows)
        else:
            rows, _ = tdp.tweet_data_working.df.shape

        if algorithm == 'serial':
            build_overlap_matrix(tdp.tweet_data_working.df, rows, filename)
        elif algorithm == 'parallel':
            build_overlap_matrix_parallel(tdp.tweet_data_working.df, rows, filename)
        else:
            build_overlap_matrix_parallel_block(tdp.tweet_data_working.df, rows, filename)

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

        logger.info("Processing: Working")
        print("Processing: Working")
        self.tweet_data_working = tdu.TweetData("working", self.config)
        self.tweet_data_working.create_dataframe(self.df, 'User-ID', 'latitude-median-working-tweets',
                                            'longitude-median-working-tweets', 'area-working-tweets',
                                            'x/y-working-tweets', 'theta-working-tweets', 'medians-distance')
        logger.info(self.tweet_data_working)
        self.write_to_json(self.tweet_data_working.df, "tweets_median_working.json")
        logger.info("Processed.")
        print("Processed.")

        logger.info("Processing: All")
        print("Processing: All")
        self.tweet_data_all = tdu.TweetData("all", self.config)
        self.tweet_data_all.create_dataframe(self.df, 'User-ID', 'latitude-mean-all-tweets',
                                                'longitude-mean-all-tweets', 'area-all-tweets', 'x/y-all-tweets',
                                                'theta-all-tweets', 'medians-distance')
        logger.info(self.tweet_data_all)
        self.write_to_json(self.tweet_data_all.df, "tweet_mean_all.json")
        logger.info("Processed.")
        print("Processed.")

        logger.info("Processing: Non-Working")
        print("Processing: Non-Working")
        self.tweet_data_non_working = tdu.TweetData("non_working", self.config)
        self.tweet_data_non_working.create_dataframe(self.df, 'User-ID', 'latitude-median-nonworking-tweets',
                                                'longitude-median-nonworking-tweets', 'area-nonworking-tweets',
                                                'x/y-nonworking-tweets', 'theta-nonworking-tweets', 'medians-distance')
        logger.info(self.tweet_data_non_working)
        self.write_to_json(self.tweet_data_non_working.df, "tweets_median_non_working.json")
        logger.info("Processed.")
        print("Processed.")

    '''
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
    '''

    def read_from_json(self, mean_all, working, non_working):
        self.tweet_data_all = tdu.TweetData("all")
        self.tweet_data_all.read_from_json(mean_all)
        #self.process_color(self.tweet_data_all.df, config)

        self.tweet_data_working = tdu.TweetData("working")
        self.tweet_data_working.read_from_json(working)
        #self.process_color(self.tweet_data_working.df, config)

        self.tweet_data_non_working = tdu.TweetData("non_working")
        self.tweet_data_non_working.read_from_json(non_working)
        #self.process_color(self.tweet_data_non_working.df, config)

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
        dist = au.calculate_distance_non_projected_havershine(
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

    tweet_spatial_analysis_config = cu.TweetSpatialAnalysisConfig("tweet_spatial_analysis.ini")
    logger.info(tweet_spatial_analysis_config)

    file_open = fu.FileOpen("data", "tweet-data.csv")
    logger.info(file_open)

    pre_processor = TweetDataPreProcessing(file_open, tweet_spatial_analysis_config)
    logger.info(pre_processor)
    pre_processor.process()


if __name__ == '__main__':
    main()


