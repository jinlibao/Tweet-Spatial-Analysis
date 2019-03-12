
import configparser
import logging
import logging.config
import os

from utils import file_utilities as fu

logger = logging.getLogger()


class TweetSpatialAnalysisConfig:

    def __init__(self, config_filename):
        config = configparser.ConfigParser()
        config.read(config_filename)

        self.filename = config_filename
        self.latitude = self.parse_range_str(config['RANGES']['latitude'])
        self.longitude = self.parse_range_str(config['RANGES']['longitude'])
        self.count = self.parse_range_str(config['RANGES']['count'])
        self.area = self.parse_range_str(config['RANGES']['area'])
        self.distance = self.parse_range_str(config['RANGES']['distance'])
        self.ratio = self.parse_range_str(config['RANGES']['ratio'])
        self.dissolve = self.parse_range_str(config['RANGES']['dissolve'])

        self.bins_count = self.parse_range_str(config['HISTOGRAMS']['bins_count'])
        self.bins_count_text = self.create_bins_text(self.bins_count)

        self.bins_ratio = self.parse_range_str(config['HISTOGRAMS']['bins_ratio'])
        self.bins_ratio_text = self.create_bins_text(self.bins_ratio)

        absolute_folder = config['SIBLING_DATA']['folder']
        folder_details = fu.FileOpen(absolute_folder)
        self.sibling_data_folder = folder_details.folder


    def parse_range_str(self, range_str):
        range_list = range_str.split()
        logger.debug("Range List: " + str(range_list))
        range_list_floats = []
        for value in range_list:
            range_list_floats.append(float(value))

        return range_list_floats

    def create_bins_text(self, bins_list):
        bins_text = []
        for idx in range(0, len(bins_list) - 1):
            #print(idx)
            #print(len(bins_list))
            if bins_list[idx + 1] - bins_list[idx] == 1:
                bins_text.append(str(int(bins_list[idx])))
            else:
                #print(str(int(self.bins_count[idx])))
                #print(str(int(self.bins_count[idx + 1])))
                if idx == len(bins_list) - 2:
                    text = str(int(bins_list[idx])) + " - " + str(int(bins_list[idx + 1]))
                else:
                    text = str(int(bins_list[idx])) + " - " + str(int(bins_list[idx + 1] - 1))

                bins_text.append(text)

        return bins_text

    def __str__(self):
        config_str = "TweetSpatialAnalysisConfig:\nFilename: " + str(self.filename)
        config_str += "\nRanges:"
        config_str += "\nLatitude: " + str(self.latitude[0]) + " to " + str(self.latitude[1])
        config_str += "\nLongitude: " + str(self.longitude[0]) + " to " + str(self.longitude[1])
        config_str += "\nCount: " + str(self.count[0]) + " to " + str(self.count[1]) + " step: " + str(self.count[2])
        config_str += "\nArea: " + str(self.area[0]) + " to " + str(self.area[1]) + " step: " + str(self.area[2])
        config_str += "\nDistance: " + str(self.distance[0]) + " to " + str(self.distance[1]) + " step: " + str(self.distance[2])
        config_str += "\nRatio: " + str(self.ratio[0]) + " to " + str(self.ratio[1]) + " step: " + str(self.ratio[2])
        config_str += "\nDissolve: " + str(self.dissolve[0]) + " to " + str(self.dissolve[1]) + " step: " + str(self.dissolve[2])
        config_str += "\nHistograms:"
        config_str += "\nBins Count: " + str(self.bins_count)
        config_str += "\nBins Count Text: " + str(self.bins_count_text)
        config_str += "\nBins Ratio: " + str(self.bins_ratio)
        config_str += "\nBins Ratio Text: " + str(self.bins_ratio_text)
        config_str += "\nSibling Data: Folder: " + str(self.sibling_data_folder)

        return config_str

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def main():
    logging.basicConfig(
        filename='logs/tweet_spatial_analysis_config.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s')

    logger.info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    logger.info("Tweet Spatial Analysis: Config:")

    tweet_spatial_analysis_config = TweetSpatialAnalysisConfig("conf/tweet_spatial_analysis.ini")

    logger.info(tweet_spatial_analysis_config)


if __name__ == '__main__':
    main()
