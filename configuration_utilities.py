
import configparser
import logging
import logging.config

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

    def parse_range_str(self, range_str):
        range_list = range_str.split()
        logger.debug("Range List: " + str(range_list))
        range_list_floats = []
        for value in range_list:
            range_list_floats.append(float(value))

        return range_list_floats

    def __str__(self):
        config_str = "TweetSpatialAnalysisConfig:\nFilename: " + str(self.filename)
        config_str += "\nLatitude: " + str(self.latitude[0]) + " to " + str(self.latitude[1])
        config_str += "\nLongitude: " + str(self.longitude[0]) + " to " + str(self.longitude[1])
        config_str += "\nCount: " + str(self.count[0]) + " to " + str(self.count[1]) + " step: " + str(self.count[2])
        config_str += "\nArea: " + str(self.area[0]) + " to " + str(self.area[1]) + " step: " + str(self.area[2])
        config_str += "\nDistance: " + str(self.distance[0]) + " to " + str(self.distance[1]) + " step: " + str(self.distance[2])
        config_str += "\nRatio: " + str(self.ratio[0]) + " to " + str(self.ratio[1]) + " step: " + str(self.ratio[2])
        return config_str

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def main():
    logging.basicConfig(
        filename='logs/tweet_spatial_analysis_config.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s')

    logger.info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    logger.info("Tweet Spatial Analysis: Config:")

    tweet_spatial_analysis_config = TweetSpatialAnalysisConfig("tweet_spatial_analysis.ini")

    logger.info(tweet_spatial_analysis_config)


if __name__ == '__main__':
    main()
