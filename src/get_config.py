import configparser
import os
def get_config(config_filename = 'config.ini'):
    config_file = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'config' + os.sep + config_filename
    # Check config file
    if not os.path.exists(config_file):
        raise FileNotFoundError('File not found: ' + config_file + \
                            '. Please configure it first by copying/editing ' +  'config.ini.dist')
    # Load config
    config = configparser.ConfigParser()
    config.read(config_file)
    return config
