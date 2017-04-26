import os
from argparse import ArgumentParser

import pytoml as toml

from attrdict import AttrDict
from toolz.dicttoolz import merge as merge_dicts


def parse_settings_toml(config_filepath, defaults=None):
    """
    Grabs key/value pairs from a TOML format config file,
    returns an argparse options-style data container object
    with the config. A dictionary can be provided with defaults,
    which any option specified in the config file will override.

    (The defaults dictionary might can be generated via passing a
    pre-configured argparse.ArgumentParser into get_parser_defaults)

    :param config_filepath: The path to a TOML config file.
    :type config_filepath: basestring
    :param defaults: A dictionary of defaults
    :type defaults: dict
    :return: A data container object of config options.
    :rtype: object
    """
    if defaults is None:
        defaults = dict()

    with open(config_filepath, 'r') as f:
        config = toml.load(f)

    config = merge_dicts(defaults, config)
    config = AttrDict(config)

    return config


def get_parser_defaults(parser):
    """
    Given an argparse.ArgumentParser pre-configured with args/options via
    add_argument, return a dictionary of {dest: default}, containing the
    options attribute names and their default value.

    :param parser: A pre-configured parser instance.
    :type parser: argparse.ArgumentParser
    :return: A dictionary of {dest: default}, containing the options attribute
             names and their default value
    :rtype: dict
    """
    defaults = dict()
    for a in parser._actions:
        defaults[a.dest] = a.default

    if 'help' in defaults:
        del defaults['help']

    return defaults


def get_config_toml(config_file=None,
                    defaults=None,
                    default_config_filename='uploader_config.toml'):

    if not defaults:
        defaults = dict()

    config_options = dict()
    if config_file:
        config_options = parse_settings_toml(config_file, defaults=defaults)
    elif os.path.isfile(default_config_filename):
        config_options = parse_settings_toml(default_config_filename,
                                             defaults=defaults)

    return config_options
