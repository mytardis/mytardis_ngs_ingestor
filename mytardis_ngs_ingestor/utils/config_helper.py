import os
from argparse import ArgumentParser
import collections

import pytoml as toml

from attrdict import AttrDict
from toolz.dicttoolz import merge as merge_dicts


# http://stackoverflow.com/a/26853961
# Deprecated by use of toolz.dicttoolz.merge_dicts
# def merge_dicts(*dict_args):
#     """
#     Given any number of dicts, shallow copy and merge into a new dict,
#     precedence goes to key value pairs in latter dicts.
#     """
#     result = {}
#     for dictionary in dict_args:
#         result.update(dictionary)
#     return result


# Based on: https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
def recursive_dict_merge(dct, merge_dct):
    """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: The merged dictionary
    :rtype: dict
    """
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                     and isinstance(merge_dct[k], collections.Mapping)):
            recursive_dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]

    return dct


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
    else:
        # Copy
        defaults = dict(defaults)

    with open(config_filepath, 'r') as f:
        config = toml.load(f)

    config = recursive_dict_merge(defaults, config)
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

    config_options = AttrDict()
    if config_file:
        config_options = parse_settings_toml(config_file, defaults=defaults)
    elif os.path.isfile(default_config_filename):
        config_options = parse_settings_toml(default_config_filename,
                                             defaults=defaults)

    return config_options
