"""MyTardis NGS Ingestor

Client for https://github.com/mytardis/mytardis
and https://github.com/pansapiens/mytardis-seqfac
"""

from setuptools import setup, find_packages
from codecs import open
from os import path
import mytardis_ngs_ingestor

here = path.abspath(path.dirname(__file__))


def get_requirements():
    from pip.req import parse_requirements
    from pip.download import PipSession

    # parse_requirements() returns generator of pip.req.InstallRequirement
    # objects
    pip_reqs = parse_requirements(path.join(here, 'requirements.txt'),
                                  session=PipSession())
    # reqs is a list of requirements
    requirements = [str(ir.req) for ir in pip_reqs]

    return requirements

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='mytardis_ngs_ingestor',
    version=mytardis_ngs_ingestor.__version__,

    description='MyTardis NGS Ingestor',
    long_description=long_description,
    url='https://github.com/pansapiens/mytardis_ngs_ingestor',

    author='Andrew Perry',
    author_email='Andrew.Perry@monash.edu',

    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: System Administrators',
        'Topic :: Scientific/Engineering :: Bio-Informatics'

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        # 'Programming Language :: Python :: 3.4',
    ],

    keywords='DNA RNA bioinformatics sequencing illumina client api data files',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=get_requirements(),

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        'dev': ['prospector',
                'check-manifest'],
        'test': ['tox',
                 'pytest'
                 'coverage',
                 'requests_mock'],
    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    # package_data={
    #      'sample': ['sample.dat'],
    # },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'illumina_uploader='
            'mytardis_ngs_ingestor.illumina_uploader:run_in_console',
            'illumina_autoprocess='
            'mytardis_ngs_ingestor.autoprocessing.autoprocess:run_in_console',
        ],
    },
    test_suite='tests',
)
