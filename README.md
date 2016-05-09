mytardis_ngs_ingestor
=====================

`mytardis_ngs_ingestor` is a data ingestion client to upload 
next-generation sequencing data (NGS) to a MyTardis data management 
server.

It runs on Linux and currently supports runs from Illumina HiSeq, 
NextSeq and MiSeq instruments.

Dependencies
------------

  * A [MyTardis](https://github.com/mytardis/mytardis) server with 
    the [mytardis-seqfac](https://github.com/pansapiens/mytardis-seqfac) 
    app installed and configured
  * (Optional) [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Python dependencies are installed from requirements.txt


Installation
------------

Create a Python virtualenv:
```sh
pip install virtualenv
virtualenv ~/.virtualenvs/mytardis_ngs_ingestor
source ~/.virtualenvs/mytardis_ngs_ingestor/bin/activate
```

Install Python dependencies:
```sh
pip install -r requirements.txt
```

While not recommended, the package can alternatively be installed 
system-wide via:
```sh
sudo python setup.py install
```

and the script invoked using `illumina_uploader`.

Configuration
-------------

The program will look for a configuration file named `uploader_config.yaml` 
in the current working directory.

The path to the configuration file can be specified using the  `--config` 
commandline option.

There is an annotated example configuration file `uploader_config_example.yaml`.
Copy this to `uploader_config.yaml` and edit it with your site requirements.

All configuration options can be overridden with a corresponding 
commandline option (eg, `fastqc_bin: /usr/local/bin/fastqc` in the 
config file becomes `--fastqc-bin=/usr/local/bin/fastqc`)


Running
-------

`illumina_uploader.py` does not run `bcl2fastq` automatically. It is a
assumed that the run has already been demultiplexed.

To ingest a run:

```sh
python mytardis_uploader/illumina_uploader.py --path=/mnt/bigdisk/160915_FHT451_0119_AC6AMWACXZ/ \
--bcl2fastq-output-path="{run_path}/{run_id}.bcl2fastq"
```

where `{run_path}/{run_id}.bcl2fastq` is the path where the bcl2fastq 
(fastq.gz files in project/sample directories) is located.

In this example, `{run_path}` is `/mnt/bigdisk/160915_FHT451_0119_AC6AMWACXZ/` 
and `{run_id}` is `160915_FHT451_0119_AC6AMWACXZ`.


How a 'run' is structured in the MyTardis data model
----------------------------------------------------

When a run is ingested, two types of Experiment are produced - 
a single *Run Experiment* and one or more *Project Experiments*.

The *Run Experiment* represents the whole run and should be only 
accessible to Facility Managers.

The *Project Experiments* correspond to distinct projects produced after 
demultiplexing (eg the 'SampleProjects' in SampleSheet.csv). These 
contain a subset of the data from the entire run and are shared with 
end-users by the Facilty Manager. A *Project Experiment* is also 
create for the reads where a barcode could not be assigned.

Datasets are created containing:

  * FASTQ files
  * FastQC reports
  
In addition, the *Run Experiment* contains a Dataset for configuration
and log files that the facility may wish to preserve (currently only 
SampleSheet.csv).

'Raw' .bcl files are not ingested.


Development
-----------

Here is a quick overview of how the project is structured. 

Instances of the `mytardis_uploader.MyTardisUploader` class handle 
(most) REST requests to the MyTardis server.  `mytardis_uploader.py` 
also contains basic config and commandline parsing functions 
(`get_config`, `add_config_args` and `validate_config`). There is no 
domain-specific code related NGS runs in `mytardis_uploader.py`.
The ingestor script `illumina_uploader.py` adds additional 
config/commandline options and instatiates MyTardisUploader instances
for making REST requests (one per MyTardis StorageBox). It handles 
parsing of metadata from Illumina sequencing runs and registering/uploading 
files to MyTardis in a structure suitable of multiplexed runs from 
multi-user facilities. `mytardis_models.py` and `models.py` provide 
data structures for holding run metadata and can serialize themselves to 
JSON for use by `MyTardisUploader`.

`illumina_uploader.py` is due for a refactor - run specific parsing 
functions should be split out into their own library, and logic around 
different filenaming and directoty structures (eg bcl2fastq v1.8.4 vs. 
2.17) should be abstracted behind an API to simplify some of the 
conditionals in `illumina_uploader.py`.
