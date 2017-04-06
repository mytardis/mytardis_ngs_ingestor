import os
import shutil
from tempfile import mkdtemp
import atexit

import logging
logger = logging.getLogger()

# a module level list of temporary directories that have been
# created, so these can be cleaned up upon premature exit
TMPDIRS = []


def create_tmp_dir(*args, **kwargs):
    tmpdir = mkdtemp(*args, **kwargs)
    global TMPDIRS
    TMPDIRS.append(tmpdir)
    return tmpdir


@atexit.register
def _cleanup_tmp():
    for tmpdir in TMPDIRS:
        if os.exists(tmpdir):
            shutil.rmtree(tmpdir)
            logger.info("Removed temp directory: %s", tmpdir)
