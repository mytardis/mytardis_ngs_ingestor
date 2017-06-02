import os


def is_in_tree(fpath, base):
    common = os.path.commonprefix([fpath, base])
    if not common or common == os.sep:
        return False
    return True
