import os


def is_in_tree(fpath, base):
    common = os.path.commonprefix([fpath, base])
    if not common or common == os.sep:
        return False
    return True


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]
