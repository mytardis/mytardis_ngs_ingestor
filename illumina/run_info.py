import csv


def parse_samplesheet(file_path, standardize_keys=True):

    # Old plain CSV format, IEM v3:
    # FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,
    # Operator,SampleProject
    #
    # last line: #_IEMVERSION_3_TruSeq LT,,,,,,,,,

    # Newer INI-style CSV format, IEM v4:
    # [Header],,,,,,,,
    # IEMFileVersion,4,,,,,,,
    # .. >snip< ..
    # Assay,TruSeq LT,,,,,,,
    # [Reads],,,,,,,,
    # .. >snip< ..
    # [Settings],,,,,,,,
    # .. >snip< ..
    # [Data],,,,,,,,
    # Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
    # .. >snip< ..
    #

    lines = []
    with open(file_path, "rU") as f:
        lines = f.readlines()

    chemistry = None

    # IEM v4 INI-style CSV
    if '[Header]' in lines[0]:
        section = None
        for i, l in enumerate(lines):
            if l[0] == '[':
                section = l[1:].split(']')[0]
            if section == 'Header' and l.startswith('Assay,'):
                chemistry = l.split(',')[1].strip()
            if section == 'Data':
                data_index = i
                break

        reader = csv.DictReader(lines[data_index+1:])

        if standardize_keys:
            samples = []
            # remove any underscores to make names consistent between
            # old and new style samplesheets (eg Sample_ID -> SampleID)
            for r in [row for row in reader]:
                r = {k.replace('_', ''): r[k] for k in r.keys()}
                samples.append(r)
        else:
            samples = [row for row in reader]

        return samples, chemistry
    else:  # Plain CSV (IEM v3 ?)
        reader = csv.DictReader(lines)
        samples = [row for row in reader]
        lastlinebit = samples[-1:][0].get('FCID', None)
        if lastlinebit is not None:
            chemistry = lastlinebit.split('_')[-1].strip()
        del samples[-1:]
        return samples, chemistry
