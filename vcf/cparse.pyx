from model import _Call

cdef int INTEGER = 0
cdef int STRING = 1
cdef int FLOAT = 2
cdef int FLAG = 3

cdef list _map(func, iterable, bad=['.', '']):
    '''``map``, but make bad values None.'''
    return [func(x) if x not in bad else None
            for x in iterable]

cdef _parse_filter(str filt_str):
    '''Parse the FILTER field of a VCF entry into a Python list

    NOTE: this method has a python equivalent and care must be taken
    to keep the two methods equivalent
    '''
    if filt_str == '.':
        return None
    elif filt_str == 'PASS':
        return []
    else:
        return filt_str.split(';')

def parse_samples(
        list names, list samples, samp_fmt,
        list samp_fmt_types, list samp_fmt_nums, site):

    sample_data = {}
    cdef char *name
    cdef char *fmt
    cdef char *sample
    cdef int entry_type
    cdef int i, j
    # cdef list samp_data = []
    cdef list sampvals
    n_samples = len(samples)
    n_formats = len(samp_fmt._fields)

    for i in range(n_samples):
        name = names[i]
        sample = samples[i]

        # parse the data for this sample
        sampdat = [None] * n_formats

        sampvals = sample.split(':')

        for j in range(n_formats):
            if j >= len(sampvals):
                break
            vals = sampvals[j]

            # short circuit the most common
            if samp_fmt._fields[j] == 'GT':
                sampdat[j] = vals
                continue
            # genotype filters are a special case
            elif samp_fmt._fields[j] == 'FT':
                sampdat[j] = _parse_filter(vals)
                continue
            elif not vals or vals == '.':
                sampdat[j] = None
                continue

            entry_type = samp_fmt_types[j]
            # TODO: entry_num is None for unbounded lists
            entry_num = samp_fmt_nums[j]

            # we don't need to split single entries
            if entry_num == 1:
                if entry_type == INTEGER:
                    try:
                        sampdat[j] = int(vals)
                    except ValueError:
                        sampdat[j] = float(vals)
                elif entry_type == FLOAT:
                    sampdat[j] = float(vals)
                else:
                    sampdat[j] = vals
                continue

            vals = vals.split(',')
            if entry_type == INTEGER:
                try:
                    sampdat[j] = _map(int, vals)
                except ValueError:
                    sampdat[j] = _map(float, vals)
            elif entry_type == FLOAT:
                sampdat[j] = _map(float, vals)
            else:
                sampdat[j] = vals

        # create a call object
        call = _Call(site, name, samp_fmt(*sampdat))
        samp_data[name] = call

    return samp_data

def parse_info(str info_str, infos, dict reserved_info_codes):
    '''Parse the INFO field of a VCF entry into a dictionary of Python
    types.

    '''
    if info_str == '.':
        return {}

    cdef list entries = info_str.split(';')
    cdef list vals
    cdef dict retdict = {}
    cdef int entry_type
    cdef str entry
    cdef str entry_key
    cdef str entry_val

    for entry in entries:
        entry_key, _, entry_val = entry.partition('=')
        try:
            entry_type = infos[entry_key].type_code
        except KeyError:
            try:
                entry_type = reserved_info_codes[entry_key]
            except KeyError:
                if entry_val:
                    entry_type = STRING
                else:
                    entry_type = FLAG

        if entry_type == INTEGER:
            vals = entry_val.split(',')
            try:
                retdict[entry_key] = _map(int, vals)
            # Allow specified integers to be flexibly parsed as floats.
            # Handles cases with incorrectly specified header types.
            except ValueError:
                retdict[entry_key] = _map(float, vals)
        elif entry_type == FLOAT:
            vals = entry_val.split(',')
            retdict[entry_key] = _map(float, vals)
        elif entry_type == FLAG:
            retdict[entry_key] = True
        elif entry_type == STRING:
            try:
                vals = entry_val.split(',') # commas are reserved characters indicating multiple values
                retdict[entry_key] = _map(str, vals)
            except AttributeError:
                entry_type = FLAG
                retdict[entry_key] = True

        try:
            if infos[entry_key].num == 1 and entry_type != FLAG:
                retdict[entry_key] = retdict[entry_key][0]
        except KeyError:
            pass

    return retdict

def format_info(dict info, info_order):
    if not info:
        return '.'
    def order_key(str field):
        # Order by header definition first, alphabetically second.
        return info_order[field], field
    return ';'.join(_stringify_pair(f, info[f]) for f in
                    sorted(info, key=order_key))

def format_sample(str fmt, sample):
    cdef str gt
    cdef list result

    if hasattr(sample.data, 'GT'):
        gt = sample.data.GT
    else:
        gt = './.' if 'GT' in fmt else ''

    result = [gt] if gt else []
    # Following the VCF spec, GT is always the first item whenever it is present.
    for field in sample.data._fields:
        value = getattr(sample.data, field)
        if field == 'GT':
            continue
        if field == 'FT':
            result.append(_format_filter(value))
        else:
            result.append(_stringify(value))
    return ':'.join(result)

cdef str _format_filter(flt):
    if flt == []:
        return 'PASS'
    return _stringify(flt, none='.', delim=';')

cdef str _stringify(x, none='.', delim=','):
    if type(x) == type([]):
        return delim.join(_write_map(str, x, none))
    return str(x) if x is not None else none

cdef str _stringify_pair(x, y, none='.', delim=','):
    if isinstance(y, bool):
        return str(x) if y else ""
    return "%s=%s" % (str(x), _stringify(y, none=none, delim=delim))

cdef list _write_map(func, iterable, none='.'):
    '''``map``, but make None values none.'''
    return [func(x) if x is not None else none
            for x in iterable]
