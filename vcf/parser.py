import codecs
import collections
import csv
import gzip
import itertools
import os
import re
import sys
from .utils import get_uncompressed_size


try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

try:
    import pysam
except ImportError:
    pysam = None

try:
    import cparse
except ImportError:
    cparse = None

from .model import _Call, _Record, make_calldata_tuple
from .model import _Substitution, _Breakend, _SingleBreakend, _SV


# Metadata parsers/constants
RESERVED_INFO = {
    "AA": "String",
    "AC": "Integer",
    "AF": "Float",
    "AN": "Integer",
    "BQ": "Float",
    "CIGAR": "String",
    "DB": "Flag",
    "DP": "Integer",
    "END": "Integer",
    "H2": "Flag",
    "H3": "Flag",
    "MQ": "Float",
    "MQ0": "Integer",
    "NS": "Integer",
    "SB": "String",
    "SOMATIC": "Flag",
    "VALIDATED": "Flag",
    "1000G": "Flag",
    # Keys used for structural variants
    "IMPRECISE": "Flag",
    "NOVEL": "Flag",
    "SVTYPE": "String",
    "SVLEN": "Integer",
    "CIPOS": "Integer",
    "CIEND": "Integer",
    "HOMLEN": "Integer",
    "HOMSEQ": "String",
    "BKPTID": "String",
    "MEINFO": "String",
    "METRANS": "String",
    "DGVID": "String",
    "DBVARID": "String",
    "DBRIPID": "String",
    "MATEID": "String",
    "PARID": "String",
    "EVENT": "String",
    "CILEN": "Integer",
    "DPADJ": "Integer",
    "CN": "Integer",
    "CNADJ": "Integer",
    "CICN": "Integer",
    "CICNADJ": "Integer",
}

RESERVED_FORMAT = {
    "GT": "String",
    "DP": "Integer",
    "FT": "String",
    "GL": "Float",
    "GLE": "String",
    "PL": "Integer",
    "GP": "Float",
    "GQ": "Integer",
    "HQ": "Integer",
    "PS": "Integer",
    "PQ": "Integer",
    "EC": "Integer",
    "MQ": "Integer",
    # Keys used for structural variants
    "CN": "Integer",
    "CNQ": "Float",
    "CNL": "Float",
    "NQ": "Integer",
    "HAP": "Integer",
    "AHAP": "Integer",
}

INTEGER = 0
STRING = 1
FLOAT = 2
FLAG = 3


def _encode_type(field_type):
    return {
        "Integer": INTEGER,
        "String": STRING,
        "Character": STRING,
        "Float": FLOAT,
        "Numeric": FLOAT,
        "Flag": FLAG,
    }[field_type]


RESERVED_INFO_CODES = {k: _encode_type(v) for k, v in RESERVED_INFO.items()}
RESERVED_FORMAT_CODES = {k: _encode_type(v) for k, v in RESERVED_FORMAT.items()}

# Spec is a bit weak on which metadata lines are singular, like fileformat
# and which can have repeats, like contig
SINGULAR_METADATA = ["fileformat", "fileDate", "reference"]

# Conversion between value in file and Python value
field_counts = {
    ".": None,  # Unknown number of values
    "A": -1,  # Equal to the number of alternate alleles in a given record
    "G": -2,  # Equal to the number of genotypes in a given record
    "R": -3,  # Equal to the number of alleles including reference in a given record
}


_Info = collections.namedtuple(
    "Info", ["id", "num", "type", "desc", "source", "version", "type_code"]
)
_Filter = collections.namedtuple("Filter", ["id", "desc"])
_Alt = collections.namedtuple("Alt", ["id", "desc"])
_Format = collections.namedtuple("Format", ["id", "num", "type", "desc", "type_code"])
_SampleInfo = collections.namedtuple(
    "SampleInfo", ["samples", "gt_bases", "gt_types", "gt_phases"]
)
_Contig = collections.namedtuple("Contig", ["id", "length"])


class _vcf_metadata_parser(object):
    """Parse the metadata in the header of a VCF file."""

    def __init__(self):
        super(_vcf_metadata_parser, self).__init__()
        self.info_pattern = re.compile(
            r"""\#\#INFO=<
            ID=(?P<id>[^,]+),\s*
            Number=(?P<number>-?\d+|\.|[AGR])?,\s*
            Type=(?P<type>Integer|Float|Flag|Character|String),\s*
            Description="(?P<desc>[^"]*)"
            (?:,\s*Source="(?P<source>[^"]*)")?
            (?:,\s*Version="?(?P<version>[^"]*)"?)?
            >""",
            re.VERBOSE,
        )
        self.filter_pattern = re.compile(
            r"""\#\#FILTER=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >""",
            re.VERBOSE,
        )
        self.alt_pattern = re.compile(
            r"""\#\#ALT=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >""",
            re.VERBOSE,
        )
        self.format_pattern = re.compile(
            r"""\#\#FORMAT=<
            ID=(?P<id>.+),\s*
            Number=(?P<number>-?\d+|\.|[AGR]),\s*
            Type=(?P<type>.+),\s*
            Description="(?P<desc>.*)"
            >""",
            re.VERBOSE,
        )
        self.contig_pattern = re.compile(
            r"""\#\#contig=<
            ID=(?P<id>[^>,]+)
            (,.*length=(?P<length>-?\d+))?
            .*
            >""",
            re.VERBOSE,
        )
        self.meta_pattern = re.compile(r"""##(?P<key>.+?)=(?P<val>.+)""")

    def vcf_field_count(self, num_str):
        """Cast vcf header numbers to integer or None"""
        if num_str is None:
            return None
        elif num_str not in field_counts:
            # Fixed, specified number
            return int(num_str)
        else:
            return field_counts[num_str]

    def read_info(self, info_string):
        """Read a meta-information INFO line."""
        match = self.info_pattern.match(info_string)
        if not match:
            raise SyntaxError("One of the INFO lines is malformed: %s" % info_string)

        num = self.vcf_field_count(match.group("number"))

        info = _Info(
            match.group("id"),
            num,
            match.group("type"),
            match.group("desc"),
            match.group("source"),
            match.group("version"),
            _encode_type(match.group("type")),
        )

        return (match.group("id"), info)

    def read_filter(self, filter_string):
        """Read a meta-information FILTER line."""
        match = self.filter_pattern.match(filter_string)
        if not match:
            raise SyntaxError(
                "One of the FILTER lines is malformed: %s" % filter_string
            )

        filt = _Filter(match.group("id"), match.group("desc"))

        return (match.group("id"), filt)

    def read_alt(self, alt_string):
        """Read a meta-information ALTline."""
        match = self.alt_pattern.match(alt_string)
        if not match:
            raise SyntaxError("One of the ALT lines is malformed: %s" % alt_string)

        alt = _Alt(match.group("id"), match.group("desc"))

        return (match.group("id"), alt)

    def read_format(self, format_string):
        """Read a meta-information FORMAT line."""
        match = self.format_pattern.match(format_string)
        if not match:
            raise SyntaxError(
                "One of the FORMAT lines is malformed: %s" % format_string
            )

        num = self.vcf_field_count(match.group("number"))

        form = _Format(
            match.group("id"),
            num,
            match.group("type"),
            match.group("desc"),
            _encode_type(match.group("type")),
        )

        return (match.group("id"), form)

    def read_contig(self, contig_string):
        """Read a meta-contigrmation INFO line."""
        match = self.contig_pattern.match(contig_string)
        if not match:
            raise SyntaxError(
                "One of the contig lines is malformed: %s" % contig_string
            )
        length = self.vcf_field_count(match.group("length"))
        contig = _Contig(match.group("id"), length)
        return (match.group("id"), contig)

    def read_meta_hash(self, meta_string):
        # assert re.match("##.+=<", meta_string)
        items = meta_string.split("=", 1)
        # Removing initial hash marks
        key = items[0].lstrip("#")
        # N.B., items can have quoted values, so cannot just split on comma
        val = OrderedDict()
        state = 0
        k = ""
        v = ""
        for c in items[1].strip("[<>]"):

            if state == 0:  # reading item key
                if c == "=":
                    state = 1  # end of key, start reading value
                else:
                    k += c  # extend key
            elif state == 1:  # reading item value
                if v == "" and c == '"':
                    v += c  # include quote mark in value
                    state = 2  # start reading quoted value
                elif c == ",":
                    val[k] = v  # store parsed item
                    state = 0  # read next key
                    k = ""
                    v = ""
                else:
                    v += c
            elif state == 2:  # reading quoted item value
                if c == '"':
                    v += c  # include quote mark in value
                    state = 1  # end quoting
                else:
                    v += c
        if k != "":
            val[k] = v
        return key, val

    def read_meta(self, meta_string):
        if re.match("##.+=<", meta_string):
            return self.read_meta_hash(meta_string)
        match = self.meta_pattern.match(meta_string)
        if not match:
            # Spec only allows key=value, but we try to be liberal and
            # interpret anything else as key=none (and all values are parsed
            # as strings).
            return meta_string.lstrip("#"), "none"
        return match.group("key"), match.group("val")


class Reader(object):
    """Reader for a VCF v 4.0 file, an iterator returning ``_Record objects``"""

    def __init__(
        self,
        fsock=None,
        filename=None,
        compressed=None,
        prepend_chr=False,
        strict_whitespace=False,
        encoding="ascii",
    ):
        """Create a new Reader for a VCF file.

        You must specify either fsock (stream) or filename.  Gzipped streams
        or files are attempted to be recogized by the file extension, or gzipped
        can be forced with ``compressed=True``

        'prepend_chr=True' will put 'chr' before all the CHROM values, useful
        for different sources.

        'strict_whitespace=True' will split records on tabs only (as with VCF
        spec) which allows you to parse files with spaces in the sample names.
        """
        super(Reader, self).__init__()

        if not (fsock or filename):
            raise Exception("You must provide at least fsock or filename")

        if fsock:
            self._reader = fsock
            if filename is None and hasattr(fsock, "name"):
                filename = fsock.name
                if compressed is None:
                    compressed = filename.endswith(".gz")
        elif filename:
            if compressed is None:
                compressed = filename.endswith(".gz")
            self._reader = open(filename, "rb" if compressed else "rt")
        self.filename = filename
        self._total_bytes = 0
        self._read_bytes = 0
        if compressed:
            self._reader = gzip.GzipFile(fileobj=self._reader)
            if sys.version > "3":
                self._reader = codecs.getreader(encoding)(self._reader)

        if strict_whitespace:
            self._separator = "\t"
        else:
            self._separator = "\t| +"

        self._row_pattern = re.compile(self._separator)
        self._alt_pattern = re.compile(r"[\[\]]")

        self.reader = (line.strip() for line in self._reader if line.strip())

        #: metadata fields from header (string or hash, depending)
        self.metadata = None
        #: INFO fields from header
        self.infos = None
        #: FILTER fields from header
        self.filters = None
        #: ALT fields from header
        self.alts = None
        #: FORMAT fields from header
        self.formats = None
        #: contig fields from header
        self.contigs = None
        self.samples = None
        self._sample_indexes = None
        self._header_lines = []
        self._column_headers = []
        self._tabix = None
        self._prepend_chr = prepend_chr
        self._parse_metainfo()
        self._format_cache = {}
        self.encoding = encoding

    def __iter__(self):
        self._read_bytes = 0
        return self

    def read_bytes(self):
        """Return read bytes from uncompress data. Usefull to have a progress bar"""
        return self._read_bytes

    def total_bytes(self):
        """Return total bytes from uncompress data. Usefull to have a progress bar"""

        if self._total_bytes == 0:
            # Compute only once. It can be time consumming
            self._total_bytes = get_uncompressed_size(self.filename)

        return self._total_bytes

    def _parse_metainfo(self):
        """Parse the information stored in the metainfo of the VCF.

        The end user shouldn't have to use this.  She can access the metainfo
        directly with ``self.metadata``."""
        for attr in ("metadata", "infos", "filters", "alts", "contigs", "formats"):
            setattr(self, attr, OrderedDict())

        parser = _vcf_metadata_parser()

        line = next(self.reader)
        while line.startswith("##"):
            self._header_lines.append(line)

            if line.startswith("##INFO"):
                key, val = parser.read_info(line)
                self.infos[key] = val

            elif line.startswith("##FILTER"):
                key, val = parser.read_filter(line)
                self.filters[key] = val

            elif line.startswith("##ALT"):
                key, val = parser.read_alt(line)
                self.alts[key] = val

            elif line.startswith("##FORMAT"):
                key, val = parser.read_format(line)
                self.formats[key] = val

            elif line.startswith("##contig"):
                key, val = parser.read_contig(line)
                self.contigs[key] = val

            else:
                key, val = parser.read_meta(line)
                if key in SINGULAR_METADATA:
                    self.metadata[key] = val
                else:
                    if key not in self.metadata:
                        self.metadata[key] = []
                    self.metadata[key].append(val)

            line = next(self.reader)

        fields = self._row_pattern.split(line[1:])
        self._column_headers = fields[:9]
        self.samples = fields[9:]
        self._sample_indexes = dict([(x, i) for (i, x) in enumerate(self.samples)])

    def _map(self, func, iterable, bad=[".", "", "NA"]):
        """``map``, but make bad values None."""
        return [func(x) if x not in bad else None for x in iterable]

    def _parse_filter(self, filt_str):
        """Parse the FILTER field of a VCF entry into a Python list

        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent
        """
        if filt_str == ".":
            return None
        elif filt_str == "PASS":
            return []
        else:
            return filt_str.split(";")

    def _parse_info(self, info_str):
        """Parse the INFO field of a VCF entry into a dictionary of Python
        types.

        """
        if info_str == ".":
            return {}

        entries = info_str.split(";")
        retdict = {}

        for entry in entries:
            entry = entry.split("=", 1)
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type_code
            except KeyError:
                try:
                    entry_type = RESERVED_INFO_CODES[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = STRING
                    else:
                        entry_type = FLAG

            if entry_type == INTEGER:
                vals = entry[1].split(",")
                try:
                    val = self._map(int, vals)
                # Allow specified integers to be flexibly parsed as floats.
                # Handles cases with incorrectly specified header types.
                except ValueError:
                    val = self._map(float, vals)
            elif entry_type == FLOAT:
                vals = entry[1].split(",")
                val = self._map(float, vals)
            elif entry_type == FLAG:
                val = True
            elif entry_type == STRING:
                try:
                    vals = entry[1].split(
                        ","
                    )  # commas are reserved characters indicating multiple values
                    val = self._map(str, vals)
                except IndexError:
                    entry_type = FLAG
                    val = True

            try:
                if self.infos[ID].num == 1 and entry_type != FLAG:
                    val = val[0]
            except KeyError:
                pass

            retdict[ID] = val

        return retdict

    def _parse_sample_format(self, samp_fmt):
        """Parse the format of the calls in this _Record"""
        samp_fmt = make_calldata_tuple(samp_fmt.split(":"))

        for fmt in samp_fmt._fields:
            try:
                entry_type = self.formats[fmt].type_code
                entry_num = self.formats[fmt].num
            except KeyError:
                entry_num = None
                try:
                    entry_type = RESERVED_FORMAT_CODES[fmt]
                except KeyError:
                    entry_type = STRING
            samp_fmt._types.append(entry_type)
            samp_fmt._nums.append(entry_num)
        return samp_fmt

    def _parse_samples(self, samples, samp_fmt, site):
        """Parse a sample entry according to the format specified in the FORMAT
        column.

        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent
        """

        # check whether we already know how to parse this format
        if samp_fmt not in self._format_cache:
            self._format_cache[samp_fmt] = self._parse_sample_format(samp_fmt)
        samp_fmt = self._format_cache[samp_fmt]

        if cparse:
            return cparse.parse_samples(
                self.samples, samples, samp_fmt, samp_fmt._types, samp_fmt._nums, site
            )

        samp_data = []
        _map = self._map

        nfields = len(samp_fmt._fields)

        for name, sample in zip(self.samples, samples):

            # parse the data for this sample
            sampdat = [None] * nfields

            for i, vals in enumerate(sample.split(":")):

                # short circuit the most common
                if samp_fmt._fields[i] == "GT":
                    sampdat[i] = vals
                    continue
                # genotype filters are a special case
                elif samp_fmt._fields[i] == "FT":
                    sampdat[i] = self._parse_filter(vals)
                    continue
                elif not vals or vals == ".":
                    sampdat[i] = None
                    continue

                entry_num = samp_fmt._nums[i]
                entry_type = samp_fmt._types[i]

                # we don't need to split single entries
                if entry_num == 1:
                    if entry_type == INTEGER:
                        try:
                            sampdat[i] = int(vals)
                        except ValueError:
                            sampdat[i] = float(vals)
                    elif entry_type == FLOAT:
                        sampdat[i] = float(vals)
                    else:
                        sampdat[i] = vals
                    continue

                vals = vals.split(",")
                if entry_type == INTEGER:
                    try:
                        sampdat[i] = _map(int, vals)
                    except ValueError:
                        sampdat[i] = _map(float, vals)
                elif entry_type == FLOAT:
                    sampdat[i] = _map(float, vals)
                else:
                    sampdat[i] = vals

            # create a call object
            call = _Call(site, name, samp_fmt(*sampdat))
            samp_data.append(call)

        return samp_data

    def _parse_alt(self, str):
        if self._alt_pattern.search(str) is not None:
            # Paired breakend
            items = self._alt_pattern.split(str)
            remoteCoords = items[1].split(":")
            chr = remoteCoords[0]
            if chr[0] == "<":
                chr = chr[1:-1]
                withinMainAssembly = False
            else:
                withinMainAssembly = True
            pos = remoteCoords[1]
            orientation = str[0] == "[" or str[0] == "]"
            remoteOrientation = re.search(r"\[", str) is not None
            if orientation:
                connectingSequence = items[2]
            else:
                connectingSequence = items[0]
            return _Breakend(
                chr,
                pos,
                orientation,
                remoteOrientation,
                connectingSequence,
                withinMainAssembly,
            )
        elif str[0] == "." and len(str) > 1:
            return _SingleBreakend(True, str[1:])
        elif str[-1] == "." and len(str) > 1:
            return _SingleBreakend(False, str[:-1])
        elif str[0] == "<" and str[-1] == ">":
            return _SV(str[1:-1])
        else:
            return _Substitution(str)

    def __next__(self):
        """Return the next record in the file."""
        line = next(self.reader)
        self._read_bytes += len(line)
        row = self._row_pattern.split(line.rstrip())
        chrom = row[0]
        if self._prepend_chr:
            chrom = "chr" + chrom
        pos = int(row[1])

        if row[2] != ".":
            ID = row[2]
        else:
            ID = None

        ref = row[3]
        alt = self._map(self._parse_alt, row[4].split(","))

        try:
            qual = int(row[5])
        except ValueError:
            try:
                qual = float(row[5])
            except ValueError:
                qual = None

        filt = self._parse_filter(row[6])
        if cparse is not None:
            info = cparse.parse_info(row[7], self.infos, RESERVED_INFO_CODES)
        else:
            info = self._parse_info(row[7])

        try:
            fmt = row[8]
        except IndexError:
            fmt = None
        else:
            if fmt == ".":
                fmt = None

        record = _Record(
            chrom, pos, ID, ref, alt, qual, filt, info, fmt, self._sample_indexes
        )

        if fmt is not None:
            samples = self._parse_samples(row[9:], fmt, record)
            record.samples = samples

        return record

    def fetch(self, chrom, start=None, end=None):
        """Fetches records from a tabix-indexed VCF file and returns an
        iterable of ``_Record`` instances

        chrom must be specified.

        The start and end coordinates are in the zero-based,
        half-open coordinate system, similar to ``_Record.start`` and
        ``_Record.end``. The very first base of a chromosome is
        index 0, and the the region includes bases up to, but not
        including the base at the end coordinate. For example
        ``fetch('4', 10, 20)`` would include all variants
        overlapping a 10 base pair region from the 11th base of
        through the 20th base (which is at index 19) of chromosome
        4. It would not include the 21st base (at index 20). See
        http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms
        for more information on the zero-based, half-open coordinate
        system.

        If end is omitted, all variants from start until the end of
        the chromosome chrom will be included.

        If start and end are omitted, all variants on chrom will be
        returned.

        requires pysam

        """
        if not pysam:
            raise Exception('pysam not available, try "pip install pysam"?')
        if not self.filename:
            raise Exception('Please provide a filename (or a "normal" fsock)')

        if not self._tabix:
            self._tabix = pysam.Tabixfile(self.filename, encoding=self.encoding)

        if self._prepend_chr and chrom[:3] == "chr":
            chrom = chrom[3:]

        self.reader = self._tabix.fetch(chrom, start, end)
        return self


class Writer(object):
    """VCF Writer. On Windows Python 2, open stream with 'wb'."""

    # Reverse keys and values in header field count dictionary
    counts = dict((v, k) for k, v in field_counts.items())

    def __init__(self, stream, template, lineterminator="\n"):
        self.writer = csv.writer(
            stream,
            delimiter="\t",
            lineterminator=lineterminator,
            quoting=csv.QUOTE_NONE,
        )
        self.template = template
        self.stream = stream

        # Order keys for INFO fields defined in the header (undefined fields
        # get a maximum key).
        self.info_order = collections.defaultdict(
            lambda: len(template.infos),
            dict(list(zip(iter(template.infos.keys()), itertools.count()))),
        )

        two = '##{key}=<ID={0},Description="{1}">\n'
        four = '##{key}=<ID={0},Number={num},Type={2},Description="{3}">\n'
        _num = self._fix_field_count
        for (key, vals) in template.metadata.items():
            if key in SINGULAR_METADATA:
                vals = [vals]
            for val in vals:
                if isinstance(val, dict):
                    values = ",".join(
                        "{0}={1}".format(key, value) for key, value in list(val.items())
                    )
                    stream.write("##{0}=<{1}>\n".format(key, values))
                else:
                    stream.write("##{0}={1}\n".format(key, val))
        for line in template.infos.values():
            stream.write(four.format(key="INFO", *line, num=_num(line.num)))
        for line in template.formats.values():
            stream.write(four.format(key="FORMAT", *line, num=_num(line.num)))
        for line in template.filters.values():
            stream.write(two.format(key="FILTER", *line))
        for line in template.alts.values():
            stream.write(two.format(key="ALT", *line))
        for line in template.contigs.values():
            if line.length:
                stream.write("##contig=<ID={0},length={1}>\n".format(*line))
            else:
                stream.write("##contig=<ID={0}>\n".format(*line))

        self._write_header()

    def _write_header(self):
        # TODO: write INFO, etc
        self.stream.write(
            "#"
            + "\t".join(self.template._column_headers + self.template.samples)
            + "\n"
        )

    def write_record(self, record):
        """write a record to the file"""
        ffs = self._map(str, [record.CHROM, record.POS, record.ID, record.REF]) + [
            self._format_alt(record.ALT),
            record.QUAL or ".",
            self._format_filter(record.FILTER),
            self._format_info(record.INFO),
        ]
        if record.FORMAT:
            ffs.append(record.FORMAT)

        samples = [
            self._format_sample(record.FORMAT, sample) for sample in record.samples
        ]
        self.writer.writerow(ffs + samples)

    def flush(self):
        """Flush the writer"""
        try:
            self.stream.flush()
        except AttributeError:
            pass

    def close(self):
        """Close the writer"""
        try:
            self.stream.close()
        except AttributeError:
            pass

    def _fix_field_count(self, num_str):
        """Restore header number to original state"""
        if num_str not in self.counts:
            return num_str
        else:
            return self.counts[num_str]

    def _format_alt(self, alt):
        return ",".join(self._map(str, alt))

    def _format_filter(self, flt):
        if flt == []:
            return "PASS"
        return self._stringify(flt, none=".", delim=";")

    def _format_info(self, info):
        if cparse:
            return cparse.format_info(info, self.info_order)

        if not info:
            return "."

        def order_key(field):
            # Order by header definition first, alphabetically second.
            return self.info_order[field], field

        return ";".join(
            self._stringify_pair(f, info[f]) for f in sorted(info, key=order_key)
        )

    def _format_sample(self, fmt, sample):
        if cparse:
            return cparse.format_sample(fmt, sample)

        if hasattr(sample.data, "GT"):
            gt = sample.data.GT
        else:
            gt = "./." if "GT" in fmt else ""

        result = [gt] if gt else []
        # Following the VCF spec, GT is always the first item whenever it is present.
        for field in sample.data._fields:
            value = getattr(sample.data, field)
            if field == "GT":
                continue
            if field == "FT":
                result.append(self._format_filter(value))
            else:
                result.append(self._stringify(value))
        return ":".join(result)

    def _stringify(self, x, none=".", delim=","):
        if type(x) == type([]):
            return delim.join(self._map(str, x, none))
        return str(x) if x is not None else none

    def _stringify_pair(self, x, y, none=".", delim=","):
        if isinstance(y, bool):
            return str(x) if y else ""
        return "%s=%s" % (str(x), self._stringify(y, none=none, delim=delim))

    def _map(self, func, iterable, none="."):
        """``map``, but make None values none."""
        return [func(x) if x is not None else none for x in iterable]


def __update_readme():
    import sys, vcf

    file("README.rst", "w").write(vcf.__doc__)


# backwards compatibility
VCFReader = Reader
VCFWriter = Writer
