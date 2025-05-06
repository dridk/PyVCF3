#!/usr/bin/env python
"""
A VCFv4.0 and 4.1 parser for Python.

Online version of PyVCF documentation is available at http://pyvcf.rtfd.org/
"""


from vcf.parser import Reader, Writer
from vcf.parser import VCFReader, VCFWriter
from vcf.filters import Base as Filter
from vcf.parser import RESERVED_INFO, RESERVED_FORMAT
from vcf.sample_filter import SampleFilter

VERSION = '1.0.4'
