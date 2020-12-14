# pylint: disable=E1101,C0301
"""
Classes for handling FASTQ formatted sequences.

Reader -> Fastq
Writer <- Fastq
"""

import os
import sys
from six import PY3, PY2, string_types
from subprocess import Popen, PIPE
from io import TextIOWrapper
import os.path
from pkg_resources import get_distribution

if PY2:
    import string


class Gzip(object):
    """ Call system gzip and maintain interface compatibility with python
    gzip module. Tries bgzip first and falls back on system gzip. """

    def __init__(self, filename, mode, gzip='bgzip'):
        self.mode = mode
        self.name = filename
        name, ext = os.path.splitext(filename)
        self._out = None
        if name == '-' and 'w' in mode:
            self.stream = sys.stdout
            self.p = None
            self.write = self.stream.write
        elif name == '-' and 'a' in mode:
            self.stream = sys.stdout
            self.p = None
            self.write = self.stream.write
        elif name == '-' and 'r' in mode:
            self.stream = sys.stdin
            self.p = None
        elif ext == '.gz':
            try:
                self.stream, self.p = self.open(filename, mode, gzip)
            except OSError:
                self.stream, self.p = self.open(filename, mode, 'gzip')
        else:
            self.stream = open(filename, mode)
            self.p = None
            self.write = self.stream.write

    def __iter__(self):
        return self

    def __next__(self):
        return next(self)

    def __next__(self):
        return next(self.stream)

    def open(self, filename, mode, gzip):
        if 'r' in mode:
            p = Popen(['gzip', '-dc', filename], stdout=PIPE)
            if 'b' in mode:
                fh = p.stdout
            elif PY3:
                fh = TextIOWrapper(p.stdout)
            elif PY2:
                fh = p.stdout
        elif 'w' in mode:
            self._out = open(filename, 'wb', 0)
            p = Popen([gzip, '-c'], stdin=PIPE, stdout=self._out)
            fh = p.stdin
            if 'b' in mode and PY3:
                self.write = self._write_str
            elif PY3:
                self.write = self._write_bytes
            elif PY2:
                self.write = self._write_str
        elif 'a' in mode:
            self._out = open(filename, 'ab', 0)
            p = Popen([gzip, '-c'], stdin=PIPE, stdout=self._out)
            fh = p.stdin
            if 'b' in mode and PY3:
                self.write = self._write_str
            elif PY3:
                self.write = self._write_bytes
            elif PY2:
                self.write = self._write_str
        return (fh, p)

    def _write_str(self, string):
        self.stream.write(string)

    def _write_bytes(self, string):
        self.stream.write(string.encode())

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self.p:
            self.p.communicate()  # might be a safer way to flush PIPE
        if self._out:
            self._out.close()


class GzipWriter(Gzip):
    def __init__(self, filename):
        super(GzipWriter, self).__init__(filename, 'w')


class GzipReader(Gzip):
    def __init__(self, filename):
        super(GzipReader, self).__init__(filename, 'r')


class Fastq(object):
    """
    A class to hold features from fastq reads.
    """

    def __init__(self, name='-', seq='', strand='+', qual=None, comment=None):
        self.name = name
        self.seq = seq
        self.strand = strand
        if qual:
            self.qual = qual
        else:
            self.qual = '#' * len(seq)
        self.comment = comment
        self.i = int()
        assert isinstance(name, string_types)
        assert isinstance(seq, string_types)

    def __iter__(self):
        return self

    def __next__(self):
        if self.i < len(self):
            value, self.i = self[self.i], self.i + 1
            return value
        else:
            raise StopIteration()

    def __getitem__(self, key):
        return self.__class__(self.name, self.seq[key], self.strand,
                              self.qual[key], self.comment)

    def __next__(self):
        return next(self)

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.name[0] != '@':
            self.name = ''.join(['@', self.name])
        if self.comment:
            return '\n'.join([
                '{0} {1}'.format(self.name, self.comment), self.seq,
                self.strand, self.qual
            ]) + '\n'
        else:
            return '\n'.join([self.name, self.seq, self.strand, self.qual
                              ]) + '\n'

    def __len__(self):
        return len(self.seq)

    def reverse(self):
        """ Returns reverse ordered self
        >>> x = Fastq(name='test', seq='TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT',qual='#######EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22')
        >>> x.reverse()
        @test
        TATAGGTTTGTAGAATGCAAGAGTTATGGAGGTATTTTT
        +
        22=CFFACDD+AAF9<CHG@BFE+<A4?<4CE#######
        <BLANKLINE>
        """
        if self.comment:
            return self.__class__(self.name, self.seq[::-1], self.strand,
                                  self.qual[::-1], self.comment[::-1])
        else:
            return self.__class__(self.name, self.seq[::-1], self.strand,
                                  self.qual[::-1])

    def complement(self):
        """ Returns the compliment of self. This only affects the sequence slot.
        >>> x = Fastq(name='test', seq='TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT',qual='#######EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22')
        >>> x.complement()
        @test
        AAAAATACCTCCATAACTCTTGCATTCTACAAACCTATA
        +
        #######EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22
        <BLANKLINE>
        """
        if PY3:
            seqtable = str.maketrans('ACTGNRYSWKM', 'TGACNYRSWKM')
            contable = str.maketrans('ACTGNRYSWKMswkmXx', 'TGACNYRswkmSWKMxX')
        elif PY2:
            seqtable = string.maketrans('ACTGNRYSWKM', 'TGACNYRSWKM')
            contable = string.maketrans('ACTGNRYSWKMswkmXx',
                                        'TGACNYRswkmSWKMxX')
        if self.comment:
            return self.__class__(self.name,
                                  str(self.seq).translate(seqtable),
                                  self.strand, self.qual,
                                  self.comment.translate(contable))
        else:
            return self.__class__(self.name,
                                  str(self.seq).translate(seqtable),
                                  self.strand, self.qual)

    def revcomplement(self):
        """ Take the reverse compliment of self.
        >>> x = Fastq(name='test', seq='TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT',qual='#######EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22')
        >>> x.revcomplement()
        @test
        ATATCCAAACATCTTACGTTCTCAATACCTCCATAAAAA
        +
        22=CFFACDD+AAF9<CHG@BFE+<A4?<4CE#######
        <BLANKLINE>
        """
        return self.reverse().complement()

    def pair(self):
        """ Return 0 1 or 2 depending on whether read is 1 or 2 in a pair, or unpaired (0).

        """
        n = self.name[-2:]
        if n[0] != '/':
            return 0
        else:
            return int(n[1])

    @property
    def safename(self):
        """Return self.name without paired-end identifier if it exists"""
        if self.name[-2] == '/':
            return self.name[:-2]
        else:
            return self.name

    def gc(self):
        """ Return the GC content of self as an int
        >>> x = Fastq(name='test', seq='TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT',qual='#######EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22')
        >>> x.gc()
        30
        """
        g = self.seq.count('G')
        c = self.seq.count('C')
        return int((g + c) / len(self) * 100)


class Reader:
    """
    A class to read the name, sequence, strand and qualities from a fastq file
    """

    def __init__(self, f):
        try:
            name, ext = os.path.splitext(f)
            if ext == '.gz':
                self.file = Gzip(''.join([name, ext]), 'r')
            else:
                self.file = open(f)
        except AttributeError:
            self.file = open(f)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self)

    def __next__(self):
        try:
            name_data = next(self.file).strip()
            try:
                name, comment = name_data.split()
            except ValueError:
                name = name_data
                comment = None
            seq = next(self.file).strip()
            strand = next(self.file).strip()
            qual = next(self.file).strip()
            if comment:
                return Fastq(
                    name=name[1:],
                    seq=seq,
                    strand=strand,
                    qual=qual,
                    comment=comment)
            else:
                return Fastq(name=name[1:], seq=seq, strand=strand, qual=qual)
        except StopIteration:
            raise StopIteration

    def subsample(self, n):
        """ Draws every nth read from self. Returns Fastq. """
        n = n * 4
        for i, line in enumerate(self.file):
            if i % n == 0:
                name = line.strip().split()[0]
            elif i % n == 1:
                seq = line.strip()
            elif i % n == 2:
                strand = line.strip()
            elif i % n == 3:
                qual = line.strip()
                yield Fastq(name=name, seq=seq, strand=strand, qual=qual)

    def fileno(self):
        return self.file.fileno()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class Writer(object):
    """ Take a Fastq class object and file name, open file and write read """

    def __init__(self, f):
        try:
            name, ext = os.path.splitext(f.name)
            if ext == '.gz':
                self.file = Gzip(''.join([name, ext]), 'w')
            else:
                self.file = f
        except AttributeError:
            self.file = f

    def write(self, read):
        """ Write fastq with conversion string appended to read name """
        self.file.write(str(read))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


if __name__ == "__main__":
    import doctest
    doctest.testmod()
