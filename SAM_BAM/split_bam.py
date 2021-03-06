## category General
## desc Splits a BAM file into smaller pieces
"""
Splits a BAM file into smaller pieces

Given a BAM file, this script will split it into smaller BAM files with a
limit on the number of reads included.

Or it will also split a BAM file into a separate BAM file for each reference
that is included.
"""

import os
import sys
import pysam
from eta import ETA


def usage():
    print __doc__
    print """
Usage: ./split_bam {-n num | -ref} in.bam out_template_name

out_template_name will be the template for the smaller BAM files.  They will
be named "out_template_name.N.bam" where out_template_name is the given
argument and N is the file number.

Options:
    -n      The number of reads to include in sub-files
            (default: 1000000)

    -ref    Split by references
"""
    sys.exit(1)


def bam_iter(bam, quiet=False, show_ref_pos=False, callback=None):
    '''
    >>> [x.qname for x in bam_iter(bam_open(os.path.join(os.path.dirname(__file__), 't', 'test.bam')), quiet=True)]
    ['A', 'B', 'E', 'C', 'D', 'F', 'Z']
    '''
    if not quiet and bam.filename:
        eta = ETA(os.stat(bam.filename).st_size)
    else:
        eta = None

    if os.path.exists('%s.bai' % bam.filename):
        # This is an indexed file, so it is ref sorted...
        # Meaning that we should show chrom:pos, instead of read names
        show_ref_pos = True

    for read in bam:
        pos = bam.tell()
        bgz_offset = pos >> 16

        if not quiet and eta:
            if callback:
                eta.print_status(bgz_offset, extra=callback(read))
            elif (show_ref_pos):
                if read.tid > -1:
                    eta.print_status(bgz_offset, extra='%s:%s %s' % (bam.getrname(read.tid), read.pos, read.qname))
                else:
                    eta.print_status(bgz_offset, extra='unmapped %s' % (read.qname))
            else:
                eta.print_status(bgz_offset, extra='%s' % read.qname)
        yield read

    if eta:
        eta.done()

def bam_split(infile, out_template, read_count=1000000, reference=False, quiet=False):
    bamfile = pysam.Samfile(infile, "rb")
    outfile = None

    file_count = 0

    count = 0
    fname = ""
    lastref = -1
    for read in bam_iter(bamfile):
        if not outfile or (not reference and count >= read_count) or (reference and lastref != read.tid):
            if outfile:
                outfile.close()
            file_count += 1
            count = 0
            if reference:
                if read.tid >= 0:
                    fname = '%s.%s.bam' % (out_template, bamfile.getrname(read.tid))
                else:
                    fname = None
            else:
                fname = '%s.%s.bam' % (out_template, file_count)

            if fname:
                outfile = pysam.Samfile(fname, "wb", template=bamfile)
            else:
                outfile = None

        if outfile:
            outfile.write(read)
            count += 1

        lastref = read.tid

    bamfile.close()
    if outfile:
        outfile.close()
    if not quiet:
        sys.stderr.write("Split into %s files" % (file_count))


if __name__ == '__main__':
    infile = None
    outfile = None
    num = 1000000
    reference = False
    last = None

    for arg in sys.argv[1:]:
        if last == '-n':
            num = int(arg)
            last = None
        elif arg == '-ref':
            reference = True
        elif arg == '-h':
                usage()
        elif arg in ['-n']:
            last = arg
        elif not infile:
            infile = arg
        elif not outfile:
            outfile = arg

    if not infile or not outfile:
        usage()
    else:
        bam_split(infile, outfile, num, reference=reference)
