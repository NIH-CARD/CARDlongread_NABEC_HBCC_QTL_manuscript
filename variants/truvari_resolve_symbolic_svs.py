
#!/usr/bin/env python3

"""
Given a VCF, fill in the <DEL>, <INV>, <DUP> ALT sequence
"""

import sys
import pysam
import truvari

MAX_SV = 100_000_000 # Filter things smaller than this
# MAX_SV = 100000000

RC = str.maketrans("ATCG", "TAGC")
def do_rc(s):
    """
    Reverse complement a sequence
    """
    return s.translate(RC)[::-1]

def resolve(entry, ref):
    """
    """
    if entry.start > ref.get_reference_length(entry.chrom):
        return entry
    if entry.alts[0] in ['<CNV>', '<INS>']:
        return entry

    seq = ref.fetch(entry.chrom, entry.start, entry.stop)
    if entry.alts[0] == '<DEL>':
        entry.ref = seq
        entry.alts = [seq[0]]
    elif entry.alts[0] == '<INV>':
        entry.ref = seq
        entry.alts = [do_rc(seq)]
    elif entry.alts[0] == '<DUP>':
        entry.info['SVTYPE'] = 'INS'
        entry.ref = seq[0]
        entry.alts = [seq]
        entry.stop = entry.start + 1

    return entry

if __name__ == '__main__':
    vcf = pysam.VariantFile(sys.argv[1])
    ref = pysam.FastaFile(sys.argv[2])
    n_header = vcf.header.copy()
    
    out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
    for entry in vcf:
        if truvari.entry_size(entry) >= MAX_SV:
            continue
        if entry.alts[0].startswith("<"):
            entry = resolve(entry, ref)
        if entry.qual is None: 
            entry.qual = 0
        try:
            out.write(entry)
        except Exception:
            sys.stderr.write(f"{entry}\n{type(entry)}\n")
