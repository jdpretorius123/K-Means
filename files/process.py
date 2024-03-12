"""
Write K-Means Output.

This module allows the user to setup
and run the K-Means algorithm.

Functions
---------
writeClusters(argv: list[str]) -> None:
    Writes clustering results.
"""

import os
from file import FastaFile
from sequence import Sequence
from kernel import Spectrum, Mismatch
from kmeans import Kmeans


def writeClusters(argv: list[str]) -> None:
    """Write clustering results."""
    infile: str = argv[0]
    ff: FastaFile = FastaFile(infile)
    seqs: list[Sequence] = ff.getSequences()
    size: int = int(argv[1])
    name: str = argv[2]
    mismatch: int = int(argv[3])
    kmers: list[str] = ff.findKmers(size)
    k: int = int(argv[4])

    if name == "SPECTRUM KERNEL":
        spec: Spectrum = Spectrum(kmers, size, name)
        kmeans: Kmeans = Kmeans(seqs, spec, k)
    else:
        mis: Mismatch = Mismatch(kmers, size, name, mismatch)
        kmeans = Kmeans(seqs, mis, k)

    limit: int = int(argv[5])
    kmeans.execute(limit)

    outfile: str = argv[6]
    kmeans.write(outfile)


def runKMeans() -> None:
    """Run K-Means algorithm."""
    infile: str = "kmeans/kmeans.fasta"
    sizes: list[str] = ["2", "6"]
    names: list[str] = ["SPECTRUM KERNEL", "MISMATCH KERNEL"]
    mismatch: str = "1"
    clusters: list[str] = ["2", "3", "5"]
    limit: str = "1000"
    outfile: str = "kmeans/kmeans.txt"

    if os.path.isfile(outfile):
        os.remove(outfile)

    argv: list[str] = list()
    for size in sizes:
        for name in names:
            for c in clusters:
                argv = [infile, size, name, mismatch, c, limit, outfile]
                writeClusters(argv)
