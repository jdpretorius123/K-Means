"""
Kmeans Class.

This module allows the user to run
the K-Means Clustering algorithm.

Classes
-------
Kmeans
"""

from typing import TextIO
from cluster import Cluster
from kernel import Mismatch, Spectrum
from sequence import Sequence

STATS = dict[str, tuple[float, int]]


class Kmeans:
    """A class to represent Kmeans clustering."""

    def __init__(
        self,
        seqs: list[Sequence],
        kernel: Spectrum | Mismatch,
        k: int
    ) -> None:
        """Construct all attributes for Kmeans."""
        self.seqs = seqs
        self.kernel = kernel
        self.k = k
        self.clusters: list[Cluster] = list()

    @property
    def seqs(self) -> list[Sequence]:
        """All Sequences in infile."""
        return self._seqs

    @seqs.setter
    def seqs(self, seqs: list[Sequence]) -> None:
        self._seqs = seqs

    @property
    def kernel(self) -> Spectrum | Mismatch:
        """Kernel for K-Means algorithm."""
        return self._kernel

    @kernel.setter
    def kernel(self, kernel: Spectrum | Mismatch) -> None:
        self._kernel = kernel

    @property
    def k(self) -> int:
        """Number of cluseters."""
        return self._k

    @k.setter
    def k(self, k: int) -> None:
        self._k = k

    def _initSequences(self) -> None:
        """Initialize Sequences."""
        size: int = self.kernel.size
        for seq in self.seqs:
            seq.buildKmers(size)
            v: list[int] = self.kernel.featureVector(seq.kmers)
            seq.features = v

    def _initClusters(self) -> None:
        """Initialize cluster centroids."""
        size: int = len(self.kernel.kmers)
        for _ in range(self.k):
            cluster: Cluster = Cluster()
            cluster._init(size)
            self.clusters.append(cluster)

    def _init(self) -> None:
        """Initialize clusters and sequences."""
        self._initSequences()
        self._initClusters()

    def _computeDot(self, seq: Sequence) -> list[float]:
        """Compute all dot products between seq and centroids."""
        products: list[float] = list()
        for c in self.clusters:
            dot: float = self.kernel.dotProduct(seq.features, c.centroid)
            products.append(dot)
        return products

    def _chooseCluster(self, products: list[float]) -> int:
        """Choose cluster."""
        target: float = max(products)
        idx: int = 0
        for i in range(len(products)):
            if products[i] == target:
                idx = i
        return idx

    def _assignSeqs(self) -> None:
        """Assign Sequences to clusters."""
        for seq in self.seqs:
            products: list[float] = self._computeDot(seq)
            idx: int = self._chooseCluster(products)
            self.clusters[idx]._add(seq)

    def _updateClusters(self) -> None:
        """Update cluster log and centroids."""
        for c in self.clusters:
            c._record()
            c._update()

    def _checkClusters(self) -> bool:
        """Check clusters for changing centroids."""
        count: int = 0
        for c in self.clusters:
            if c._check():
                count += 1
        if count == self.k:
            return True
        return False

    def _clearClusters(self) -> None:
        """Clear clusters."""
        for c in self.clusters:
            c._clear()

    def execute(self, limit: int) -> None:
        """Execute K-Means Clustering algorithm."""
        run: int = 1
        self._init()
        while (not self._checkClusters()) and (run < limit):
            self._clearClusters()
            self._assignSeqs()
            self._updateClusters()
            run += 1

    def _printCluster(self, stats: STATS) -> None:
        """Print cluster."""
        keys: list[str] = list(stats.keys())
        for key in keys:
            percent: float = stats[key][0]
            count: int = stats[key][1]
            print(f"\t{key} = {percent} ({count})")

    def _printStats(self) -> None:
        """Print statistics."""
        for i in range(len(self.clusters)):
            print(f"Cluster {i + 1}")
            stats: STATS = self.clusters[i]._collectStats()
            self._printCluster(stats)

    def _printKernel(self) -> None:
        """Print kernel."""
        name: str = self.kernel.name
        size: int = self.kernel.size
        print(f"{name} (KMER={size}, CLUSTERS={self.k}):")

    def print(self) -> None:
        """Print K-Means algorithm results."""
        self._printKernel()
        self._printStats()

    def _writeCluster(self, path: str, i: int, stats: STATS) -> None:
        """Write cluster."""
        file: TextIO = open(path, "a")
        file.write(f"Cluster {i + 1}\n")
        keys: list[str] = list(stats.keys())
        for key in keys:
            percent: float = stats[key][0]
            count: int = stats[key][1]
            file.write(f"\t{key} = {percent} ({count})\n")
        file.close()

    def _writeStats(self, path: str) -> None:
        """Print statistics."""
        for i in range(len(self.clusters)):
            stats: STATS = self.clusters[i]._collectStats()
            self._writeCluster(path, i, stats)

    def _writeKernel(self, path: str) -> None:
        """Print kernel."""
        file: TextIO = open(path, "a")
        name: str = self.kernel.name
        size: int = self.kernel.size
        file.write(f"{name} (KMER={size}, CLUSTERS={self.k}):\n")
        file.close()

    def write(self, path: str) -> None:
        """Write K-Means algorithm results."""
        file: TextIO = open(path, "a")
        self._writeKernel(path)
        self._writeStats(path)
        file.write("\n")
        file.close()
