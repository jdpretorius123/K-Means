"""
Kernel Class.

This module allows the user to apply a
Spectrum and/or Mismatch kernel to Kmers.

Classes
-------
Kernel
Spectrum
Mismatch
"""

class Kernel:
    """A class to represent Kernel for Sequence comparison."""

    def __init__(self, kmers: list[str], size: int, name: str) -> None:
        """Construct all attributes for Kernel."""
        self.kmers = kmers
        self.size = size
        self.name = name

    @property
    def kmers(self) -> list[str]:
        """All unique kmers."""
        return self._kmers

    @kmers.setter
    def kmers(self, kmers: list[str]) -> None:
        self._kmers = kmers

    @property
    def size(self) -> int:
        """Size of kmer."""
        return self._size

    @size.setter
    def size(self, size: int) -> None:
        self._size = size

    @property
    def name(self) -> str:
        """Name of Kernel."""
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        self._name = name

    def featureVector(self, kmers: list[str]) -> list[int]:
        """Compute feature vector."""
        raise NotImplementedError(
            "featureVector not defined for parent class Kernel."
        )

    def featureDict(self, v: list[int]) -> dict[str, int]:
        """Create feature dictionary."""
        featureDict: dict[str, int] = dict(zip(self.kmers, v))
        return featureDict

    def dotProduct(self, v1: list[int], v2: list[float]) -> float:
        """Compute dot product between feature vectors."""
        dot: float = 0
        for i in range(len(v1)):
            dot += v1[i] * v2[i]
        return dot


class Spectrum(Kernel):
    """A class to represent a Spectrum kernel."""

    def __init__(self, kmers: list[str], size: int, name: str) -> None:
        """Construct all attributes for Spectrum."""
        super().__init__(kmers, size, name)

    def featureVector(self, kmers: list[str]) -> list[int]:
        """Compute feature vector."""
        count: int = 0
        v: list[int] = list()
        for i in range(len(self.kmers)):
            for j in range(len(kmers)):
                if kmers[j] == self.kmers[i]:
                    count += 1
            v.append(count)
            count = 0
        return v


class Mismatch(Kernel):
    """A class to represent a Mismatch kernel."""

    def __init__(
        self, kmers: list[str], size: int, name: str, mismatch: int
    ) -> None:
        """Construct all attributes for Mismatch."""
        super().__init__(kmers, size, name)
        self.mismatch = mismatch

    @property
    def mismatch(self) -> int:
        """Number of allowed mismatches."""
        return self._mismatch

    @mismatch.setter
    def mismatch(self, mismatch: int) -> None:
        self._mismatch = mismatch

    def _numMismatch(self, ref: str, cand: str) -> int:
        """Count number of mismatches between reference and candidate."""
        count: int = 0
        for i in range(len(ref)):
            if ref[i] != cand[i]:
                count += 1
        return count

    def featureVector(self, kmers: list[str]) -> list[int]:
        """Compute feature vectors for Sequences."""
        count: int = 0
        v: list[int] = list()
        for i in range(len(self.kmers)):
            for j in range(len(kmers)):
                if self._numMismatch(kmers[j], self.kmers[i]) <= self.mismatch:
                    count += 1
            v.append(count)
            count = 0
        return v
