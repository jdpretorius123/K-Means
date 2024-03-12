"""
Cluster Class.

This module allows the user to create a
cluster for the K-Means Clustering Algorithm.

Classes
-------
Cluster
"""

import random
from sequence import Sequence


class Cluster:
    """A class to represent a K-Means Cluster."""

    def __init__(self) -> None:
        """Construct all attributes for Cluster."""
        self.centroid: list[float] = list()
        self.log: list[str] = list()
        self.members: list[Sequence] = list()

    @property
    def centroid(self) -> list[float]:
        """Cluster centroid."""
        return self._centroid

    @centroid.setter
    def centroid(self, centroid: list[float]) -> None:
        self._centroid = centroid

    @property
    def log(self) -> list[str]:
        """Past centroids."""
        return self._log

    @log.setter
    def log(self, log: list[str]) -> None:
        self._log = log

    @property
    def members(self) -> list[Sequence]:
        """Cluster members."""
        return self._members

    @members.setter
    def members(self, members: list[Sequence]) -> None:
        self._members = members

    def _init(self, size: int) -> None:
        """Initialize centroid."""
        choices: list[int] = list(range(size))
        self.centroid = [random.choice(choices) for _ in range(size)]

    def _add(self, member: Sequence) -> None:
        """Add member to cluster."""
        self.members.append(member)

    def _update(self) -> None:
        """Update cluster centroid."""
        count: int = len(self.members)
        if count != 0:
            total: float = 0
            for i in range(len(self.centroid)):
                for j in range(count):
                    total += self.members[j].features[i]
                self.centroid[i] = total / count
                total = 0

    def _record(self) -> None:
        """Record past centroid."""
        record: str = str(self.centroid)
        self.log.append(record)

    def _clear(self) -> None:
        """Clear membership."""
        self.members.clear()

    def _check(self) -> bool:
        """Check if past equals current."""
        if len(self.log) > 0:
            previous: str = self.log[-1]
            current: str = str(self.centroid)
            return current == previous
        return False

    def _calcStats(self, id: str) -> tuple[float, int]:
        """Calculate cluster membership statistics."""
        count: int = 0
        total: int = len(self.members)
        for member in self.members:
            if member.id == id:
                count += 1
        percent: float = 0
        if total != 0:
            percent = count / total
        percent = round(percent, 2)
        stats: tuple[float, int] = percent, count
        return stats

    def _collectStats(self) -> dict[str, tuple[float, int]]:
        """Collect cluster membership statistics."""
        idStats: list[tuple[float, int]] = list()
        ids: list[str] = ["intergenic", "intron", "exon"]
        for id in ids:
            stats: tuple[float, int] = self._calcStats(id)
            idStats.append(stats)
        return dict(zip(ids, idStats))

    def _printCentroid(self) -> None:
        """Print centroid."""
        centroid: list[str] = [str(round(val, 2)) for val in self.centroid]
        limit: int = 15
        for i in range(0, len(centroid), limit):
            stop: int = i + limit
            print(" ".join(centroid[i:stop]) + "\n")

    def _printMembers(self) -> None:
        """Print cluster members."""
        count: int = len(self.members)
        limit: int = 60
        if count != 0:
            for i in range(count):
                for j in range(0, len(self.members[i].seqStr), limit):
                    stop: int = j + limit
                    print(self.members[i].seqStr[j:stop])
                print()

    def print(self) -> None:
        """Print cluster."""
        self._printCentroid()
        print()
        self._printMembers()
