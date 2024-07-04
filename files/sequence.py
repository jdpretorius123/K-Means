"""
Sequence Class.

This module allows the user to store one DNA
sequence (exon, intron, intergenic) from a fasta 
file in a Sequence instance.

Classes
-------
Sequence

Attributes
----------
- seqStr: Sequence as string of nucleotides
- seqLst: list of nucleotides in Sequence
- id: Sequence id which can be one of exon, 
    intron, or intergenic
- kmers: list of kmer nucleotide subsequences
    in Sequence
- features: list of counts of all possible 
    kmers in Sequence

Methods
------
- seqStr(seq: str) -> None:
- seqList(seq: list[str]) -> None:
- id(id: str) -> None:
- kmers(kmers: list[str]) -> None:
- features(features: list[int]) -> None:
- toList() -> None:
- reverse() -> None:
- getBase(pos: int) -> str:
- getLength() -> int:
- buildKmers(size: int) -> None:
- print(rep: str) -> None:
"""


class Sequence:
    """A class to represent a sequence."""

    def __init__(self, seq: str, id: str) -> None:
        """
        Construct instance of Sequence.
        
        Parameters
        ----------
        seq: nucleotide sequence used to create Sequence
        id: Sequence id which can be one of exon, intron,
            or intergenic
        """
        self.seqStr = seq
        self.seqLst: list[str] = list(seq)
        self.id = id
        self.kmers: list[str] = list()
        self.features: list[int] = list()

    @property
    def seqStr(self) -> str:
        """
        Get or set nucelotides of Sequence.
        
        Changing the str sequence requires changing
        all other attributes.
        """
        return self._seqStr

    @seqStr.setter
    def seqStr(self, seq: str) -> None:
        if isinstance(seq, str):
            self._seqStr = seq
        else:
            raise ValueError('"seq" must be a str')

    @property
    def seqLst(self) -> list[str]:
        """
        Get of set list of nucleotides of Sequence.
        
        ''.join(self.seqList) must equal self.seqStr.
        Changing only one or the other is not acceptable.
        """
        return self._seqLst

    @seqLst.setter
    def seqLst(self, seq: list[str]) -> None:
        self._seqLst = seq

    @property
    def id(self) -> str:
        """
        Get or set Sequence id.
        
        Altering the id will make for incorrect
        clustering statistics.
        """
        return self._id

    @id.setter
    def id(self, id: str) -> None:
        self._id = id

    @property
    def kmers(self) -> list[str]:
        """
        Get or set Sequence kmers.
        
        Kmers are built from the original Sequence, so
        changing one requires changing both and possibly
        all attributes.
        """
        return self._kmers

    @kmers.setter
    def kmers(self, kmers: list[str]) -> None:
        self._kmers = kmers

    @property
    def features(self) -> list[int]:
        """
        Get or set Sequence feature vector.
        
        The feature vector is built from the Sequence
        kmers. Changing the feature vector requires
        changing both and possibly all attributes.
        """
        return self._features

    @features.setter
    def features(self, features: list[int]) -> None:
        self._features = features

    def _buildKmers(self, size: int) -> None:
        """
        Helper method for buildKmers(self, size: int) -> None.
        
        Parameters
        ----------
        size: kmer size
        """
        kmers: list[str] = list()
        length: int = self.getLength()
        stop: int = length - (size - 1)
        for i in range(stop):
            end: int = i + size
            kmer: str = self.seqStr[i:end]
            kmers.append(kmer)
        self.kmers = kmers

    def toList(self) -> None:
        """Convert Sequence to nucelotide list."""
        self.seqLst = list(self.seqStr)

    def reverse(self) -> None:
        """Reverse Sequence nucleotide list."""
        self.seqLst.reverse()

    def getBase(self, pos: int) -> str:
        """
        Return nucleotide at pos in Sequence.
        
        Parameters
        ----------
        pos: position of nucleotide to return

        Return
        ------
        base: nucleotide at pos in Sequence
        """
        length: int = self.getLength()
        if pos >= length:
            raise ValueError(f"{pos} greater than length {length}")
        base: str = self.seqStr[pos]
        return base

    def getLength(self) -> int:
        """
        Return Sequence length.
        
        Return
        ------
        length: number of nucleotides in Sequence
        """
        length: int = len(self.seqStr)
        return length

    def buildKmers(self, size: int) -> None:
        """
        Build Sequence kmers with length size. 
        
        Kmers are subsequences of Sequence. For example,
        consider the sequence ATTAG. Kmers of size 2 for 
        this sequence are AT, TT, TA, AG.
        
        Parameters
        ----------
        size: kmer size
        """
        self._buildKmers(size)

    def print(self, rep: str) -> None:
        """
        Print Sequence as string, list, or kmers.
        
        Parameters
        ----------
        rep: version of Sequence to print which can be
            one of string, list, and kmers
        """
        match rep:
            case "string":
                print(self.seqStr)
            case "list":
                print(" ".join(self.seqLst))
            case "kmers":
                print(" ".join(self.kmers))
            case _:
                raise ValueError(
                    f"{rep} must be one of 'string', 'list', or 'kmers'"
                    )
