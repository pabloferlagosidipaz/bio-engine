"""
Sequence Data Utilities
=======================

This module contains fast, low-level algorithms for DNA sequence manipulation,
including complement generation and IUPAC consensus merging.
"""

from itertools import zip_longest

# IUPAC DNA nucleotide codes and their complements
IUPAC_COMPLEMENT_MAP: dict[str, str] = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K',
    'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
    'N': 'N', '-': '-', ' ': ' ', '.': '-'
}

# Precompute translation table for string translation.
# Includes lowercase to uppercase complement mappings.
_TRANS_TABLE = str.maketrans(IUPAC_COMPLEMENT_MAP)
_TRANS_TABLE.update({ord(k.lower()): v for k, v in IUPAC_COMPLEMENT_MAP.items()})

IUPAC_TO_BASES: dict[str, frozenset[str]] = {
    'A': frozenset({'A'}), 'C': frozenset({'C'}), 'G': frozenset({'G'}), 'T': frozenset({'T'}),
    'R': frozenset({'A', 'G'}), 'Y': frozenset({'C', 'T'}), 'S': frozenset({'C', 'G'}), 'W': frozenset({'A', 'T'}),
    'K': frozenset({'G', 'T'}), 'M': frozenset({'A', 'C'}),
    'B': frozenset({'C', 'G', 'T'}), 'D': frozenset({'A', 'G', 'T'}), 'H': frozenset({'A', 'C', 'T'}), 'V': frozenset({'A', 'C', 'G'}),
    'N': frozenset({'A', 'C', 'G', 'T'}), '-': frozenset({'-'}), '.': frozenset({'-'})
}

# Reverse mapping: frozenset of bases -> IUPAC symbol
# Variations without gaps/Ns are included for robust lookup.
BASES_TO_IUPAC: dict[frozenset[str], str] = {}
for symbol, bases in IUPAC_TO_BASES.items():
    BASES_TO_IUPAC[bases] = symbol

def get_complement(sequence: str) -> str:
    """
    Returns the complementary sequence of a DNA string supporting IUPAC codes.
    Maps A-T, T-A, C-G, G-C, R-Y, etc.
    Returns uppercase sequence.
    """
    return sequence.translate(_TRANS_TABLE)

def get_reverse_complement(sequence: str) -> str:
    """
    Returns the reverse complementary sequence of a DNA string.
    """
    return sequence.translate(_TRANS_TABLE)[::-1]

def get_iupac_symbol(bases: set[str] | list[str] | str | frozenset[str]) -> str:
    """
    Returns the IUPAC code for a set of bases.
    Gaps are ignored if other bases are present.
    """
    if isinstance(bases, str):
        return bases.upper()

    # Standardize to uppercase and unique bases
    bases_set = {b.upper() for b in bases}

    # Remove gaps/N if others exist to find a more specific consensus
    if len(bases_set) > 1:
        bases_set.discard('-')
        bases_set.discard('.')

    if len(bases_set) > 1:
        bases_set.discard('N')

    if not bases_set:
        return '-'

    fs = frozenset(bases_set)
    return BASES_TO_IUPAC.get(fs, 'N')

def iupac_to_bases(symbol: str) -> frozenset[str]:
    """
    Returns the set of bases represented by an IUPAC symbol.
    """
    if not symbol:
        return frozenset()
    return IUPAC_TO_BASES.get(symbol.upper(), IUPAC_TO_BASES['N'])

def get_iupac_consensus(seq1: str, seq2: str) -> list[str]:
    """
    Generates an IUPAC consensus sequence from two aligned sequences.
    """
    if not seq1: return list(seq2.upper())
    if not seq2: return list(seq1.upper())


    return [
        b1.upper() if b1.upper() == b2.upper() else get_iupac_symbol({b1, b2})
        for b1, b2 in zip_longest(seq1, seq2, fillvalue='-')
    ]

def two_iupac_consensus(seq1: str, seq2: str) -> list[str]:
    """
    Generates an IUPAC consensus sequence from two aligned sequences, 
    properly merging IUPAC symbols.
    """
    if not seq1: return list(seq2.upper())
    if not seq2: return list(seq1.upper())

    consensus = []

    _to_bases = IUPAC_TO_BASES.get
    _n_set = IUPAC_TO_BASES['N']
    _get_symbol = BASES_TO_IUPAC.get

    for b1, b2 in zip_longest(seq1, seq2, fillvalue='-'):
        b1u, b2u = b1.upper(), b2.upper()
        if b1u == b2u:
            consensus.append(b1u)
        else:
            # Union of bases represented by each IUPAC symbol
            bases1 = _to_bases(b1u, _n_set)
            bases2 = _to_bases(b2u, _n_set)
            merged_bases = bases1 | bases2

            # Check precomputed mapping or fallback to resolving the symbol
            res = _get_symbol(merged_bases)
            if res:
                consensus.append(res)
            else:
                 consensus.append(get_iupac_symbol(merged_bases))

    return consensus
