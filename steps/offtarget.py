"""Shared off-target search, reused from the crispr-guide-design project.

This is the same 2-bit packed comparison described and tested there:
https://github.com/mhommii/crispr-guide-design (see steps/04_off_target_search.py
and tests/test_mismatch_counting.py). It is copied here rather than imported
so this repository runs on its own.

Each base becomes 2 bits (A=00, C=01, G=10, T=11), so a 20-base protospacer
is one 40-bit integer. Comparing a guide to a site is an XOR; folding each
2-bit pair to one bit and counting set bits gives the mismatch count for
every site at once.
"""

import gzip
import shutil
from pathlib import Path

import numpy as np

PROTOSPACER_LENGTH = 20
SEED_LENGTH = 12
CHUNK = 4_000_000

WEIGHTS = (4 ** np.arange(PROTOSPACER_LENGTH - 1, -1, -1)).astype(np.uint64)
SEED_MASK = np.uint64((1 << (2 * SEED_LENGTH)) - 1)
PAIR_MASK = np.uint64(0x5555555555555555)


def load_chromosome(archive: Path, unpacked: Path) -> np.ndarray:
    """Read a UCSC .fa.gz chromosome into an array of base codes (N = 255)."""
    if not unpacked.exists():
        unpacked.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(archive, "rb") as source, unpacked.open("wb") as target:
            shutil.copyfileobj(source, target)
    with unpacked.open("rb") as handle:
        handle.readline()  # drop the ">chrN" header
        sequence = handle.read().replace(b"\n", b"").upper()
    lookup = np.full(256, 255, dtype=np.uint8)
    for base, code in zip(b"ACGT", range(4)):
        lookup[base] = code
    return lookup[np.frombuffer(sequence, dtype=np.uint8)]


def _pack(starts: np.ndarray, source: np.ndarray) -> np.ndarray:
    packed = np.empty(len(starts), dtype=np.uint64)
    for begin in range(0, len(starts), CHUNK):
        block = starts[begin:begin + CHUNK]
        window = source[block[:, None] + np.arange(PROTOSPACER_LENGTH)]
        packed[begin:begin + CHUNK] = window.astype(np.uint64) @ WEIGHTS
    return packed


def _sites_on_strand(source: np.ndarray):
    # Element k of (is_g[1:-1] & is_g[2:]) is True when source[k+1] and
    # source[k+2] are both G, so the NGG PAM starts at k.
    is_g = source == 2
    pam_starts = np.flatnonzero(is_g[1:-1] & is_g[2:])
    pam_starts = pam_starts[pam_starts >= PROTOSPACER_LENGTH]
    starts = pam_starts - PROTOSPACER_LENGTH

    valid = np.ones(len(starts), dtype=bool)
    for begin in range(0, len(starts), CHUNK):
        block = starts[begin:begin + CHUNK]
        window = source[block[:, None] + np.arange(PROTOSPACER_LENGTH + 3)]
        valid[begin:begin + CHUNK] = ~(window == 255).any(axis=1)
    starts = starts[valid]
    return _pack(starts, source), starts


def index_chromosome(codes: np.ndarray):
    """Every NGG site on both strands, as packed 20-mers."""
    forward_packed, forward_starts = _sites_on_strand(codes)

    reverse_codes = 3 - codes[::-1]
    reverse_codes[codes[::-1] == 255] = 255
    reverse_packed, reverse_starts = _sites_on_strand(reverse_codes)

    packed = np.concatenate([forward_packed, reverse_packed])
    starts = np.concatenate([forward_starts, reverse_starts])
    strands = np.concatenate([
        np.zeros(len(forward_packed), dtype=np.uint8),
        np.ones(len(reverse_packed), dtype=np.uint8),
    ])
    return packed, starts, strands


def encode(protospacer: str) -> np.uint64:
    value = 0
    for base in protospacer:
        value = (value << 2) | "ACGT".index(base)
    return np.uint64(value)


def decode(value: int) -> str:
    return "".join("ACGT"[(value >> (2 * shift)) & 0b11]
                   for shift in range(PROTOSPACER_LENGTH - 1, -1, -1))


def count_mismatches(packed_sites: np.ndarray, protospacer: str):
    """Total and seed-region mismatches of one guide against every site."""
    difference = packed_sites ^ encode(protospacer)
    folded = (difference | (difference >> np.uint64(1))) & PAIR_MASK
    return np.bitwise_count(folded), folded
