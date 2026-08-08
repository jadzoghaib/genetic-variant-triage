"""Targeted ClinVar region fetch via the tabix index + HTTP range requests.

NCBI serves the 193 MB ClinVar VCF at ~100 KB/s (per-IP throttled; parallel
connections do not help). A full download is ~18-32 minutes.

The VCF is bgzip-compressed with a tabix (.tbi) index, and the server honours
byte ranges. So instead of downloading the bulk file we:

  1. fetch the 610 KB .tbi once,
  2. read its linear index to find the compressed byte offset of a genomic
     region,
  3. range-request only those bytes,
  4. decompress the BGZF members and parse the records.

That turns a 193 MB download into a few hundred KB per gene, and removes the
bulk download from the product architecture entirely — ClinVar becomes an
on-demand per-target fetch, exactly like AlphaFold.

Tabix format: https://samtools.github.io/hts-specs/tabix.pdf
"""

from __future__ import annotations

import gzip
import struct
import zlib
from dataclasses import dataclass

import httpx

CLINVAR_VCF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_TBI = CLINVAR_VCF + ".tbi"

LINEAR_SHIFT = 14  # tabix linear index window = 16 kb


@dataclass
class TabixIndex:
    """Linear index only — sufficient for region lookup, and far simpler than
    the binning index. Maps reference name -> tuple of virtual offsets, one per
    16 kb window."""

    linear: dict[str, tuple[int, ...]]
    file_size: int

    def byte_range(self, chrom: str, beg: int, end: int, pad_windows: int = 2) -> tuple[int, int]:
        """Compressed byte range [start, stop) covering the region.

        Virtual offsets pack the compressed block offset in the high 48 bits
        and the offset within the decompressed block in the low 16.
        """
        if chrom not in self.linear:
            raise KeyError(f"{chrom!r} not in index; have e.g. {list(self.linear)[:5]}")
        ioffs = self.linear[chrom]
        if not ioffs:
            raise KeyError(f"{chrom!r} has an empty linear index")

        lo = min(beg >> LINEAR_SHIFT, len(ioffs) - 1)
        # Windows with no overlapping record are stored as 0; walk back to the
        # last populated window so we start before the region, never after it.
        while lo > 0 and ioffs[lo] == 0:
            lo -= 1
        start = ioffs[lo] >> 16

        hi = (end >> LINEAR_SHIFT) + pad_windows
        if hi < len(ioffs) and ioffs[hi] != 0:
            stop = (ioffs[hi] >> 16) + 65536  # +1 block so the last record is whole
        else:
            stop = self.file_size
        return start, min(stop, self.file_size)


def load_index(url: str = CLINVAR_TBI, vcf_url: str = CLINVAR_VCF) -> TabixIndex:
    raw = gzip.decompress(httpx.get(url, timeout=180).content)
    if raw[:4] != b"TBI\x01":
        raise ValueError(f"not a tabix index: magic={raw[:4]!r}")

    n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = struct.unpack_from("<8i", raw, 4)
    off = 36
    names = [n.decode() for n in raw[off:off + l_nm].split(b"\x00") if n]
    off += l_nm

    linear: dict[str, tuple[int, ...]] = {}
    for i in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", raw, off)
        off += 4
        for _ in range(n_bin):
            _bin_id, n_chunk = struct.unpack_from("<Ii", raw, off)
            off += 8 + n_chunk * 16  # skip the binning index entirely
        (n_intv,) = struct.unpack_from("<i", raw, off)
        off += 4
        linear[names[i]] = struct.unpack_from(f"<{n_intv}Q", raw, off)
        off += n_intv * 8

    size = int(httpx.head(vcf_url, timeout=60).headers["content-length"])
    return TabixIndex(linear=linear, file_size=size)


def _bgzf_decompress(data: bytes) -> bytes:
    """Decompress concatenated gzip members, tolerating a truncated tail.

    A byte range will almost always end mid-block; that final partial member is
    expected and is simply dropped.
    """
    out = bytearray()
    pos = 0
    while pos < len(data):
        d = zlib.decompressobj(31)
        try:
            out += d.decompress(data[pos:])
        except zlib.error:
            break
        consumed = len(data) - pos - len(d.unused_data)
        if consumed <= 0:
            break
        pos += consumed
    return bytes(out)


def fetch_region(index: TabixIndex, chrom: str, beg: int, end: int,
                 vcf_url: str = CLINVAR_VCF) -> list[str]:
    """VCF data lines overlapping chrom:beg-end."""
    start, stop = index.byte_range(chrom, beg, end)
    r = httpx.get(vcf_url, headers={"Range": f"bytes={start}-{stop - 1}"}, timeout=300)
    r.raise_for_status()
    text = _bgzf_decompress(r.content).decode("utf-8", errors="replace")

    lines = text.split("\n")
    kept = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        f = line.split("\t", 2)
        if len(f) < 3 or f[0] != chrom:
            continue
        try:
            pos = int(f[1])
        except ValueError:
            continue  # truncated first/last line
        if beg <= pos <= end:
            kept.append(line)
    return kept
