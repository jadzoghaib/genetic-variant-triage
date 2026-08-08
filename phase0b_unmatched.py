# /// script
# requires-python = ">=3.11"
# dependencies = ["duckdb", "pandas"]
# ///
"""
Phase 0b — diagnose the 74 ClinVar missense SNVs that did not match AlphaMissense.

The schema decision this settles: does Locus need per-transcript variant effects,
or can it stand on MANE Select alone?

Discriminator: for each unmatched record, is its genomic POS present in the
AlphaMissense file at all?
  - POS absent entirely  -> not coding on the MANE Select transcript. ClinVar is
    asserting 'missense' against some other transcript. => multi-transcript reality.
  - POS present, ALT absent -> that substitution is synonymous/nonsense on MANE,
    or AlphaMissense excluded it. => single-transcript is fine, edge case only.
  - POS present, REF disagrees -> coordinate/normalisation problem. => serious.
"""

from pathlib import Path

import duckdb
import pandas as pd

DATA = Path(__file__).parent / "data"
DB = Path(__file__).parent / "phase0_scratch.duckdb"

con = duckdb.connect(str(DB))

am = pd.concat(
    [pd.read_csv(DATA / f"am_{g}_hg38.csv").assign(symbol=g)
     for g in ("BRCA1", "TP53", "PTEN")],
    ignore_index=True,
)
am["chrom"] = am["CHROM"].str.removeprefix("chr")
con.register("am_raw", am)

unmatched = con.execute("""
    SELECT symbol, chrom, pos, ref, alt, clnsig, clnsig_raw, clinvar_id
    FROM joined WHERE is_missense AND is_snv AND NOT matched
""").df()
print(f"unmatched ClinVar missense SNVs: {len(unmatched)}\n")

con.register("um", unmatched)

print("=" * 74)
print("A. Is the genomic position present in AlphaMissense at all?")
print("=" * 74)
print(con.execute("""
    SELECT um.symbol,
           COUNT(*) AS unmatched,
           COUNT(*) FILTER (WHERE p.pos IS NULL) AS pos_absent_from_AM,
           COUNT(*) FILTER (WHERE p.pos IS NOT NULL) AS pos_present_alt_missing
    FROM um
    LEFT JOIN (SELECT DISTINCT symbol, POS AS pos FROM am_raw) p
      ON p.symbol = um.symbol AND p.pos = um.pos
    GROUP BY um.symbol ORDER BY um.symbol
""").df().to_string(index=False))

print("\n" + "=" * 74)
print("B. Where the position IS present, does the reference base agree?")
print("   (disagreement would mean a coordinate/normalisation defect)")
print("=" * 74)
ref_check = con.execute("""
    SELECT um.symbol, um.pos, um.ref AS clinvar_ref, a.REF AS am_ref, um.alt, um.clinvar_id
    FROM um JOIN (SELECT DISTINCT symbol, POS, REF FROM am_raw) a
      ON a.symbol = um.symbol AND a.POS = um.pos
    WHERE a.REF <> um.ref
""").df()
print(f"reference-base disagreements: {len(ref_check)}")
if len(ref_check):
    print(ref_check.to_string(index=False))

print("\n" + "=" * 74)
print("C. For positions present in AM, which ALTs does AM carry vs what ClinVar asked for?")
print("=" * 74)
sample = con.execute("""
    SELECT um.symbol, um.pos, um.ref, um.alt AS clinvar_alt,
           string_agg(DISTINCT a.ALT, '/') AS am_alts_at_pos,
           string_agg(DISTINCT a.protein_variant, ',') AS am_protein_variants
    FROM um JOIN am_raw a ON a.symbol = um.symbol AND a.POS = um.pos
    GROUP BY um.symbol, um.pos, um.ref, um.alt
    LIMIT 10
""").df()
print(sample.to_string(index=False) if len(sample) else "  (none — every unmatched position is absent from AM)")

print("\n" + "=" * 74)
print("D. What transcripts is ClinVar asserting 'missense' against?")
print("=" * 74)
mc = con.execute("""
    SELECT mc, COUNT(*) AS n
    FROM joined WHERE is_missense AND is_snv AND NOT matched
    GROUP BY mc ORDER BY n DESC LIMIT 6
""").df()
print(mc.to_string(index=False))

print("\n" + "=" * 74)
print("E. Where do the unmatched positions sit relative to the AM coding span?")
print("=" * 74)
print(con.execute("""
    SELECT um.symbol,
           MIN(um.pos) AS unmatched_min, MAX(um.pos) AS unmatched_max,
           (SELECT MIN(POS) FROM am_raw a WHERE a.symbol = um.symbol) AS am_min,
           (SELECT MAX(POS) FROM am_raw a WHERE a.symbol = um.symbol) AS am_max
    FROM um GROUP BY um.symbol ORDER BY um.symbol
""").df().to_string(index=False))

con.close()
