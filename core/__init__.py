"""Locus core logic.

Pure functions over dataframes and scalars. No I/O, no database, no UI imports
— that separation is what makes the rules unit-testable without a database and
what makes swapping Streamlit for React a re-skin rather than a rewrite.

All database access lives in queries.py; all rendering lives above this layer.
"""
