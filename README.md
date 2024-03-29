Macromolecular Data Frames
==========================

[![Last release](https://img.shields.io/pypi/v/macromol_dataframe.svg)](https://pypi.python.org/pypi/macromol_dataframe)
[![Python version](https://img.shields.io/pypi/pyversions/macromol_dataframe.svg)](https://pypi.python.org/pypi/macromol_dataframe)
[![Documentation](https://img.shields.io/readthedocs/macromol_dataframe.svg)](https://macromol-dataframe.readthedocs.io/en/latest/)
[![Test status](https://img.shields.io/github/actions/workflow/status/kalekundert/macromol_dataframe/test.yml?branch=master)](https://github.com/kalekundert/macromol_dataframe/actions)
[![Test coverage](https://img.shields.io/codecov/c/github/kalekundert/macromol_dataframe.svg)](https://app.codecov.io/github/kalekundert/macromol_dataframe)
[![Last commit](https://img.shields.io/github/last-commit/kalekundert/macromol_dataframe?logo=github)](https://github.com/kalekundert/macromol_dataframe)

*Macromol Dataframe* is a library meant to help with processing macromolecular 
coordinate data, e.g. mmCIF files downloaded from the Protein Data Bank (PDB). 
The key idea behind this library is that the best way to work with such data is 
by using data frames, specifically [`polars.DataFrame`].  The advantages of 
this approach are:

- Flexibility: Data frames are general-purpose data processing tools, and are 
  more than capable of accommodating any kind of analysis.

- Performance: Data frames are meant for processing huge quantities of data, 
  and are accordingly well-optimized.  Polars in particular achieves very good 
  performance by using techniques such as execution planning, SIMD 
  instructions, and multi-threading.

- Familiarity: Data scientists work with data frames all the time, so using 
  them here lowers the learning curve and makes this library easy to get 
  started with.  There's not much to learn!

Here's an example showing how to load a specific biological assembly from an 
mmCIF file:

```pycon
>>> import macromol_dataframe as mmdf
>>> df = mmdf.read_biological_assembly('6uad.cif.gz', model_id='1', assembly_id='1')
>>> df.select('seq_id', 'comp_id', 'atom_id', 'x', 'y', 'z')
shape: (2_312, 6)
┌────────┬─────────┬─────────┬───────────┬──────────┬──────────┐
│ seq_id ┆ comp_id ┆ atom_id ┆ x         ┆ y        ┆ z        │
│ ---    ┆ ---     ┆ ---     ┆ ---       ┆ ---      ┆ ---      │
│ i64    ┆ str     ┆ str     ┆ f64       ┆ f64      ┆ f64      │
╞════════╪═════════╪═════════╪═══════════╪══════════╪══════════╡
│ 2      ┆ ASN     ┆ N       ┆ -9.89268  ┆ 25.4788  ┆ 9.32073  │
│ 2      ┆ ASN     ┆ CA      ┆ -11.30656 ┆ 25.42029 ┆ 8.91019  │
│ 2      ┆ ASN     ┆ C       ┆ -12.19303 ┆ 26.2788  ┆ 9.79681  │
│ 2      ┆ ASN     ┆ O       ┆ -12.48258 ┆ 25.8771  ┆ 10.91766 │
│ 2      ┆ ASN     ┆ CB      ┆ -11.82931 ┆ 23.99427 ┆ 8.9393   │
│ …      ┆ …       ┆ …       ┆ …         ┆ …        ┆ …        │
│ null   ┆ HOH     ┆ O       ┆ -41.101   ┆ 23.389   ┆ 7.03     │
│ null   ┆ HOH     ┆ O       ┆ -4.60757  ┆ 22.48844 ┆ 9.93407  │
│ null   ┆ HOH     ┆ O       ┆ -22.48104 ┆ 27.68223 ┆ -4.26327 │
│ null   ┆ HOH     ┆ O       ┆ -38.8232  ┆ 17.99957 ┆ 9.24767  │
│ null   ┆ HOH     ┆ O       ┆ -40.22527 ┆ 15.63538 ┆ 7.88049  │
└────────┴─────────┴─────────┴───────────┴──────────┴──────────┘
```

[`polars.DataFrame`]: https://docs.pola.rs/py-polars/html/reference/dataframe/index.html#dataframe
