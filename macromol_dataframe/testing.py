import polars as pl
import numpy as np
import math

from .coords import make_coord_frame
from io import StringIO

def matrix(params):
    io = StringIO(params)
    return np.loadtxt(io, dtype=float)

def vector(params):
    return np.array([eval(x, math.__dict__) for x in params.split()])

def frame(params):
    origin = coord(params['origin'])
    rot_vec_rad = vector(params['rot_vec_rad'])
    return make_coord_frame(origin, rot_vec_rad)

def frames(params):
    return [frame(x) for x in params]

def coord(params):
    return matrix(params)

def coords(params):
    coords = matrix(params)
    coords.shape = (1, *coords.shape)[-2:]
    return coords

def atoms_fwf(params):
    dtypes = {
            'chain_id': str,
            'subchain_id': str,
            'seq_id': int,
            'comp_id': str,
            'atom_id': str,
            'element': str,
            'x': float,
            'y': float,
            'z': float,
            'occupancy': float,
            'b_factor': float,
    }
    col_aliases = {
            'chain': 'chain_id',
            'subchain': 'subchain_id',
            'resn': 'comp_id',
            'resi': 'seq_id',
            'e': 'element',
            'q': 'occupancy',
            'b': 'b_factor',
    }

    io = StringIO(params)

    # If there isn't a header row, assume that only the element and coordinate 
    # columns were given.  The rest of the columns are given default values.
    # 
    # If there is a header row, create whichever columns it specifies.  Note 
    # that the 'x', 'y', and 'z' columns are always required.

    header = io.readline().split()
    if {'x', 'y', 'z'} <= set(header):
        pass
    else:
        header = ['e', 'x', 'y', 'z']
        io.seek(0)

    rows = []
    for line in io.readlines():
        rows.append(line.split())

    df = (
            pl.DataFrame(rows, header, orient='row')
            .rename(lambda x: col_aliases.get(x, x))
            .with_columns(
                pl.col('*').replace(['.', '?'], None)
            )
            .cast({
                col_aliases.get(k, k): dtypes.get(k, str)
                for k in header
            })
    )

    if 'comp_id' not in df.columns:
        df = df.with_columns(comp_id=pl.lit('ALA'))
    if 'element' not in df.columns:
        df = df.with_columns(element=pl.lit('C'))
    if 'occupancy' not in df.columns:
        df = df.with_columns(occupancy=pl.lit(1.0))

    return df

def atoms_csv(params):
    dtypes = {
            'symmetry_mate': pl.Int32,
            'chain_id': str,
            'subchain_id': str,
            'alt_id': str,
            'seq_id': int,
            'comp_id': str,
            'atom_id': str,
            'element': str,
            'x': float,
            'y': float,
            'z': float,
            'occupancy': float,
            'b_factor': float,
    }
    df = pl.read_csv(StringIO(params))
    return (
            df
            .with_columns(
                pl.col(pl.String).str.strip_chars()
            )
            .cast({k: dtypes.get(k, str) for k in df.columns})
    )


