import parametrize_from_file as pff
import macromol_dataframe as mmdf
import polars as pl
import polars.testing
import numpy as np

from test_coords import frame, coords
from io import StringIO

with_py = pff.Namespace()

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

def atoms_csv(csv):
    return (
            pl.read_csv(StringIO(csv))
            .with_columns(
                pl.col(pl.String).str.strip_chars()
            )
            .cast({
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
            })
    )


@pff.parametrize(
        schema=pff.cast(
            atoms=atoms_fwf,
            homogeneous=with_py.eval,
            expected=coords,
        ),
)
def test_get_atom_coords(atoms, homogeneous, expected):
    coords = mmdf.get_atom_coords(atoms, homogeneous=homogeneous)
    np.testing.assert_array_equal(coords, expected)

@pff.parametrize(
        schema=pff.cast(
            atoms=atoms_fwf,
            coords=coords,
            expected=atoms_fwf,
        ),
)
def test_replace_atom_coords(atoms, coords, expected):
    actual = mmdf.replace_atom_coords(atoms, coords)
    pl.testing.assert_frame_equal(actual, expected)

@pff.parametrize(
        schema=pff.cast(
            atoms_x=atoms_fwf,
            frame_xy=frame,
            expected_y=atoms_fwf,
        ),
)
def test_transform_atom_coords(atoms_x, frame_xy, expected_y):
    atoms_y = mmdf.transform_atom_coords(atoms_x, frame_xy)
    pl.testing.assert_frame_equal(atoms_y, expected_y)



