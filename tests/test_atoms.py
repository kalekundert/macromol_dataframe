import parametrize_from_file as pff
import macromol_dataframe as mmdf
import polars as pl
import polars.testing
import numpy as np

from macromol_dataframe.testing import atoms_fwf
from test_coords import frame, coords

with_py = pff.Namespace()

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



