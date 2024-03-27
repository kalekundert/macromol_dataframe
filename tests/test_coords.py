import numpy as np
import macromol_dataframe as mmdf
import parametrize_from_file as pff

from pytest import approx
from hypothesis import given, settings
from hypothesis.extra.numpy import arrays
from scipy.spatial.transform import Rotation
from macromol_dataframe.testing import *
from math import pi

with_math = pff.Namespace('from math import *')

def float_bounds(max=100):
    return dict(
            min_value=-max,
            max_value=max,
            allow_nan=False,
            allow_infinity=False,
    )


def test_homogenize_coords_1d():
    coords = np.array([1, 2, 3])
    expected = np.array([1, 2, 3, 1])
    np.testing.assert_array_equal(
            mmdf.homogenize_coords(coords),
            expected,
    )

def test_homogenize_coords_2d():
    coords = np.array([
        [1, 2, 3],
        [2, 3, 4],
    ])
    expected = np.array([
        [1, 2, 3, 1],
        [2, 3, 4, 1],
    ])
    np.testing.assert_array_equal(
            mmdf.homogenize_coords(coords),
            expected,
    )

@pff.parametrize(
        schema=pff.cast(frame_xy=matrix, coords_x=coords, expected_y=coords),
)
def test_transform_coords(frame_xy, coords_x, expected_y):
    coords_x = mmdf.homogenize_coords(coords_x)
    coords_y = mmdf.transform_coords(coords_x, frame_xy)
    expected_y = mmdf.homogenize_coords(expected_y)
    assert coords_y == approx(expected_y)

@pff.parametrize(
        schema=pff.cast(
            origin=coord,
            rot_vec_rad=vector,
            coords_x=coords,
            expected_y=coords,
        ),
)
def test_make_coord_frame(origin, rot_vec_rad, coords_x, expected_y):
    # It's not enough to test the 16 numbers making up the matrix, I need to 
    # test that it transforms things in the way it should.
    frame_xy = mmdf.make_coord_frame(origin, rot_vec_rad)
    test_transform_coords(frame_xy, coords_x, expected_y)

    np.testing.assert_allclose(
            mmdf.get_origin(frame_xy),
            origin,
    )
    np.testing.assert_allclose(
            mmdf.get_rotation_matrix(frame_xy),
            Rotation.from_rotvec(rot_vec_rad).as_matrix(),
    )

@settings(deadline=None)
@given(
        arrays(float, 3, elements=float_bounds()),
        arrays(float, 3, elements=float_bounds(2*pi)),
        arrays(float, 3, elements=float_bounds()),
)
def test_invert_coord_frame(origin, rot_vec_rad, coords_x):
    frame_xy = mmdf.make_coord_frame(origin, rot_vec_rad)
    frame_yx = mmdf.invert_coord_frame(frame_xy)

    coords_x = mmdf.homogenize_coords(coords_x)
    coords_y = mmdf.transform_coords(coords_x, frame_xy)
    coords_x2 = mmdf.transform_coords(coords_y, frame_yx)

    assert coords_x2 == approx(coords_x)


