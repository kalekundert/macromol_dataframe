import polars as pl
import numpy as np
import math

from .coords import make_coord_frame_from_rotation_vector
from .error import TidyError
from io import StringIO
from functools import partial

NOT_GIVEN = object()

def matrix(params):
    io = StringIO(params)
    return np.loadtxt(io, dtype=float)

def vector(params):
    return np.array([eval(x, math.__dict__) for x in params.split()])

def frame(params):
    origin = coord(params['origin'])
    rot_vec_rad = vector(params['rot_vec_rad'])
    return make_coord_frame_from_rotation_vector(origin, rot_vec_rad)

def frames(params):
    return [frame(x) for x in params]

def coord(params):
    return matrix(params)

def coords(params):
    coords = matrix(params)
    coords.shape = (1, *coords.shape)[-2:]
    return coords

def atoms_fwf(params):
    return dataframe(
            params,
            dtypes={
                'seq_id': int,
                'x': float,
                'y': float,
                'z': float,
                'occupancy': float,
                'b_factor': float,
            },
            col_aliases={
                'chain': 'chain_id',
                'subchain': 'subchain_id',
                'segi': 'subchain_id',
                'resn': 'comp_id',
                'resi': 'seq_id',
                'atom': 'atom_id',
                'elem': 'element',
                'e': 'element',
                'q': 'occupancy',
                'b': 'b_factor',
            },
            default_cols={
                'comp_id': 'ALA',
                'element': 'C',
                'occupancy': 1.0,
            },
            default_header=['e', 'x', 'y', 'z'],
            required_cols=['x', 'y', 'z'],
    )

def atoms_csv(params):
    dtypes = {
            'symmetry_mate': pl.Int32,
            'chain_id': str,
            'subchain_id': str,
            'alt_id': str,
            'seq_id': int,
            'seq_label': str,
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

def dataframe(
        params=NOT_GIVEN,
        *,
        exprs={},
        dtypes={},
        col_aliases={},
        default_header=[],
        default_cols={},
        required_cols=[],
):
    """
    Instantiate a data frame from a string containing a fixed-width formatted 
    table.

    Arguments:
        exprs:
            A dictionary mapping column names (aliases ok) to polars 
            expressions.  These expressions provide the means to parse more 
            complicated field values.  Only those expressions corresponding to 
            columns in the table will be used; the rest will be silently 
            ignored.  Expressions are applied before data types (see the 
            *dtypes* argument below), and it is allowed to specify both for the 
            same column.

        dtypes:
            A dictionary mapping column names (aliases ok) to data types.  The 
            string values parsed from the data will be converted to these 
            types.  Any columns not included in this mapping will be left as 
            strings.

        col_aliases:
            A dictionary mapping "input" column names to "output" column names.  
            Because the data frames used in test cases are typed by hand, it is 
            often convenient to allow shorthand column headers.  However, the 
            output column names must match those expected by the APIs being 
            tested, which are typically more verbose.

        default_header:
            A list of strings giving the names to use for each column (aliases 
            ok), if no column names are specified.  Column names are considered 
            to be specified if the first row of the table contains all of the 
            required column names, or all of the default column names if 
            *required_cols* isn't specified.

        default_cols:
            A dictionary mapping column names (aliases ok) to default values.  
            If any of these columns are not present in the table, they will be 
            added to the dataframe with the associated value in every row.  
            Note that the default columns are added before the *dtypes* are 
            evaluated.  This means that even if you specify a default value of 
            the wrong type (e.g. null, int instead of float, etc.), the column 
            will still end up with the correct data type.

        required_cols:
            A list of column names (aliases ok) that must be present in the 
            output dataframe.  If a required column is missing, and exception 
            will be raised.
    """

    if params is NOT_GIVEN:
        kwargs = locals()
        del kwargs['params']
        return partial(dataframe, **kwargs)

    rows = [line.split() for line in params.splitlines()]
    header, body = rows[0], rows[1:]

    def resolve_aliases_dict(d):
        return {
                col_aliases.get(k, k): v
                for k, v in d.items()
        }

    def resolve_aliases_list(l):
        return [col_aliases.get(k, k) for k in l]

    exprs = resolve_aliases_dict(exprs)
    dtypes = resolve_aliases_dict(dtypes)
    default_header = resolve_aliases_list(default_header)
    default_cols = resolve_aliases_dict({
        k: pl.lit(v)
        for k, v in default_cols.items()
    })
    required_cols = resolve_aliases_list(required_cols)
    header = resolve_aliases_list(header)

    if default_header:
        if set(required_cols or default_header) - set(header):
            header = default_header
            body = rows

    for row in body:
        if len(row) != len(header):
            raise InputError(
                    "rows must have same number of column as header",
                    info=[f"header: {header!r}"],
                    blame=[f"offending row: {row!r}"],
            )

    def drop_missing_cols(d):
        return {k: v for k, v in d.items() if k in header}

    def drop_existing_cols(d):
        return {k: v for k, v in d.items() if k not in header}

    schema = {k: str for k in header}
    df = (
            pl.DataFrame(body, schema, orient='row')
    )
    df = (
            df
            .with_columns(
                pl.col('*').replace(['.', '?'], None)
            )
            .with_columns(**drop_missing_cols(exprs))
            .with_columns(**drop_existing_cols(default_cols))
            .cast(drop_missing_cols(dtypes))
    )

    if set(required_cols) - set(df.columns):
        raise InputError(
                "missing required columns",
                info=[f"required columns: {required_cols!r}"],
                blame=[f"specified columns: {df.columns!r}"],
        )

    return df

class InputError(TidyError):
    pass
