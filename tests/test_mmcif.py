import pytest
import parametrize_from_file as pff
import macromol_dataframe as mmdf
import macromol_dataframe.mmcif as _mmdf
import polars as pl
import polars.testing
import numpy as np
import gemmi.cif

from macromol_dataframe.testing import dataframe, atoms_csv, matrix
from pathlib import Path

with_pl = pff.Namespace('import polars as pl')
with_mmdf = pff.Namespace(mmdf=mmdf, Column=_mmdf.Column)

def assembly_gen(params):
    return dataframe(
            params,
            required_cols=['assembly_id', 'subchain_ids', 'oper_expr'],
            exprs={
                'subchain_ids': pl.col('subchain_ids').str.split(',')
            },
    )

def oper_map(params):
    return {k: matrix(v) for k, v in params.items()}

def entities(params):
    return dataframe(
            params,
            col_aliases={
                'mw': 'formula_weight_Da',
            },
            dtypes={
                'mw': float
            },
    )


@pff.parametrize(
        schema=[
            pff.cast(schema=with_mmdf.eval, expected=with_pl.eval),
            with_mmdf.error_or('expected'),
        ],
)
def test_extract_dataframe(mmcif, prefix, schema, expected, error):
    cif = gemmi.cif.read_string(mmcif).sole_block()
    with error:
        df = _mmdf._extract_dataframe(cif, prefix, schema)
        pl.testing.assert_frame_equal(df, expected, check_exact=False)

@pff.parametrize(
        schema=[
            pff.cast(
                asym_atoms=atoms_csv,
                assembly_gen=assembly_gen,
                oper_map=oper_map,
                entities=entities,
            ),
            with_mmdf.error_or('expected'),
        ],
)
def test_read_mmcif(
        tmp_path,
        mmcif,
        pdb_id,
        asym_atoms,
        assembly_gen,
        oper_map,
        entities,
        error,
):
    cif_path = tmp_path / 'mock.cif'
    cif_path.write_text(mmcif)

    with error:
        struct = mmdf.read_mmcif(cif_path)

        assert struct.id == pdb_id
        pl.testing.assert_frame_equal(
                struct.asym_atoms, asym_atoms,
                check_exact=False,
                check_column_order=False,
        )
        pl.testing.assert_frame_equal(
                struct.assembly_gen, assembly_gen,
                check_exact=False,
                check_column_order=False,
        )
        for k in oper_map:
            np.testing.assert_allclose(
                    struct.oper_map[k], oper_map[k]
            )
        pl.testing.assert_frame_equal(
                struct.entities, entities,
                check_exact=False,
                check_column_order=False,
        )

@pff.parametrize(
        schema=[
            pff.cast(expected=atoms_csv),
            with_mmdf.error_or('expected'),
        ],
)
def test_read_biological_assembly(tmp_path, mmcif, model_id, assembly_id, expected, error):
    cif_path = tmp_path / 'mock.cif'
    cif_path.write_text(mmcif)

    with error:
        atoms = mmdf.read_biological_assembly(
                cif_path,
                model_id=model_id,
                assembly_id=assembly_id,
        )
        pl.testing.assert_frame_equal(
                atoms, expected,
                check_exact=False,
                check_column_order=False,
        )

@pff.parametrize(
        schema=pff.cast(expected=atoms_csv),
)
def test_read_asymmetric_unit(tmp_path, mmcif, expected):
    cif_path = tmp_path / 'mock.cif'
    cif_path.write_text(mmcif)

    atoms = mmdf.read_asymmetric_unit(cif_path)
    pl.testing.assert_frame_equal(
            atoms, expected,
            check_exact=False,
            check_column_order=False,
    )

def test_write_mmcif(tmp_path):
    test_dir = Path(__file__).parent 
    in_path = test_dir / 'pdb' / '1fav.cif.gz'
    out_path = tmp_path / '1fav.cif'
    ref_path = test_dir / 'ref' / '1fav.cif'

    atoms = mmdf.read_biological_assembly(in_path, model_id='1', assembly_id='1')
    mmdf.write_mmcif(out_path, atoms)

    if not ref_path.exists():
        pytest.fail(f"""\
No reference path found.

Output from this test run was written to:
{out_path}

Check this output manually.  If it looks good, copy it to:
{ref_path}
""")

    # If the two files aren't the same, pytest takes forever to calculate a 
    # diff, so we go out of our way to avoid that.
    if ref_path.read_text() != out_path.read_text():
        pytest.fail(f"""\
The current mmCIF output doesn't match the reference output.

Output from this test run was written to:
{out_path}

Check this output manually.  If it looks good, copy it to:
{ref_path}
""")

