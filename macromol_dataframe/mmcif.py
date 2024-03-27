import polars as pl
import polars.selectors as cs
import numpy as np
import gemmi.cif

from .atoms import transform_atom_coords
from .coords import Frame
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from typing import Any
from collections.abc import Callable

def read_mmcif(in_path: Path, model_id: str, assembly_id: str):
    cif = gemmi.cif.read(str(in_path)).sole_block()

    with _add_path_to_mmcif_error(in_path):
        asym_atoms = _extract_atom_site(cif)
        struct_assembly_gen, struct_oper_map = _extract_struct_assembly_gen(
                cif, asym_atoms,
        )
        return _make_biological_assembly(
                _select_model(asym_atoms, model_id), 
                struct_assembly_gen,
                struct_oper_map,
                assembly_id,
        )

def read_mmcif_asymmetric_unit(in_path: Path):
    cif = gemmi.cif.read(str(in_path)).sole_block()
    return _extract_atom_site(cif)

def write_mmcif(out_path: Path, atoms: pl.DataFrame, name: str = None):
    col_map = {
            'chain_id': 'auth_asym_id',
            'subchain_id': 'label_asym_id',
            'alt_id': 'label_alt_id',
            'seq_id': 'label_seq_id',
            'comp_id': 'label_comp_id',
            'atom_id': 'label_atom_id',
            'element': 'type_symbol',
            'x': 'Cartn_x',
            'y': 'Cartn_y',
            'z': 'Cartn_z',
            'occupancy': 'occupancy',
            'b_factor': 'B_iso_or_equiv',
    }

    block = gemmi.cif.Block(name or out_path.stem)
    loop = block.init_loop('_atom_site.', list(col_map.values()))

    # Give each chain a unique name, otherwise symmetry mates will appear to be 
    # duplicate atoms, and this can confuse downstream programs (e.g. pymol).

    atoms_str = (
            atoms
            .with_columns(
                pl.concat_str(
                    pl.col('chain_id'),
                    pl.col('symmetry_mate') + 1,
                )
            )
            .with_columns(
                cs.float().round(3),
            )
            .with_columns(
                pl.col('*').cast(str).replace(None, '?')
            )
    )

    for row in atoms_str.iter_rows(named=True):
        loop.add_row([row[k] for k in col_map])

    options = gemmi.cif.WriteOptions()
    options.align_loops = 30

    block.write_file(str(out_path), options)

@contextmanager
def _add_path_to_mmcif_error(path: Path):
    try:
        yield

    except MmcifError as err:
        err.info = [f'path: {path}', *err.info]
        raise

def _extract_dataframe(cif, key_prefix, schema):
    # Gemmi automatically interprets `?` and `.`, but this leads to a few 
    # problems.  First is that it makes column dtypes dependent on the data; if 
    # a column doesn't have any non-null values, polars won't know that it 
    # should be a string.  Second is that gemmi distinguishes between `?` 
    # (null) and `.` (false).  This is a particularly unhelpful distinction 
    # when the column in question is supposed to contain float data, because 
    # the latter then becomes 0 rather than null.
    #
    # To avoid these problems, when initially loading the data frame, we 
    # explicitly specify a schema where each column is a string.  Doing this 
    # happens to convert any booleans present in the data to null, thereby 
    # solving both of the above problems at once.

    loop = cif.get_mmcif_category(f'_{key_prefix}.')
    df = pl.DataFrame(loop, {k: str for k in loop})

    if df.is_empty():
        schema = {k: v.dtype for k, v in schema.items()}
        return pl.DataFrame([], schema)

    # Check for missing required columns:
    missing_cols = [
            v.name
            for v in schema.values()
            if v.required and v.name not in df.columns
    ]
    if missing_cols:
        err = MmcifError("missing required column(s)")
        err.info = [f"category: _{key_prefix}.*"]
        err.blame = [f"missing column(s): {missing_cols}"]
        raise err

    return (
            df

            # Fill in missing optional columns:
            .with_columns([
                pl.lit(None, dtype=str).alias(v.name)
                for v in schema.values()
                if not v.required and v.name not in df.columns
            ])

            # Cast, rename, and sort desired columns:
            .select([
                pl.col(v.name).cast(v.dtype).alias(k)
                for k, v in schema.items()
            ])

            # Remove all-null rows:
            .filter(~pl.all_horizontal(pl.all().is_null()))
    )

def _extract_atom_site(cif):
    return _extract_dataframe(
            cif, 'atom_site',
            schema=dict(
                model_id=Column('pdbx_PDB_model_num'), 
                chain_id=Column('auth_asym_id'),
                subchain_id=Column('label_asym_id'),
                alt_id=Column('label_alt_id'),
                seq_id=Column('label_seq_id', dtype=int),
                comp_id=Column('label_comp_id'),
                atom_id=Column('label_atom_id'),
                element=Column('type_symbol', required=True),
                x=Column('Cartn_x', dtype=float, required=True),
                y=Column('Cartn_y', dtype=float, required=True),
                z=Column('Cartn_z', dtype=float, required=True),
                occupancy=Column('occupancy', dtype=float),
                b_factor=Column('B_iso_or_equiv', dtype=float),
            ),
    )

def _extract_struct_assembly_gen(cif, asym_atoms):
    """
    Construct the rules needed to build biological assemblies from the 
    asymmetric unit.  If the mmCIF file doesn't contain this information, 
    assume that the asymmetric unit is a biological assembly.
    """
    struct_oper_list = (
            _extract_dataframe(
                cif, 'pdbx_struct_oper_list',
                schema=dict(
                    id=Column('id', required=True),
                    matrix_11=Column('matrix[1][1]', dtype=float, required=True),
                    matrix_12=Column('matrix[1][2]', dtype=float, required=True),
                    matrix_13=Column('matrix[1][3]', dtype=float, required=True),
                    vector_1=Column('vector[1]', dtype=float, required=True),
                    matrix_21=Column('matrix[2][1]', dtype=float, required=True),
                    matrix_22=Column('matrix[2][2]', dtype=float, required=True),
                    matrix_23=Column('matrix[2][3]', dtype=float, required=True),
                    vector_2=Column('vector[2]', dtype=float, required=True),
                    matrix_31=Column('matrix[3][1]', dtype=float, required=True),
                    matrix_32=Column('matrix[3][2]', dtype=float, required=True),
                    matrix_33=Column('matrix[3][3]', dtype=float, required=True),
                    vector_3=Column('vector[3]', dtype=float, required=True),
                ),
            )
    )

    if struct_oper_list.is_empty():
        struct_oper_map = {
                '1': np.eye(4)
        }
        struct_assembly_gen = pl.DataFrame([
            dict(
                assembly_id='1',
                subchain_ids=asym_atoms['subchain_id'].unique(),
                oper_expr='1',
            ),
        ])

    else:
        struct_oper_map = {
                x['id']: np.array([
                    [x['matrix_11'], x['matrix_12'], x['matrix_13'], x['vector_1']],
                    [x['matrix_21'], x['matrix_22'], x['matrix_23'], x['vector_2']],
                    [x['matrix_31'], x['matrix_32'], x['matrix_33'], x['vector_3']],
                    [             0,              0,              0,             1],
                ])
                for x in struct_oper_list.iter_rows(named=True)
        }
        struct_assembly_gen = (
                _extract_dataframe(
                    cif, 'pdbx_struct_assembly_gen',
                    schema=dict(
                        assembly_id=Column('assembly_id', required=True),
                        subchain_ids=Column('asym_id_list', required=True),
                        oper_expr=Column('oper_expression', required=True),
                    ),
                )
                .with_columns(
                    pl.col('subchain_ids').str.split(',')
                )
        )

    return struct_assembly_gen, struct_oper_map

def _select_model(asym_atoms, model_id: str):
    assert isinstance(model_id, str)

    no_models_specified = (
            asym_atoms
            .select(pl.col('model_id').is_null().all())
            .item()
    )
    if no_models_specified:
        pass
    else:
        asym_atoms = (
                asym_atoms
                .filter(
                    pl.col('model_id') == model_id
                )
        )

    return asym_atoms.drop('model_id')

def _make_biological_assembly(
        asym_atoms,
        struct_assembly_gen,
        struct_oper_map,
        assembly_id,
):
    oper_exprs = (
            struct_assembly_gen
            .filter(pl.col('assembly_id') == assembly_id)
    )

    if oper_exprs.is_empty():
        known_assemblies = \
                struct_assembly_gen['assembly_id'].unique().to_list()

        err = MmcifError("can't find biological assembly")
        err.info += [f"known assemblies: {known_assemblies}"]
        err.blame = [f"unknown assembly: {assembly_id!r}"]
        raise err

    bio_atoms = []

    for row in oper_exprs.iter_rows(named=True):
        subchain_ids = row['subchain_ids']
        frames = _parse_oper_expression(row['oper_expr'], struct_oper_map)

        for i, frame in enumerate(frames):
            sym_atoms = (
                    transform_atom_coords(
                        asym_atoms.filter(
                            pl.col('subchain_id').is_in(subchain_ids)
                        ),
                        frame,
                    )
                    .with_columns(
                        symmetry_mate=pl.lit(i),
                    )
            )
            bio_atoms.append(sym_atoms)

    return pl.concat(bio_atoms).rechunk()

def _parse_oper_expression(expr: str, oper_map: dict[str, Frame]):
    # According to the PDBx/mmCIF specification [1], it's possible for the 
    # operation expression to contain parenthetical expressions.  This would 
    # indicate that each transformation in one set of parentheses should be 
    # combined separately with each transformation in the next.  It's also 
    # possible for ranges of numbers to be specified with dashes.
    #
    # Handling the above cases would add a lot of complexity to this function.  
    # For now, I haven't found any structures that use this advanced syntax, so 
    # I didn't take the time to implement it.  But this is something I might 
    # have to come back to later on.

    if '(' in expr:
        err = MmcifError("unsupported expression")
        err.info = [
                "parenthetical expressions in biological assembly transformations are not currently supported",
                "handling these cases properly is not trivial, and at the time this code was written, there were no examples of such expressions in the PDB",
        ]
        err.blame = [f"expression: {expr}"]
        raise err

    return [
            oper_map[k]
            for k in expr.split(',')
    ]

@dataclass
class Column:
    name: str
    dtype: Callable[[str], Any] = str
    required: bool = False

class MmcifError(Exception):
    # This class mimics `tidyexc`.  I didn't want to use `tidyexc` directly, 
    # because it has weird incompatibilities with multiprocessing, and this 
    # code is meant to be used in pytorch data loaders.  Although I don't know 
    # the exact cause of the incompatibilities, I assume they have to do either 
    # with class-level state or maintaining a data structure of arbitrary 
    # objects.  This code gets rid of those aspects of `tidyexc` and just keeps 
    # track of a few strings that will be used to make an error message.

    def __init__(self, brief=None):
        self.brief = brief
        self.info = []
        self.blame = []

    def __str__(self):
        info_strs = ['• ' + x for x in self.info]
        blame_strs = ['✖ ' + x for x in self.blame]
        msg_strs = [self.brief, *info_strs, *blame_strs]
        return '\n'.join(msg_strs)

