import polars as pl
import polars.selectors as cs
import numpy as np
import gemmi.cif
import functools
import operator as op

from .atoms import Atoms, transform_atom_coords
from .coords import Frame
from .error import TidyError
from parsy import ParseError
from functools import reduce
from itertools import product
from more_itertools import flatten
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from typing import Dict, Union

class Structure:
    asym_atoms: Atoms
    assembly_gen: pl.DataFrame
    oper_map: Dict[str, Frame]
    entities: pl.DataFrame

    def __init__(self, id):
        self.id = id
        # Attributes set in `read_mmcif()`.

    def __repr__(self):
        return f'<Structure {self.id}>'

def read_mmcif(cif_path: Path) -> Structure:
    """
    Parse the information in an mmCIF file into a number of data frames.

    Arguments:
        cif_path:
            The path to the mmCIF file to read.

    This function should be used when neither `read_biological_assembly()` nor 
    `read_asymmetric_unit()` provide all of the information you want.  This 
    function returns more information, but in a less convenient format.

    In principle, this function should return *all* of the information encoded 
    in the mmCIF file.  In practice, it only returns the information that 
    someone has had some use for, and taken the time to implement.  If there's 
    any information that you feel is missing, let me know by opening a new 
    issue and/or pull request.  I'm very open to returning more information 
    from this function.
    """
    cif = gemmi.cif.read(str(cif_path)).sole_block()

    with _add_path_to_mmcif_error(cif_path):
        struct = Structure(cif.name)
        struct.asym_atoms = _extract_atom_site(cif)
        struct.assemblies = _extract_struct_assembly(cif)
        struct.assembly_gen, struct.oper_map = \
                _extract_struct_assembly_gen(cif, struct.asym_atoms)
        struct.entities = _extract_entities(cif)
        struct.polymers = _extract_polymers(cif)
    
    return struct

def read_biological_assembly(
        cif_path: Path,
        *,
        model_id: str,
        assembly_id: str,
) -> Atoms:
    """
    Parse a single biological assembly from the given mmCIF file.

    Arguments:
        cif_path:
            The path containing the mmCIF file to read.

        model_id:
            The id number of the model to read.  Valid ids are given by the 
            `_atom_site.pdbx_PDB_model_num` field in the mmCIF file.  If the 
            file doesn't specify any model numbers, this argument will be 
            ignored.

        assembly_id:
            The id string of the assembly to generate.  Valid ids are given by 
            the `_pdbx_struct_assembly` loop in the mmCIF file.

    Returns:
        A dataframe containing a row for each atom in the biological assembly.  
        See `make_biological_assembly()` for a more detailed description of 
        this dataframe.

    It is assumed that the caller already knows something about the file in 
    question, namely (i) which assembly is of interest and (ii) which model is 
    of interest.  For example, this information might come from a separate 
    program that scans the PDB and identifies assemblies of interest.  If you 
    don't have this information in advance, and you want to work it out from 
    the contents of the file itself, you need to parse the mmCIF file and 
    generate the biological assembly in separate steps.  See `read_mmcif()` and 
    `make_biological_assembly()`.  This is exactly what this function does 
    under the hood.

    The reason for returning only one assembly (e.g. as opposed to returning 
    coordinates for every assembly) is that generating assemblies is relatively 
    expensive.  Each atom in the assembly has to undergo a coordinate 
    transformation.  For performance-critical applications, like machine 
    learning, it is much better to transform only those coordinates that are 
    actually needed.
    """
    struct = read_mmcif(cif_path)

    with _add_path_to_mmcif_error(cif_path):
        return make_biological_assembly(
                select_model(struct.asym_atoms, model_id), 
                struct.assembly_gen,
                struct.oper_map,
                assembly_id,
        )

def read_asymmetric_unit(cif_path: Path) -> Atoms:
    """
    Parse coordinates for every atom in the asymmetric unit.

    Arguments:
        cif_path:
            The path containing the mmCIF file to read.

    This is basically a simplified version of `read_mmcif()` that only returns 
    atomic coordinates and not any of the other relationships encoded in the 
    mmCIF file.
    """
    struct = read_mmcif(cif_path)
    return struct.asym_atoms

def write_mmcif(cif_path: Path, atoms: Atoms, name: str = None) -> None:
    """
    Write the given atoms to a new mmCIF file.

    Arguments:
        cif_path:
            The path of the file to write.  If this path already exists, it 
            will be overwritten.

        atoms:
            A dataframe containing the atoms to include in the output file.  
            This dataframe must have the following columns:

            - ``chain_id``
            - ``subchain_id``
            - ``alt_id``
            - ``seq_id``
            - ``comp_id``
            - ``atom_id``
            - ``element``
            - ``x``
            - ``y``
            - ``z``
            - ``occupancy``
            - ``b_factor``

            If present, the ``symmetry_mate`` column will be concatenated to 
            the ``chain_id`` column to create a unique id for each symmetric 
            copy of a chain.  Any other columns will be ignored.

        name:
            An identifier to include in the mmCIF file.  By default, this is 
            taken from the last component of the given path.

    Note that this function is meant to help with debugging, by providing a way 
    to visualize *atoms* data frames in programs like PyMOL or Chimera.  It's 
    not meant to export all of the information imported by `read_mmcif()`.
    """
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

    block = gemmi.cif.Block(name or cif_path.name.split('.')[0])
    loop = block.init_loop('_atom_site.', list(col_map.values()))

    # Give each chain a unique name, otherwise symmetry mates will appear to be 
    # duplicate atoms, and this can confuse downstream programs (e.g. pymol).

    if 'symmetry_mate' in atoms.columns:
        atoms = (
            atoms
            .with_columns(
                chain_id=pl.concat_str(
                    pl.col('chain_id'),
                    pl.col('symmetry_mate') + 1,
                ),
                subchain_id=pl.concat_str(
                    pl.col('subchain_id'),
                    pl.col('symmetry_mate') + 1,
                ),
            )
        )

    atoms_str = (
            atoms
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

    block.write_file(str(cif_path), options)

def select_model(asym_atoms: Atoms, model_id: str) -> Atoms:
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

def make_biological_assembly(
        asym_atoms: Atoms,
        struct_assembly_gen: pl.DataFrame,
        struct_oper_map: Dict[str, Frame],
        assembly_id: str,
) -> Atoms:
    bio_opers = (
            struct_assembly_gen
            .filter(pl.col('assembly_id') == assembly_id)
    )

    if bio_opers.is_empty():
        known_assemblies = \
                struct_assembly_gen['assembly_id'].unique().to_list()

        err = MmcifError("can't find biological assembly")
        err.info = [f"known assemblies: {known_assemblies}"]
        err.blame = [f"unknown assembly: {assembly_id!r}"]
        raise err

    bio_atoms = []
    symmetry_mate = 0

    for row in bio_opers.iter_rows(named=True):
        subchain_ids = row['subchain_ids']

        for oper_ids in row['oper_ids']:
            frame = reduce(op.matmul, (struct_oper_map[x] for x in oper_ids))
            sym_atoms = (
                    transform_atom_coords(
                        asym_atoms.filter(
                            pl.col('subchain_id').is_in(subchain_ids)
                        ),
                        frame,
                    )
                    .with_columns(
                        symmetry_mate=pl.lit(symmetry_mate),
                    )
            )
            symmetry_mate += 1

            bio_atoms.append(sym_atoms)

    return pl.concat(bio_atoms).rechunk()

def get_pdb_path(pdb_dir: Union[Path, str], pdb_id: str, suffix: str = '.cif.gz'):
    """
    Return the path to a mmCIF file identified by a PDB id, assuming that the 
    file is stored in a directory matching the structure of the PDB itself.
    """
    return Path(pdb_dir) / pdb_id[1:3] / (pdb_id + suffix)

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

    loop = {
            k: [v if isinstance(v, str) else None for v in vs]
            for k, vs in cif.get_mmcif_category(f'_{key_prefix}.').items()
    }
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
    return (
            _extract_dataframe(
                cif, 'atom_site',
                schema=dict(
                    model_id=Column('pdbx_PDB_model_num'), 
                    chain_id=Column('auth_asym_id'),
                    subchain_id=Column('label_asym_id'),
                    entity_id=Column('label_entity_id'),
                    alt_id=Column('label_alt_id'),
                    seq_id=Column('label_seq_id', dtype=int),
                    seq_label_1=Column('auth_seq_id'),
                    seq_label_2=Column('pdbx_PDB_ins_code'),
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
            .with_columns(
                pl.concat_str(
                    'seq_label_1',
                    'seq_label_2',
                    ignore_nulls=True,
                ).alias('seq_label').replace({'': None}),

                # All of the elements in the PDB are uppercase anyways, but it 
                # doesn't hurt to make sure.
                pl.col('element').str.to_uppercase(),

                # Some structures (e.g. 1mno) have atoms with negative 
                # occupancies.  I'm not aware of any structures with 
                # occupancies greater than 1, but if they exist, such values 
                # also wouldn't make any sense.  
                #
                # While there's some argument for leaving these nonsensical 
                # values so the user can deal with them how they want, I think 
                # that most users will simply not realize that this could be a 
                # problem at all.  Clipping these values may not be exactly 
                # what the user wants, but it will never be a crazy thing to 
                # do, and it has the potential to avoid subtle bugs.  Overall, 
                # I think it's worth doing.
                pl.col('occupancy').clip(0, 1),
            )
            .drop(
                'seq_label_1',
                'seq_label_2',
            )
    )

def _extract_struct_assembly(cif):
        return _extract_dataframe(
                cif, 'pdbx_struct_assembly',
                schema=dict(
                    id=Column('id', required=True),
                    details=Column('details'),
                    polymer_count=Column('oligomeric_count', dtype=int),
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
                # The `to_list()` call is necessary because old versions of 
                # polars fail to successfully infer that this column should 
                # have dtype `List[String]`, not `Object`.
                subchain_ids=asym_atoms['subchain_id'].unique().to_list(),
                oper_ids=[['1']],
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
                        oper_ids=Column('oper_expression', required=True),
                    ),
                )
                .with_columns(
                    pl.col('subchain_ids').str.split(','),
                    pl.col('oper_ids').map_elements(
                        _parse_oper_expression,
                        return_dtype=pl.List(pl.List(str)),
                    ),
                )
        )

        for row in struct_assembly_gen.iter_rows(named=True):
            oper_ids = set(flatten(row['oper_ids']))
            if unknown_oper_ids := oper_ids - set(struct_oper_map):
                err = MmcifError("missing assembly operation")
                err.info = [f"assembly: {row['assembly_id']}"]
                err.blame = [f"operations referenced but not defined: {','.join(x for x in sorted(unknown_oper_ids))}"]
                raise err

    return struct_assembly_gen, struct_oper_map

def _extract_entities(cif):
    return _extract_dataframe(
            cif, 'entity',
            schema=dict(
                id=Column('id', dtype=str, required=True),
                type=Column('type', dtype=str),
                formula_weight_Da=Column('formula_weight', dtype=float),
            ),
    )

def _extract_polymers(cif):
    return (
            _extract_dataframe(
                cif, 'entity_poly',
                schema=dict(
                    entity_id=Column('entity_id', dtype=str, required=True),
                    strand_ids=Column('pdbx_strand_id', dtype=str),
                    type=Column('type', dtype=str),
                    sequence=Column('pdbx_seq_one_letter_code', dtype=str),
                    sequence_canonical=Column('pdbx_seq_one_letter_code_can', dtype=str),
                    has_nonstandard_linkage=Column('nstd_linkage', dtype=str),
                    has_nonstandard_monomer=Column('nstd_monomer', dtype=str),
                ),
            )
            .with_columns(
                pl.col(
                    'has_nonstandard_linkage',
                    'has_nonstandard_monomer',
                ).replace_strict({
                    'n': False,
                    'no': False,
                    'y': True,
                    'yes': True,
                }),
                pl.col('sequence', 'sequence_canonical').str.replace('\n', ''),
                pl.col('strand_ids').str.split(','),
            )
    )

def _parse_oper_expression(expr: str):
    parser = _make_oper_expression_parser()

    try:
        return parser.parse(expr)

    except ParseError as err1:
        err2 = MmcifError("failed to parse assembly operation expression")
        err2.info = [f"expression: {expr}"]
        err2.blame = [str(err1)]
        raise err2

@functools.lru_cache(maxsize=None)
def _make_oper_expression_parser():
    from parsy import generate, string, regex

    # This is more restrictive than the regular expression given by the 
    # PDBx/mmCIF dictionary for `_pdbx_struct_oper_list.id` [1].  However, that 
    # regular expression allows hyphens and commas, which would be ambiguous in 
    # the context of these expressions.  Alpha-numeric identifiers are enough 
    # for every structure in the PDB.
    #
    # [1]: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_pdbx_struct_oper_list.id.html
    oper_id = regex(r'\w+').desc('oper id')

    # It's not totally clear how non-numeric operation id should be treated by 
    # ranges, so we only allow numeric ids in ranges.
    oper_id_int = regex(r'[0-9]+').desc('oper id')

    @generate
    def oper_id_range():
        a = yield oper_id_int
        yield string('-')
        b = yield oper_id_int
        return [str(x) for x in range(int(a), int(b)+1)]

    @generate
    def oper_list():
        oper_id_term = oper_id_range | oper_id.map(lambda x: [x])
        oper_ids = yield oper_id_term.sep_by(string(','))
        return [[x] for x in flatten(oper_ids)]

    oper_group = string('(') >> oper_list << string(')')

    @generate
    def oper_cartesian_product():
        term = yield oper_group
        terms = []

        while term:
            terms.append(term)
            term = yield oper_group.optional()

        return [
                sum(x, start=[])
                for x in product(*terms)
        ]

    return oper_cartesian_product | oper_list

@dataclass
class Column:
    name: str
    dtype: pl._typing.PolarsDataType = str
    required: bool = False

class MmcifError(TidyError):
    pass
