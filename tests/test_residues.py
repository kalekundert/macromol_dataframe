import macromol_dataframe as mmdf
import parametrize_from_file as pff
import polars as pl

from macromol_dataframe.testing import dataframe
from polars.testing import assert_frame_equal
from pathlib import Path

@pff.parametrize(
        schema=pff.cast(
            atoms=dataframe(
                dtypes={
                    'resi': int,
                    'expected': int,
                    'symm': int,
                },
                col_aliases={
                    'model': 'model_id',
                    'symm': 'symmetry_mate',
                    'chain': 'subchain_id',
                    'resi': 'seq_id',
                    'atom': 'atom_id',
                },
            ),
        ),
)
def test_assign_residue_ids(atoms):
    # The residue ids aren't guaranteed to be any particular value, so we can't 
    # check that each residue got the "right" id.  Instead, we check that each 
    # atoms is grouped with all of the other atoms from its residue.

    atoms = (
            mmdf.assign_residue_ids(atoms)
            .with_row_index()
    )

    def extract_residues(key):
        residues = (
                atoms
                .group_by(key)
                .agg('index')
                .get_column('index')
                .to_list()
        )
        return {frozenset(x) for x in residues}

    actual = extract_residues('residue_id')
    expected = extract_residues('expected')

    assert actual == expected

def test_assign_residue_ids_4ous():
    # I included this test case because it broke an old version of my code by 
    # virtue of having three symmetry mates, each with a calcium atom.  My old 
    # code mistook these calcium atoms for a "residue" with 3 Cα atoms.  The 
    # new code is not really susceptible to this kind of mistake, but I still 
    # like having this test case.

    test_dir = Path(__file__).parent 
    cif_path = test_dir / 'pdb' / '4ous.cif.gz'

    atoms = mmdf.read_biological_assembly(cif_path, model_id='1', assembly_id='1')
    atoms = mmdf.assign_residue_ids(atoms)

    n = atoms.select(pl.col('residue_id').n_unique()).item()

    # I counted the actual number of residues by hand, in PyMOL:
    #
    # Monomer: 132
    # Assembly: 3*132 = 396
    #
    # Note that the calcium and the waters don't have sequence ids, so they are 
    # excluded.

    assert n == 396

explode_atoms = dataframe(
        exprs={
            'alt_id': pl.col('alt_id').replace({'.': None}),
        },
        dtypes={
            'res_id': int,
            'x': float,
            'y': float,
            'z': float,
        },
        col_aliases={
            'res_id': 'residue_id',
            'atom': 'atom_id',
        },
)
@pff.parametrize(
        schema=pff.cast(
            atoms=explode_atoms,
            expected=explode_atoms,
        ),
)
def test_explode_residue_conformations(atoms, expected):
    assert_frame_equal(
        mmdf.explode_residue_conformations(atoms),
        expected,
        check_row_order=False,
    )

def test_explode_residue_conformations_4rek():
    # `4rek` is a very high resolution structure, so it has a lot of alternate 
    # conformations, which makes it a good test case:

    test_dir = Path(__file__).parent 
    cif_path = test_dir / 'pdb' / '4rek.cif.gz'

    atoms = mmdf.read_asymmetric_unit(cif_path)
    atoms = mmdf.prune_hydrogen(atoms)
    atoms = mmdf.assign_residue_ids(atoms)
    atoms = mmdf.explode_residue_conformations(atoms).drop('residue_id')

    def check_residue(resi, expected_str):
        expected = dataframe(
                expected_str,
                exprs={
                    'alt_id': pl.col('alt_id').replace({'.': None}),
                },
                dtypes={
                    'seq_id': int,
                    'x': float,
                    'y': float,
                    'z': float,
                    'occupancy': float,
                },
        )
        cols = expected_str.splitlines()[0].split()
        actual = atoms.filter(seq_id=resi).select(cols)

        assert_frame_equal(
            actual,
            expected,
            check_row_order=False,
            check_column_order=False,
        )

    # G1 has no alternate conformations:

    check_residue(1, '''\
chain_id  seq_id  comp_id  alt_id  atom_id  element       x       y        z  occupancy
       A       1  GLY           .  N        N        24.622   1.567   47.504       1.00
       A       1  GLY           .  CA       C        25.726   0.861   48.124       1.00
       A       1  GLY           .  C        C        26.844   0.511   47.157       1.00
       A       1  GLY           .  O        O        26.866   0.966   46.016       1.00
''')

    # V3 has alternate conformations for Cα and all the sidechain atoms, but 
    # none of the other backbone atoms:

    check_residue(3, '''\
chain_id  seq_id  comp_id  alt_id  atom_id  element       x       y       z  occupancy
       A       3  VAL           A  N        N        29.166  -2.179  44.934       1.00
       A       3  VAL           B  N        N        29.166  -2.179  44.934       1.00
       A       3  VAL           A  CA       C        29.165  -3.431  44.195       0.44
       A       3  VAL           B  CA       C        29.160  -3.428  44.177       0.56
       A       3  VAL           A  C        C        30.472  -3.459  43.393       1.00
       A       3  VAL           B  C        C        30.472  -3.459  43.393       1.00
       A       3  VAL           A  O        O        30.896  -2.443  42.875       1.00
       A       3  VAL           B  O        O        30.896  -2.443  42.875       1.00
       A       3  VAL           A  CB       C        27.918  -3.512  43.306       0.44
       A       3  VAL           B  CB       C        27.971  -3.523  43.168       0.56
       A       3  VAL           A  CG1      C        27.915  -2.363  42.376       0.44
       A       3  VAL           B  CG1      C        27.924  -4.910  42.521       0.56
       A       3  VAL           A  CG2      C        27.833  -4.841  42.550       0.44
       A       3  VAL           B  CG2      C        26.623  -3.206  43.835       0.56
''')

    # P40 (and its downstream neighbors) have alternate conformations for every 
    # atom:
    
    check_residue(40, '''\
chain_id  seq_id  alt_id  comp_id  atom_id  element       x       y       z  occupancy
       A      40       A  PRO      N        N         2.207  -7.917  30.839       0.69
       A      40       B  PRO      N        N         1.982  -7.459  30.890       0.31
       A      40       A  PRO      CA       C         1.120  -7.727  29.885       0.69
       A      40       B  PRO      CA       C         0.949  -7.142  29.887       0.31
       A      40       A  PRO      C        C         0.419  -6.395  30.124       0.69
       A      40       B  PRO      C        C         0.399  -5.713  29.944       0.31
       A      40       A  PRO      O        O         0.082  -6.028  31.254       0.69
       A      40       B  PRO      O        O         0.147  -5.189  31.032       0.31
       A      40       A  PRO      CB       C         0.186  -8.904  30.181       0.69
       A      40       B  PRO      CB       C        -0.189  -8.110  30.234       0.31
       A      40       A  PRO      CG       C         1.086  -9.954  30.763       0.69
       A      40       B  PRO      CG       C         0.441  -9.227  30.944       0.31
       A      40       A  PRO      CD       C         2.105  -9.202  31.555       0.69
       A      40       B  PRO      CD       C         1.673  -8.702  31.624       0.31
''')

    # M50 has one backbone conformation, but two sidechain conformations:

    check_residue(50, '''\
chain_id  seq_id  alt_id  comp_id  atom_id  element       x       y       z  occupancy
       A      50       A  MET      N        N         2.551  -1.402  15.571       1.00 
       A      50       B  MET      N        N         2.551  -1.402  15.571       1.00 
       A      50       A  MET      CA       C         2.652  -0.157  14.847       1.00 
       A      50       B  MET      CA       C         2.652  -0.157  14.847       1.00 
       A      50       A  MET      C        C         1.305   0.416  14.422       1.00 
       A      50       B  MET      C        C         1.305   0.416  14.422       1.00 
       A      50       A  MET      O        O         1.213   1.604  14.187       1.00 
       A      50       B  MET      O        O         1.213   1.604  14.187       1.00 
       A      50       A  MET      CB       C         3.510  -0.294  13.590       0.19 
       A      50       B  MET      CB       C         3.643  -0.310  13.689       0.81 
       A      50       A  MET      CG       C         4.991  -0.081  13.809       0.19 
       A      50       B  MET      CG       C         5.032  -0.624  14.206       0.81 
       A      50       A  MET      SD       S         5.810  -1.653  14.052       0.19 
       A      50       B  MET      SD       S         6.276  -0.736  12.909       0.81 
       A      50       A  MET      CE       C         5.857  -2.274  12.370       0.19 
       A      50       B  MET      CE       C         5.837  -2.281  12.114       0.81 
''')
