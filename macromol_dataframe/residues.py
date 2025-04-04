import polars as pl

def assign_residue_ids(atoms, *, drop_null_ids=True, maintain_order=False):
    """
    Assign a unique numeric identifier to each residue.

    Arguments:
        atoms:
            A dataframe of atom coordinates.  Typically this dataframe would be 
            loaded by :func:`read_mmcif()` or similar, but only the following 
            columns are used:

            - ``model_id`` (optional)
            - ``symmetry_mate`` (optional)
            - ``subchain_id``
            - ``seq_id``

            Two atoms are considered to belong to the same residue if and only 
            if they have matching values for all of the above identifiers.

        drop_null_ids:
            If True, atoms that have a null ``seq_id`` will be dropped from the 
            dataframe.  These atoms (e.g. solvent molecules) don't belong to 
            any residue.  If False, these atoms may be grouped into a single 
            "residue".

    Returns:
        A dataframe with the same columns as `atoms`, plus a new column called 
        `residue_id`.

    Note that the residues identified by this function are not necessarily 
    *protein* residues, and are not necessarily complete, so they are not 
    guaranteed to have Cα atoms or anything like that.
    """
    
    id_cols = []

    if 'model_id' in atoms.columns:
        id_cols += ['model_id']
    if 'symmetry_mate' in atoms.columns:
        id_cols += ['symmetry_mate']

    id_cols += ['subchain_id', 'seq_id']

    atoms = atoms.lazy()

    if drop_null_ids:
        atoms = atoms.drop_nulls('seq_id')

    return (
            atoms
            .group_by(id_cols, maintain_order=maintain_order)
            .agg(atom_cols=pl.struct(pl.col('*')))
            .with_row_index('residue_id')
            .explode('atom_cols')
            .select('residue_id', pl.col('atom_cols').struct.field('*'))
            .collect()
    )

def explode_residue_conformations(atoms, id_name='alt_id'):
    """
    Explode any single-conformation atoms belonging to residues with multiple 
    alternate conformations.

    Arguments:
        atoms:
            A dataframe of atom coordinates.  This dataframe must have 
            the following columns:

            - ``residue_id``, e.g. created by :func:`assign_residue_ids()`.
            - ``alt_id``

        id_name:
            The name to use for the column containing the exploded alternate 
            location ids.  By default, the exploded ids will overwrite the 
            original ids, but sometime it can be useful to have both.

    Structures can use alternate location identifiers to specify multiple 
    conformations for a single atom.  Within a residue (and sometime between 
    nearby residues), atoms with the same alternate id are meant to be part of 
    the same conformation, and atoms with no alternate id are part of all 
    conformations.  

    The purpose of this function is to allow you to group by each discrete 
    residue conformation.  This requires making one copy of each "no id" atom 
    for each different conformation the residue can adopt, then giving the 
    appropriate alternate id to each copy.
    """

    # This logic fails for empty inputs, due to pola-rs/polars#22006.  
    # Fortunately, it's easy to work around this.
    if atoms.is_empty():
        return atoms.with_columns(pl.col('alt_id').alias(id_name))

    return (
            atoms
            .lazy()

            # Make a list of all the alternate conformations within each 
            # residue:
            .with_columns(
                    alt_ids=pl.col('alt_id')
                        .over('residue_id', mapping_strategy='join')
                        .list.drop_nulls()
                        .list.unique()
            )

            # Specially handle residues with no alternate conformations:
            .with_columns(
                    alt_ids=pl.when(pl.col('alt_ids').list.len() == 0)
                      .then([None])
                      .otherwise(pl.col('alt_ids'))
            )

            # Make a list of all the alternate conformations that each atom 
            # belongs to: 
            .with_columns(
                    pl.when(pl.col('alt_id').is_null())
                      .then(pl.col('alt_ids'))
                      .otherwise(pl.concat_list('alt_id'))
                      .alias(id_name)
            )

            .drop('alt_ids')
            .explode(id_name)
            .collect()
    )

