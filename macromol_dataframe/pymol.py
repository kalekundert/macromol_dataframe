import polars as pl
from textwrap import dedent

def from_pymol(sele='all', state=0):
    from pymol import cmd
    rows = []
    cmd.iterate_state(
            state, sele, dedent('''\
            rows.append(
                dict(
                    model_id=state,
                    chain_id=chain,
                    subchain_id=segi,
                    alt_id=alt,
                    seq_id=resi,
                    comp_id=resn,
                    atom_id=name,
                    element=elem,
                    x=x,
                    y=y,
                    z=z,
                    occupancy=q,
                    b_factor=b,
                )
            )'''),
            space=locals(),
    )
    return pl.DataFrame(rows)

def set_ascii_dataframe_format():
    """
    Instruct polars to render dataframes in a format that works nicely with 
    PyMOL.
    """
    pl.Config.set_tbl_formatting('ASCII_FULL_CONDENSED')
    pl.Config.set_tbl_width_chars(240)
    pl.Config.set_tbl_cols(15)
