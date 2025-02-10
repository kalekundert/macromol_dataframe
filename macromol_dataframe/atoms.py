import polars as pl

from .coords import (
        Coords3, Coords4, Frame, transform_coords, homogenize_coords,
)

from typing import Union
from typing_extensions import TypeAlias

Atoms: TypeAlias = pl.DataFrame

def prune_hydrogen(atoms):
    """
    Remove all hydrogen atoms from the structure.

    There are two reasons for doing this: 

    - The position of each hydrogen atom is strongly determined by the 
      positions of the neighboring heavy atoms, so including hydrogens doesn't 
      add much information.

    - High resolution structures are much more likely to include hydrogens, but 
      we don't want our model to be able to distinguish between high- and
      low-resolution structures.  That means either adding hydrogens to low- 
      resolution structures, or removing them from high-resolution structures.  
      The latter is easier and doesn't involve adding data that's not in the 
      actual dataset.
    """
    return atoms.filter(~pl.col('element').is_in(['H', 'D']))

def prune_water(atoms):
    """
    Remove all water molecules from the structure.

    There are two reasons for doing this:

    - High-resolution structures model many more water molecules than 
      low-resolution structures, and we don't want our model to be able to 
      discriminate on that basis.

    - Even in high-resolution structures, many of the waters that end up in 
      the model are not very highly structured, and may not be informative.  
      Removing all the water makes the input less noisy.  Plus, the presence of 
      highly-structured waters can still be inferred from the macromolecule 
      itself.
    """
    return atoms.filter(~pl.col('comp_id').is_in(['HOH', 'DOD']))


def get_atom_coords(
        atoms: Atoms,
        *,
        homogeneous: bool=False,
) -> Union[Coords3, Coords4]:
    coords = atoms.select('x', 'y', 'z').to_numpy()
    return homogenize_coords(coords) if homogeneous else coords

def replace_atom_coords(atoms: Atoms, coords: Union[Coords3, Coords4]):
    return (
            atoms.with_columns(
                x=coords[:,0],
                y=coords[:,1],
                z=coords[:,2],
            )
    )

def transform_atom_coords(
        atoms_x: Atoms,
        frame_xy: Frame,
) -> Atoms:
    coords_x = get_atom_coords(atoms_x, homogeneous=True)
    coords_y = transform_coords(coords_x, frame_xy)
    return replace_atom_coords(atoms_x, coords_y)
