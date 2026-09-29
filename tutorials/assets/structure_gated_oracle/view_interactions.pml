# Usage: pymol structure.cif tutorials/assets/structure_gated_oracle/view_interactions.pml
select ligand, chain B and organic
select designated_receptor_atoms, chain A and ((resi 71 and name OE1+OE2) or (resi 109+168 and name N))
show sticks, ligand or designated_receptor_atoms
color marine, ligand
color salmon, designated_receptor_atoms
distance glu71_hbond, ligand and donors, chain A and resi 71 and name OE1+OE2, mode=2
distance met109_hbond, ligand and acceptors, chain A and resi 109 and name N, mode=2
distance asp168_hbond, ligand and acceptors, chain A and resi 168 and name N, mode=2
select dfg_marker, chain A and (resi 71+155+169 and name CA)
show spheres, dfg_marker
set sphere_scale, 0.35, dfg_marker
zoom ligand or designated_receptor_atoms, 8
