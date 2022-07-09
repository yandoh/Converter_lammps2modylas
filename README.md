# Converter_lammps2modylas
Converter from lammps input (*.data) to modylas input (.mdff, .mdxyz)

This Pthon script enables users to convert from a LAMMPS's .data file created by the Winmostar to .mdff and .mdxyz files to be read by the MODYLAS.

Use of Python 3 is presumed.

About the Winmostar, see https://winmostar.com/jp/.

[IMPORTANT NOTICE]

Present version of the script supports one-component system only.

Multi-component system will be supported later.

---
Step 1) Convsertion from .data to .mdff and .mdxyz

>python Convert_gaff_nspecies.py sessionname.data 

If distance constraints are introduced for X-H bonds, enter the following command,

>python Convert_gaff_nspecies.py sessionname.data shake

Then, sessionname.mdff and sessionname.mdxyz will be created.

---
Step 2) Convsertion from .mdff to .mdff.rearranged (also, .mdxyz to .mdxyz.rearranged)

>python Convert_mdff_nspecies.py sessionname.mdff   

`Note 1`: This .mdff file is the resultant in step 1).

`Note 2`: The command is common regardless of distance constraints or not.

Then, sessionname.mdff.rearranged and sessionname.mdxyz.rearranged will be created.

`Note 3`: In these .rearranged files, order of atoms is rearranged so that atom's in the same segment is serially packed in <segment> tags.
This is essential to run MD calcualion with distance constraints by modylas software. 
  
Further, information files to be read by Fortran programs to analyze MD trajectories, managed by another repository "MD_analysis" will be automatically created.
