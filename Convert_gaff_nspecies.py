#Copyright (c) 2022 Yoshimichi ANDOH
#Released under the MIT license
#https://opensource.org/licenses/mit-license.php

#!/usr/bin/env python3
import sys

#
# read arguments for .py
#
arguments=sys.argv
#print(len(arguments))
if len(arguments) == 1:
  print('You do NOT specified .data file.')
  sys.exit()

#
# Enter filename
#
filename=arguments[1]
temp=filename.split(".")
sessionname=temp[0]

#
# Enter shake_switch
#
shake_switch=0
if len(arguments) == 3:
  check=arguments[2]
  if check=="shake":
    shake_switch=1
  else:
    print('Second argument must be \"shake\"')
    sys.exit()
#print(shake_switch)
#sys.exit()

with open(filename, mode="r", encoding="utf-8") as f:
#with open("./thf_liquid.data", mode="r", encoding="utf-8") as f:
  read_data = list(f)

length_data=len(read_data)
#print(length_data)

#debug: check input2
#for i in range(0,length_data):
# print(read_data[i])

#
### open file for output
#
f_mdff=open(sessionname+'.mdff','w')
f_mdxyz=open(sessionname+'.mdxyz','w')


### input atoms, bonds, angles, dihedrals, and impropers
for i in range(1,6):
  parts = read_data[i].split()
  if i==1:
    natoms=parts[0]
  elif i==2:
    nbonds=parts[0]
  elif i==3:
    nangles=parts[0]
  elif i==4:
    ndihedrals=parts[0]
  elif i==5:
    nimpropers=parts[0]

#debug info
#print(natoms,nbonds,nangles,ndihedrals,nimpropers)
#sys.exit()
  
### input types of atoms, bonds, angles, dihedrals, and impropers
for i in range(7,12):
  parts = read_data[i].split()
  if i==7:
    natom_types=parts[0]
  if i==8:
    nbond_types=parts[0]
  if i==9:
    nangle_types=parts[0]
  if i==10:
    ndihedral_types=parts[0]
  if i==11:
    nimproper_types=parts[0]
# print(parts)

#debug info
#print(natom_types,nbond_types,nangle_types,ndihedral_types,nimproper_types)
#sys.exit()

### input boxsize
for i in range(13,16):
  parts = read_data[i].split()
  if i==13:
    lx=parts[1]
  if i==14:
    ly=parts[1]
  if i==15:
    lz=parts[1]

#debug info
#print(lx,ly,lz)
#sys.exit()

line_id_mass=0
line_id_atom=0
line_id_bond=0
line_id_angle=0
line_id_dihed=0
line_id_improper=0
line_id_paircoef=0
line_id_bondcoef=0
line_id_anglecoef=0
line_id_dihedcoef=0
line_id_impropercoef=0
### input mass ###
for i in range(17,length_data):
  parts = read_data[i].split()
  if "Masses" in parts:
#   print(i,"Masses")
    line_id_mass=i
  elif "Coeffs" in parts:
    if "Pair" in parts:
#     print(i,"P-Ceff")
      line_id_paircoef=i
    elif "Bond" in parts:
#     print(i,"Bnd-Ceff")
      line_id_bondcoef=i
    elif "Angle" in parts:
#     print(i,"Ang-Ceff")
      line_id_anglecoef=i
    elif "Dihedral" in parts:
#     print(i,"Dih-Ceff")
      line_id_dihedcoef=i
    elif "Improper" in parts:
#     print(i,"Imp-Ceff")
      line_id_impropercoef=i
  elif "Atoms" in parts:
#   print(i,"Atoms")
    line_id_atom=i
  elif "Bonds" in parts:
#   print(i,"Bonds")
    line_id_bond=i
  elif "Angles" in parts:
#   print(i,"Angles")
    line_id_angle=i
  elif "Dihedrals" in parts:
#   print(i,"Diheds")
    line_id_dihed=i
  elif "Impropers" in parts:
#   print(i,"Imps")
    line_id_improper=i

#debug info#
#print(line_id_mass,line_id_atom,line_id_bond,line_id_angle,line_id_dihed,line_id_improper)
#print(line_id_paircoef,line_id_bondcoef,line_id_anglecoef,line_id_dihedcoef,line_id_impropercoef)
#sys.exit()

### input Masses information ###
start_lno=line_id_mass+2
end_lno=start_lno+int(natom_types)
mass_list=[]
typelno_to_typeid=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[0]),float(parts[1]),parts[3])
  mass_list.append(temp)
  temp=(i-start_lno,int(parts[0]),parts[3])
  typelno_to_typeid.append(temp)
  
#print(mass_list)
#print(typelno_to_typeid)
#sys.exit()

# list to convert from type_id to type_line_no
typeid_to_typelno=[]
end_lno=int(natom_types)
for itypeid in range(1,end_lno+1):
  for j in range(0,end_lno):
    temp=typelno_to_typeid[j]
    jlineno=temp[0]
    jtypeid=temp[1]
    jatomname=temp[2]
    if itypeid==jtypeid:
      temp=(itypeid,jlineno,jatomname)
      typeid_to_typelno.append(temp)
      break

#debug info
#for i in range(0,end_lno):
# print(i,typeid_to_typelno[i])
#sys.exit()

### input pair coeffs information ###
start_lno=line_id_paircoef+2
end_lno=start_lno+int(natom_types)
paircoef_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[0]),float(parts[1]),float(parts[2]))
  paircoef_list.append(temp)
  
#print(paircoef_list)


### input bond coeffs information ###
start_lno=line_id_bondcoef+2
end_lno=start_lno+int(nbond_types)
bondcoef_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[0]),float(parts[1]),float(parts[2]))
  bondcoef_list.append(temp)
  
#print(bondcoef_list)


### input angle coeffs information ###
start_lno=line_id_anglecoef+2
end_lno=start_lno+int(nangle_types)
anglecoef_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[0]),float(parts[1]),float(parts[2]))
  anglecoef_list.append(temp)
  
#print(anglecoef_list)


### input dihedral coeffs information ###
start_lno=line_id_dihedcoef+2
end_lno=start_lno+int(ndihedral_types)
dihedcoef_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[0]),float(parts[1]),int(parts[2]),int(parts[3]),float(parts[4]))
  dihedcoef_list.append(temp)
  
#print(dihedcoef_list)

### input improper coeffs information ###
start_lno=line_id_impropercoef+2
end_lno=start_lno+int(nimproper_types)
impropercoef_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[0]),float(parts[1]),int(parts[2]),int(parts[3]))
  impropercoef_list.append(temp)
  
#print(impropercoef_list)

### input Atoms information ###
start_lno=line_id_atom+2
end_lno=start_lno+int(natoms)
charge_list=[]
atomid_to_typeid=[]
xyz_list=[]
#
molid_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
#(1)
  temp=(int(parts[0]),float(parts[3]))
  charge_list.append(temp)
#(2)
  temp=(int(parts[0]),int(parts[2]))
  atomid_to_typeid.append(temp)
#(3)
  temp=(float(parts[4]),float(parts[5]),float(parts[6]))
  xyz_list.append(temp)
#(4) molid
  temp=(int(parts[0]),int(parts[1]),int(parts[2]),float(parts[3]))
  molid_list.append(temp)
  
## note: this eq. is correct for one-component system only! ##
nmolecules=parts[1]
#print(nmolecules)
#sys.exit()
#print(charge_list)
#print(atomid_to_typeid)
#print(xyz_list)

## count number of molecule for each molecule kind
natomlines=len(molid_list)
#print(natomlines)
atoms_in_same_mole=0
#atom_ids_list=[]
charge_list_in_mole=[]
molnumber=1
atoms_in_each_mole=[]
hatom_in_each_mole=[]
#atomids_in_each_mole=[]
charges_in_each_mole=[]
for i in range(0,natomlines):
  temp1=molid_list[i]
  molid1=temp1[1]
# atmid1=temp1[2]
  charge1=temp1[3]
  if atoms_in_same_mole == 0:
    atoms_in_same_mole = 1
    hatom_in_each_mole.append(0)
#   atom_ids_list.append(atmid1)
    charge_list_in_mole.append(charge1)
    continue
  else:
    temp0=molid_list[i-1]
    molid0=temp0[1]
#   atmid0=temp0[2]
#   print(i,i-1,molid0,molid1)
    charge0=temp0[3]
    if molid0 == molid1:
      atoms_in_same_mole = atoms_in_same_mole + 1
#     atom_ids_list.append(atmid1)
      charge_list_in_mole.append(charge1)
    else:
      atoms_in_each_mole.append(atoms_in_same_mole)
#     atomids_in_each_mole.append(atom_ids_list)
      charges_in_each_mole.append(charge_list_in_mole)
      molnumber = molnumber + 1
      atoms_in_same_mole = 1
#     atom_ids_list=[]
      charge_list_in_mole=[charge1]
      hatom_in_each_mole.append(i)
atoms_in_each_mole.append(atoms_in_same_mole)
#atomids_in_each_mole.append(atom_ids_list)
charges_in_each_mole.append(charge_list_in_mole)

#print(molnumber)
#print(atoms_in_each_mole)
#print(hatom_in_each_mole)
#
#print(len(charges_in_each_mole))
#for i in range(0,len(charges_in_each_mole)):
#  temp=charges_in_each_mole[i]
#  print(i,temp)
#sys.exit()

## count number of molecule species
species_number=0
molnumber_in_species=0
atoms_in_each_species=[]
molecules_in_each_species=[]
#
global_hatom_in_each_species=[0]
for i in range(0,molnumber):
   hatom=hatom_in_each_mole[i]
   natoms_per_mole1=atoms_in_each_mole[i]
   charges_list1=charges_in_each_mole[i]
#  print(hatom,natoms_per_mole1)
   if species_number == 0:
     species_number = 1
     molnumber_in_species= 1
     continue
   else:
     natoms_per_mole0=atoms_in_each_mole[i-1]
     charges_list0=charges_in_each_mole[i-1]
     hatom0=hatom_in_each_mole[i-1]
#    print(i,natoms_per_mole0,natoms_per_mole1)
#    if natoms_per_mole0 == natoms_per_mole1:  ## same species
     if charges_list0 == charges_list1:  ## same species
       molnumber_in_species=molnumber_in_species+1
       continue
     else:   ## different species
       atoms_in_each_species.append(natoms_per_mole0)
       molecules_in_each_species.append(molnumber_in_species)
#      print("hatom,hatom0=",hatom,hatom0)
       global_hatom_in_each_species.append(hatom)
#      molecules_in_each_species[species_number]=1
       species_number=species_number+1
       molnumber_in_species=1
atoms_in_each_species.append(natoms_per_mole1)
molecules_in_each_species.append(molnumber_in_species)

#debug
print("nspecies=",species_number)
print("nav=",atoms_in_each_species)
print("nmv=",molecules_in_each_species)
#for i in range(0,len(global_hatom_in_each_species)):
#  print(i,global_hatom_in_each_species[i],atoms_in_each_species[i])
#sys.exit()

### input Bonds information ###
start_lno=line_id_bond+2
end_lno=start_lno+int(nbonds)
bonds_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[1]),int(parts[2]),int(parts[3]))
  bonds_list.append(temp)
  
#print(bonds_list)
#print("len(bonds_list)=",len(bonds_list))


### input Angles information ###
start_lno=line_id_angle+2
end_lno=start_lno+int(nangles)
angles_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[1]),int(parts[2]),int(parts[3]),int(parts[4]))
  angles_list.append(temp)
  
#print(angles_list)


### input Dihedrals information ###
start_lno=line_id_dihed+2
end_lno=start_lno+int(ndihedrals)
dihedrals_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[1]),int(parts[2]),int(parts[3]),int(parts[4]),int(parts[5]))
  dihedrals_list.append(temp)
  
#print(dihedrals_list)


### input impropers information ###
start_lno=line_id_improper+2
end_lno=start_lno+int(nimpropers)
impropers_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# print(parts)
  temp=(int(parts[1]),int(parts[2]),int(parts[3]),int(parts[4]),int(parts[5]))
  impropers_list.append(temp)
  
#print(impropers_list)
#sys.exit()

###
### output 1st block in .mdff ###
###

#f_mdff.write("<forcefield>")

f_mdff.write("<forcefield>\n")
f_mdff.write("  type=gaff\n")
f_mdff.write("  special_divide_lj=2.0\n")
f_mdff.write("  special_divide_coulomb=1.2\n")
f_mdff.write("</forcefield>\n")
f_mdff.write("\n")
f_mdff.write("<system>\n")
for i in range(0,species_number):
  f_mdff.write("  " +str(i)+" "+str(molecules_in_each_species[i])+"\n")
#f_mdff.write("  " +str(0)+" "+str(nmolecules)+"\n")
f_mdff.write("</system>\n")
f_mdff.write("\n")

#sys.exit()

###
### output 2nd block in .mdff ###
###
natoms_per_mole=int(int(natoms)/int(nmolecules))

f_mdff.write("<topology and parameters>\n")
f_mdff.write(" nspecies ="+str(species_number)+"\n")
#f_mdff.write(" nspecies = 1\n")

for i in range(0,species_number):
  f_mdff.write("  <species>\n")
  f_mdff.write("    id ="+str(i)+"\n")
# f_mdff.write("    id = 0\n")
  f_mdff.write("    natom ="+str(atoms_in_each_species[i])+"\n")
# f_mdff.write("    natom ="+str(natoms_per_mole)+"\n")
#####
#mass
#####
  f_mdff.write("    <mass>\n")
  for j in range(0,atoms_in_each_species[i]):
    hatom=global_hatom_in_each_species[i]
    atomid=hatom+j
    temp=atomid_to_typeid[atomid]
    typeid=temp[1]
    for k in range(0,len(mass_list)):
      temp2=mass_list[k]
      if(typeid==temp2[0]):
        value=temp2[1]
        break
    f_mdff.write("    "+str(value)+" # "+str(temp2[2])+"\n")
  f_mdff.write("    </mass>\n")
########
#epsilon
########
  f_mdff.write("    <epsilon>\n")
  for j in range(0,atoms_in_each_species[i]):
    hatom=global_hatom_in_each_species[i]
    atomid=hatom+j
    temp=atomid_to_typeid[atomid]
    typeid=temp[1]
    for k in range(0,len(paircoef_list)):
      temp2=paircoef_list[k]
      if(typeid==temp2[0]):
        value=temp2[1]
        break
    f_mdff.write("    "+str(value)+"\n")
  f_mdff.write("    </epsilon>\n")
####
# r
####
  f_mdff.write("    <r>\n")
  for j in range(0,atoms_in_each_species[i]):
    hatom=global_hatom_in_each_species[i]
    atomid=hatom+j
    temp=atomid_to_typeid[atomid]
    typeid=temp[1]
    for k in range(0,len(paircoef_list)):
      temp2=paircoef_list[k]
      if(typeid==temp2[0]):
        value=temp2[2]
        break
    f_mdff.write("    "+str(value)+"\n")
  f_mdff.write("    </r>\n")
#########
# charge
#########
  f_mdff.write("    <charge>\n")
  for j in range(0,atoms_in_each_species[i]):
    hatom=global_hatom_in_each_species[i]
    atomid=hatom+j
    temp=charge_list[atomid]
    value=temp[1]
    f_mdff.write("    "+str(value)+"\n")
  f_mdff.write("    </charge>\n")
###############
# lj void pair
###############
  nvoid=0
  f_mdff.write("    <lj void pair>\n")
  end_lno=len(bonds_list)
  hatom=global_hatom_in_each_species[i]
  eatom=global_hatom_in_each_species[i]+atoms_in_each_species[i]
# print("hatm,eatom=",hatom,eatom)
# print("bonds_list=",len(bonds_list))
  for m in range(0,end_lno):
    temp=bonds_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
#   iacheck=temp[1]+hatom-1
#   jacheck=temp[2]+hatom-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
#   print("iacheck,jacheck=",iacheck,jacheck)
    ia=temp[1]-1 -hatom
    ja=temp[2]-1 -hatom
    f_mdff.write("     "+str(ia)+" "+str(ja)+"\n")
    nvoid+=1
#
# print("angles_list=",len(angles_list))
  end_lno=len(angles_list)
  for m in range(0,end_lno):
    temp=angles_list[m]
    iacheck=temp[1]-1
    kacheck=temp[3]-1
#   iacheck=temp[1]+hatom-1
#   kacheck=temp[3]+hatom-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if kacheck < hatom or kacheck >= eatom:
      continue
#   print("iacheck,kacheck=",iacheck,kacheck)
    ia=temp[1]-1 -hatom
    ka=temp[3]-1 -hatom
    f_mdff.write("     "+str(ia)+" "+str(ka)+"\n")
    nvoid+=1
#
  f_mdff.write("    </lj void pair>\n")
  f_mdff.write("    nvoidpair_lj ="+str(nvoid)+"\n")
##################
# lj special pair
##################
  f_mdff.write("    <lj special pair>\n")
  f_mdff.write("    </lj special pair>\n")
  f_mdff.write("    nspecialpair_lj =0\n")
####################
# coulomb void pair
####################
  nvoid=0
  f_mdff.write("    <coulomb void pair>\n")
  hatom=global_hatom_in_each_species[i]
  eatom=global_hatom_in_each_species[i]+atoms_in_each_species[i]
# print("hatm,eatom=",hatom,eatom)
  end_lno=len(bonds_list)
  for m in range(0,end_lno):
    temp=bonds_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
#   iacheck=temp[1]+hatom-1
#   jacheck=temp[2]+hatom-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
#   print("iacheck,jacheck=",iacheck,jacheck)
    ia=temp[1]-1 -hatom
    ja=temp[2]-1 -hatom
    f_mdff.write("     "+str(ia)+" "+str(ja)+"\n")
    nvoid+=1
#
  end_lno=len(angles_list)
  for m in range(0,end_lno):
    temp=angles_list[m]
    iacheck=temp[1]-1
    kacheck=temp[3]-1
#   iacheck=temp[1]+hatom-1
#   kacheck=temp[3]+hatom-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if kacheck < hatom or kacheck >= eatom:
      continue
    ia=temp[1]-1 -hatom
    ka=temp[3]-1 -hatom
    f_mdff.write("     "+str(ia)+" "+str(ka)+"\n")
    nvoid+=1
#
  f_mdff.write("    </coulomb void pair>\n")
  f_mdff.write("    nvoidpair_coulomb ="+str(nvoid)+"\n")
#############
# shake pair
#############
  shake_on_list=[]
  shake_off_list=[]
  nshake_on=0
  nshake_off=0
  hatom=global_hatom_in_each_species[i]
  eatom=global_hatom_in_each_species[i]+atoms_in_each_species[i]
# print("hatom,eatom=",hatom,eatom)
  end_lno=len(bonds_list)
  for m in range(0,end_lno):
    temp=bonds_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
#   iacheck=temp[1]+hatom-1
#   jacheck=temp[2]+hatom-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
    btype=temp[0]-1
    ia=temp[1]-1
    ja=temp[2]-1
    tempbond=bondcoef_list[btype]
    kbond=tempbond[1] ## K_bond
    lbond=tempbond[2] ## b0
#   tempi=atomid_to_typeid[ia+hatom]
#   tempj=atomid_to_typeid[ja+hatom]
    tempi=atomid_to_typeid[ia]
    tempj=atomid_to_typeid[ja]
    itype=tempi[1]-1
    jtype=tempj[1]-1
#
    tempi=typeid_to_typelno[itype]
    tempj=typeid_to_typelno[jtype]
    iatomname=tempi[2]
    jatomname=tempj[2]
#   print(i,m,iatomname,jatomname)
    temp=(int(ia),int(ja),float(kbond),float(lbond),iatomname,jatomname)
#
    hxbond=0
    if iatomname.startswith("h"):
      hxbond=1
    elif jatomname.startswith("h"):
      hxbond=1
    elif iatomname.startswith("H"):
      hxbond=1
    elif jatomname.startswith("H"):
      hxbond=1

    if hxbond==1:
      shake_on_list.append(temp)
      nshake_on+=1
    else:
      shake_off_list.append(temp)
      nshake_off+=1

# print(nshake_on,nshake_off)
# print(i,shake_on_list)
# print(i,shake_off_list)

  end_lno=int(atoms_in_each_species[i])
  flg=[]
  for m in range(0,end_lno):
    flg.append(0)
  shake_group_global=[]
  shake_group_local=[]
  for m in range(0,end_lno):
    if flg[m]==1:
      continue
    shake_group_local.append(m)
    flg[m]==1
    for j in range(0,nshake_on):
      temp=shake_on_list[j]
      ia=temp[0]-hatom
      ja=temp[1]-hatom
      if ia==m:
        shake_group_local.append(ja)
        flg[ja]=1
    shake_group_global.append(shake_group_local)
    shake_group_local=[]
  nshakegroup=len(shake_group_global)
#
  print(nshakegroup)
  print(shake_group_global)
# print("nshake_on,nshake_off=",nshake_on,nshake_off)
# print(len(shake_on_list),len(shake_off_list))
# sys.exit()

  if shake_switch == 1:
  # <shake>
    f_mdff.write("    <shake pair>\n")
    for m in range(0,nshake_on):
      temp=shake_on_list[m]
      ia=temp[0]-hatom
      ja=temp[1]-hatom
      f_mdff.write("     "+str(ia)+" "+str(ja)                 +" "+str(temp[3])+" # "+str(temp[4])+" "+str(temp[5])+"\n")
#     f_mdff.write("     "+str(temp[0])+" "+str(temp[1])                 +" "+str(temp[3])+" # "+str(temp[4])+" "+str(temp[5])+"\n")
    f_mdff.write("    </shake pair>\n")
  # <bond>
    f_mdff.write("    nbond="+str(nshake_off)+"\n")
    f_mdff.write("    <bond>\n")
    for m in range(0,nshake_off):
      temp=shake_off_list[m]
      ia=temp[0]-hatom
      ja=temp[1]-hatom
      f_mdff.write("     "+str(ia)+" "+str(ja)+" "+str(temp[2])+" "+str(temp[3])+" # "+str(temp[4])+" "+str(temp[5])+"\n")
#     f_mdff.write("     "+str(temp[0])+" "+str(temp[1])+" "+str(temp[2])+" "+str(temp[3])+" # "+str(temp[4])+" "+str(temp[5])+"\n")
    f_mdff.write("    </bond>\n")
  else:
  # <shake> (empty)
    f_mdff.write("    <shake pair>\n")
    f_mdff.write("    </shake pair>\n")
  # <bond>
#   end_lno=int(atoms_in_each_species[i])
#   f_mdff.write("    nbond="+str(end_lno)+"\n")
    f_mdff.write("    nbond="+str(nshake_on+nshake_off)+"\n")
    f_mdff.write("    <bond>\n")
#   for m in range(0,end_lno):
#     temp=bonds_list[m]
#     btype=temp[0]-1
#     ia=temp[1]-1 #-hatom
#     ja=temp[2]-1 #-hatom
#     temp=bondcoef_list[btype]
#     kbond=temp[1]
#     lbond=temp[2]
#     f_mdff.write("     "+str(ia)+" "+str(ja)+" "+str(kbond)+" "+str(lbond)+"\n")
    for m in range(0,nshake_on):
      temp=shake_on_list[m]
      ia=temp[0]-hatom
      ja=temp[1]-hatom
      f_mdff.write("     "+str(ia)+" "+str(ja)+" "+str(temp[2])+" "+str(temp[3])+" # "+str(temp[4])+" "+str(temp[5])+"\n")
#     f_mdff.write("     "+str(temp[0])+" "+str(temp[1])+" "+str(temp[2])+" "+str(temp[3])+" # "+str(temp[4])+" "+str(temp[5])+"\n")
    for m in range(0,nshake_off):
      temp=shake_off_list[m]
      ia=temp[0]-hatom
      ja=temp[1]-hatom
      f_mdff.write("     "+str(ia)+" "+str(ja)+" "+str(temp[2])+" "+str(temp[3])+" # "+str(temp[4])+" "+str(temp[5])+"\n")
#     f_mdff.write("     "+str(temp[0])+" "+str(temp[1])+" "+str(temp[2])+" "+str(temp[3])+" # "+str(temp[4])+" "+str(temp[5])+"\n")
    f_mdff.write("    </bond>\n")
########
# angle
########
  nangles_in_each_species=0
  end_lno=len(angles_list)
  for m in range(0,end_lno):
    temp=angles_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
    kacheck=temp[3]-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
    if kacheck < hatom or kacheck >= eatom:
      continue
    nangles_in_each_species+=1
#
  f_mdff.write("    nangle="+str(nangles_in_each_species)+"\n")
  f_mdff.write("    <angle>\n")
  for m in range(0,end_lno):
    temp=angles_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
    kacheck=temp[3]-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
    if kacheck < hatom or kacheck >= eatom:
      continue
    btype=temp[0]-1
    ia=temp[1]-1
    ja=temp[2]-1
    ka=temp[3]-1
    temp=anglecoef_list[btype]
    kangle=temp[1]
    langle=temp[2]
    f_mdff.write("     "+str(ia-hatom)+" "+str(ja-hatom)+" "+str(ka-hatom)+" "+str(kangle)+" "+str(langle)+"\n")
  f_mdff.write("    </angle>\n")
############
# dihedrals 
############
  ndihedrals_in_each_species=0
  end_lno=len(dihedrals_list)
  for m in range(0,end_lno):
    temp=dihedrals_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
    kacheck=temp[3]-1
    lacheck=temp[4]-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
    if kacheck < hatom or kacheck >= eatom:
      continue
    if lacheck < hatom or lacheck >= eatom:
      continue
    ndihedrals_in_each_species+=1

  f_mdff.write("    ndihedral="+str(ndihedrals_in_each_species)+"\n")
  f_mdff.write("    <dihedral>\n")
  for m in range(0,end_lno):
    temp=dihedrals_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
    kacheck=temp[3]-1
    lacheck=temp[4]-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
    if kacheck < hatom or kacheck >= eatom:
      continue
    if lacheck < hatom or lacheck >= eatom:
      continue
    btype=temp[0]-1
    ia=temp[1]-1
    ja=temp[2]-1
    ka=temp[3]-1
    la=temp[4]-1
    temp=dihedcoef_list[btype]
    kdihed=temp[1]
    ndihed=temp[2]
    pdihed=temp[3]
    f_mdff.write("     "+str(ia-hatom)+" "+str(ja-hatom)+" "+str(ka-hatom)+" "+str(la-hatom)+" "+str(kdihed)+" "+str(ndihed)+" "+str(pdihed)+"\n")
  f_mdff.write("    </dihedral>\n")
###################
# improper torsion
###################
  nimpropers_in_each_species=0
  end_lno=len(impropers_list)
  for m in range(0,end_lno):
    temp=impropers_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
    kacheck=temp[3]-1
    lacheck=temp[4]-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
    if kacheck < hatom or kacheck >= eatom:
      continue
    if lacheck < hatom or lacheck >= eatom:
      continue
    nimpropers_in_each_species+=1

  f_mdff.write("    nitorsion="+str(nimpropers_in_each_species)+"\n")
  f_mdff.write("    <itorsion>\n")
  for m in range(0,end_lno):
    temp=impropers_list[m]
    iacheck=temp[1]-1
    jacheck=temp[2]-1
    kacheck=temp[3]-1
    lacheck=temp[4]-1
    if iacheck < hatom or iacheck >= eatom:
      continue
    if jacheck < hatom or jacheck >= eatom:
      continue
    if kacheck < hatom or kacheck >= eatom:
      continue
    if lacheck < hatom or lacheck >= eatom:
      continue
    btype=temp[0]-1
    ia=temp[1]-1
    ja=temp[2]-1
    ka=temp[3]-1
    la=temp[4]-1
    temp=impropercoef_list[btype]
    kimprp=temp[1]
    dimprp=temp[2]
    nimprp=temp[3]
    f_mdff.write("     "+str(ia-hatom)+" "+str(ja-hatom)+" "+str(ka-hatom)+" "+str(la-hatom)+" "+str(kimprp)+" "+str(dimprp)+" "+str(nimprp)+"\n")
  f_mdff.write("    </itorsion>\n")
###########
# segments
###########
  f_mdff.write("    <segments>\n")
  f_mdff.write("      nsegment= "+str(nshakegroup)+"\n")
  for m in range(0,nshakegroup):
    temp=shake_group_global[m]
    natoms_in_segment=len(temp)
    f_mdff.write("      <segment>\n")
    f_mdff.write("        ID= "+str(m)+"\n")
    f_mdff.write("        natom="+str(natoms_in_segment)+"\n")
    f_mdff.write("        <atom>\n")
    for j in range(0,natoms_in_segment):
      f_mdff.write("           "+str(temp[j])+"\n")
    f_mdff.write("        </atom>\n")
    f_mdff.write("      </segment>\n")
  f_mdff.write("    </segments>\n")
#
  f_mdff.write("  </species>\n")

f_mdff.write("</topology and parameters>\n")
f_mdff.close()

###
### output 1st block in .mdxyz ###
###
end_lno=int(natoms)
f_mdxyz.write("<atom>\n")
f_mdxyz.write("  natom="+str(natoms)+"\n")
f_mdxyz.write("  <positions>\n")
for i in range(0,end_lno):
  temp=xyz_list[i]
  x=temp[0]
  y=temp[1]
  z=temp[2]
  f_mdxyz.write("  "+str(x)+" "+str(y)+" "+str(z)+"\n")
f_mdxyz.write("  </positions>\n")
f_mdxyz.write("  <velocities>\n")
f_mdxyz.write("  </velocities>\n")
f_mdxyz.write("</atom>\n")

###
### output 2nd block in .mdxyz ###
###
end_lno=5
f_mdxyz.write("<thermostat>\n")
f_mdxyz.write("  nthermostat=5\n")
f_mdxyz.write("  <positions>\n")
for i in range(0,end_lno):
  f_mdxyz.write("    "+str(0.0)+"\n")
f_mdxyz.write("  </positions>\n")
f_mdxyz.write("  <velocities>\n")
for i in range(0,end_lno):
  f_mdxyz.write("    "+str(0.0)+"\n")
f_mdxyz.write("  </velocities>\n")
f_mdxyz.write("</thermostat>\n")
f_mdxyz.write("<barostat>\n")
f_mdxyz.write("  nbarostat=5\n")
f_mdxyz.write("  <positions>\n")
for i in range(0,end_lno):
  f_mdxyz.write("    "+str(0.0)+"\n")
f_mdxyz.write("  </positions>\n")
f_mdxyz.write("  <velocities>\n")
for i in range(0,end_lno):
  f_mdxyz.write("    "+str(0.0)+"\n")
f_mdxyz.write("  </velocities>\n")
f_mdxyz.write("</barostat>\n")

###
### output 3rd block in .mdxyz ###
###
f_mdxyz.write("<periodic cell>\n")
f_mdxyz.write("  <length>\n")
f_mdxyz.write("    "+str(lx)+" "+str(ly)+" "+str(lz)+"\n")
f_mdxyz.write("  </length>\n")
f_mdxyz.write("  <angle>\n")
f_mdxyz.write("    "+str(90.0)+" "+str(90.0)+" "+str(90.0)+"\n")
f_mdxyz.write("  </angle>\n")
f_mdxyz.write("  <vboxg>\n")
f_mdxyz.write("    "+str(0.0)+" "+str(0.0)+" "+str(0.0)+"\n")
f_mdxyz.write("    "+str(0.0)+" "+str(0.0)+" "+str(0.0)+"\n")
f_mdxyz.write("    "+str(0.0)+" "+str(0.0)+" "+str(0.0)+"\n")
f_mdxyz.write("  </vboxg>\n")
f_mdxyz.write("</periodic cell>\n")

f_mdxyz.close()

sys.exit()
