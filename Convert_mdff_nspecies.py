#Copyright (c) 2022 Yoshimichi ANDOH
#Released under the MIT license
#https://opensource.org/licenses/mit-license.php

#!/usr/bin/env python3
import sys

#
# read arguments for .py
#
arguments=sys.argv
if len(arguments) == 1:
  print('You do NOT specified .data file.')
  sys.exit()

#
# Enter filename
#
filename=arguments[1]
temp=filename.split(".")
sessionname=temp[0]
if temp[1] != 'mdff':
  print('ERROR: Input file must be .mdff.')
  sys.exit()

#
# Enter shake_switch
#
shake_switch=0
if len(arguments) == 3:
  check=arguments[2]
  if check=="shake":
    shake_switch=1
  else:
    print('ERROR: Second argument must be \"shake\"')
    sys.exit()

with open(filename, mode="r", encoding="utf-8") as f:
  read_data = list(f)

length_data=len(read_data)

#debug: check input2
#for i in range(0,length_data):
# print(read_data[i])
#sys.exit()

#
### open file for output
#
f_mdff=open(sessionname+'.mdff.rearranged','w')

# parse parameters
nspecies=0
for i in range(0,length_data):
  parts = read_data[i].split("=")
  if "nspecies" in parts[0]:
    nspecies=int(parts[1])
#   print(nspecies)
    break
#if nspecies != 1:
#  print("ERROR: nspecies != 1")
#  sys.exit()


nvoidpair_lj=[]
nspecialpair_lj=[]
nvoidpair_coulomb=[]
nbond=[]
nangle=[]
ndihedral=[]
nitorsion=[]
nsegment=[]
for i in range(0,length_data):
  parts = read_data[i].split("=")
  if "nvoidpair_lj" in parts[0]:
    nvoidpair_lj_in=int(parts[1])
    nvoidpair_lj.append(nvoidpair_lj_in)
  if "nspecialpair_lj" in parts[0]:
    nspecialpair_lj_in=int(parts[1])
    nspecialpair_lj.append(nspecialpair_lj_in)
  if "nvoidpair_coulomb" in parts[0]:
    nvoidpair_coulomb_in=int(parts[1])
    nvoidpair_coulomb.append(nvoidpair_coulomb_in)
  if "nbond" in parts[0]:
    nbond_in=int(parts[1])
    nbond.append(nbond_in)
  if "nangle" in parts[0]:
    nangle_in=int(parts[1])
    nangle.append(nangle_in)
  if "ndihedral" in parts[0]:
    ndihedral_in=int(parts[1])
    ndihedral.append(ndihedral_in)
  if "nitorsion" in parts[0]:
    nitorsion_in=int(parts[1])
    nitorsion.append(nitorsion_in)
  if "nsegment" in parts[0]:
    nsegment_in=int(parts[1])
    nsegment.append(nsegment_in)

##debug
print(nbond,nangle,ndihedral,nitorsion)
print(nsegment)
#sys.exit()
  
# parse line numbers
line_id_system_start=0
line_id_system_end=0
line_id_nspecies=0
#
line_id_species_start=[]
line_id_species_end=[]
line_id_mass_start=[]
line_id_mass_end=[]
line_id_eps_start=[]
line_id_eps_end=[]
line_id_r_start=[]
line_id_r_end=[]
line_id_charge_start=[]
line_id_charge_end=[]
line_id_ljvoid_start=[]
line_id_ljvoid_end=[]
line_id_ljspecial_start=[]
line_id_ljspecial_end=[]
line_id_clvoid_start=[]
line_id_clvoid_end=[]
line_id_shake_start=[]
line_id_shake_end=[]
line_id_bond_start=[]
line_id_bond_end=[]
line_id_angle_start=[]
line_id_angle_end=[]
line_id_dihed_start=[]
line_id_dihed_end=[]
line_id_improper_start=[]
line_id_improper_end=[]
line_id_segment_start=[]
line_id_segment_end=[]

### input line numbers ###
for i in range(0,length_data):
  parts = read_data[i].split()
#
  if "<system>" in parts:
    line_id_system_start=i
  if "</system>" in parts:
    line_id_system_end=i
  if "nspecies" in parts:
    line_id_nspecies=i
#
  if "<species>" in parts:
    line_id_species_start.append(i)
  if "</species>" in parts:
    line_id_species_end.append(i)
  if "<mass>" in parts:
    line_id_mass_start.append(i)
  if "</mass>" in parts:
    line_id_mass_end.append(i)
  if "<epsilon>" in parts:
    line_id_eps_start.append(i)
  if "</epsilon>" in parts:
    line_id_eps_end.append(i)
  if "<r>" in parts:
    line_id_r_start.append(i)
  if "</r>" in parts:
    line_id_r_end.append(i)
  if "<charge>" in parts:
    line_id_charge_start.append(i)
  if "</charge>" in parts:
    line_id_charge_end.append(i)
  if "void" in parts:
    if "<lj" in parts:
      line_id_ljvoid_start.append(i)
    if "</lj" in parts:
      line_id_ljvoid_end.append(i)
    if "<coulomb" in parts:
      line_id_clvoid_start.append(i)
    if "</coulomb" in parts:
      line_id_clvoid_end.append(i)
  if "special" in parts:
    if "<lj" in parts:
      line_id_ljspecial_start.append(i)
    if "</lj" in parts:
      line_id_ljspecial_end.append(i)
  if "<shake" in parts:
    line_id_shake_start.append(i)
  if "</shake" in parts:
    line_id_shake_end.append(i)
  if "<bond>" in parts:
    line_id_bond_start.append(i)
  if "</bond>" in parts:
    line_id_bond_end.append(i)
  if "<angle>" in parts:
    line_id_angle_start.append(i)
  if "</angle>" in parts:
    line_id_angle_end.append(i)
  if "<dihedral>" in parts:
    line_id_dihed_start.append(i)
  if "</dihedral>" in parts:
    line_id_dihed_end.append(i)
  if "<itorsion>" in parts:
    line_id_improper_start.append(i)
  if "</itorsion>" in parts:
    line_id_improper_end.append(i)
  if "<segments>" in parts:
    line_id_segment_start.append(i)
  if "</segments>" in parts:
    line_id_segment_end.append(i)

###debug
#print(line_id_system_start,line_id_system_end)
#print(line_id_mass_start,line_id_mass_end)
#print(line_id_eps_start,line_id_eps_end)
#print(line_id_r_start,line_id_r_end)
#print(line_id_charge_start,line_id_charge_end)
#print(line_id_ljvoid_start,line_id_ljvoid_end)
#print(line_id_ljspecial_start,line_id_ljspecial_end)
#print(line_id_clvoid_start,line_id_clvoid_end)
#print(line_id_shake_start,line_id_shake_end)
#print(line_id_bond_start,line_id_bond_end)
#print(line_id_angle_start,line_id_angle_end)
#print(line_id_dihed_start,line_id_dihed_end)
#print(line_id_improper_start,line_id_improper_end)
#print(line_id_segment_start,line_id_segment_end)
#sys.exit()

#
# Output header of .mdff
#
for i in range(0,line_id_nspecies+1):
  line=read_data[i].split("\n")
# print(line)
  f_mdff.write(line[0]+str("\n"))
#sys.exit()


#
# input nmole (sys_info) #
#
start_lno=line_id_system_start+1
end_lno  =line_id_system_end
sysmol_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# temp  =(int(parts[0]),parts[1])
  temp  =int(parts[1])
  sysmol_list.append(temp)

#
# input natom (sys_info) #
#
#natom=[]
sysatom_list=[]
for i in range(0,nspecies):
  start_lno=line_id_species_start[i]+2
  end_lno=line_id_species_start[i]+3
  for m in range(start_lno,end_lno):
    parts = read_data[m].split("=")
    natomin=int(parts[1])
    sysatom_list.append(natomin)

print(sysmol_list)
print(sysatom_list)

#
# record <species>
#
for i in range(0,nspecies):
# input segments
  natoms_in_segment=[]
  line_id_atomtag_segment=[]
  start_lno=line_id_segment_start[i]
  end_lno=line_id_segment_end[i]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split("=")
    if "natom" in parts[0]:
      value=int(parts[1])
      natoms_in_segment.append(value)
      temp=[m+2,m+2+value]
      line_id_atomtag_segment.append(temp)
# atomid_origignal
  atomid_org=[]
  for m in range(0,nsegment[i]):
    temp=line_id_atomtag_segment[m]
    start_lno=temp[0]
    end_lno=temp[1]
    for j in range(start_lno,end_lno):
      value=int(read_data[j])
      atomid_org.append(value)
#
  if len(atomid_org) != sysatom_list[i]:
    print("Error len(atomid_org) != sysatom_list")
    sys.exit()
  id_old2new=[]
  for m in range(0,sysatom_list[i]):
    id_old2new.append(0)
  for m in range(0,sysatom_list[i]):
    newid=atomid_org[m]
    id_old2new[newid]=m
#
  f_mdff.write("  <species>\n")
  f_mdff.write("    id ="+str(i)+"\n")
  f_mdff.write("    natom ="+str(sysatom_list[i])+"\n")
#######
# mass
#######
  start_lno=line_id_mass_start[i]+1
  end_lno  =line_id_mass_end[i]
  mass_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(float(parts[0]),parts[2])
    mass_list.append(temp)
  new_mass_list=[]
  for m in range(0,len(mass_list)):
    new_mass_list.append(0)
  for oldid in range(0,len(mass_list)):
    newid=id_old2new[oldid]
    new_mass_list[newid]=mass_list[oldid]
  f_mdff.write("    <mass>\n")
  for m in range(0,len(new_mass_list)):
    value=new_mass_list[m]
    f_mdff.write("    "+str(value[0])+" # "+str(value[1])+"\n")
  f_mdff.write("    </mass>\n")
##########
# epsilon
##########
  start_lno=line_id_eps_start[i]+1
  end_lno  =line_id_eps_end[i]
  eps_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(float(parts[0]))
    eps_list.append(temp)
  new_eps_list=[]
  for m in range(0,len(eps_list)):
    new_eps_list.append(0)
  for oldid in range(0,len(eps_list)):
    newid=id_old2new[oldid]
    new_eps_list[newid]=eps_list[oldid]
  f_mdff.write("    <epsilon>\n")
  for m in range(0,len(new_eps_list)):
    value=new_eps_list[m]
    f_mdff.write("    "+str(value)+"\n")
  f_mdff.write("    </epsilon>\n")
####
# r
####
  start_lno=line_id_r_start[i]+1
  end_lno  =line_id_r_end[i]
  r_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(float(parts[0]))
    r_list.append(temp)
  new_r_list=[]
  for m in range(0,len(r_list)):
    new_r_list.append(0)
  for oldid in range(0,len(r_list)):
    newid=id_old2new[oldid]
    new_r_list[newid]=r_list[oldid]
  f_mdff.write("    <r>\n")
  for m in range(0,len(new_r_list)):
    value=new_r_list[m]
    f_mdff.write("    "+str(value)+"\n")
  f_mdff.write("    </r>\n")
#########
# charge
#########
  start_lno=line_id_charge_start[i]+1
  end_lno  =line_id_charge_end[i]
  charge_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(float(parts[0]))
    charge_list.append(temp)
  new_charge_list=[]
  for m in range(0,len(charge_list)):
    new_charge_list.append(0)
  for oldid in range(0,len(charge_list)):
    newid=id_old2new[oldid]
    new_charge_list[newid]=charge_list[oldid]
  f_mdff.write("    <charge>\n")
  for m in range(0,len(new_charge_list)):
    value=new_charge_list[m]
    f_mdff.write("    "+str(value)+"\n")
  f_mdff.write("    </charge>\n")
##########
# lj void
##########
  start_lno=line_id_ljvoid_start[i]+1
  end_lno  =line_id_ljvoid_end[i]
  ljvoid_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(int(parts[0]),int(parts[1]))
    ljvoid_list.append(temp)
  new_ljvoid_list=[]
  for m in range(0,len(ljvoid_list)):
    temp=ljvoid_list[m]
    oldi=temp[0]
    oldj=temp[1]
    newi=id_old2new[oldi]
    newj=id_old2new[oldj]
    new_ljvoid_list.append([newi,newj])
  f_mdff.write("    <lj void pair>\n")
  for m in range(0,len(new_ljvoid_list)):
    value=new_ljvoid_list[m]
    f_mdff.write("    "+str(value[0])+" "+str(value[1])+"\n")
  f_mdff.write("    </lj void pair>\n")
  f_mdff.write("    nvoidpair_lj="+str(nvoidpair_lj[i])+"\n")
#
  if len(new_ljvoid_list) != nvoidpair_lj[i]:
    print("ERROR: new_ljvoid_list do not match nvoidpair_lj.")
    sys.exit()
#####################################
# lj special pair (always 0 in gaff)
#####################################
  f_mdff.write("    <lj special pair>\n")
  f_mdff.write("    </lj special pair>\n")
  f_mdff.write("    nspecialpair_lj =0\n")
###############
# coulomb void
###############
  start_lno=line_id_clvoid_start[i]+1
  end_lno  =line_id_clvoid_end[i]
  clvoid_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(int(parts[0]),int(parts[1]))
    clvoid_list.append(temp)
  new_clvoid_list=[]
  for m in range(0,len(clvoid_list)):
    temp=clvoid_list[m]
    oldi=temp[0]
    oldj=temp[1]
    newi=id_old2new[oldi]
    newj=id_old2new[oldj]
    new_clvoid_list.append([newi,newj])
  f_mdff.write("    <coulomb void pair>\n")
  for m in range(0,len(new_clvoid_list)):
    value=new_clvoid_list[m]
    f_mdff.write("    "+str(value[0])+" "+str(value[1])+"\n")
  f_mdff.write("    </coulomb void pair>\n")
  f_mdff.write("    nvoidpair_coulomb="+str(nvoidpair_coulomb[i])+"\n")
#
  if len(new_clvoid_list) != nvoidpair_lj[i]:
    print("ERROR: new_clvoid_list do not match nvoidpair_coulomb.")
    sys.exit()
########
# shake
########
  start_lno=line_id_shake_start[i]+1
  end_lno  =line_id_shake_end[i]
  nshakepair=end_lno-start_lno
  shake_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(int(parts[0]),int(parts[1]),float(parts[2]),parts[4],parts[5])
    shake_list.append(temp)
  new_shake_list=[]
  for m in range(0,len(shake_list)):
    temp=shake_list[m]
    oldi=temp[0]
    oldj=temp[1]
    length=temp[2]
    namei=temp[3]
    namej=temp[4]
    newi=id_old2new[oldi]
    newj=id_old2new[oldj]
    new_shake_list.append([newi,newj,length,namei,namej])
  f_mdff.write("    <shake pair>\n")
  for m in range(0,len(new_shake_list)):
    value=new_shake_list[m]
    f_mdff.write("    "+str(value[0])+" "+str(value[1])+" "+str(value[2])+" # "+str(value[3])+" "+str(value[4])+"\n")
  f_mdff.write("    </shake pair>\n")
#
  if len(new_shake_list) != nshakepair:
    print("ERROR: new_shake_list do not match nshakepair.")
    print(len(new_shake_list), nshakepair)
    sys.exit()
#######
# bond
#######
  start_lno=line_id_bond_start[i]+1
  end_lno  =line_id_bond_end[i]
  bond_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    if len(parts) == 7:
      temp  =(int(parts[0]),int(parts[1]),float(parts[2]), \
              float(parts[3]),parts[5],parts[6])
    else:
      temp  =(int(parts[0]),int(parts[1]),float(parts[2]), \
              float(parts[3]))
    bond_list.append(temp)
  new_bond_list=[]
  for m in range(0,len(bond_list)):
    temp=bond_list[m]
    oldi=temp[0]
    oldj=temp[1]
    kbond=temp[2]
    length=temp[3]
    if len(temp)==7:
      namei=temp[4]
      namej=temp[5]
    newi=id_old2new[oldi]
    newj=id_old2new[oldj]
    if len(temp)==7:
      new_bond_list.append([newi,newj,kbond,length,namei,namej])
    else:
      new_bond_list.append([newi,newj,kbond,length])
  f_mdff.write("    nbond="+str(nbond[i])+"\n")
  f_mdff.write("    <bond>\n")
  for m in range(0,len(new_bond_list)):
    value=new_bond_list[m]
    if len(value)==7:
      f_mdff.write("    "+str(value[0])+" "+str(value[1])+" " \
                   +str(value[2])+str(value[3])+" # " \
                   +str(value[4])+" "+str(value[5])+"\n")
    else:
      f_mdff.write("    "+str(value[0])+" "+str(value[1])+" " \
                   +str(value[2])+" "+str(value[3])+"\n")
  f_mdff.write("    </bond>\n")
#
  if len(new_bond_list) != nbond[i]:
    print("ERROR: new_bond_list do not match nbond.")
    sys.exit()
########
# angle
########
  start_lno=line_id_angle_start[i]+1
  end_lno  =line_id_angle_end[i]
  angle_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(int(parts[0]),int(parts[1]),int(parts[2]),float(parts[3]),float(parts[4]))
    angle_list.append(temp)
  new_angle_list=[]
  for m in range(0,len(angle_list)):
    temp=angle_list[m]
    oldi=temp[0]
    oldj=temp[1]
    oldk=temp[2]
    kangle=temp[3]
    angle0=temp[4]
    newi=id_old2new[oldi]
    newj=id_old2new[oldj]
    newk=id_old2new[oldk]
    new_angle_list.append([newi,newj,newk,kangle,angle0])
  f_mdff.write("    nangle="+str(nangle[i])+"\n")
  f_mdff.write("    <angle>\n")
  for m in range(0,len(new_angle_list)):
    value=new_angle_list[m]
    f_mdff.write("    "+str(value[0])+" "+str(value[1])+" "+str(value[2])+" "+str(value[3])+" "+str(value[4])+"\n")
  f_mdff.write("    </angle>\n")
#
  if len(new_angle_list) != nangle[i]:
    print("ERROR: new_angle_list do not match nangle.")
    sys.exit()
###########
# dihedral
###########
  start_lno=line_id_dihed_start[i]+1
  end_lno  =line_id_dihed_end[i]
  dihedral_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3]),float(parts[4]),int(parts[5]),int(parts[6]))
    dihedral_list.append(temp)
  new_dihedral_list=[]
  for m in range(0,len(dihedral_list)):
    temp=dihedral_list[m]
    oldi=temp[0]
    oldj=temp[1]
    oldk=temp[2]
    oldl=temp[3]
    kdihedral=temp[4]
    ntheta=temp[5]
    phi=temp[6]
    newi=id_old2new[oldi]
    newj=id_old2new[oldj]
    newk=id_old2new[oldk]
    newl=id_old2new[oldl]
    new_dihedral_list.append([newi,newj,newk,newl,kdihedral,ntheta,phi])
  f_mdff.write("    ndihedral="+str(ndihedral[i])+"\n")
  f_mdff.write("    <dihedral>\n")
  for m in range(0,len(new_dihedral_list)):
    value=new_dihedral_list[m]
    f_mdff.write("    "+str(value[0])+" "+str(value[1])+" " \
         +str(value[2])+" "+str(value[3])+" "+str(value[4])+" " \
         +str(value[5])+" "+str(value[6])+"\n")
  f_mdff.write("    </dihedral>\n")
#
  if len(new_dihedral_list) != ndihedral[i]:
    print("ERROR: new_dihedral_list do not match ndihedral.")
#   print(len(new_dihedral_list),ndihedral[i])
    sys.exit()
###########
# itorsion
###########
  start_lno=line_id_improper_start[i]+1
  end_lno  =line_id_improper_end[i]
  improper_list=[]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split()
    temp  =(int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3]),float(parts[4]),int(parts[5]),int(parts[6]))
    improper_list.append(temp)
  new_improper_list=[]
  for m in range(0,len(improper_list)):
    temp=improper_list[m]
    oldi=temp[0]
    oldj=temp[1]
    oldk=temp[2]
    oldl=temp[3]
    kimp=temp[4]
    dimp=temp[5]
    nimp=temp[6]
    newi=id_old2new[oldi]
    newj=id_old2new[oldj]
    newk=id_old2new[oldk]
    newl=id_old2new[oldl]
    new_improper_list.append([newi,newj,newk,newl,kimp,dimp,nimp])
  f_mdff.write("    nitorsion="+str(nitorsion[i])+"\n")
  f_mdff.write("    <itorsion>\n")
  for m in range(0,len(new_improper_list)):
    value=new_improper_list[m]
    f_mdff.write("    "+str(value[0])+" "+str(value[1])+" " \
         +str(value[2])+" "+str(value[3])+" "+str(value[4])+" " \
         +str(value[5])+" "+str(value[6])+"\n")
  f_mdff.write("    </itorsion>\n")
#
  if len(new_improper_list) != nitorsion[i]:
    print("ERROR: new_improper_list do not match nitorsion.")
#   print(len(new_improper_list),nitorsion[i])
    sys.exit()
##########
# segment
##########
  f_mdff.write("    <segments>\n")
  f_mdff.write("      nsegment="+str(nsegment[i])+"\n")
  idatom=0
  for m in range(0,nsegment[i]):
    natominsegment=natoms_in_segment[m]
    f_mdff.write("      <segment>\n")
    f_mdff.write("        ID= "+str(m)+"\n")
    f_mdff.write("        natom= "+str(natominsegment)+"\n")
    f_mdff.write("        <atom>\n")
    for j in range(0,natominsegment):
      f_mdff.write("          "+str(idatom)+"\n")
      idatom=idatom+1
    f_mdff.write("        </atom>\n")
    f_mdff.write("      </segment>\n")
  f_mdff.write("    </segments>\n")
#######
  f_mdff.write("  </species>\n")
f_mdff.write("</topology and parameters>\n")

print("convert of .mdff ends.")

sys.exit()

### input segments ###
natoms_in_segment=[]
line_id_atomtag_segment=[]
for i in range(0,nspecies):
  start_lno=line_id_segment_start[i]
  end_lno=line_id_segment_end[i]
  for m in range(start_lno,end_lno):
    parts = read_data[m].split("=")
    if "natom" in parts[0]:
      value=int(parts[1])
      natoms_in_segment.append(value)
      temp=[i+2,i+2+value]
      line_id_atomtag_segment.append(temp)
sys.exit()

atomid_org=[]
for i in range(0,nsegment):
  temp=line_id_atomtag_segment[i]
  start_lno=temp[0]
  end_lno=temp[1]
  for j in range(start_lno,end_lno):
    value=int(read_data[j])
    atomid_org.append(value)

###debug
#print(natoms_in_segment)
#print(line_id_atomtag_segment)
#print(atomid_org)

if len(atomid_org) != natom:
  print("Error len(atomid_org) != natom")
  sys.exit()

id_old2new=[]
for i in range(0,natom):
  id_old2new.append(0)

for i in range(0,natom):
  newid=atomid_org[i]
  id_old2new[newid]=i


### input/output masse information ###
start_lno=line_id_mass_start+1
end_lno  =line_id_mass_end
mass_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
# temp  =(float(parts[0]))
  temp  =(float(parts[0]),parts[2])
  mass_list.append(temp)
new_mass_list=[]
for i in range(0,len(mass_list)):
  new_mass_list.append(0)
for oldid in range(0,len(mass_list)):
  newid=id_old2new[oldid]
  new_mass_list[newid]=mass_list[oldid]
#print(new_mass_list)  
f_mdff.write("    <mass>\n")
for i in range(0,len(new_mass_list)):
  value=new_mass_list[i]
# f_mdff.write("    "+str(value)+"\n")
  f_mdff.write("    "+str(value[0])+" # "+str(value[1])+"\n")
f_mdff.write("    </mass>\n")

### input/output epsilon information ###
start_lno=line_id_eps_start+1
end_lno  =line_id_eps_end
eps_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(float(parts[0]))
  eps_list.append(temp)
new_eps_list=[]
for i in range(0,len(eps_list)):
  new_eps_list.append(0)
for oldid in range(0,len(eps_list)):
  newid=id_old2new[oldid]
  new_eps_list[newid]=eps_list[oldid]
#print(new_eps_list)  
f_mdff.write("    <epsilon>\n")
for i in range(0,len(new_eps_list)):
  value=new_eps_list[i]
  f_mdff.write("    "+str(value)+"\n")
f_mdff.write("    </epsilon>\n")

### input/output r information ###
start_lno=line_id_r_start+1
end_lno  =line_id_r_end
r_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(float(parts[0]))
  r_list.append(temp)
new_r_list=[]
for i in range(0,len(r_list)):
  new_r_list.append(0)
for oldid in range(0,len(r_list)):
  newid=id_old2new[oldid]
  new_r_list[newid]=r_list[oldid]
#print(new_r_list)  
f_mdff.write("    <r>\n")
for i in range(0,len(new_r_list)):
  value=new_r_list[i]
  f_mdff.write("    "+str(value)+"\n")
f_mdff.write("    </r>\n")

### input/output charge information ###
start_lno=line_id_charge_start+1
end_lno  =line_id_charge_end
charge_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(float(parts[0]))
  charge_list.append(temp)
new_charge_list=[]
for i in range(0,len(charge_list)):
  new_charge_list.append(0)
for oldid in range(0,len(charge_list)):
  newid=id_old2new[oldid]
  new_charge_list[newid]=charge_list[oldid]
#print(new_charge_list)  
f_mdff.write("    <charge>\n")
for i in range(0,len(new_charge_list)):
  value=new_charge_list[i]
  f_mdff.write("    "+str(value)+"\n")
f_mdff.write("    </charge>\n")

#sys.exit()

### input/output ljvoid information ###
start_lno=line_id_ljvoid_start+1
end_lno  =line_id_ljvoid_end
ljvoid_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(int(parts[0]),int(parts[1]))
  ljvoid_list.append(temp)
#
new_ljvoid_list=[]
for i in range(0,len(ljvoid_list)):
  temp=ljvoid_list[i]
  oldi=temp[0]
  oldj=temp[1]
  newi=id_old2new[oldi]
  newj=id_old2new[oldj]
  new_ljvoid_list.append([newi,newj])
#
f_mdff.write("    <lj void pair>\n")
for i in range(0,len(new_ljvoid_list)):
  value=new_ljvoid_list[i]
  f_mdff.write("    "+str(value[0])+" "+str(value[1])+"\n")
f_mdff.write("    </lj void pair>\n")
f_mdff.write("    nvoidpair_lj="+str(nvoidpair_lj)+"\n")

if len(new_ljvoid_list) != nvoidpair_lj:
  print("ERROR: new_ljvoid_list do not match nvoidpair_lj.")
  sys.exit()

# lj special pair (always 0 in gaff)
f_mdff.write("    <lj special pair>\n")
f_mdff.write("    </lj special pair>\n")
f_mdff.write("    nspecialpair_lj =0\n")

### input/output clvoid information ###
start_lno=line_id_clvoid_start+1
end_lno  =line_id_clvoid_end
clvoid_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(int(parts[0]),int(parts[1]))
  clvoid_list.append(temp)
#
new_clvoid_list=[]
for i in range(0,len(clvoid_list)):
  temp=clvoid_list[i]
  oldi=temp[0]
  oldj=temp[1]
  newi=id_old2new[oldi]
  newj=id_old2new[oldj]
  new_clvoid_list.append([newi,newj])
#
f_mdff.write("    <coulomb void pair>\n")
for i in range(0,len(new_clvoid_list)):
  value=new_clvoid_list[i]
  f_mdff.write("    "+str(value[0])+" "+str(value[1])+"\n")
f_mdff.write("    </coulomb void pair>\n")
f_mdff.write("    nvoidpair_coulomb="+str(nvoidpair_coulomb)+"\n")

if len(new_clvoid_list) != nvoidpair_lj:
  print("ERROR: new_clvoid_list do not match nvoidpair_coulomb.")
  sys.exit()

### input/output shake information ###
start_lno=line_id_shake_start+1
end_lno  =line_id_shake_end
nshakepair=end_lno-start_lno
shake_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(int(parts[0]),int(parts[1]),float(parts[2]),parts[4],parts[5])
  shake_list.append(temp)
#
new_shake_list=[]
for i in range(0,len(shake_list)):
  temp=shake_list[i]
  oldi=temp[0]
  oldj=temp[1]
  length=temp[2]
  namei=temp[3]
  namej=temp[4]
  newi=id_old2new[oldi]
  newj=id_old2new[oldj]
  new_shake_list.append([newi,newj,length,namei,namej])
#
f_mdff.write("    <shake pair>\n")
for i in range(0,len(new_shake_list)):
  value=new_shake_list[i]
  f_mdff.write("    "+str(value[0])+" "+str(value[1])+" "+str(value[2])+" # "+str(value[3])+" "+str(value[4])+"\n")
f_mdff.write("    </shake pair>\n")

if len(new_shake_list) != nshakepair:
  print("ERROR: new_shake_list do not match nshakepair.")
  print(len(new_shake_list), nshakepair)
  sys.exit()

### input/output bond information ###
start_lno=line_id_bond_start+1
end_lno  =line_id_bond_end
bond_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  if len(parts) == 7:
    temp  =(int(parts[0]),int(parts[1]),float(parts[2]), \
            float(parts[3]),parts[5],parts[6])
  else:
    temp  =(int(parts[0]),int(parts[1]),float(parts[2]), \
            float(parts[3]))
  bond_list.append(temp)
#
new_bond_list=[]
for i in range(0,len(bond_list)):
  temp=bond_list[i]
  oldi=temp[0]
  oldj=temp[1]
  kbond=temp[2]
  length=temp[3]
  if len(temp)==7:
    namei=temp[4]
    namej=temp[5]
  newi=id_old2new[oldi]
  newj=id_old2new[oldj]
  if len(temp)==7:
    new_bond_list.append([newi,newj,kbond,length,namei,namej])
  else:
    new_bond_list.append([newi,newj,kbond,length])
#
f_mdff.write("    nbond="+str(nbond)+"\n")
f_mdff.write("    <bond>\n")
for i in range(0,len(new_bond_list)):
  value=new_bond_list[i]
  if len(value)==7:
    f_mdff.write("    "+str(value[0])+" "+str(value[1])+" " \
                 +str(value[2])+str(value[3])+" # " \
                 +str(value[4])+" "+str(value[5])+"\n")
  else:
    f_mdff.write("    "+str(value[0])+" "+str(value[1])+" " \
                 +str(value[2])+" "+str(value[3])+"\n")
f_mdff.write("    </bond>\n")

if len(new_bond_list) != nbond:
  print("ERROR: new_bond_list do not match nbond.")
  sys.exit()

### input/output angle information ###
start_lno=line_id_angle_start+1
end_lno  =line_id_angle_end
angle_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(int(parts[0]),int(parts[1]),int(parts[2]),float(parts[3]),float(parts[4]))
  angle_list.append(temp)
#
new_angle_list=[]
for i in range(0,len(angle_list)):
  temp=angle_list[i]
  oldi=temp[0]
  oldj=temp[1]
  oldk=temp[2]
  kangle=temp[3]
  angle0=temp[4]
  newi=id_old2new[oldi]
  newj=id_old2new[oldj]
  newk=id_old2new[oldk]
  new_angle_list.append([newi,newj,newk,kangle,angle0])
#
f_mdff.write("    nangle="+str(nangle)+"\n")
f_mdff.write("    <angle>\n")
for i in range(0,len(new_angle_list)):
  value=new_angle_list[i]
  f_mdff.write("    "+str(value[0])+" "+str(value[1])+" "+str(value[2])+" "+str(value[3])+" "+str(value[4])+"\n")
f_mdff.write("    </angle>\n")

if len(new_angle_list) != nangle:
  print("ERROR: new_angle_list do not match nangle.")
  sys.exit()

### input/output dihedral information ###
start_lno=line_id_dihedral_start+1
end_lno  =line_id_dihedral_end
dihedral_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3]),float(parts[4]),int(parts[5]),int(parts[6]))
  dihedral_list.append(temp)
#
new_dihedral_list=[]
for i in range(0,len(dihedral_list)):
  temp=dihedral_list[i]
  oldi=temp[0]
  oldj=temp[1]
  oldk=temp[2]
  oldl=temp[3]
  kdihedral=temp[4]
  ntheta=temp[5]
  phi=temp[6]
  newi=id_old2new[oldi]
  newj=id_old2new[oldj]
  newk=id_old2new[oldk]
  newl=id_old2new[oldl]
  new_dihedral_list.append([newi,newj,newk,newl,kdihedral,ntheta,phi])
#
f_mdff.write("    ndihedral="+str(ndihedral)+"\n")
f_mdff.write("    <dihedral>\n")
for i in range(0,len(new_dihedral_list)):
  value=new_dihedral_list[i]
  f_mdff.write("    "+str(value[0])+" "+str(value[1])+" " \
       +str(value[2])+" "+str(value[3])+" "+str(value[4])+" " \
       +str(value[5])+" "+str(value[6])+"\n")
f_mdff.write("    </dihedral>\n")

if len(new_dihedral_list) != ndihedral:
  print("ERROR: new_dihedral_list do not match ndihedral.")
  print(len(new_dihedral_list),ndihedral)
  sys.exit()

### input/output itorsion information ###
start_lno=line_id_improper_start+1
end_lno  =line_id_improper_end
improper_list=[]
for i in range(start_lno,end_lno):
  parts = read_data[i].split()
  temp  =(int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3]),float(parts[4]),int(parts[5]),int(parts[6]))
  improper_list.append(temp)
#
new_improper_list=[]
for i in range(0,len(improper_list)):
  temp=improper_list[i]
  oldi=temp[0]
  oldj=temp[1]
  oldk=temp[2]
  oldl=temp[3]
  kimp=temp[4]
  dimp=temp[5]
  nimp=temp[6]
  newi=id_old2new[oldi]
  newj=id_old2new[oldj]
  newk=id_old2new[oldk]
  newl=id_old2new[oldl]
  new_improper_list.append([newi,newj,newk,newl,kimp,dimp,nimp])
#
f_mdff.write("    nitorsion="+str(nitorsion)+"\n")
f_mdff.write("    <itorsion>\n")
for i in range(0,len(new_improper_list)):
  value=new_improper_list[i]
  f_mdff.write("    "+str(value[0])+" "+str(value[1])+" " \
       +str(value[2])+" "+str(value[3])+" "+str(value[4])+" " \
       +str(value[5])+" "+str(value[6])+"\n")
f_mdff.write("    </itorsion>\n")

if len(new_improper_list) != nitorsion:
  print("ERROR: new_improper_list do not match nitorsion.")
  print(len(new_improper_list),nitorsion)
  sys.exit()

### input/output segment information ###
f_mdff.write("    <segments>\n")
f_mdff.write("      nsegment="+str(nsegment)+"\n")
idatom=0
for i in range(0,nsegment):
  natominsegment=natoms_in_segment[i]
  f_mdff.write("      <segment>\n")
  f_mdff.write("        ID= "+str(i)+"\n")
  f_mdff.write("        natom= "+str(natominsegment)+"\n")
  f_mdff.write("        <atom>\n")
  for j in range(0,natominsegment):
    f_mdff.write("          "+str(idatom)+"\n")
    idatom=idatom+1
  f_mdff.write("        </atom>\n")
  f_mdff.write("      </segment>\n")
f_mdff.write("    </segments>\n")
f_mdff.write("  </species>\n")
f_mdff.write("</topology and parameters>\n")

print("convert of .mdff ends.")

#
### open file for output
#
f_mdxyz=open(sessionname+'.mdxyz.rearranged','w')

#
# read old .mdxyz
#
filename=sessionname+".mdxyz"
with open(filename, mode="r", encoding="utf-8") as fxyz:
  read_data_xyz = list(fxyz)

length_data_xyz=len(read_data_xyz)

##print(length_data_xyz)

line_id_position_start=0
line_id_position_end=0
natomall=0
for i in range(0,length_data_xyz):
  parts = read_data_xyz[i].split()
# print(parts)
  if "natom" in parts[0]:
    temp=parts[0].split("=")
    natomall=int(temp[1])
  if "<positions>" in parts:
    line_id_position_start=i
  if "</positions>" in parts:
    line_id_position_end=i
    break

#print(line_id_position_start,line_id_position_end)
#print(natomall)
#sys.exit()

position_new=[]
for i in range(0,natomall):
  position_new.append(0)

start_lno=line_id_position_start+1
end_lno=line_id_position_end
idold=0
idmol=0
for i in range(start_lno,end_lno):
  parts = read_data_xyz[i].split()
# print(parts)
  x=parts[0]
  y=parts[1]
  z=parts[2]
  idoldmol=idold-idmol*natom
# print(idoldmol,idmol)
  newidmol=id_old2new[idoldmol]
  newid=newidmol+idmol*natom
# print(idoldmol,idmol,idold,newid)
  position_new[newid]=[x,y,z]
  idold=idold+1
  if idoldmol==natom-1:
    idmol=idmol+1

#sys.exit()

###
### output new .mdxyz ###
###
for i in range(0,start_lno):
  line=read_data_xyz[i].split("\n")
  f_mdxyz.write(line[0]+"\n")

for i in range(0,natomall):
  xyz=position_new[i]
  x=xyz[0]
  y=xyz[1]
  z=xyz[2]
  f_mdxyz.write("  "+x+" "+y+" "+z+"\n")

for i in range(end_lno,length_data_xyz):
  line=read_data_xyz[i].split("\n")
  f_mdxyz.write(line[0]+"\n")
  
print("convert of .mdxyz ends.")



#
### open file for sys_info
#
f_sysinfo=open('sys_info','w')
f_sysinfo.write("species name \n")
f_sysinfo.write(str(nspecies)+"\n")
for i in range(0,nspecies):
  f_sysinfo.write(str(sysatom_list[i])+"  "+str(sysmol_list[i])+"\n")
f_sysinfo.close()
print("./sys_info created.")

#
### open file for massinfo.mdff
#
f_massinfo=open('massinfo.mdff','w')

f_massinfo.write("<mass>\n")
for i in range(0,len(new_mass_list)):
  value=new_mass_list[i]
  f_massinfo.write("    "+str(value[0])+" # "+str(value[1])+"\n")
f_massinfo.write("</mass>\n")

f_massinfo.close()
print("./massinfo created.")

#
### open file for seginfo.mdff
#
f_seginfo=open('seginfo.mdff','w')

f_seginfo.write("    <segments>\n")
f_seginfo.write("      "+str(nsegment)+"\n")
idatom=0
for i in range(0,nsegment):
  natominsegment=natoms_in_segment[i]
  f_seginfo.write("      <segment>\n")
  f_seginfo.write("        "+str(i)+"\n")
  f_seginfo.write("        "+str(natominsegment)+"\n")
  f_seginfo.write("        <atom>\n")
  for j in range(0,natominsegment):
    f_seginfo.write("          "+str(idatom)+"\n")
    idatom=idatom+1
  f_seginfo.write("        </atom>\n")
  f_seginfo.write("      </segment>\n")
f_seginfo.write("    </segments>\n")

f_seginfo.close()
print("./seginfo created.")
