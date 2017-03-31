#!/usr/bin/python

import os
import sys

def read_atomnames(line_idx,lines,num_items,pfile):
    i, tmp_data = 0, []
    while i < num_items:
        idx = i % 20
        if idx == 0:
            pfile.write(lines[line_idx])
            line_idx += 1
        tmp_data.append(lines[line_idx][idx*4:idx*4+4])
        i += 1
    if num_items == 0: line_idx += 1
    return tmp_data, line_idx

#-----------------------------

def read_atomtypes(line_idx,lines,num_items,pfile):
        i, tmp_data = 0, []
        while i < num_items:
            idx = i % 10
            if idx == 0:
                pfile.write(lines[line_idx])
                line_idx += 1
            tmp_data.append(lines[line_idx][idx*8:idx*8+8])
            i += 1
        if num_items == 0: line_idx += 1
        return tmp_data, line_idx

#-----------------------------

def read_pointers(line_idx,lines,num_items,pfile):
        i, tmp_data = 0, []
        while i < num_items:
            idx = i % 10
            if idx == 0:
                pfile.write(lines[line_idx])
                line_idx += 1
            tmp_data.append(lines[line_idx][idx*8:idx*8+8])
            i += 1
        if num_items == 0: line_idx += 1
        return tmp_data, line_idx

#-----------------------------

def read_bcoef(line_idx,lines,num_items):
        i, tmp_data = 0, []
        while i < num_items:
            idx = i % 5
            if idx == 0:
                line_idx += 1
            tmp_data.append(lines[line_idx][idx*16:idx*16+16])
            i += 1
        if num_items == 0: line_idx += 1
        return tmp_data, line_idx

#=============================

line_idx = 0
temp_data = []
bindex = []
npiindex = []
bzero = '  0.00000000E+00'

print("\nThis is a simple script designed to zero out the Lennard-Jones")
print("B values for dummy atoms in our multisite ion model.\n")

pname = input('Please input topology file: ')

pfile_lines = open(pname,"r").readlines()

new_pname = pname + '.mod'

print('Modified parameter file will be %s' % new_pname)

new_pfile = open(new_pname,"w")

while line_idx < len(pfile_lines):

    line = pfile_lines[line_idx]

    if line[0:5] == '%FLAG':
        this_flag = line[6:].strip()
        
        if this_flag == 'POINTERS':
            new_pfile.write(pfile_lines[line_idx])
            line_idx += 1
            temp_data, line_idx = read_pointers(line_idx,pfile_lines,30,new_pfile)
            natoms = int(temp_data[0])
            ntypes = int(temp_data[1])
            print('Topology file has %s atoms and %s atom types' % (natoms,ntypes))
                
        elif this_flag == 'ATOM_NAME':
            new_pfile.write(pfile_lines[line_idx])
            line_idx += 1
            dummy_atom = 0
            temp_data, line_idx = read_atomnames(line_idx,pfile_lines,natoms,new_pfile)
            while temp_data[dummy_atom] != 'D1  ': dummy_atom +=1
            print('D1 atom at position %d in atom list'  % dummy_atom)

        elif this_flag == 'ATOM_TYPE_INDEX':
            new_pfile.write(pfile_lines[line_idx])
            line_idx += 1
            temp_data, line_idx = read_atomtypes(line_idx,pfile_lines,natoms,new_pfile)
            dummy_type = int(temp_data[dummy_atom])
            print('D1 atom has type %d' % dummy_type)
            for i in range(1,ntypes+1): npiindex.append(ntypes*(i-1)+dummy_type)

        elif this_flag == 'NONBONDED_PARM_INDEX':
            new_pfile.write(pfile_lines[line_idx])
            line_idx += 1
            temp_data, line_idx = read_atomtypes(line_idx,pfile_lines,ntypes*ntypes,new_pfile)
            for i in npiindex: 
                bindex.append(int(temp_data[i-1]))
            print('Making substitution in LENNARD_JONES_BCOEF at positions: %s' % bindex)

        elif this_flag == 'LENNARD_JONES_BCOEF':
            new_pfile.write(pfile_lines[line_idx])
            line_idx += 1
            new_pfile.write(pfile_lines[line_idx])
            temp_data, line_idx = read_bcoef(line_idx, pfile_lines,int(ntypes*(ntypes+1)/2))
            for i in range(1,len(temp_data)+1):
                idx = i % 5
                if i in bindex:
                    new_pfile.write(bzero)
                else:
                    new_pfile.write(temp_data[i-1])
                if idx==0:
                    new_pfile.write('\n')
            new_pfile.write('\n')
            line_idx +=1
            
        else:
            new_pfile.write(line)
            line_idx += 1

    else:
        new_pfile.write(line)
        line_idx += 1

new_pfile.close()
sys.exit(1)
