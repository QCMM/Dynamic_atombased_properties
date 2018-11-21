import os, os.path
import openbabel as ob
import psi4
import numpy as np
from openmoltools import amber
import mdtraj
import shutil
from openmoltools.utils import getoutput
from openmoltools.amber import check_for_errors

"""
dynamic_qm_properties.py
A tool to obtain dynamic qm descriptors

Handles the primary functions
"""


def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())

def input2xyz(input_file):
    """
    Converts input_file to xyz format using the OpenBabel module. 
    More info about this module in http://openbabel.org/docs/dev/UseTheLibrary/PythonDoc.html
    """
    inputName=input_file.split('.')[0]
    inputFormat=input_file.split('.')[1]
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats(inputFormat, "xyz")
    mol=ob.OBMol()
    obConversion.ReadFile(mol,input_file)
    obConversion.WriteFile(mol,inputName+'.xyz')
    return(inputName+'.xyz')#deberia retornar un objeto mas complejo, pero por ahora solamente el nombre


def input2pdb(input_file):
    """
    Converts input_file to pdb format using the OpenBabel module. 
    More info about this module in http://openbabel.org/docs/dev/UseTheLibrary/PythonDoc.html
    """
    inputName=input_file.split('.')[0]
    inputFormat=input_file.split('.')[1]
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats(inputFormat, "pdb")
    mol=ob.OBMol()
    obConversion.ReadFile(mol,input_file)
    obConversion.WriteFile(mol,inputName+'.pdb')
    return(inputName+'.pdb')#deberia retornar un objeto mas complejo, pero por ahora solamente el nombre

def input2mol2(input_file):
    """
    Converts input_file to mol2 format using the OpenBabel module. 
    More info about this module in http://openbabel.org/docs/dev/UseTheLibrary/PythonDoc.html
    """
    inputName=input_file.split('.')[0]
    inputFormat=input_file.split('.')[1]
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats(inputFormat, "mol2")
    mol=ob.OBMol()
    obConversion.ReadFile(mol,input_file)
    obConversion.WriteFile(mol,inputName+'.mol2')
    return(inputName+'.mol2')#deberia retornar un objeto mas complejo, pero por ahora solamente el nombre

def xyz2gcrt(xyz_file,charge,mult):#Necesita mult y carga, datos que no estan en XYZ, pero si en MOL2 creo.
    """
    Converts XYZ file to gcrt format. The latter can be read by Antechamber.
    The different variables "line" are meant to satisfy the requirements of gcrt
    """
    system_name=xyz_file.split('.')[0]
    output_format='.gcrt'
    output_file=system_name+output_format
    line1='--Link1--\n'
    line2='%chk=molecule\n'
    line3='#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt\n'#default method
    line4=' \n'
    line5='remark goes here\n'
    line6=line4
    line7=str(charge)+'   '+str(mult)+'\n'
    line8=line4
    with open(output_file,'w') as temp1:
        f=open(xyz_file,'r')
        f_lines=f.readlines()
        temp1.writelines([line1,line2,line3,line4, line5, line6, line7])
        temp1.writelines(f_lines[2:])
        temp1.writelines([line8])
        f.close()
    return output_file

def psi4calc(mol):
    """
    Creates a molecule object for Psi4 from  XYZ file
    """
    geom=psi4.qcdb.Molecule.init_with_xyz(mol)
    mol=psi4.geometry(geom.create_psi4_string_from_molecule())
    return mol

def gcrt2prmtop(gcrt_file):
    system_name=gcrt_file.split('.')[0]
    amber.run_antechamber(system_name, gcrt_file, charge_method=None, input_format='gcrt')
    amber.run_tleap(system_name, system_name+'.gaff.mol2',system_name+'.frcmod',system_name+'.prmtop',system_name+'.crd')
    return None


def run_tleap_solv(molecule_name, gaff_mol2_filename, frcmod_filename, prmtop_filename=None, inpcrd_filename=None, log_debug_output=False, leaprc='leaprc.gaff'):
    """Run AmberTools tleap to create simulation files for AMBER
    Parameters
    ----------
    molecule_name : str
        The name of the molecule
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    prmtop_filename : str, optional, default=None
        Amber prmtop file produced by tleap, defaults to molecule_name
    inpcrd_filename : str, optional, default=None
        Amber inpcrd file produced by tleap, defaults to molecule_name
    log_debug_output : bool, optional, default=False
        If true, will send output of tleap to logger.
    leaprc : str, optional, default = 'leaprc.gaff'
        Optionally, specify alternate leaprc to use, such as `leaprc.gaff2`
    Returns
    -------
    prmtop_filename : str
        Amber prmtop file produced by tleap
    inpcrd_filename : str
        Amber inpcrd file produced by tleap
    """
    if prmtop_filename is None:
        prmtop_filename = "%s.prmtop" % molecule_name
    if inpcrd_filename is None:
        inpcrd_filename = "%s.inpcrd" % molecule_name

    #Get absolute paths for input/output
    gaff_mol2_filename = os.path.abspath( gaff_mol2_filename )
    frcmod_filename = os.path.abspath( frcmod_filename )
    prmtop_filename = os.path.abspath( prmtop_filename )
    inpcrd_filename = os.path.abspath( inpcrd_filename )

    #Work in a temporary directory, on hard coded filenames, to avoid any issues AMBER may have with spaces and other special characters in filenames
    with mdtraj.utils.enter_temp_directory():
        shutil.copy( gaff_mol2_filename, 'file.mol2' )
        shutil.copy( frcmod_filename, 'file.frcmod' )

        tleap_input = """
    source oldff/leaprc.ff99SB
    source %s
    source leaprc.water.tip3p
    LIG = loadmol2 file.mol2
    check LIG
    loadamberparams file.frcmod
    solvatebox LIG TIP3PBOX 8 
    addions LIG Cl- 0 
    saveamberparm LIG out.prmtop out.inpcrd
    quit
""" % leaprc

        file_handle = open('tleap_commands', 'w')
        file_handle.writelines(tleap_input)
        file_handle.close()

        cmd = "tleap -f %s " % file_handle.name
        if log_debug_output: logger.debug(cmd)

        output = getoutput(cmd)
        if log_debug_output: logger.debug(output)

        check_for_errors( output, other_errors = ['Improper number of arguments'] )

        #Copy back target files
        shutil.copy( 'out.prmtop', prmtop_filename )
        shutil.copy( 'out.inpcrd', inpcrd_filename )

    return prmtop_filename, inpcrd_filename

