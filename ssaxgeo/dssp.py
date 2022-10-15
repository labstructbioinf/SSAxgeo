import numpy as np
import pandas as pd
import subprocess
import math
import re

### FUNCTIONS ##################################################################
def writeDSSPcsv(dssp_array, filepath, prefix):
    '''Write DSSP csv output'''

    with open(filepath+prefix+"_dssp.csv", 'w') as dssp_csv:
        dssp_csv.write("#index,dssp_ss\n")
        for index, row in dssp_array.iterrows():
            #print(row["SS"])
            aa_index = str(index)
            ss = row["ss"]
            dssp_csv.write(aa_index+','+ss+"\n")
    dssp_csv.close()

def os_cmd(cmd):
    subprocess.call([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def runDSSP(coord_flpath, outprfx):
    ''' Run DSSP on specific of entry on a PDB List array.'''
    # get dir
    flnme = coord_flpath.split('/')[-1]
    dir = coord_flpath[0:len(coord_flpath) - len(flnme)]

    try:
        cmd = "dssp -i {} -o {}".format(coord_flpath, dir+outprfx+'.dssp')
        os_cmd(cmd)
        #print(pdbcode+pdbchain)
        #writeDSSPcsv(dssp_array, filepath=flpath,prefix=pdbcode+pdbchain)
        return 0

    except:
        print("!!WARNING : FAIL for {}".format(flpath))
        # -- TODO --
        # get list of failed files to check in the future
        # ----------
        return 1
#-------------------------------------------------------------------------------

class DSSPData:
  '''
  This class is designed to deal with DSSP output.
  This code is a modified version of the one available at:
  https://gist.github.com/jlhg/5181883
  '''

  def __init__(self, file):
      self.num    = []
      self.resnum = []
      self.moltyp = []
      self.aa     = []
      self.struct = []
      self.bp1    = []
      self.bp2    = []
      self.acc    = []
      self.h_nho1 = []
      self.h_ohn1 = []
      self.h_nho2 = []
      self.h_ohn2 = []
      self.tco    = []
      self.kappa  = []
      self.alpha  = []
      self.phi    = []
      self.psi    = []
      self.xca    = []
      self.yca    = []
      self.zca    = []

      input_handle = open(file, 'r')

      line_num = 0
      start=False
      for line in input_handle:

        if( re.search('#', line) ):
          start=True
          continue

        if( start ):
          self.num.append(    line[0:5].strip() )
          self.resnum.append( line[5:10].strip() )
          self.moltyp.append( line[10:12].strip() )
          self.aa.append(     line[12:14].strip() )
          self.struct.append( line[14:25] )
          self.bp1.append(    line[25:29].strip() )
          self.bp2.append(    line[29:34].strip() )
          self.acc.append(    line[34:38].strip() )
          self.h_nho1.append( line[38:50].strip() )
          self.h_ohn1.append( line[50:61].strip() )
          self.h_nho2.append( line[61:72].strip() )
          self.h_ohn2.append( line[72:83].strip() )
          self.tco.append(    line[83:91].strip() )
          self.kappa.append(  line[91:97].strip() )
          self.alpha.append(  line[97:103].strip() )
          self.phi.append(    line[103:109].strip() )
          self.psi.append(    line[109:115].strip() )
          self.xca.append(    line[115:122].strip() )
          self.yca.append(    line[122:129].strip() )
          self.zca.append(    line[129:136].strip() )

  def getDSSP_data(self):
      # 1 - Get only SS assignment, original PDB resnum and AA seq
      self.data = []
      for i, each in enumerate(self.struct):
          ss = each[0:3].replace(' ', '')
          #h_nho1
          h_nho1_list =  self.h_nho1[i].split(',')
          h_nho1_A = int(h_nho1_list[0])
          h_nho1_en = float(h_nho1_list[1])

          #h_nho2
          h_nho2_list =  self.h_nho2[i].split(',')
          h_nho2_A = int(h_nho2_list[0])
          h_nho2_en = float(h_nho2_list[1])

          #h_ohn1
          h_ohn1_list =  self.h_ohn1[i].split(',')
          h_ohn1_A = int(h_ohn1_list[0])
          h_ohn1_en = float(h_ohn1_list[1])

          #h_ohn2
          h_ohn2_list =  self.h_ohn2[i].split(',')
          h_ohn2_A = int(h_ohn2_list[0])
          h_ohn2_en = float(h_ohn2_list[1])

          self.data.append([self.resnum[i], self.aa[i], ss,
                            h_nho1_A, h_nho1_en,
                            h_nho2_A, h_nho2_en,
                            h_ohn1_A, h_ohn1_en,
                            h_ohn2_A, h_ohn2_en,
                            self.phi[i], self.psi[i]])

      Data = pd.DataFrame(self.data)
      Data.columns = ['resnum', 'aa', 'ss',
                      'h_nho1_A','h_nho1_en',
                      'h_nho2_A','h_nho2_en',
                      'h_ohn1_A','h_ohn1_en',
                      'h_ohn2_A','h_ohn2_en',
                      'phi', 'psi']
      #self.data = Data.drop_duplicates('resnum')
      self.data = Data
      return self.data


if __name__ == "__name__":
    #flnm = '/home/amarinho/PosDocBox/PDBxgeo/localPDB/pdbs/zg/4zgnB.dssp'
    #flnm = '/media/amarinho/wHD_Seagate/PosDocBox/PDBxgeo/localPDB/pdbs/27/127lA.dssp'

    # --| DEBUG |--
    # multiconf
    #flnm = '/media/amarinho/wHD_Seagate/PosDocBox/PDBxgeo/localPDB/pdbs/ak/1akpA.dssp'

    # non-continuous chain
    #flnm = '/home/amarinho/PosDocBox/PDBxgeo/localPDB/pdbs/if/1ifgA.dssp'
    # -------------
    DSSPdt = DSSPData(flnm)
    with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
        print(DSSPdt.getDSSP_data())
