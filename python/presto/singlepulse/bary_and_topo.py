#! /usr/bin/env python

"""
   Original code found in presto. Written by Scott M. Ransom.
   Modified by Chitrang Patel to read information from PSRFITs file.
   Modified to return topocentric and corresponding barycentric
   times. 
"""
from __future__ import print_function
from presto.presto.prestoswig import *
import numpy as Num
from presto import psr_utils
from presto import psr_constants
from presto import psrfits
# add filterbank module
from presto import filterbank, sigproc
import warnings


def read_inffile(filename):
   """
   read_inffile(filename):
       Return an infodata 'C' structure containing the data from the
       'inf' file in 'filename'.  'filename' should not include the
       '.inf' suffix.
   """
   id = infodata()
   print("Reading information from", "\""+filename+".inf\"")
   readinf(id, filename)
   return id


def bary_to_topo(infofilenm, rawdatafile=False, ephem="DE200"):
   """
   bary_to_topo(infofilenm, ephem="DE200"):
      Returns the barycentric and topocentric times evert 10 seconds.
      The data for the observation must be found in the info file.
   """
   if infofilenm[-4:]==".inf":
       infofilenm = infofilenm[:-4]
       filetype = 'inf'
   elif infofilenm[-5:]==".fits":
       infofilenm = infofilenm
       filetype = 'PSRFITS'
    # add the filterbank file type so code doesn't crash
   elif infofilenm.endswith(".fil"):
       infofilenm = infofilenm
       filetype = 'FILTERBANK'
   else:
       raise ValueError("file type not recognized. Must be .inf, or .fits")
   if filetype=="inf": 
       obs = read_inffile(infofilenm)
       T = obs.N * obs.dt
       dt = 10.0
       tto = obs.mjd_i + obs.mjd_f
       tts = Num.arange(tto, tto + (T + dt) / psr_constants.SECPERDAY, dt / psr_constants.SECPERDAY)
       nn = len(tts)
       bts = Num.zeros(nn, 'd')
       vel = Num.zeros(nn, 'd')
       ra = psr_utils.coord_to_string(obs.ra_h, obs.ra_m, obs.ra_s)
       dec = psr_utils.coord_to_string(obs.dec_d, obs.dec_m, obs.dec_s)
       if (obs.telescope == 'Parkes'):  tel = 'PK'
       elif (obs.telescope == 'Effelsberg'):  tel = 'EB'
       elif (obs.telescope == 'Arecibo'):  tel = 'AO'
       elif (obs.telescope == 'MMT'):  tel = 'MT'
       elif (obs.telescope == 'GBT'):  tel = 'GB'
       elif (obs.telescope == 'CHIME'): tel = 'CH'
       else:
          print("Telescope not recognized.")
          return 0
   elif filetype=="PSRFITS": 
       if not rawdatafile:
           rawdatafile = psrfits.PsrfitsFile(infofilenm)
       T = rawdatafile.specinfo.T
       dt = 10.0
       tto = rawdatafile.specinfo.start_MJD[0]
       # Actually psr_utils.SECPERDAY throws an error hence I have hardcoded that in the
       # filterbank case, maybe it should be updated here as well???
       tts = Num.arange(tto, tto + (T + dt) / psr_constants.SECPERDAY, dt / psr_constants.SECPERDAY)
       nn = len(tts)
       bts = Num.zeros(nn, 'd')
       vel = Num.zeros(nn, 'd')
       ra = rawdatafile.specinfo.ra_str
       dec = rawdatafile.specinfo.dec_str
       if (rawdatafile.specinfo.telescope == 'Parkes'):  tel = 'PK'
       elif (rawdatafile.specinfo.telescope == 'Effelsberg'):  tel = 'EB'
       elif (rawdatafile.specinfo.telescope == 'Arecibo'):  tel = 'AO'
       elif (rawdatafile.specinfo.telescope == 'MMT'):  tel = 'MT'
       elif (rawdatafile.specinfo.telescope == 'GBT'):  tel = 'GB'
       elif (rawdatafile.specinfo.telescope == 'CHIME'): tel = 'CH'
       else:
          print("Telescope not recognized.")
          return 0
      
   elif filetype=="FILTERBANK":
       # Get parameters needed for barycentric correction in case of 
       # filterbank file. Here telescope name is hardcoded to GBT because the 
       # filterbank module doesn't seem to read the telescope name from the 
       # file header. 
       if not rawdatafile:
           rawdatafile = filterbank.FilterbankFile(infofilenm)
       T = rawdatafile.nspec * rawdatafile.dt
       dt = 10.0
       tto = rawdatafile.tstart
       tts = Num.arange(tto, tto + (T + dt) / psr_constants.SECPERDAY, dt / psr_constants.SECPERDAY)
       nn = len(tts)
       bts = Num.zeros(nn, 'd')
       vel = Num.zeros(nn, 'd')
       ra_j = rawdatafile.src_raj
       dec_j = rawdatafile.src_dej
       ra = psr_utils.coord_to_string(ra_j // 10000, ra_j % 10000 // 100, ra_j % 10000 % 100)
       if dec_j > 0:
           dec = psr_utils.coord_to_string(dec_j // 10000, dec_j % 10000 // 100, dec_j % 10000 % 100)
       elif dec_j < 0:
           dec_j = -1*dec_j  # take the magnitude to convet properly
           dec = psr_utils.coord_to_string(-1*( dec_j // 10000), dec_j % 10000 // 100, dec_j % 10000 % 100)
       
       # Get telescope name from the telescope ID in the header        
       telID = rawdatafile.telescope_id
       telescope = sigproc.ids_to_telescope[telID]
       if (telescope == 'Parkes'):  tel = 'PK'
       elif (telescope == 'Effelsberg'):  tel = 'EB'
       elif (telescope == 'Arecibo'):  tel = 'AO'
       elif (telescope == 'MMT'):  tel = 'MT'
       elif (telescope == 'GBT'):  tel = 'GB'
       elif (telescope == 'CHIME'): tel = 'CH'
       else:
          print("Telescope not recognized.")
          return 0
   
   barycenter(tts, bts, vel, ra, dec, tel, ephem)
   avgvel = Num.add.reduce(vel) / nn
   tts = Num.arange(nn, dtype='d') * dt
   bts = (bts - bts[0]) * psr_constants.SECPERDAY
   return tts, bts
