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
from presto import psrfits
# add filterbank module
from presto import filterbank
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
       tts = Num.arange(tto, tto + (T + dt) / psr_utils.SECPERDAY, dt / psr_utils.SECPERDAY)
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
       tts = Num.arange(tto, tto + (T + dt) / psr_utils.SECPERDAY, dt / psr_utils.SECPERDAY)
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
       tts = Num.arange(tto, tto + (T + dt) / 86400, dt / 86400)  # 86400 are seconds per day
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
       tel = 'GB'
       warnings.warn("Barycenter time will be wrong because telescope name header not read correctly.\
                     Assuming GBT as the telescope location")

   # There was a mismatch in number of parameters passed to barycenter function. 
   # I think this must be crashing for psrfits as well !!
   barycenter(tts, bts, vel, ra, dec, tel, ephem)
   avgvel = Num.add.reduce(vel) / nn
   tts = Num.arange(nn, dtype='d') * dt
   bts = (bts - bts[0]) * 86400  # 86400 are seconds per day
   return tts, bts
