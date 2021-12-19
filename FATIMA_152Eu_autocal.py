#!/usr/bin/python3
import os.path as path
import sys
import datetime
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit
from operator import itemgetter #for sorting lists of lists
#functions can change mutable objects (e.g. lists) but not immutable ones (e.g. int)

def test (array):
   array.append(9)

def read_histogram_twoColumn (infile, xarray, yarray, rebin=0):
   if (rebin < 2):
     with open(infile, "r") as infile:
       for line in infile:
          content = line.split()
          if len(content) > 1:
             xarray.append(float(content[0]))
             yarray.append(float(content[1]))
          else:
             return
   else:
     nrebin = int(rebin)
     rebinContent = 0.
     xpos = 0
     with open(infile, "r") as infile:
       for i, line in enumerate(infile):
          content = line.split()
          if len(content) > 1:
             xcont = float(content[0])
             ycont = float(content[1])
          else:
             return
          if (i > 0 and (i%nrebin) == 0):
             xarray.append(xpos)
             yarray.append(rebinContent)
             xpos = xcont
             rebinContent = ycont
          elif (i == 0):
             xpos = xcont
             rebinContent = ycont
          else:
             rebinContent += ycont

def write_calibration_to_file(oname, detname, calarray):
  with open(oname, "w") as ofile:
    detname += "  "
    ofile.write(detname)
    for c in reversed(calarray):
    #  value = c + "  "
    #  ofile.write(value)
      ofile.write("%.5e  " % c)
    ofile.write("\n")

def total_int (yarray):
   tot_int = 0;
   for p in yarray:
      tot_int += p
   return tot_int

def cut_spec_to_interesting (xarray, yarray, divisions = 20, cutoff = 0.003):
   #First adjust the length to match multiple of divisions
   rest = len(yarray)%divisions
   if (rest != 0):
      for i in range(rest):
         yarray.pop()
         xarray.pop()
   #Integrate each division
   t_int = total_int(yarray)
   part_integral = []
   for part in range(divisions):
      tmp_int = 0
      for i in range(int(part*len(yarray)/divisions), int((part+1)*len(yarray)/divisions)):
         tmp_int += yarray[i]
      part_integral.append(tmp_int/t_int)
   #Throw out divisions that don't make the cut
   bins_per_div = len(xarray)/divisions
   for i, part in enumerate(reversed(part_integral)):
      if part < cutoff:
         del xarray[int(-bins_per_div):]
         del yarray[int(-bins_per_div):]
      else:
         break

#return in [bin, prominence[bin]]
#bincut = 10 is approximately 3 sigma in a naive way
def get_prominence (xarray, yarray, bincut = 10):
   prominence = [0]*len(xarray)
   max_peak = max(yarray)
   for i, part in enumerate(xarray):
      left_min = 0
      right_min = 0
#      if xarray[i] != 199:
#        prominence[i] = 0
#        continue
      #Cut off bottom 5 bins to suppress possible ADC/digitiser noise
      if i < 5:
         prominence[i] = 0
         continue

      if yarray[i] > bincut: 
         j=i#+3
         right_min = yarray[i]
         while (j < len(xarray)):
            if (yarray[i] < yarray[j]):
               break
            right_min = min(right_min, yarray[j])
            j+=1
   
         j=i#-3
         left_min = yarray[i]
         while (j > -1):
            #if (yarray[i] < yarray[j] and abs(yarray[i] - yarray[j]) > np.sqrt(yarray[i])*3):
            if (yarray[i] < yarray[j]):
               break
            left_min = min(left_min, yarray[j])
            j-=1
         prominence[i] = yarray[i] - 1.0*max(left_min, right_min) - np.sqrt(yarray[i])
      else:
         prominence[i] = 0

   return prominence


#peaklist is in format [bin, y[bin]] with bin number, to get channel use x[bin]
                #min relative hight (rel. to the maximum peak (usually x ray peak))
def get_peak_list (xarray, yarray, this_prom, min_rel_hight = 0.005, bin_crit = 5):
   peaklist = []
   max_prom = max(this_prom)
   max_prom_index = this_prom.index(max_prom)
   for i,p in enumerate(this_prom):
      if p > max_prom*min_rel_hight:
         peaklist.append([i,yarray[i]])
         #peaklist.append([xarray[i],yarray[i]])
   #This is now probably obsolete as the peak-finding function has been updated to
   #generally give only one point per peak. But doesn't hurt to have it here in case
   #it still happens sometimes 
   clean_list = []
   i=0
   while(i<len(peaklist)):
      steps = 0
      #change the 8 at one point to an approximate width of a peak
      for j in range(8):
         if i+j >= len(peaklist):
            break
         if xarray[peaklist[i+j][0]] - xarray[peaklist[i][0]] < xarray[-1]*bin_crit/len(xarray) :
            steps += 1
         else:
            break
      clean_list.append(peaklist[int(i + steps/2)])
      i += steps
   return clean_list
   #return peaklist

#returns a factor bins_per_keV based on estimating the 1408 as the last "prominent" peak
def estimate_calibration_factor (peaklist):
   return (peaklist[-1][0]/1408.0)
   

#returns the element of peaklist which is closest to match the energy based on the coarse claibration.
#returns [energy, -1,-1] if there was no peak in peak list that was acceptably close
#return array is [peakname, bin, y[bin], index_in_peaklist]
#(fail_crit is distance in keV)
def assign_peak (peaklist, bins_per_keV, energy, fail_crit = 30):
   distance = 1000000
   best_peak = [distance,-1,-1,-1]
   e_chn = energy*bins_per_keV
   #could be made more efficient by assuming a sorted peak list -- break after distance increases
   #and possibly start in the middle of the peaklist and check which direction to go
   #but will probably not give much of a speed increase

   #print("\n", energy)
   while fail_crit < 70:
     for i,p in enumerate(peaklist):
        distance = abs(e_chn - p[0])
        if  distance < best_peak[0]:
           best_peak = [distance, p[0], p[1], i]
           #print (best_peak, "--", p[0]/bins_per_keV)
     if best_peak[0]/bins_per_keV > fail_crit:
        fail_crit += 10 
     else:
        return [energy, best_peak[1], best_peak[2], best_peak[3]] 
   return [energy, -1, -1, -1]
#returns list of positions in peak list that best match the main peak energies
#returns [energy, -1,-1] if assignment was not possible
def assign_main_peaks (peaklist, per_keV, this_prom):
   main_peak_energies = [40, 121, 244, 344, 779, 1408]
   #main_peak_energies = [244]
   main_peak_list = []
   fcrit = 40
   max_peak = 0
   max_peak_chn = 0
   max_peak_index = -1
   for p in peaklist:
      max_peak = max(max_peak, p[1])
   for i,p in enumerate(peaklist):
      if p[1] == max_peak:
         max_peak_index = i
         max_peak_chn = p[0]
         break
   goon = 1

   if len(peaklist) < 1:
      return 0

   while (goon):
      del main_peak_list[:]
      for p in main_peak_energies:
         main_peak_list.append(assign_peak(peaklist, per_keV, p, fcrit))
      if len(main_peak_list) < 1:
         return 0   
      #Check if the list is complete
      #The assignment assumes an approximately linear calibration taking 0 and 1408
      #as references. This most easily fails at low energies, so the xray and the 121
      #should be checked individually if the peak height of the assigned peaks makes
      #sense. The xray is also often cut by a threshold. If it is not found, the program
      #should just silently go on without using it.
      xray_not_found = 0
      g121_not_found = 0
      badness = 0
      start_i = 0
      if main_peak_energies[0] == 40:
         if main_peak_list[0][1] == -1:
            xray_not_found = 1
            del main_peak_list[0:1]
            del main_peak_energies[0:1]
            print("X-ray peak could not be identified.")
            badness += 1
         elif main_peak_list[0][1] != max_peak:
           if max_peak_chn < main_peak_list[1][1]:
              main_peak_list.append([40, peaklist[max_peak_index][0], peaklist[max_peak_index][1], max_peak_index])
              del main_peak_list[0:1]
              main_peak_list = sorted(main_peak_list, key = itemgetter(0))
           else:
              del main_peak_list[0:1]
              del main_peak_energies[0:1]
              print("X-ray peak could not be identified.")
              xray_not_found = 1

      p121_i = -1
      p244_i = -1
      p121 = [121,-1,-1,-1]
      p244 = [244,-100,-1,-1]
      energy_list_i_121 = -1
      for li,l in enumerate(main_peak_energies):
         if l == 121:
           for i,p in enumerate(main_peak_list):
              if p[0] == 121:
                 p121 = p
                 p121_i = i
              if p[0] == 244:
                 p244 = p
                 p244_i = i
           if p121[2] < p244[2]:
              if p121[3] + 1 < p244[3]:
                 maxininterval = 0
                 for j in range(p121[3] + 1, p244[3]):
                    maxininterval = max(maxininterval, peaklist[j][1])
                 for i,p in enumerate(peaklist):
                    if p[1] == maxininterval:
                      main_peak_list.append([121, p[0], p[1], i])
                      del main_peak_list[p121_i:p121_i + 1]
                      main_peak_list = sorted(main_peak_list, key = itemgetter(0))
           else:
              if p121[2] == -1:
                  g121_not_found = 1
                  del main_peak_list[p121_i:p121_i + 1]
                  del main_peak_energies[li:li+1]
                  print("121 keV peak could not be identified.")

      for i in range(0,len(main_peak_list)):
         if main_peak_list[i][1] == -1:
            badness += 1

      if badness == 0:
         for e,i,p,ind in main_peak_list:
            plt.plot(Cutx[i],p,'C8o')
         return main_peak_list
      #If badness is there, try again with larger fcrit to still get a match.
      elif fcrit < 60:
         fcrit = 60
         continue
      else:
         print("Could not assign peaks. Is this a 152Eu decay spectrum?")
         exit(1)
   print("should never get here -- in assign_main_peaks")
   exit(1)
   

def assign_side_peaks (peaklist, per_keV, this_prom):
   side_peak_energies = [411444, 867, 964, 10861112]
   #side_peak_energies = [411444, 867, 964, 10861112]
   #main_peak_energies = [244]
   side_peak_list = []
   fcrit = 40
   goon = 1
   while (goon):
      for p in side_peak_energies:
         if (p == 411444):
            side_peak_list.append(assign_peak(peaklist, per_keV, 444., 35))
            index_in_peaklist = side_peak_list[-1][3]
            #the 411 is easily wrongly assigned as the 444 because they are close together
            #the fit algorithm later expects to have the 444, so we must make sure that
            #this one is not wrongly assigned. Check the peak in peaklist just before the
            #one assigned to 444. It should be 344, and therefore much higher. If not, just
            #take the next peak and hope for the best.
            print(peaklist[index_in_peaklist][0], peaklist[index_in_peaklist - 1][0],\
             (peaklist[index_in_peaklist][0] - peaklist[index_in_peaklist - 1][0])/per_keV)
            if peaklist[index_in_peaklist - 1][1] < 2*side_peak_list[-1][2]:
               side_peak_list[-1][0] = 411444
            else:
               if (peaklist[index_in_peaklist + 1][0]/per_keV - 444.) < 60:
                  side_peak_list.pop()
                  tmparray = []
                  tmparray.append(411444)
                  #tmparray.append(0)
                  #tmparray.append(0)
                  tmparray.append(peaklist[index_in_peaklist + 1][0])
                  tmparray.append(peaklist[index_in_peaklist + 1][1])
                  tmparray.append(index_in_peaklist + 1)
                  side_peak_list.append(tmparray)
               else:
                  side_peak_list[-1][0] = 411444
         elif (p == 10861112):
            side_peak_list.append(assign_peak(peaklist, per_keV, 1112., 25))
            side_peak_list[-1][0] = 10861112
         else:
            side_peak_list.append(assign_peak(peaklist, per_keV, p, fcrit))
   
      #Check if the list is complete
      #   121_not_found = 0
      badness = 0
      start_i = 0
      clean_list = []
      for i in range(start_i,len(side_peak_list)):
         if side_peak_list[i][1] == -1:
            print("Could not assign", side_peak_list[i][0], " peak(s).")
            badness += 1
#            del side_peak_list[i:i+1]
         else:
           clean_list.append(side_peak_list[i])
      print(clean_list)
      if badness == 0:
         for e,i,p,ind in clean_list:
            print(Cutx[i],p)
            plt.plot(Cutx[i],p,'C8o')
         return clean_list
      #If badness is there, try different things to still get a match.
      #first check if only the xray is not found. This could be due to a high threshold.
      #remove the xray and check again. If it works, high threshold is probably the issue
      elif badness > 0 and badness < len(side_peak_energies):
         print ("Could not assign all minor peaks. Check the peak finder sensitivity to improve calibration.")
         for e,i,p,ind in clean_list:
            print(Cutx[i],p)
            plt.plot(Cutx[i],p,'C1o')
         return clean_list
      else:
         print ("Could not assign any minor peaks. Check peak finder sensitivity to improve calibration.")
         return []
         

#returns width in channels based on FATIMA mean fwhm measured at S4, GSI (see NIM A paper)
def get_fwhm_bins (bins_per_keV, energy):
   return bins_per_keV*(117.3*energy**(-0.53172)*energy/100)

def get_fwhm_keV (energy):
   return (117.3*energy**(-0.53172)*energy/100)


def fit_peak (xarray, yarray, peakarray, xl, xr, init_par):
   try:
      popt, pcov = curve_fit(gauss, xarray[xl:xr], yarray[xl:xr], init_par)
   except RuntimeError:
      print ("Fit did not converge for", peakarray[0])
      return [peakarray[0], -1]
   #popt = init_par
   plt.plot(xarray[xl], yarray[xl], 'C2<')
   plt.plot(xarray[xr], yarray[xr], 'C2>')
   xrng = np.linspace(0, xarray[-1], len(xarray))
   #print (xl, xr, xrng[xl], xrng[xr])
   plt.plot(xrng[xl:xr], gauss(xrng[xl:xr], popt[0], popt[1], popt[2], popt[3], popt[4]), 'r-')
   gtop =  gauss(popt[2], popt[0], popt[1], popt[2], popt[3], popt[4])
   #print (peakarray[0], popt[0], popt[1], popt[2], gtop)
   plt.plot(popt[2], gtop, 'C2^')
   return [peakarray[0], popt[0], popt[1], popt[2], gtop] 

def fit_double_peak (xarray, yarray, peakarray, xl, xr, init_par):
   #print (init_par)
   try:
     popt, pcov = curve_fit(doublegauss, xarray[xl:xr], yarray[xl:xr], init_par)
   except RuntimeError:
     print ("Fit did not converge for", peakarray[0])
     popt = init_par
     plt.plot(xarray[xl], yarray[xl], 'r<')
     plt.plot(xarray[xr], yarray[xr], 'r>')
     xrng = np.linspace(0, xarray[-1], len(xarray))
     plt.plot(xrng[xl:xr],\
      doublegauss(xrng[xl:xr], popt[0], popt[1], popt[2],\
      popt[3], popt[4], popt[5], popt[6], popt[7]), 'r-')
     dgtop1 = doublegauss(popt[2], popt[0], popt[1], popt[2], popt[3],\
      popt[4], popt[5], popt[6], popt[7])
     dgtop2 = doublegauss(popt[5], popt[0], popt[1], popt[2], popt[3],\
      popt[4], popt[5], popt[6], popt[7]) 
     plt.plot(popt[2], dgtop1, 'r^')
     plt.plot(popt[5], dgtop2, 'r^')
     return [peakarray[0], -1]

   #popt = init_par
   #print (peakarray[0], popt)
   plt.plot(xarray[xl], yarray[xl], 'C2<')
   plt.plot(xarray[xr], yarray[xr], 'C2>')
   xrng = np.linspace(0, xarray[-1], len(xarray))
   plt.plot(xrng[xl:xr],\
    doublegauss(xrng[xl:xr], popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7]), 'r-')
   dgtop1 = doublegauss(popt[2], popt[0], popt[1], popt[2], popt[3],\
    popt[4], popt[5], popt[6], popt[7])
   dgtop2 = doublegauss(popt[5], popt[0], popt[1], popt[2], popt[3],\
    popt[4], popt[5], popt[6], popt[7]) 
   #print (peakarray[0], popt[0], popt[1], popt[2], dgtop1,\
   # popt[3], popt[4], popt[5], dgtop2)
   plt.plot(popt[2], dgtop1, 'C2^')
   plt.plot(popt[5], dgtop2, 'C2^')
   return [peakarray[0], popt[0], popt[1], popt[2], dgtop1,\
    popt[3], popt[4], popt[5], dgtop2]


#returns fit results for the main 152Eu peaks
def fit_main_peaks (xarray, yarray, main_peak_list, per_keV, binw):
   fit_results = []
   fit_result = [] #peakname, area, sigma, position, yvalue
   init_par   = [0,0,0,0,0] #area, sigma, position, slope, absz

   
   for p in main_peak_list:
      #sigma
      fwhm_bins = get_fwhm_bins(per_keV, p[0])
      fwhm_chn  = (xarray[int(fwhm_bins)] - xarray[0])
      init_par[1] = fwhm_chn/2.355
      #position
      pos_chn = xarray[p[1]]
      init_par[2] = pos_chn
      #calculate fit region limits
      #also used to get linear background
      if (p[0] == 40):
        xtmp_l_chn = pos_chn - 2.*fwhm_chn
        xtmp_r_chn = pos_chn + 3.*fwhm_chn
        xtmp_l_bins = int(p[1] - 2.*fwhm_bins)
        xtmp_r_bins = int(p[1] + 3.*fwhm_bins)
      else:
        xtmp_l_chn = pos_chn - 2.*fwhm_chn
        xtmp_r_chn = pos_chn + 2.*fwhm_chn
        xtmp_l_bins = int(p[1] - 2.*fwhm_bins)
        xtmp_r_bins = int(p[1] + 2.*fwhm_bins)
      #slope
      init_par[3] = (yarray[xtmp_r_bins] - yarray[xtmp_l_bins])/(xtmp_r_chn - xtmp_l_chn) 
      #absz
      init_par[4] = yarray[xtmp_r_bins] - init_par[3]*xtmp_r_chn
      #area
      tmp_area = 0
      bg_estimate = yarray[xtmp_l_bins]
      for t in range(xtmp_l_bins, xtmp_r_bins):
         tmp_area += yarray[t] - bg_estimate
      init_par[0] = binw*tmp_area
      #print (xtmp_l_bins, xtmp_r_bins)
      #print (init_par)
      fit_result = fit_peak(xarray, yarray, p, xtmp_l_bins, xtmp_r_bins, init_par)
      if fit_result[1] != -1:
         fit_results.append(fit_result)
   return fit_results 

def fit_side_peaks (xarray, yarray, side_peak_list, per_keV, binw):
   fit_results = [] #peakname, area, sigma, position, yvalue
   tr = []
   init_par          = [0,0,0,0,0] #area, sigma, position, slope, absz
   init_par_double   = [0,0,0,0,0,0,0,0] #area1, sigma1, position1, area2, sigma2, position2, slope, absz
   for p in side_peak_list:
      del tr[:]
      #sigma
      if p[0] == 411444:
        fwhm_bins = get_fwhm_bins(per_keV, 444)
      elif p[0] == 10861112:
        fwhm_bins = get_fwhm_bins(per_keV, 1112)
      else:
        fwhm_bins = get_fwhm_bins(per_keV, p[0])
      fwhm_chn  = (xarray[int(fwhm_bins)] - xarray[0])
      init_par[1] = fwhm_chn/2.355
      #position
      pos_chn = xarray[p[1]]
      init_par[2] = pos_chn
      #calculate fit region limits
      #also used to get linear background
      xtmp_l_chn = pos_chn - 2.*fwhm_chn
      xtmp_r_chn = pos_chn + 2.*fwhm_chn
      xtmp_l_bins = int(p[1] - 2.*fwhm_bins)
      xtmp_r_bins = int(p[1] + 2.*fwhm_bins)
      #slope
      init_par[3] = (yarray[xtmp_r_bins] - yarray[xtmp_l_bins])/(xtmp_r_chn - xtmp_l_chn) 
      #absz
      init_par[4] = yarray[xtmp_r_bins] - init_par[3]*xtmp_r_chn
      #area
      tmp_area = 0
      bg_estimate = yarray[xtmp_r_bins]
      for t in range(xtmp_l_bins, xtmp_r_bins):
         tmp_area += yarray[t] - bg_estimate
      init_par[0] = tmp_area*binw
      if (p[0] == 411444):
         xtmp_l_bins = int(p[1] - 3.5*fwhm_bins)
         init_par_double[0] = init_par[0]/2
         init_par_double[3] = init_par[0]/2
         init_par_double[1] = fwhm_chn/2.355
         init_par_double[4] = fwhm_chn/2.355
         #slope
         init_par_double[6] = (yarray[xtmp_r_bins] - yarray[xtmp_l_bins])/(xtmp_r_chn - xtmp_l_chn) 
         #absz
         init_par_double[7] = yarray[xtmp_r_bins] - init_par_double[6]*xtmp_r_chn
         init_par_double[5] = init_par[2]
         init_par_double[2] = init_par[2] - 2*fwhm_chn
         tr = fit_double_peak(xarray, yarray, p, xtmp_l_bins, xtmp_r_bins, init_par_double)
         if tr[3] < tr[7]:
           fit_results.append([411,tr[1],tr[2],tr[3],tr[4]])
           fit_results.append([444,tr[5],tr[6],tr[7],tr[8]])
         else:
           fit_results.append([444,tr[1],tr[2],tr[3],tr[4]])
           fit_results.append([411,tr[5],tr[6],tr[7],tr[8]])
      elif (p[0] == 10861112):
         xtmp_l_bins = int(p[1] - 3*fwhm_bins)
         init_par_double[0] = init_par[0]/2
         init_par_double[3] = init_par[0]/2
         init_par_double[1] = fwhm_chn/2.355
         init_par_double[4] = fwhm_chn/2.355
         #slope
         init_par_double[6] = (yarray[xtmp_r_bins] - yarray[xtmp_l_bins])/(xtmp_r_chn - xtmp_l_chn) 
         #absz
         init_par_double[7] = yarray[xtmp_r_bins] - init_par_double[6]*xtmp_r_chn
         init_par_double[5] = init_par[2]
         init_par_double[2] = init_par[2] - fwhm_chn/2
         tr = fit_double_peak(xarray, yarray, p, xtmp_l_bins, xtmp_r_bins, init_par_double)
         if tr[3] < tr[7]:
           fit_results.append([1086,tr[1],tr[2],tr[3],tr[4]])
           fit_results.append([1112,tr[5],tr[6],tr[7],tr[8]])
         else:
           fit_results.append([1086,tr[1],tr[2],tr[3],tr[4]])
           fit_results.append([1112,tr[5],tr[6],tr[7],tr[8]])
      else:
         fit_results.append(fit_peak(xarray, yarray, p, xtmp_l_bins, xtmp_r_bins, init_par))
   return fit_results

def check_peaks (peaklist, original):
  peakheight = {
    40:   66.,
    121:  30.,
    244:  7.8,
    344:  10.0,
    411:  1.8,
    444:  1.9,
    779:  1.9,
    867:  1.2,
    964:  1.5,
    1086: 0.9,
    1112: 1.2,
    1408: 1
  }
  Iref = -1
  for p in peaklist:
    if p[0] == 1408:
      Iref = p[2]
  if Iref < 0:
    print("1408 keV line not identified.")
    return 0
  bad = []
  for i,p in enumerate(peaklist):
    if p[0] == 411444:
      e = 444
    elif p[0] == 10861112:
      e = 1112
    else:
      e = p[0]

    #if the height is very off, try to re-assign the peak. Test the two or three
    #previous ones, but avoid double assignments! Doesn't work very well for 244
    #when having doubles as well as singles histograms...
    if e == 244:
      continue
    if abs(p[2]/Iref - peakheight[e]) > peakheight[e]*0.3:
      print ("The following peak seems to be wrongly assigned:")
      print(p, abs(p[2]/Iref), peakheight[e], peakheight[e]*0.3, Iref)
      bad.append(i)
    else:
      print("good ", p, abs(p[2]/Iref), peakheight[e], peakheight[e]*0.3, Iref)

  if len(bad) > int(len(peaklist)/2):
    print ("Too many peaks are off. Will not try re-assignment.")
  else:
    for i in bad:
      already_in_list = 0
      if peaklist[i][0] == 411444:
        e = 444
      elif peaklist[i][0] == 10861112:
        e = 1112
      else:
        e = peaklist[i][0]
      for j in range(peaklist[i][3] - 3, peaklist[i][3] + 3):
        if abs(original[j][1]/Iref - peakheight[e]) < peakheight[e]*0.3:
          print("candidate", original[j][0])
          for k, chli in enumerate(peaklist):
            if k != i and original[j][0] == chli[1]:
              already_in_list = 1
          if already_in_list == 0:
            print("Found a candidate to re-assign")
            tmparray = []
            if e == 444:
              tmparray.append(411444)
            elif e == 1112:
              tmparray.append(10861112)
            else:
              tmparray.append(e)
            tmparray.append(original[j][0])
            tmparray.append(original[j][1])
            tmparray.append(j)
            peaklist.append(tmparray)
            print(tmparray)
            break
          else:
            already_in_list = 0
    for i in bad:
      peaklist.pop(i)
          
          

def calibrate (peaklist):
   return 1

def gauss(xx,a,s,m,c,b):
   return a/(s*2.506631)*np.exp(-0.5*(xx-m)**2/s**2) + xx*c + b

def doublegauss(xx,a1,s1,m1,a2,s2,m2,c,b):
   return a1/(s1*2.506631)*np.exp(-0.5*(xx-m1)**2/s1**2) + a2/(s2*2.506631)*np.exp(-0.5*(xx-m2)**2/s2**2) +xx*c + b

def straight(xx,b,c):
   return b*xx + c

def poly3(xx, a,b,c,d):
   return a*xx*xx*xx + b*xx*xx + c*xx + d

def poly4(xx, a,b,c,d,e):
   return a*xx*xx*xx*xx + b*xx*xx*xx + c*xx*xx + d*xx + e


##############
#MAIN:
##############

Ax = []
Ay = []
Cutx = []
Cuty = []
#infile = "first_test/152Eu_det19.gnu"
#infile = "cpar001.spe.gnu"
#infile = "LaBr3.tv.gnu"
#infile = "LaBr3.xraycut.gnu"

if(len(sys.argv) < 2):
   print ("Usage: autocal.py <histogram file name> (<detnumber>)")
   exit(0)
elif(len(sys.argv) == 2):
   DETNAME = ""
else:
   DETNAME = str(sys.argv[2])
infile = str(sys.argv[1])

plt.subplot(2,1,1)

print ("Calibrating for file", infile)
#Getting the base of the file name for use as output file names later and for setting the figure title
basename = path.basename(infile)
iod = basename.index('.')
pdfname = basename[:iod] + ".pdf"
outname = basename[:iod] + ".cal"
barename = basename[:iod]
now = datetime.datetime.now()
titlestr = now.strftime("%Y-%m-%d  %H:%M") + "     " + basename

plt.title(titlestr)

##############
#-----Preparing the histogram
##############
avg_counts_per_bin = 0
binning_factor = 1
while avg_counts_per_bin < 200:
  del Cutx[:]
  del Cuty[:]
  #print ("binning", binning_factor)
  read_histogram_twoColumn(infile, Ax, Ay, binning_factor)
  Cutx = Ax
  Cuty = Ay
  cut_spec_to_interesting(Cutx, Cuty)
  if len(Cutx) < 300:
    print ("Looks like the spectrum", infile, "doesn't have enough statistics.")
    break
  avg_counts_per_bin = total_int(Cuty)/len(Cuty)
  binning_factor += 1    

binwidth = Cutx[1] - Cutx[0] 
 
plt.plot(Cutx, Cuty, 'b-')
print ("Did rebin with factor", binning_factor)

##############
#-----Finding and assigning the peaks and fitting them
##############
prominence = get_prominence(Cutx, Cuty)
#print (prominence)
#plt.step(Cutx, prominence, where='mid') 

peaks = get_peak_list(Cutx, Cuty, prominence)

for i,p in peaks:
  plt.plot(Cutx[i],p,'C4o')

linCal = estimate_calibration_factor(peaks)

mainPeaks = assign_main_peaks(peaks, linCal, prominence)
sidePeaks  = assign_side_peaks(peaks, linCal, prominence)
allPeaks  = mainPeaks + sidePeaks

check_peaks(mainPeaks, peaks)
#peak check doesn't work yet for sidepeaks... it wants the 1408 in the list...
#check_peaks(sidePeaks, peaks)

results = []
results += fit_main_peaks(Cutx, Cuty, mainPeaks, linCal, binwidth)
results += fit_side_peaks(Cutx, Cuty, sidePeaks, linCal, binwidth)

print ("RESULTS OF THE PEAKFINDER:")
print ("========")
print ("peak, area, sigma, position, high-point")

results = sorted(results, key = itemgetter(0))
CalE   = []
CalChn = []

energyvalue = {
  40:   41.3122,
  121:  121.7817,
  244:  244.6975,
  344:  344.2785,
  411:  411.1,
  444:  443.965,
  779:  778.9040,
  867:  867.378,
  964:  964.079,
  1086: 1085.869,
  1112: 1112.074,
  1408: 1408.006
}

for r in results:
  CalE.append(energyvalue[r[0]])
  CalChn.append(r[3])
  print (r)


plt.yscale('log')
plt.ylim(3,)

##############
#-----Doing the calibration fit
##############
plt.subplot(2,1,2)

#Attempt a fit with a 3rd order polynomial
init_par3 = [1e-15, 1e-9, linCal, 0]
params3, pcov3 = curve_fit(poly3, CalChn, CalE, init_par3)

var = 0
for i,r in enumerate(results):
  var += (CalE[i] -  poly3(r[3], params3[0], params3[1], params3[2], params3[3]))**2
mean_dev = np.sqrt(var)/len(r)
print("----Fit result for 3rd order polynomial:")
print(params3)
print("----mean deviation with 3rd order poly= %.4f" % mean_dev, "keV\n")

used_poly4 = 0
#If 3rd degree fit was not good enough, try with 4th order
if mean_dev > 0.5:
  print("Mean deviation is too large with polynomial of order 3. Trying 4th order...")
  used_poly4 = 1
  init_par4 = [1e-16, 1e-11, 1e-8, linCal, 0]
  params4, pcov4 = curve_fit(poly4, CalChn, CalE, init_par4)
  var = 0
  for i,r in enumerate(results):
    var += (CalE[i] -  poly4(r[3], params4[0], params4[1], params4[2], params4[3], params4[4]))**2
  mean_dev = np.sqrt(var)/len(r)
  print("----Fit result for 4th order polynomial:")
  print(params4)
  print("----mean deviation with 4th order poly= %.4f" % mean_dev, "keV\n")

#plot the result and write parameters to file

if used_poly4:
  for i,r in enumerate(results):
#    print (r[0], r[3], poly4(r[3], params4[0], params4[1], params4[2], params4[3], params4[4]))
    plt.plot(CalE[i], CalE[i] -  poly4(r[3], params4[0], params4[1], params4[2],\
     params4[3], params4[4]), 'C3x')
    write_calibration_to_file(outname, barename, params4)
else:
  for i,r in enumerate(results):
#    print (r[0], r[3], poly3(r[3], params3[0], params3[1], params3[2], params3[3]))
    plt.plot(CalE[i], CalE[i] -  poly3(r[3], params3[0], params3[1], params3[2], params3[3]), 'C3x')
    write_calibration_to_file(outname, barename, params3)

plt.savefig("test.pdf")

plt.show()

