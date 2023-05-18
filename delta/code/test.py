import os
import numpy as np
#import pandas as pd
from aicsimageio import AICSImage
#import pickle
#from pathlib import Path
#import glob

from skimage import filters, io, morphology, feature, img_as_ubyte
from skimage.measure import label, regionprops
from tifffile import imsave, imread
#import sys
import csv

dirpath = 'tiff'
fieldnames = ['date','genotype','name', 'sampleID', 'shape', 'gfp','brDisc', 'brDisc_gfp','brDisc_gfpNEG','br_ratio','dl','dl_gfp','dl_gfpNEG','dl_sd','dl_neg_sd','dl_ratio', 'dl_bkg'] 

with open('output.csv', 'w') as out:
	out.writer = csv.DictWriter(out, delimiter = ',', fieldnames = fieldnames)
	out.writer.writeheader()
	for file in os.listdir(dirpath):

		date = file.split('_')[0] 
		genotype = file.split('_')[1]
		name = file.split('_')[2].split('.')[0]
	#	print(name)
	#	print(genotype)
	
		file_path = os.path.join(dirpath, file)
#		print(file_path)
		image = AICSImage(file_path)
		shape = image.shape
		dapi = image.data[:,0,:,:,:]
		gfp = image.data[:,1,:,:,:]
		brd = image.data[:,2,:,:,:]
		dl = image.data[:,3,:,:,:]
	
#		print(name)
#		print(shape)
		
		dapi_gauss = filters.gaussian(dapi, sigma = 5)	
		dapi_threshold = filters.threshold_mean(dapi_gauss)
		dapi_mask = dapi_gauss > 0.1

			
		gfp_gauss = filters.gaussian(gfp, sigma = 10)
		gfp_threshold = filters.threshold_minimum(gfp_gauss)
		gfp_cut = np.where(dapi_mask == False, 0, gfp_gauss)
		#gfp_mask = gfp_cut >= gfp_threshold
		gfp_mask = gfp_gauss >= 0.09
		
		gfpNEG_dapi = np.where(gfp_mask == True, 0, dapi_gauss)
		gfpNEG_mask = gfpNEG_dapi > 0.1
#		gfpNEG_mask = np.array(np.subtract(gfp_mask, dapi_mask, dtype = float), dtype = bool)
	
		brDisc_gfp = np.zeros_like(brd)
		brDisc_gfp[gfp_mask] = brd[gfp_mask]
		brDisc_gfp_signal = np.mean(brDisc_gfp[gfp_mask])
	
		brDisc_gfpNEG = np.zeros_like(brd)
		brDisc_gfpNEG[gfpNEG_mask] = brd[gfpNEG_mask]
		brDisc_gfpNEG_signal = np.mean(brDisc_gfpNEG[gfpNEG_mask])
	
		br_bkg = np.mean(brd[~dapi_mask])	

		br_ratio = (brDisc_gfp_signal-br_bkg) / (brDisc_gfpNEG_signal-br_bkg)
	###
		dl_gfp = np.zeros_like(dl)
		dl_gfp[gfp_mask] = dl[gfp_mask]
		dl_gfp_signal = np.mean(dl_gfp[gfp_mask])
		dl_gfp_sd = np.std(dl_gfp[gfp_mask])
	
		dl_gfpNEG = np.zeros_like(dl)
		dl_gfpNEG[gfpNEG_mask] = dl[gfpNEG_mask]
		dl_gfpNEG_signal = np.mean(dl_gfpNEG[gfpNEG_mask])
		dl_gfpNEG_sd = np.std(dl_gfpNEG[gfpNEG_mask])

		dl_bkg = np.mean(dl[~dapi_mask])	
			
		dl_ratio = (dl_gfp_signal-dl_bkg) / (dl_gfpNEG_signal-dl_bkg)
	
		output = { 
			'date' : date,
			'genotype' : genotype,
			'name' : name,
			'sampleID' : f'{genotype}_{name}',
			'shape' : '.'.join(map(str,shape)),
			#'gfp' : np.mean(gfp),
			#'brDisc' : np.mean(brd),
			'brDisc_gfp' : brDisc_gfp_signal,
			'brDisc_gfpNEG' : brDisc_gfpNEG_signal,
			'br_ratio' : br_ratio,
			#'dl' : np.mean(dl),
			'dl_gfp' : dl_gfp_signal,
			'dl_gfpNEG' : dl_gfpNEG_signal,
			'dl_sd' : dl_gfp_sd,
			'dl_neg_sd' : dl_gfpNEG_sd,
			'dl_ratio' : dl_ratio,
			'dl_bkg' : dl_bkg
		}
		
		out.writer.writerow(output)	
	
		io.imsave(f'output/{date}_{genotype}_{name}-gfp.tiff', gfp)
		io.imsave(f'output/{date}_{genotype}_{name}-dapi.tiff', dapi)
		io.imsave(f'output/{date}_{genotype}_{name}-brD.tiff', brd)
		io.imsave(f'output/{date}_{genotype}_{name}-gfp_mask.tiff', gfp_mask)
		io.imsave(f'output/{date}_{genotype}_{name}-gfpNEG_mask.tiff', gfpNEG_mask)
		io.imsave(f'output/{date}_{genotype}_{name}-dapi_mask.tiff', dapi_mask)
		io.imsave(f'output/{date}_{genotype}_{name}-dl_gfp.tiff', dl_gfp)
		io.imsave(f'output/{date}_{genotype}_{name}-dl_gfpNEG.tiff', dl_gfpNEG)
#for i in img.get_iter_image():
#	for z in i.get_iter_z():
#		print(z)
#		io.imsave(f'test-{i}-{z}.tiff',np.array(z))
#print(image.shape)
#print(image)

#
#print(dapi)
#print(dapi.shape)
#io.imsave('dapi.tiff', dapi)
#io.imsave('gfp.tiff', gfp)
#io.imsave('br.tiff', brd)
#io.imsave('dl.tiff', dl)

#	print(stack)
#	AICSImage.save(stack, 'test.tiff')
#	io.imsave('test-gfp.tiff', gfp)
