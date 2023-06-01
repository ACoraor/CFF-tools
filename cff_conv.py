#Module cff_conv
#Created by Aria Coraor

import numpy as np
import argparse
import matplotlib
import jtools
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from subtool import psub
import pandas as pd
import os

def main():
	""" Plot convergence of CFF error over epochs.
	"""
	#Load the file
	txt = load_stdout(args.file)
	
	#Process to rmses
	rmses = get_rmse(txt)
	
	#Plot mse, grad
	plot_data(rmses)

def load_stdout(fn):
	"""Load the stdout as plaintext.

	Parameters:
		fn: *str*
			Path to SSAGES stdout log.
	Returns:
		txt: *list of str*
			Each line of the text file.
	"""
	with open(fn,'r') as f:
		txt = f.readlines()
	return txt

def get_rmse(txt):
	"""Process the raw logs into columns of iter, mse, gamma, mu, grad.

	Parameters:
		txt: *list of str*
			Each line of the text file.
	Returns:
		data: *pd.DataFrame*, shape: (n_iters, 5)
			Iteration data from CFF. Columns are:
			<Iter, MSE, gamma, mu, grad>
	"""
	#Get iterlines and gamma lines exclusively.
	iters = [ele for ele in txt if ele[:5] == "iter:" or ele[:6] == "gamma1")]
	
	#Get indices of last iterations before gamma1, which are before biasing.
	gamma_inds = [ind for ind in range(len(iters)) if iters[ind][:6] == "gamma1"]
	gamma_inds = np.array(gamma_inds,dtype=int)
	last_iters = [iters[i] for i in gamma_inds - 1]
	
	#Turn last_iters into numerical data
	data = np.zeros((len(last_iters),5))
	labels = ["iter", "mse", "gamma", "mu", "grad"]
	df = pd.DataFrame()	

	#Iterate dp by dp, getting values column by column
	for i,it in enumerate(last_iters):
		mse_ind = it.find('mse:')
		data[i,0] = int(it[6:mse_ind-1])
		
		gamma_ind = it.find("gamma:")
		data[i,1] = float(it[mse_ind+5:gamma_ind-1])

		mu_ind = it.find("mu:")
		data[i,2] = float(it[gamma_ind+7:mu_ind-1])
		
		grad_ind = it.find('grad:')
		data[i,3] = int(it[mu_ind+4:grad_ind-1])

		data[i,4] = float(it[grad_ind+6:-1])
		
	
	return txt

def plot_data(data):

def plot_data(data):
	"""Plot the data on the first 2 cvs, averaging over the 3rd.

	Parameters:
		data: *np.array*, shape (n,8)
			Raw output from the CFF method, cut off at 7 lines to be numpy-readable.
			First 3 columns are cvs, and last column is estimated free energy.
		json: *dict*
			Loaded SSAGES input json used to create the datafile.
		fh: *bool*
			If true, cut first 10% of data and divide by t. Result should appear
			flat.
	"""
	#Verify data feng integrity
	#print("data[:,-1].shape:",data[:,-1].shape)
	#print("Minimum free engs:",sorted(data[:,-1])[:10])
	if fh:
		data = data[len(data)/10:]
	
	#First: Construct 2D datastructure
	grid = np.asarray(json['methods'][0]['grid']['number_points'])
	feng = np.zeros(grid)
	#Now, calculate the binmap from positions in data to bin in feng
	lower = np.asarray(json['methods'][0]['grid']['lower'])
	upper = np.asarray(json['methods'][0]['grid']['upper'])
	n_cvs = len(json['CVs'])
	

	delta = ( upper - lower)/grid
	#Assume grid is the same
	xvals = np.linspace(lower+delta*0.5,upper-delta*0.5,grid[0])
	#real_grid = (xvals[1:] + xvals[:-1])/2
	#xvals = real_grid
	deltas = xvals[1,:] - xvals[0,:]
	#print("xvals:",xvals)
	#print("deltas:",deltas)
	#input()
	#print("xvals.shape:",xvals.shape)
	for d in data:
		#print("pos:",d[:n_cvs])
		inds = tuple(np.round((d[:n_cvs] - xvals[0,:])/deltas).astype(int))
		#print("inds:",inds)
		#print("gridval:",xvals[0,:] + inds*deltas)
		#print("")
		#Check that inds haven't been added already
		feng[inds] = d[-1]

	#If n_cvs > 2, average over last cvs
	#print("n_cvs:",n_cvs)
	if n_cvs > 2:
		axes = tuple(np.arange(2,n_cvs))
		ps = np.exp(-feng)
		sum_ps = np.sum(ps,axis=axes)
		feng = -np.log(sum_ps)
	#print("feng shape:",feng.shape)

	#Actually plot
	plt.clf()
	#Flip rb to axis 1
	YY,XX = np.meshgrid(xvals[:,0],xvals[:,1])
	
	b = plt.contourf(XX,YY,np.transpose(feng),levels=np.arange(0,np.max(feng)),cmap='Spectral_r')
	a = plt.contour(XX,YY,np.transpose(feng),levels=np.arange(0,np.max(feng)),
		colors='black',linewidths=0.5)
	#print("feng:",feng)
	cbar = plt.colorbar(b)
	cbar.ax.set_ylabel("Free Energy (kT)")
	#cbar.add_lines(a)
	cv2 = "$r_\\alpha (\\AA)$"
	cv1 = "$r_\\beta (\\AA)$"
	plt.xlabel(cv1)
	plt.ylabel(cv2)

	if not fh:
		fn = "fes.pdf"
	else:
		fn = "fh.pdf"
	dn = os.path.basename(os.getcwd())
	fn = dn + "_" + fn

	plt.tight_layout()
	plt.savefig(fn)
	print("Saved FES at %s" % fn)



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-f','--file', type=str, default='stdout',help='Path to CFF. logfile')
	parser.add_argument('-i','--input', type=str, default='input.json',help='Input json file')
	#parser.add_argument('-a','--all',default=False,action='store_const',const=True, help="Calculate helical parameters for all datafiles in NRL range. Output to angles, radii files.")
	args = parser.parse_args()
	main()
