#Module plot_cff
#Created by Aria Coraor

import numpy as np
import argparse
import matplotlib
import jtools
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from subtool import psub
import os

def main():
	""" Plot cv1 vs. cv2 FES.
	"""
	try:
		data = np.loadtxt(args.file)
	except:
		print("Data file can't be numpy-read. Trying again with truncation.")
		cmd = 'sed -n "7,\\$p" %s > tmp' % args.file
		psub(cmd)
		data = np.loadtxt("tmp")
		os.remove("tmp")
		print("Manual truncation succeeded.")

	plot_data(data)
	plot_data(data,fh=True)

def gaussian(x,y, height,sigmax, sigmay):
    return height*np.exp( -( (x**2 / (2*sigmax**2)) + (y**2 / (2*sigmay**2) ) ) )

def plot_data(data,fh=False):
	"""Plot the data on the first 2 cvs, averaging over the 3rd.

	Parameters:
		data: *np.array*, shape (n,8)
			Raw output from the CFF method, cut off at 7 lines to be numpy-readable.
			First 3 columns are cvs, and last column is estimated free energy.
		fh: *bool*
            If true, cut first 10% of data and divide by t. Result should appear
            flat.
	"""
	#Verify data feng integrity
	#print("data[:,-1].shape:",data[:,-1].shape)
	#print("Minimum free engs:",sorted(data[:,-1])[:10])
	if fh:
		data = data[len(data)//2:]
	
	#First: Construct 2D datastructure
	lower = np.array([50,50])
	upper = np.array([250,250])
	ngrid = 100
	X = np.linspace(lower[0],upper[0],ngrid)
	Y = np.linspace(lower[1],upper[1],ngrid)

	XX,YY = np.meshgrid(X,Y)
	feng = np.zeros(XX.shape)

	for hill in data:
		curr_x = XX - hill[2]
		curr_y = YY - hill[1]
		feng -= gaussian(curr_x,curr_y,hill[-1],hill[4],hill[3])

	#Renormalize data
	if not fh:
		feng = feng - np.min(feng)
	else:
		#feng /= len(data)//2
		feng -= np.average(feng)
		plt.title("Second-half G variation")

	#Actually plot
	plt.clf()
	#Flip rb to axis 1
	#YY,XX = np.meshgrid(xvals[:,0],xvals[:,1])
	
	if not fh:
		b = plt.contourf(XX,YY,feng,levels=np.arange(np.min(feng),np.max(feng)),
			cmap='Spectral_r')
		plt.contour(XX,YY,feng,levels=np.arange(np.min(feng),np.max(feng)),colors='black',
			linewidths=0.5)
	else:
		#Special plotting for fh
		width = np.max(feng) - np.min(feng)
		levels = np.arange(np.floor(np.min(feng)),np.ceil(np.max(feng)))
		b = plt.contourf(XX,YY,feng,levels=levels,	cmap='viridis')
		plt.contour(XX,YY,feng,levels=levels,colors='black',
			linewidths=0.5)
		
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
	parser.add_argument('-f','--file', type=str, default=None,help='Path to hills.out file')
	#parser.add_argument('-a','--all',default=False,action='store_const',const=True, help="Calculate helical parameters for all datafiles in NRL range. Output to angles, radii files.")
	args = parser.parse_args()
	main()
