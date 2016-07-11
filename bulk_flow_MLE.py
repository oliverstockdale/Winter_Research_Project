import numpy as np

def convertAngles():
	raDeg = []
	decDeg = []
	
	# open files with position of galaxies in equatorial coordinates
	openRA = open('/Users/oliverstockdale/Dropbox/University/Summer_Research/Code/Data/raAngle.txt')
	openDEC = open('/Users/oliverstockdale/Dropbox/University/Summer_Research/Code/Data/decAngle.txt')
	
	for line in openRA.read().split('\n'):
		raDeg = np.append(raDeg, float(line))
		
	for line in openDEC.read().split('\n'):
		decDeg = np.append(decDeg, float(line))
	

	# converts degrees to radians
	raRad = np.radians(raDeg)
	decRad = np.radians(decDeg)
	
	# converts from equatorial coordinates to cartesian
	x = np.cos(raRad) * np.cos(decRad)
	y = np.sin(raRad) * np.cos(decRad)
	z = np.sin(decRad)

	
	global posVec
	# creates a vector to store coordinates
	posVec = np.zeros((len(raRad), 3))
	
	for i in range(0, len(raRad)):
		posVec[i] = (x[i], y[i], z[i])
	
	#closes all files
	openRA.close()
	openDEC.close()
	return()

def MLE():
	peculiarVel = []
	v_p = open('/Users/oliverstockdale/Dropbox/University/Summer_Research/Code/Data/peculiarVel.txt')
	
	for line in v_p.read().split('\n'):
		peculiarVel = np.append(peculiarVel, float(line))
	
	sigma_star = 300
	sigma_n =  0.2 * peculiarVel
	
	A = np.zeros((3,3))
	
	n = 0
	#while loop populates the A matrix, calculating each element and adding it to the 
	#previous value to sum over n galaxies
	
	for n in range(0, len(peculiarVel)):
		for j in range(0,3):
			for i in range(0,3):
				A[i,j] += (posVec[n, i] * posVec[n, j])/(sigma_star**2 + sigma_n[n]**2)
	
	
	A_inv = np.linalg.inv(A)
	
	def weight(i,n):
		#calculates the weight of a given galaxy with component i
		w_in = 0
		for j in range(0,3):
			w_in += A_inv[i,j] * (posVec[n,j]/(sigma_star**2 + sigma_n[n]**2))
		return w_in
	
	
	def bulkFlow(i):
		#calculates the bulk flow based off the weight 
		u = 0
		for n in range(0, len(peculiarVel)):
			u += weight(i,n) * peculiarVel[n]
		
		return u
    
	global u_x, u_y, u_z
	u_x = bulkFlow(0)
	u_y = bulkFlow(1)
	u_z = bulkFlow(2)
	
	print 'The bulk flow is ', np.sqrt(u_x**2 + u_y**2 + u_z**2)
	v_p.close()
	return()
	
def convertCart():
	decRad = np.arcsin(u_z/(np.sqrt(u_x**2 + u_y**2 + u_z**2)))
	raRad = np.arcsin((u_y/(np.sqrt(u_x**2 + u_y**2 + u_z**2)))/np.cos(decRad))
	print np.degrees(decRad)
	print np.degrees(raRad)
	return()


convertAngles()
MLE()
convertCart()	

	
