import numpy as np
import pandas as pd

#Operating conditions
T = 300.0 #K
f_H2S = 1.0 #bar
f_H2 = 1.0 #bar

#Settings
excel_path = 'mos2_data.xlsx'
n_S_100 = 6 #number of S atoms at 100% S-coverage
n_H_100 = 3 #number of H atoms at 100% H-coverage

#=====================================================================
#=====================================================================

#Fundamental constants
k_B = 8.6173303e-5 #Boltzmann constant in eV/K
h = 6.626070040e-34 #planck constant in J*s
e = 1.6021766208e-19 #elementary charge in C
c = 2.99792458e8 #speed of light in m/s
kJmoltoeV = 0.01036427 #kJ/mol to eV conversion

def read_excel_data():
#Read in data from Excel sheet
	E_Moedge = pd.read_excel(excel_path,'E_Mo_edge',header=1, index_col=0).values
	E_Sedge = pd.read_excel(excel_path,'E_S_edge',header=1, index_col=0).values
	S_vib_Moedge = pd.read_excel(excel_path,'S_vibrations_Mo_edge',
		header=None,index_col=0).values
	S_vib_Sedge = pd.read_excel(excel_path,'S_vibrations_S_edge',
		header=None,index_col=0).values
	ref_species = pd.read_excel(excel_path,'Reference_species')
	return E_Moedge, E_Sedge, S_vib_Moedge, S_vib_Sedge, ref_species

def freq_to_energy(freqs):
#Converts frequencies to energies
	vib_energies = [(100*h*c/e)*freq for freq in freqs]
	return vib_energies

def get_mu_ref(species,E_DFT,T):
#Calculates mu^0 from tabulated data
	tabulated_data_sheet = 'Tabulated_data'
	if T < 0:
		ValueError('Absolute temperature in K must be used')
	if species == 'H2':
		if T <= 1000:
			shomate_params = pd.read_excel(excel_path,
				tabulated_data_sheet).values[0:8]
		elif T > 1000 and T <= 2500:
			shomate_params = pd.read_excel(excel_path,
				tabulated_data_sheet).values[10:18]
		elif T > 2500:
			shomate_params = pd.read_excel(excel_path,
				tabulated_data_sheet).values[20:28]
		h0 = pd.read_excel(excel_path,
			tabulated_data_sheet).values[50][0]
		E_ZPE = pd.read_excel(excel_path,
			tabulated_data_sheet).values[56][0]
	elif species == 'H2S':
		if T <= 1400:
			shomate_params = pd.read_excel(excel_path,
				tabulated_data_sheet).values[30:38]
		elif T > 1400:
			shomate_params = pd.read_excel(excel_path,
				tabulated_data_sheet).values[40:48]
		h0 = pd.read_excel(excel_path,
			tabulated_data_sheet).values[53][0]
		E_ZPE = pd.read_excel(excel_path,
			tabulated_data_sheet).values[59][0]
	else:
		ValueError('Wrong species string')
	A = shomate_params[0][1]
	B = shomate_params[1][1]
	C = shomate_params[2][1]
	D = shomate_params[3][1]
	E = shomate_params[4][1]
	F = shomate_params[5][1]
	G = shomate_params[6][1]
	H = shomate_params[7][1]
	t = T/1000.0
	h = (A*t+B*t**2/2+C*t**3/3+D*t**4/4-E/t+F-H)*kJmoltoeV
	dh = h - h0
	s = (A*np.log(t)+B*t+C*t**2/2+D*t**3/3-E/(2*t**2)+G)*kJmoltoeV/1000.0
	mu_0 = dh - T*s + E_ZPE + E_DFT	
	return mu_0

def get_zpe(vib_energies):
#Calculate ZPE
	ZPE = 0.5*np.sum(vib_energies)
	return ZPE

def get_internal_energy(vib_energies,E_DFT):
#Calculate U
	ZPE = get_zpe(vib_energies)
	vib_sum = 0
	for vib_energy in vib_energies:
		vib_sum += vib_energy/(np.exp(vib_energy/(k_B*T))-1.0)
	U = E_DFT + ZPE + vib_sum
	return U

def get_entropy(vib_energies):
#Calculate S
	vib_sum = 0
	for vib_energy in vib_energies:
		vib_sum += (vib_energy/(k_B*T*(np.exp(vib_energy/(k_B*T))-1.0))
			- np.log(1.0-np.exp(-vib_energy/(k_B*T))))
	S = k_B*vib_sum
	return S

def get_helmholtz(vib_energies,E_DFT):
#Calculate F (if E_DFT = 0, this is F^vib)
	U = get_internal_energy(vib_energies,E_DFT)
	S = get_entropy(vib_energies)
	F = U-T*S
	return F

def freq_to_F_vib_H(H_vib):
#Gets vibrational frequencies of H*
	vib_energies = freq_to_energy(H_vib)
	F_vib_H = get_helmholtz(vib_energies,0)
	return F_vib_H

def freq_to_F_vib_S(S_vib,n_S):
#Gets vibrational frequencies of S*
	F_vib_S = np.zeros(len(n_S))
	for i in range(len(n_S)-1):
		row_vals_raw = S_vib[i,:]
		row_vals = row_vals_raw[~np.isnan(row_vals_raw)]
		vib_freqs = list(filter(None,row_vals))
		vib_energies = freq_to_energy(vib_freqs)
		F_vib = get_helmholtz(vib_energies,0)
		F_vib_S[i+1] = F_vib
	return F_vib_S

def run_thermo():
#Run the thermodynamic analysis
	E_Moedge, E_Sedge, S_vib_Moedge, S_vib_Sedge, ref_species = read_excel_data()
	H_vibrations_sheet = 'H_vibrations'
	vib_H_Mo_bridge = pd.read_excel(excel_path,
		H_vibrations_sheet).values[:,0]
	vib_H_Mo_inter = pd.read_excel(excel_path,
		H_vibrations_sheet).values[:,1]
	vib_H_S_basal = pd.read_excel(excel_path,
		H_vibrations_sheet).values[:,2]
	vib_H_S_edge = pd.read_excel(excel_path,
		H_vibrations_sheet).values[:,3]

	#Calculate F^vib terms
	n_S_Moedge = np.arange(0,np.max(np.shape(E_Moedge)[0]))
	n_S_Sedge = np.arange(0,np.max(np.shape(E_Sedge)[0]))
	n_H_Moedge = np.arange(0,np.max(np.shape(E_Moedge)[1]))
	n_H_Sedge = np.arange(0,np.max(np.shape(E_Sedge)[1]))

	F_vib_H_Mo_bridge = freq_to_F_vib_H(vib_H_Mo_bridge)
	F_vib_H_Mo_inter = freq_to_F_vib_H(vib_H_Mo_inter)
	F_vib_H_S_basal = freq_to_F_vib_H(vib_H_S_basal)
	F_vib_H_S_edge = freq_to_F_vib_H(vib_H_S_edge)

	F_vib_S_Moedge = freq_to_F_vib_S(S_vib_Moedge,n_S_Moedge)
	F_vib_S_Sedge = freq_to_F_vib_S(S_vib_Sedge,n_S_Sedge)

	#Calculate F_MoSxHy
	F_vib_H_Moedge = np.zeros(np.shape(E_Moedge))
	F_vib_H_Sedge = np.zeros(np.shape(E_Sedge))

	F_vib_H_Moedge[0,:] = F_vib_H_Mo_bridge*n_H_Moedge

	F_vib_H_Moedge[1,1] = F_vib_H_Mo_bridge
	F_vib_H_Moedge[1,2] = F_vib_H_Mo_bridge*2
	F_vib_H_Moedge[1,3] = F_vib_H_Mo_bridge*2+F_vib_H_S_edge

	F_vib_H_Moedge[2,1] = F_vib_H_Mo_bridge
	F_vib_H_Moedge[2,2] = F_vib_H_Mo_bridge+F_vib_H_S_edge
	F_vib_H_Moedge[2,3] = F_vib_H_Mo_bridge+F_vib_H_S_edge*2

	F_vib_H_Moedge[3,:] = F_vib_H_S_edge*n_H_Moedge
	F_vib_H_Moedge[4,:] = F_vib_H_S_edge*n_H_Moedge
	F_vib_H_Moedge[5,:] = F_vib_H_S_edge*n_H_Moedge
	F_vib_H_Moedge[6,:] = F_vib_H_Mo_bridge*n_H_Moedge

	F_vib_H_Sedge[0,:] = F_vib_H_Mo_bridge*n_H_Sedge
	
	F_vib_H_Sedge[1,1] = F_vib_H_Mo_bridge
	F_vib_H_Sedge[1,2] = F_vib_H_Mo_bridge*2
	F_vib_H_Sedge[1,3] = F_vib_H_Mo_bridge*2+F_vib_H_S_edge

	F_vib_H_Sedge[2,1] = F_vib_H_Mo_bridge
	F_vib_H_Sedge[2,2] = F_vib_H_Mo_bridge + F_vib_H_Mo_inter
	F_vib_H_Sedge[2,3] = F_vib_H_Mo_bridge + F_vib_H_Mo_inter*2

	F_vib_H_Sedge[3,:] = F_vib_H_Mo_inter*n_H_Sedge

	F_vib_H_Sedge[4,1] = F_vib_H_S_basal
	F_vib_H_Sedge[4,2] = F_vib_H_S_basal+F_vib_H_Mo_inter
	F_vib_H_Sedge[4,3] = F_vib_H_S_basal+F_vib_H_Mo_inter*2

	F_vib_H_Sedge[5,1] = F_vib_H_Mo_inter
	F_vib_H_Sedge[5,2] = F_vib_H_Mo_inter+F_vib_H_S_basal
	F_vib_H_Sedge[5,3] = F_vib_H_Mo_inter+F_vib_H_S_basal*2

	F_vib_H_Sedge[6,:] = F_vib_H_S_basal*n_H_Sedge
	
	F_Moedge = E_Moedge+F_vib_S_Moedge[:,np.newaxis]+F_vib_H_Moedge
	F_Sedge = E_Sedge+F_vib_S_Sedge[:,np.newaxis]+F_vib_H_Sedge
	F_Moedge = F_Moedge - F_Moedge[0,0]
	F_Sedge = F_Sedge - F_Sedge[-1,0]

	#Calculate chemical potentials
	if  T < 298 or T > 6000:
		print('Warning: Shomate fit parameters for H2 and H2S were not '
			+'fit for T < 298 K or T > 6000 K')
	E_H2 = ref_species.loc[0][1]
	E_H2S = ref_species.loc[1][1]
	E_S_bulk = ref_species.loc[2][1]
	E_Mo_bulk = ref_species.loc[3][1]
	E_MoS2_bulk = ref_species.loc[4][1]
	mu_H2_0 = get_mu_ref('H2',E_H2,T)
	mu_H2S_0 = get_mu_ref('H2S',E_H2S,T)
	mu_S = np.log(f_H2S/f_H2)*(k_B*T) - mu_H2_0 + mu_H2S_0
	mu_H = 0.5*(np.log(f_H2)*(k_B*T) + mu_H2_0)
	delta_mu_S = mu_S - (E_H2S - E_H2)
	delta_mu_H = mu_H - 0.5*E_H2
	Gf_MoS2 = E_MoS2_bulk - E_Mo_bulk - 2*E_S_bulk
	ref_diff = E_S_bulk - (E_H2S - E_H2)
	if (delta_mu_S < Gf_MoS2 + ref_diff) or (delta_mu_S > ref_diff):
		print('Note: MoS2 is potentially unstable for this delta mu_S')
	if delta_mu_H > 0:
		print('Note: MoS2 is potentially unstable for this delta mu_H')

	#Calculate grand potentials
	phi_Moedge = F_Moedge - n_S_Moedge[:,np.newaxis]*mu_S - n_H_Moedge*mu_H
	phi_Moedge = phi_Moedge - phi_Moedge[0,0]
	phi_Sedge = F_Sedge - n_S_Sedge[:,np.newaxis]*mu_S - n_H_Sedge*mu_H
	phi_Sedge = phi_Sedge - phi_Sedge[-1,0]

	#Determine most stable coverages and save results
	rc_Moedge = np.unravel_index(np.nanargmin(phi_Moedge),phi_Moedge.shape)
	rc_Sedge = np.unravel_index(np.nanargmin(phi_Sedge),phi_Sedge.shape)
	print('T = '+str(T)+' K')
	print('f_HS/f^0 = '+str(f_H2))
	print('f_H2S/f_H2 = '+str(f_H2S/f_H2))
	print('delta mu_S = '+str(delta_mu_S))
	print('delta mu_H = '+str(delta_mu_H))
	print('Most stable coverage on Mo-edge: theta_S = '
		+str(np.round(n_S_Moedge[rc_Moedge[0]]/n_S_100,2))+', theta_H = '
		+str(np.round(n_H_Moedge[rc_Moedge[1]]/n_H_100,2)))
	print('Most stable coverage on S-edge: theta_S = '
		+str(np.round(n_S_Sedge[rc_Sedge[0]]/n_S_100,2))
		+', theta_H = '+str(np.round(n_H_Sedge[rc_Sedge[1]]/n_H_100,2)))
	return F_Moedge, F_Sedge, phi_Moedge, phi_Sedge, delta_mu_S, delta_mu_H

def write_results(F_Moedge, F_Sedge, phi_Moedge, phi_Sedge):
#Write results to csv file
	np.savetxt('F_Mo_edge.csv',F_Moedge,delimiter=',')
	np.savetxt('F_S_edge.csv',F_Sedge,delimiter=',')
	np.savetxt('phi_Mo_edge.csv',phi_Moedge,delimiter=',')
	np.savetxt('phi_S_edge.csv',phi_Sedge,delimiter=',')

F_Moedge, F_Sedge, phi_Moedge, phi_Sedge, delta_mu_S, delta_mu_H = run_thermo()
write_results(F_Moedge, F_Sedge, phi_Moedge, phi_Sedge)
