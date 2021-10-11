import numpy as np

def Abs6(signal, time, energy):
    T = time[-1:][0]
    kappa = -np.log(1e-6)/T
    mask = np.exp(-kappa*time)
    signal = np.multiply(signal - signal[0],mask)
    sig = np.zeros(np.shape(energy)).astype(np.complex128)
    for Ei, E in enumerate(energy):
        sig[Ei] = np.sum(np.multiply(np.exp(1j*E*time),signal))
    sig *= (time[1]-time[0])
    return (-1./137.)*4*np.pi*(1/0.0001)*np.multiply(np.imag(sig),energy)

def Abs5(signal, time, energy):
    T = time[-1:][0]
    mask = 1 - 3*(time/T)**2 + 2*(time/T)**3
    signal = np.multiply(signal - signal[0],mask)
    sig = np.zeros(np.shape(energy)).astype(np.complex128)
    for Ei, E in enumerate(energy):
        sig[Ei] = np.sum(np.multiply(np.exp(1j*E*time),signal))
    sig *= (time[1]-time[0])
    return (-1./137.)*4*np.pi*(1/0.0001)*np.multiply(np.imag(sig),energy)

data = np.loadtxt('./output_icwf/dipole.dat')

energy=np.arange(0,1,0.0005)

sig = Abs6(data[:,1],data[:,0],energy)

np.savetxt('Abs.txt',np.array([energy, sig]).T)
