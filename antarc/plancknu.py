import numpy as np

  

class Blackbody():
    # function f = plancknu(nu_icm,T);
    #
    # spectral Planck function as function of wavenumbers (cm-1)
    #
    # [h]    = J*s
    # [c]    = m/s
    # [cbar] = cm/s
    # [k]    = J*K-1
    #
    #
    #    Note: LBLRTM uses ...
    #c    Constants from NIST 01/11/2002
    #
    #      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
    #     *     CLIGHT / 2.99792458E+10 /, 
    # 
    #h    = 6.62606896e-34				# J s;  CODATA 2006
    #c    = 2.99792458e8				# m/s;  NIST
    #k    = 1.3806504e-23				# J K-1; CODATA 2006
    #cbar = 2.99792458e10			    # cm/s
    #nu = nu_icm+0.
    #top = h * cbar**3 * 2 * nu**3
    #bottom = c**2 *  ( np.exp(h*cbar*nu/(k*T))-1 )
    #f = cbar * top/bottom

    def __init__(self, nu_icm):
        h_2_cbar4_over_c2 = 1.1910427584934559e-08
        h_cbar_over_k = 1.4387751601679204
        
        self.top = h_2_cbar4_over_c2 * nu_icm*nu_icm*nu_icm
        self.exp_arg = h_cbar_over_k * nu_icm
        self.Nnu = len(nu_icm)


    def set_B(self, T):
        self.T = T
        self.B = self.top / ( np.exp(self.exp_arg/T) - 1 )
    
        #[B]= cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2
        #[B]= W m-2 cm

    def get_B(self, T):
        self.set_B(T)
        return self.B
    


def plancknu(nu_icm,T):

    # function f = plancknu(nu_icm,T);
    #
    # spectral Planck function as function of wavenumbers (cm-1)
    #
    # [h]    = J*s
    # [c]    = m/s
    # [cbar] = cm/s
    # [k]    = J*K-1
    #
    #
    #    Note: LBLRTM uses ...
    #c    Constants from NIST 01/11/2002
    #
    #      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
    #     *     CLIGHT / 2.99792458E+10 /, 
    # 
    #h    = 6.62606896e-34				# J s;  CODATA 2006
    #c    = 2.99792458e8				# m/s;  NIST
    #k    = 1.3806504e-23				# J K-1; CODATA 2006
    #cbar = 2.99792458e10			    # cm/s
    h_2_cbar4_over_c2 = 1.1910427584934559e-08
    h_cbar_over_k = 1.4387751601679204
    bottom = np.exp(h_cbar_over_k*nu_icm/T) - 1
    f = h_2_cbar4_over_c2 * nu_icm**3 / bottom

    #[B]= cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2
    #[B]= W m-2 cm

    return f

	

def plancknu_Ts(nu_icm, Ts):
    # function f = plancknu(nu_icm,T);
    #
    # spectral Planck function as function of wavenumbers (cm-1)
    #
    # [h]    = J*s
    # [c]    = m/s
    # [cbar] = cm/s
    # [k]    = J*K-1
    #
    #
    #    Note: LBLRTM uses ...
    #c    Constants from NIST 01/11/2002
    #
    #      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
    #     *     CLIGHT / 2.99792458E+10 /, 
    # 
    #h    = 6.62606896e-34				# J s;  CODATA 2006
    #c    = 2.99792458e8				# m/s;  NIST
    #k    = 1.3806504e-23				# J K-1; CODATA 2006
    #cbar = 2.99792458e10			    # cm/s
    h_2_cbar4_over_c2 = 1.1910427584934559e-08
    h_cbar_over_k = 1.4387751601679204
    
    top = h_2_cbar4_over_c2 * nu_icm**3
    exp_arg = h_cbar_over_k * nu_icm
    Nnu = len(nu_icm)
    NT = len(Ts)
    f = np.zeros((Nnu,NT))
    for i in range(NT):
        bottom = np.exp(exp_arg/Ts[i]) - 1
        f[:,i] = top / bottom

    #[B]= cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2
    #[B]= W m-2 cm

    return f


def plancknu_slow(nu_icm, T):
    h    = 6.62606896e-34				# J s;  CODATA 2006
    c    = 2.99792458e8				# m/s;  NIST
    k    = 1.3806504e-23				# J K-1; CODATA 2006
    cbar = 100*c			    # cm/s

    nu = nu_icm+0.
    top = h * cbar**3 * 2 * nu**3
    bottom = c**2 *  ( np.exp(h*cbar*nu/(k*T))-1 )
    f = cbar * top/bottom
    return f
