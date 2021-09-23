"""
Methods for fitting results of Monte Carlo simulations.
"""

import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import newton
from scipy.optimize import curve_fit

class BlochFit:
    """Bloch and Multi-Bloch function as a weighted sum over N Bloch modes, F(T) = w1 * (1 - (T/Tc1)^b1) + w2 * (1 - (T/Tc2)^b2) + ..."""
    
    #number of Bloch modes to fit
    num_Bloch_modes = 1
    
    #bounds on Tc for fitting
    Tc_min = 1.0
    Tc_max = 1000.0
    
    #bounds on beta value for fitting
    b_min = 0.2
    b_max = 0.8
    
    #bounds on weights for fitting
    w_min = 0.0
    w_max = 1.0
    
    #bounds on exponent fitting (power of Bloch law)
    p_min = -20.0
    p_max = 20.0
    
    def __init__(self, 
                 Tc_min = 1.0, Tc_max = 1000.0, 
                 b_min = 0.2, b_max = 0.8,
                 p_min = -20.0, p_max = 20.0,
                 num_Bloch_modes = 1):
    
        self.Tc_min = Tc_min
        self.Tc_max = Tc_max
        
        self.b_min = b_min
        self.b_max = b_max
        
        self.p_min = p_min
        self.p_max = p_max
        
        self.num_Bloch_modes = num_Bloch_modes

    ######################################################################################################
    # Multi-Bloch fit
        
    #Multi-Bloch function as a weighted sum over N Bloch modes, F(T) = w1 * (1 - (T/Tc1)^b1) + w2 * (1 - (T/Tc2)^b2) + ...
    #Sum includes the number of Bloch modes specified here (N = num_Bloch_modes):
    #Note w1 + w2 + ... + wN = 1
    def __mBloch(self, params, T):
    
        value = 0.0
        for idx in range(self.num_Bloch_modes):
            #params is a list as [Tc1, beta1, weight1, Tc2, beta2, weight2, ..., TcN, betaN, weightN]
            Tc = params[idx * 3]
            beta = params[idx * 3 + 1]
            weight = params[idx * 3 + 2]
            value += weight * np.where(T < Tc, np.power(1.0 - T / Tc, beta), 0.0)
            
        return value            
    
    def fit_mBloch(self, TRange, mRange):
        """fit using configured multi-Bloch fitting"""
        
        #Least squares fitting with set bounds; the lambda specifies the residuals function
        bounds_min, bounds_max = [], []
        for idx in range(self.num_Bloch_modes):
            bounds_min.append(self.Tc_min)
            bounds_min.append(self.b_min)
            bounds_min.append(self.w_min)
            bounds_max.append(self.Tc_max)
            bounds_max.append(self.b_max)
            bounds_max.append(self.w_max)
        res = least_squares(lambda params: self.__mBloch(params, TRange) - mRange, bounds_max, bounds = (bounds_min, bounds_max))
        
        #find average Tc value as a weighted average
        avTc, Tcs, weights, betas = 0.0, [], [], []
        for idx in range(self.num_Bloch_modes): 
            Tcs.append(res.x[idx * 3])
            weights.append(res.x[idx * 3 + 2])
            betas.append(res.x[idx * 3 + 1])
            avTc += res.x[idx * 3 + 2] * res.x[idx * 3]
        
        return avTc, Tcs, betas, weights
    
    def mBloch_fitted(self, Tcs, beta, weights, TRange):
        """get the multi-Bloch fitted function with this"""
        
        params = []
        for idx in range(len(Tcs)):
            params.append(Tcs[idx])
            params.append(beta[idx])
            params.append(weights[idx])
            
        return self.__mBloch(params, np.asarray(TRange))
    
    def fit_mBloch_range(self, TRange, mRange, num_Bloch_modes_range):        
        """find average Tc as a function of number of Bloch modes from 1 to max_num_Bloch_modes in order to check for convergence"""
        
        Tc_vs_modes = []        
        for num_modes in range(num_Bloch_modes_range[0], num_Bloch_modes_range[1] + 1):
            print('Fitting with %d modes.' % num_modes)
            self.set_num_Bloch_modes(num_modes)
            avTc, Tcs, betas, weights = self.fit_mBloch(TRange, mRange)
            Tc_vs_modes.append(avTc)
            
        Tc = np.sum(Tc_vs_modes) / len(Tc_vs_modes)
        Tc_std = np.std(Tc_vs_modes)
            
        return Tc_vs_modes, Tc, Tc_std
    
    ######################################################################################################
    # Simple Bloch fit
    
    def __Bloch(self, T, Tc, beta):
        return np.where(T < Tc, np.power(1.0 - T / Tc, beta), 0.0)
    
    def fit_Bloch(self, TRange, mRange, mRange_std = []):
        """fit using configured Bloch fitting"""
        
        if len(mRange_std):
            for idx in range(len(mRange_std)): 
                if mRange_std[idx] < 1e-5: mRange_std[idx] = 1e-5
            popt, pcov = curve_fit(
                    self.__Bloch, 
                    np.asarray(TRange), np.asarray(mRange), sigma = np.asarray(mRange_std),
                    bounds = ((self.Tc_min, self.b_min), (self.Tc_max, self.b_max)))
        else:
            popt, pcov = curve_fit(self.__Bloch, np.asarray(TRange), np.asarray(mRange), bounds = ((self.Tc_min, self.b_min), (self.Tc_max, self.b_max)))
            
        perr = np.sqrt(np.diag(pcov))
        Tc = popt[0]
        Tc_std = perr[0]
        beta = popt[1]
        beta_std = perr[1]
        
        return Tc, beta, Tc_std, beta_std
    
    def Bloch_fitted(self, Tc, beta, TRange):
        """get the Bloch fitted function with this"""
        return self.__Bloch(TRange, Tc, beta)

    ######################################################################################################
    # Power of Bloch fit with fixed Tc and beta
    
    def __Bloch_power(self, T, Tc, beta, p):
        return np.where(T < Tc, np.power(1.0 - T / Tc, beta * p), 0.0)
    
    def fit_Bloch_power(self, TRange, kRange, Tc, beta, kRange_std = []):
        """fit using power of Bloch law with fixed Tc and beta"""
        
        if len(kRange_std):
            for idx in range(len(kRange_std)): 
                if kRange_std[idx] < 1e-5: kRange_std[idx] = 1e-5        
                popt, pcov = curve_fit(
                        lambda TRange, p: self.__Bloch_power(np.asarray(TRange), Tc, beta, p), 
                        np.asarray(TRange), np.asarray(kRange), sigma = np.asarray(kRange_std),
                        bounds = ((self.p_min,), (self.p_max,)))
        else:
            popt, pcov = curve_fit(
                    lambda TRange, p: self.__Bloch_power(np.asarray(TRange), Tc, beta, p), 
                    np.asarray(TRange), np.asarray(kRange), 
                    bounds = ((self.p_min,), (self.p_max,)))
        
        perr = np.sqrt(np.diag(pcov))
        return popt[0], perr[0]
    
    def Bloch_power_fitted(self, Tc, beta, p, TRange):
        """get the Bloch power fitted function with this"""
        return np.power(1.0 - np.asarray(TRange) / Tc, beta * p)
    
    ######################################################################################################
    # m scaling exponent fit
    
    def __m_power(self, TRange, mRange, p):
        
        return mRange ** p
    
    def fit_m_power(self, TRange, kRange, mRange, kRange_std = []):
        """fit k for exponent of m scaling"""
        
        if len(kRange_std):
            for idx in range(len(kRange_std)): 
                if kRange_std[idx] < 1e-5: kRange_std[idx] = 1e-5
            popt, pcov = curve_fit(
                    lambda TRange, p: self.__m_power(np.asarray(TRange), np.asarray(mRange), p), 
                    np.asarray(TRange), np.asarray(kRange), sigma = np.asarray(kRange_std),
                    bounds = ((self.p_min,), (self.p_max,)))
        else:
            popt, pcov = curve_fit(
                    lambda TRange, p: self.__m_power(np.asarray(TRange), np.asarray(mRange), p), 
                    np.asarray(TRange), np.asarray(kRange), 
                    bounds = ((self.p_min,), (self.p_max,)))
        
        perr = np.sqrt(np.diag(pcov))
        return popt[0], perr[0]
    
    def m_power_fitted(self, mRange, p):
        """get m scaling law to power p"""
        return np.asarray(mRange) ** p            
            
######################################################################################################
######################################################################################################
######################################################################################################
        
class CurieWeissFit:
    
    #bounds on Tc for fitting
    Tc_min = 1.0
    Tc_max = 1000.0
    
    #Boltzmann constant
    kB = 1.38064852e-23
    #Bohr magneton
    muB = 9.27400968e-24
    #free space permeability
    mu0 = 4 * np.pi * 1e-7
    
    #atomistic moment in units of muB
    mu_s = 1.0
    #external field
    Hext = 0.0
    
    def __init__(self, 
                 mu_s = 1.0, Hext = 0.0, 
                 Tc_min = 1.0, Tc_max = 1000.0):
        
        self.mu_s = mu_s
        self.Hext = Hext
        
        self.Tc_min = Tc_min
        self.Tc_max = Tc_max
    
    #solve B(x * c) - x = 0 for x
    #e.g. c = 3*Tc/T
    def __F(self, x, c, d):
        
        return np.where(c*x + d != 0.0, 1.0 / np.tanh(c*x + d) - 1.0 / (c*x + d) - x, 0.0)
    
    #Curie-Weiss law for given Tc (TCurie) and temperature range (TRange)
    #me = B(me * 3*Tc/T) for zero applied field
    #B(x) = coth(x) - 1/x
    def __Curie_Weiss(self, TRange, Tc):
        
        data, root = [], 1.0
        for Temperature in TRange:
            if Temperature > 0:
                c = 3*Tc / Temperature
                d = self.mu_s * self.muB * self.mu0 * self.Hext / (self.kB * Temperature)
                try:
                    if root != 0.0:
                        root = newton(self.__F, root, args = (c,d), maxiter = 10000, rtol = 1e-6)
                except:
                    root = 0.0
                    
                data.append(root)
            else:
                data.append(1.0)
            
        return data
    
    #Find Curie temperature by fitting Curie_Weiss to computed me vs T data.
    #Fit using given number of steps
    def Curie_Weiss_Fit(self, TRange, mRange, mRange_std = []):
        """fit using configured Curie-Weiss fitting"""
        
        if len(mRange_std):
            for idx in range(len(mRange_std)): 
                if mRange_std[idx] < 1e-5: mRange_std[idx] = 1e-5
            popt, pcov = curve_fit(
                    self.__Curie_Weiss, 
                    np.asarray(TRange), np.asarray(mRange), sigma = np.asarray(mRange_std),
                    bounds = ((self.Tc_min), (self.Tc_max)))
        else:
            popt, pcov = curve_fit(self.__Curie_Weiss, np.asarray(TRange), np.asarray(mRange), bounds = ((self.Tc_min), (self.Tc_max)))
        
        perr = np.sqrt(np.diag(pcov))
        Tc = popt[0]
        Tc_std = perr[0]
        
        return Tc, Tc_std
    
    def Curie_Weiss_fitted(self, TCurie, TRange):
        """get the Curie-Weiss fitted function with this"""
        return self.__Curie_Weiss(np.asarray(TRange), TCurie)
    
    def susrel(self, TRange, meRange, Tc, mu_s = 0.0):
        """
        compute relative longitudinal susceptibility based on Curie-Weiss law
        takes a temperature range (TRange), corresponding me values (meRange), atomic moment (mu_s) and Tc value
        """
        
        if mu_s == 0.0: mu_s = self.mu_s
        
        susrelRange = [0.0] * len(TRange)
        
        for idx in range(len(TRange)):
            
            Temperature = TRange[idx]
            if Temperature > 0.0:
                x = meRange[idx] * 3 * Tc / Temperature
                Bdash = -1 / np.sinh(x)**2 + 1 / x**2
                susrelRange[idx] = (mu_s * self.muB / (self.kB * Temperature)) * Bdash / (1 - Bdash * 3 * Tc / Temperature)
            
        return susrelRange
    
    def susrel_tran_uni(self, TRange, meRange, Ms0, K0, Tc):
        """
        compute relative transverse susceptibility for uniaxial anisotropy
        takes a temperature range (TRange), corresponding me values (meRange), Ms0 (A/m) and K0 (J/m^3) values (at 0K) and Tc value
        """
        
        Tc_idx = len(TRange)
        for idx in range(len(TRange)):
            if TRange[idx] >= Tc:
                Tc_idx = idx - 1
                break
        
        return TRange[:Tc_idx], [(Ms0 / (2*K0)) / me for me in meRange[:Tc_idx]]
    
    def Curie_Weiss(self, Temperature, Tc, mu_s = 0.0):
        """
        Curie-Weiss single temperature me value : me = B(me * 3*Tc/T)
        """
        
        if mu_s == 0.0: mu_s = self.mu_s
        
        me = 1.0
        
        if Temperature > 0:
            c = 3*Tc / Temperature
            d = mu_s * self.muB * self.mu0 * self.Hext / (self.kB * Temperature)
            try:
                me = newton(self.__F, me, args = (c,d), maxiter = 10000, rtol = 1e-6)
            except:
                me = 0.0
                
        return me
    
    def Curie_Weiss_susrel(self, Temperature, Tc, mu_s = 0.0):
        """
        Single temperature value compute relative longitudinal susceptibility based on Curie-Weiss law
        """
        
        if mu_s == 0.0: mu_s = self.mu_s
        
        me = self.Curie_Weiss(Temperature, Tc, mu_s)

        susrel_value = 0.0

        if Temperature > 0.0:
            x = me * 3 * Tc / Temperature
            Bdash = -1 / np.sinh(x)**2 + 1 / x**2
            susrel_value = (mu_s * self.muB / (self.kB * Temperature)) * Bdash / (1 - Bdash * 3 * Tc / Temperature)
            
        return susrel_value

######################################################################################################
######################################################################################################
######################################################################################################
        
    