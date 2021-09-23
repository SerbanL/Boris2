"""
Methods for computing scaling laws using Monte Carlo algorithm in atomistic meshes.
Available pre-configured methods (import as needed):
    
simulate_m_scaling(ns, fileName, Tmax = 700, Hext = 0.0)
simulate_m_scaling_constrained(ns, fileName, Tmax = 700, Hext = 0.0, direction = [90, 0])
simulate_susceptibilities(ns, fileName, Tmax = 700, Hext = 0.0)
simulate_k_uniaxial_scaling(ns, fileName, Tmax = 700, hard_axis = [90, 90], easy_axis = [90, 0])
simulate_k_cubic_scaling(ns, fileName, Tmax = 700, hard_axis = [np.degrees(np.arctan(np.sqrt(2))), 45], easy_axis1 = [90, 0], easy_axis2 = [90, 90])
simulate_m_k_uniaxial_scaling(ns, fileName, Tmax = 700, hard_axis = [90, 90], easy_axis = [90, 0])
simulate_m_k_cubic_scaling(ns, fileName, Tmax = 700, hard_axis = [np.degrees(np.arctan(np.sqrt(2))), 45], easy_axis1 = [90, 0], easy_axis2 = [90, 90])
simulate_exchange_stiffness(ns, fileName, Tmax = 700, Hext = 0.0)
"""

from NetSocks import NSClient
from AtomisticFits import BlochFit
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

class MonteCarlo:
    
    ns = NSClient()
    
    #m stage setpoints
    m_stage = [1.0, 0.6, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.02, 0.005]
    
    #thermalization and averaging number of iterations to apply when below corresponding m setpoint
    iter_thrm = [1000, 2000, 2000, 2000, 10000, 10000, 4000, 3000, 3000, 3000]
    iter_avrg = [1000, 1000, 2000, 3000, 15000, 30000, 10000, 10000, 3000, 3000]
    iter_avrg_ref = [1000, 1000, 2000, 3000, 15000, 30000, 10000, 10000, 3000, 3000]
    
    #as above but for the constrained MC algorithm
    iter_thrm_c = [1000, 2000, 4000, 10000, 10000, 10000, 10000, 5000, 3000, 3000]
    iter_avrg_c = [1000, 1000, 4000, 10000, 20000, 30000, 20000, 10000, 3000, 3000]
    iter_avrg_c_ref = [1000, 1000, 4000, 10000, 20000, 30000, 20000, 10000, 3000, 3000]
    
    #as above but for a domain wall
    iter_thrm_dw = [2000, 2000, 2000, 2000, 10000, 10000, 4000, 3000, 3000, 3000]
    iter_avrg_dw = [12500, 12500, 12500, 12500, 15000, 15000, 10000, 10000, 3000, 3000]
    iter_avrg_dw_ref = [12500, 12500, 12500, 12500, 15000, 15000, 10000, 10000, 3000, 3000]
    
    #as above but for use with SLLG, not Monte Carlo
    iter_thrm_sllg = [10000, 10000, 10000, 10000, 20000, 30000, 10000, 10000, 5000, 5000]
    iter_avrg_sllg = [5000, 5000, 10000, 10000, 15000, 20000, 10000, 5000, 3000, 3000]
    iter_avrg_sllg_ref = [5000, 5000, 10000, 10000, 15000, 20000, 10000, 5000, 3000, 3000]
    use_sllg = False
    
     #thermalization and averaging number of iterations to apply when below corresponding m setpoint
    iter_thrm_RKKY_H = [1000, 2000, 2000, 2000, 10000, 10000, 4000, 3000, 3000, 3000]
    iter_avrg_RKKY_H = [12500, 12500, 12500, 12500, 15000, 15000, 10000, 10000, 3000, 3000]
    iter_avrg_RKKY_H_ref = [12500, 12500, 12500, 12500, 15000, 15000, 10000, 10000, 3000, 3000]
    
    #save every iter_avrg_spacing
    iter_avrg_spacing = 1
    
    #if simulating susceptibilities need to used a large fixed number of averages to capture data, otherwise the variance is not correct
    #in most cases the absolute bare minimum is ~25000, but more may be needed to reduce statistical uncertainty (noise), especially for the transverse susceptibility with low K constants
    fixed_avrg_iter = 50000
    
    #only fit for domain wall width for m above this value
    dw_width_m_setpoint = 0.1
    
    #temperature step to use when below corresponding m setpoint
    step_stage = [10.0, 10.0, 5.0, 4.0, 2.0, 1.0, 2.0, 4.0, 6.0, 10.0]
    
    #if an expected Tc value is passed in then adjust the temperature steps accordingly
    #the values in step_stage are set for the reference Tc below
    step_reference_Tc = 628.0
    
    #stage counter use as index in the above lists
    stage = 0
    
    #Boltzmann constant
    kB = 1.38064852e-23
    #Bohr magneton
    muB = 9.27400968e-24
    
    ##########################################################################
    
    def __init__(self, ns = NSClient(), expected_Tc = 628.0, fixed_average = False, averaging_level = 2, use_sllg = False):
        """
        Setup Monte Carlo simulation. If expected_Tc value passed in, the default temperature scheme will be adjusted accordingly.
        If running with fixed number of averaging iterations, use fixed_average = True: this is useful for susceptibility simulations where a large number of averages are required to produce a correct output.
        Can also adjust the number of averages used by setting averaging_level as:
        0 : Very fast (quarter the number of default averages : use for testing only)
        1 : Fast (half the number of averages : bare minimum, be careful with susceptibilities as these may be incorrect)
        2 : Default (accurate/high quality)
        3 : Very high quality (double number of default averages)
        """
        
        self.use_sllg = use_sllg
                
        self.ns = ns
        self.step_stage = [value * expected_Tc / self.step_reference_Tc for value in self.step_stage]
        
        divisor = 1.0
        if averaging_level == 0: divisor = 4.0
        elif averaging_level == 1: divisor = 2.0
        elif averaging_level == 2: divisor = 1.0
        elif averaging_level == 3: divisor = 0.5
        elif averaging_level == 4: divisor = 0.25
        elif averaging_level == 5: divisor = 0.125
        elif averaging_level == 6: divisor = 0.0625
        elif averaging_level == 7: divisor = 0.03125
            
        if fixed_average:
            self.iter_avrg = [self.fixed_avrg_iter / divisor] * len(self.m_stage)
            self.iter_avrg_c = [self.fixed_avrg_iter / divisor] * len(self.m_stage)
            self.iter_avrg_dw = [self.fixed_avrg_iter / divisor] * len(self.m_stage)
            self.iter_avrg_sllg = [self.fixed_avrg_iter / divisor] * len(self.m_stage)
            self.iter_avrg_RKKY_H = [self.fixed_avrg_iter / divisor] * len(self.m_stage)
        else:
            self.iter_avrg = [value / divisor for value in self.iter_avrg_ref]
            self.iter_avrg_c = [value / divisor for value in self.iter_avrg_c_ref]
            self.iter_avrg_dw = [value / divisor for value in self.iter_avrg_dw_ref]
            self.iter_avrg_sllg = [value / divisor for value in self.iter_avrg_sllg_ref]
            self.iter_avrg_RKKY_H = [value / divisor for value in self.iter_avrg_RKKY_H_ref]
         
    ##########################################################################
        
    #thermalization routine
    def __thermalize_stage(self, Ms0, Temperature, get_domainwall = False):
        
        if Temperature == 0.0 and not self.use_sllg: 
            self.ns.computefields()
            e_anis = self.ns.showdata('e_anis')
            e_exch = self.ns.showdata('e_exch')
            dw_width = self.ns.showdata('dwpos_x')[1]
            return 1.0, 0.0, 0.0, 0.0, 0.0, e_anis, e_exch, dw_width
        
        self.ns.reset()
        self.ns.Run()
        
        #get data saved in temporary.txt : Mx, My, Mz, e_anis e_exch dw_width
        if get_domainwall:
            self.ns.dp_load('temporary.txt', [0, 1, 2, 3, 4, 6, 0, 1, 2, 3, 4, 5])
        else:
            self.ns.dp_load('temporary.txt', [0, 1, 2, 3, 4, 0, 1, 2, 3, 4])
        
        self.ns.dp_div(0, Ms0)
        mxstd = self.ns.dp_mean(0)[1]
        self.ns.dp_pow(0, 2)
        
        self.ns.dp_div(1, Ms0)
        mystd = self.ns.dp_mean(1)[1]
        self.ns.dp_pow(1, 2)
        
        self.ns.dp_div(2, Ms0)
        mzstd = self.ns.dp_mean(2)[1]
        self.ns.dp_pow(2, 2)
        
        self.ns.dp_adddp(0, 1, 0)
        self.ns.dp_adddp(0, 2, 0)
        self.ns.dp_pow(0, 0.5)
        m = self.ns.dp_mean(0)
        
        e_anis = self.ns.dp_mean(3)[0]
        e_exch = self.ns.dp_mean(4)[0]
        
        dw_width = 0.0
        #for domain wall average exclude points which are too far from the mean (ratio of 0.9 is good - not ratio of 1 means include points from 0 to 2 * mean)
        #this is enough to exlude points unreasonably far away from mean, and importantly exclude points with 0 value: these are failed dw fits.
        if get_domainwall: dw_width = self.ns.dp_mean(5, 0.9)[0]
            
        return m[0], m[1], mxstd, mystd, mzstd, e_anis, e_exch, dw_width
    
    #thermalization routine for RKKY simulation
    def __thermalize_stage_RKKY(self, Ms0, Temperature):
        
        if Temperature == 0.0: 
            self.ns.computefields()
            e_anis = self.ns.showdata('e_anis')
            t_surfexch = self.ns.showdata('t_surfexch')
            mag_t_surfexch = -np.sqrt(t_surfexch[0]**2 + t_surfexch[1]**2 + t_surfexch[2]**2)
            return 1.0, e_anis, mag_t_surfexch, 1.0
        
        self.ns.reset()
        self.ns.Run()
        
        #get data saved in temporary.txt : Mx, My, Mz, Mx2, My2, Mz2, e_anis, e_anis2, t_surexch, t_surexch2
        self.ns.dp_load('temporary.txt', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
        
        self.ns.dp_div(0, Ms0)
        self.ns.dp_div(1, Ms0)
        self.ns.dp_div(2, Ms0)
        self.ns.dp_div(3, Ms0)
        self.ns.dp_div(4, Ms0)
        self.ns.dp_div(5, Ms0)
        
        mx = self.ns.dp_mean(0)[0]
        my = self.ns.dp_mean(1)[0]
        mz = self.ns.dp_mean(2)[0]
        mx2 = self.ns.dp_mean(3)[0]
        my2 = self.ns.dp_mean(4)[0]
        mz2 = self.ns.dp_mean(5)[0]        
        #average over the bilayer - used for H field method
        mav = np.sqrt((mx + mx2)**2 + (my + my2)**2 + (mz + mz2)**2) / 2
        
        #average m in each layer
        self.ns.dp_pow(0, 2)
        self.ns.dp_pow(1, 2)
        self.ns.dp_pow(2, 2)
        self.ns.dp_adddp(0, 1, 0)
        self.ns.dp_adddp(0, 2, 0)
        self.ns.dp_pow(0, 0.5)
        m = self.ns.dp_mean(0)
        
        self.ns.dp_pow(3, 2)
        self.ns.dp_pow(4, 2)
        self.ns.dp_pow(5, 2)
        self.ns.dp_adddp(3, 4, 3)
        self.ns.dp_adddp(3, 5, 3)
        self.ns.dp_pow(3, 0.5)
        m2 = self.ns.dp_mean(3)
            
        #average energies
        e_anis = self.ns.dp_mean(6)[0]
        e_anis2 = self.ns.dp_mean(7)[0]
        
        t_surfexch_x = self.ns.dp_mean(8)[0]
        t_surfexch_y = self.ns.dp_mean(9)[0]
        t_surfexch_z = self.ns.dp_mean(10)[0]
        
        t_surfexch2_x = self.ns.dp_mean(11)[0]
        t_surfexch2_y = self.ns.dp_mean(12)[0]
        t_surfexch2_z = self.ns.dp_mean(13)[0]
            
        mag_t_surfexch = np.sqrt(t_surfexch_x**2 + t_surfexch_y**2 + t_surfexch_z**2)
        mag_t_surfexch2 = np.sqrt(t_surfexch2_x**2 + t_surfexch2_y**2 + t_surfexch2_z**2)
        
        return (m[0] + m2[0]) / 2, (e_anis + e_anis2) / 2, -(mag_t_surfexch + mag_t_surfexch2) / 2, mav
        
    ##########################################################################

    #sweep temperature up
    def __sweep_temperature_up(self, 
                               #Temperature sweep range (min max if use_default_step, else the temperature values to use in list)
                               Tsweep_range = [0, 700], use_default_step = True,
                               #thermalization, and averaging number of iterations (and where to save result of averaging)
                               iter_thrm_range = [], iter_avrg_range = [], outputFile = '',
                               #if get_domainwall, then expecting to have domain wall data to average in last field of temporary.txt file (dwpos_x data)
                               #in this case only get domain wall width if Temperature <= Tmax_dw
                               Tmax_dw = 0.0, get_domainwall = False, do_rkky = False):
        
        #find 0 temperature magnetisation
        Mvec = self.ns.showdata('<M>')
        Ms0 = np.sqrt(Mvec[0]**2 + Mvec[1]**2 + Mvec[2]**2)
        
        #find atomic moment set (units of muB)
        mu_s = self.ns.setparam(self.ns.meshfocus(), 'mu_s')
        
        n = self.ns.showmcells()
        N = n[0] * n[1] * n[2]
        
        self.ns.dataprecision(12)
        
        Trange, mrange, Krange, susrelrange = [], [], [], []
        
        Temperature, m, temp_iter = Tsweep_range[0], 1.0, 0
        while Temperature <= Tsweep_range[-1]:
            
            self.ns.temperature(Temperature)
            
            #stop getting domain wall width if m is too low (for T > Tmax_dw)
            if get_domainwall and Temperature > Tmax_dw:
                get_domainwall = False
                #delete last data field : this will be the domain wall width (dwpos_x)
                self.ns.deldata(self.ns.data() - 1)
            
            if do_rkky:
                
                m, e_anis, e_surfexch, mav = self.__thermalize_stage_RKKY(Ms0, Temperature)
                
                if len(outputFile): 
                    self.ns.SaveDataToFile(outputFile, [Temperature, m, e_anis, e_surfexch])
            else:
            
                m, mstd, mxstd, mystd, mzstd, e_anis, e_exch, dw_width = self.__thermalize_stage(Ms0, Temperature, get_domainwall)
                if Temperature > 0.0: 
                    susrel_l = (mu_s * self.muB * N / (self.kB * Temperature)) * mstd**2
                    susrel_x = (mu_s * self.muB * N / (self.kB * Temperature)) * mxstd**2
                    susrel_y = (mu_s * self.muB * N / (self.kB * Temperature)) * mystd**2
                    susrel_z = (mu_s * self.muB * N / (self.kB * Temperature)) * mzstd**2
                else: 
                    susrel_l = 0.0; susrel_x = 0.0; susrel_y = 0.0; susrel_z = 0.0
            
                if len(outputFile): 
                    self.ns.SaveDataToFile(outputFile, [Temperature, m, mstd, susrel_l, susrel_x, susrel_y, susrel_z, e_anis, e_exch, dw_width])
                    
                susrelrange.append(susrel_l)
            
            Trange.append(Temperature)
            mrange.append(m)
            Krange.append(e_anis)
            
            #increment stage?
            if self.stage != len(self.m_stage) - 1 and m < self.m_stage[self.stage + 1]:
                #m point reached, so next stage : adjust number of thermalization and averaging iterations
                #first check exactly which stage we are at since it's possible we need to increment stage by more than 1
                for self.stage in range(len(self.m_stage) - 1): 
                    if m > self.m_stage[self.stage + 1]: break
                else: self.stage = len(self.m_stage) - 1

                self.ns.editstagestop(0, 'iter', iter_thrm_range[self.stage])
                self.ns.editstagestop(1, 'iter', iter_avrg_range[self.stage])
                
            #set value of Temperature for next step
            if use_default_step:   
                #use default stepping mechanism
                Temperature += self.step_stage[self.stage]   
            else:
                #not default stepping, so expecting Tsweep_range to contain exact temperature values to simulate
                temp_iter += 1
                if temp_iter < len(Tsweep_range): Temperature = Tsweep_range[temp_iter]
                else: break
                
        return Trange, mrange, susrelrange, Krange
    
    def __sweep_temperature_up_RKKY_withH(self, 
                               bMesh, tMesh,
                               Trange, mrange, krange,
                               #thermalization, and averaging number of iterations (and where to save result of averaging)
                               iter_thrm_range = [], iter_avrg_range = [], outputFile = ''):
        
        self.ns.dataprecision(12)
        
        a = self.ns.cellsize()[0]
        Js = self.ns.setparam(tMesh, 'Js')
        J1 = Js / a**2
        Ms0 = self.ns.showms()
        K0 = self.ns.showk()
        mu0 = 4*np.pi*1e-7
        
        self.ns.meshfocus(bMesh)
        bRect = self.ns.meshrect()
        self.ns.meshfocus(tMesh)
        tRect = self.ns.meshrect()
        
        self.ns.temperature(0)
        self.ns.setangle(90, 90, bMesh)
        self.ns.setangle(90, 270, tMesh)
        
        t1 = bRect[5] - bRect[2]
        t2 = tRect[5] - tRect[2]
        
        Hsat_K = 2 * K0 / (mu0 * Ms0)
        Hsat_RKKY = (-J1 / (mu0*Ms0)) * (t1 + t2) / (t1*t2)
        Hsat = Hsat_RKKY + Hsat_K
        H = Hsat * 0.75
        self.ns.setfield(H, 90, 0)
        
        self.ns.setode('LLGStatic', 'SDesc')
        self.ns.addstage('Relax')
        self.ns.editstagestop(2, 'mxh', 1e-4)
        self.ns.savedataflag(0)
        self.ns.Run()
        self.ns.delstage(2)
        self.ns.savedataflag(1)
        
        for (Temperature, msat, k) in zip(Trange, mrange, krange):
            
            self.ns.temperature(Temperature)
            H = Hsat * 0.75
            self.ns.setfield(H, 90, 0)
            
            m, e_anis, e_surfexch, mav = self.__thermalize_stage_RKKY(Ms0, Temperature)
            
            if Temperature != 0.0:
                #based on m and msat values, and currently set field, calculate actual saturation field
                Hsat = msat*H / mav
                Hsat_K = 2 * K0 * k / (mu0 * Ms0 * m)
                Hsat_RKKY = Hsat - Hsat_K
                #now solve for J1(T)
                J1 = -Hsat_RKKY * mu0*Ms0*msat*t1*t2 / (t1 + t2)
                
            if len(outputFile): 
                self.ns.SaveDataToFile(outputFile, [Temperature, m, 0.0, J1])
            
            #increment stage?
            if self.stage != len(self.m_stage) - 1 and msat < self.m_stage[self.stage + 1]:
                #m point reached, so next stage : adjust number of thermalization and averaging iterations
                #first check exactly which stage we are at since it's possible we need to increment stage by more than 1
                for self.stage in range(len(self.m_stage) - 1): 
                    if msat > self.m_stage[self.stage + 1]: break
                else: self.stage = len(self.m_stage) - 1

                self.ns.editstagestop(0, 'iter', iter_thrm_range[self.stage])
                self.ns.editstagestop(1, 'iter', iter_avrg_range[self.stage])
    
    ##########################################################################
    
    #simple classical (standard) MC sweep
    def __run_MC_classical(self,
                           #Temperature sweep range (min max if use_default_step, else the temperature values to use in list)
                           Tsweep_range = [0, 700], use_default_step = True,
                           #file to save results to
                           outputFile = '',
                           #applied field and direction (normally easy axis direction)
                           Hext = 0.0, direction_deg = [90, 0], 
                           #if get_energies is not enabled then energy terms are not computed (not cheap, so only enable if needed)
                           get_energies = False):
        
        self.stage = 0
        
        #starting state : uniform magnetization
        self.ns.setangle(direction_deg[0], direction_deg[1], self.ns.meshfocus())
        self.ns.setfield(Hext, direction_deg[0], direction_deg[1])
        
        #first simulation stage: thermalize, no data save
        if self.use_sllg:
            self.ns.setstage('Relax')
            self.ns.editstagestop(0, 'iter', self.iter_thrm_sllg[self.stage])
        else:
            self.ns.setstage('MonteCarlo')
            self.ns.editstagestop(0, 'iter', self.iter_thrm[self.stage])
        
        #second simulation stage: save enough points to average
        if self.use_sllg:
            self.ns.addstage('Relax')
            self.ns.editstagestop(1, 'iter', self.iter_avrg_sllg[self.stage])
        else:
            self.ns.addstage('MonteCarlo')
            self.ns.editstagestop(1, 'iter', self.iter_avrg[self.stage])
            
        self.ns.editdatasave(1, 'iter', self.iter_avrg_spacing)
        
        self.ns.setdata('<M>')
        self.ns.adddata('e_anis')
        self.ns.adddata('e_exch')
        self.ns.savedatafile('temporary.txt')
        
        self.ns.mcconstrain(0)
        self.ns.mcserial(0, self.ns.meshfocus())
        if get_energies: self.ns.mccomputefields(1)
        else: self.ns.mccomputefields(0)
        self.ns.cuda(1)
        
        if self.use_sllg:
            return self.__sweep_temperature_up(Tsweep_range, use_default_step, self.iter_thrm_sllg, self.iter_avrg_sllg, outputFile)
        else:
            return self.__sweep_temperature_up(Tsweep_range, use_default_step, self.iter_thrm, self.iter_avrg, outputFile)
        
    #simple classical (standard) MC sweep
    def __run_MC_classical_RKKY(self,
                           #bottom and top mesh names
                           bMesh, tMesh,
                           #Temperature sweep range (min max if use_default_step, else the temperature values to use in list)
                           Tsweep_range = [0, 700], use_default_step = True,
                           #file to save results to
                           outputFile = '',
                           use_H_method = False, mrange = [], krange = []):
        
        self.stage = 0
        
        self.ns.setstage('MonteCarlo')
        if use_H_method: self.ns.editstagestop(0, 'iter', self.iter_thrm_RKKY_H[self.stage])
        else: self.ns.editstagestop(0, 'iter', self.iter_thrm[self.stage])
        
        self.ns.addstage('MonteCarlo')
        if use_H_method: self.ns.editstagestop(1, 'iter', self.iter_avrg_RKKY_H[self.stage])
        else: self.ns.editstagestop(1, 'iter', self.iter_avrg[self.stage])
            
        self.ns.editdatasave(1, 'iter', self.iter_avrg_spacing)
        
        self.ns.setdata('<M>', bMesh)
        self.ns.adddata('<M>', tMesh)
        self.ns.adddata('e_anis', bMesh)
        self.ns.adddata('e_anis', tMesh)
        self.ns.adddata('t_surfexch', bMesh)
        self.ns.adddata('t_surfexch', tMesh)
        self.ns.savedatafile('temporary.txt')

        self.ns.cuda(1)
        
        if use_H_method:
            self.ns.mccomputefields(0)
            return self.__sweep_temperature_up_RKKY_withH(bMesh, tMesh, Tsweep_range, mrange, krange, self.iter_thrm_RKKY_H, self.iter_avrg_RKKY_H, outputFile)
        else:
            self.ns.mccomputefields(1)
            return self.__sweep_temperature_up(Tsweep_range, use_default_step, self.iter_thrm_c, self.iter_avrg_c, outputFile, do_rkky = True)
    
    #constrained MC sweep, with energy computation enabled
    def __run_MC_constrained(self,
                             #Temperature sweep range (min max if use_default_step, else the temperature values to use in list)
                             Tsweep_range = [0, 700], use_default_step = True,
                             #file to save results to
                             outputFile = '', 
                             #applied field and constraining direction
                             Hext = 0.0, direction_deg = [90, 0]):
        
        self.stage = 0
        
        #starting state : uniform magnetization
        self.ns.setangle(direction_deg[0], direction_deg[1], self.ns.meshfocus())
        self.ns.setfield(Hext, direction_deg[0], direction_deg[1])
        
        #first simulation stage: thermalize, no data save
        self.ns.setstage('MonteCarlo')
        self.ns.editstagestop(0, 'iter', self.iter_thrm_c[self.stage])
        
        #second simulation stage: save enough points to average
        self.ns.addstage('MonteCarlo')
        self.ns.editstagestop(1, 'iter', self.iter_avrg_c[self.stage])
        self.ns.editdatasave(1, 'iter', self.iter_avrg_spacing)
        
        self.ns.setdata('<M>')
        self.ns.adddata('e_anis')
        self.ns.adddata('e_exch')
        self.ns.savedatafile('temporary.txt')
        
        theta, phi = np.radians(direction_deg[0]), np.radians(direction_deg[1])
        self.ns.mcconstrain([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)], self.ns.meshfocus())
        self.ns.mcserial(0, self.ns.meshfocus())
        
        self.ns.mccomputefields(1)
        self.ns.cuda(1)
        
        return self.__sweep_temperature_up(Tsweep_range, use_default_step, self.iter_thrm_c, self.iter_avrg_c, outputFile)
    
    #MC sweep starting from a Bloch wall state, where we also obtain domain wall width if Temperature <= Tmax_dw
    def __run_MC_BlochWall(self,
                           #Temperature sweep range (min max if use_default_step, else the temperature values to use in list)
                           Tsweep_range = [0, 700], use_default_step = True,
                           #file to save results to
                           outputFile = '', 
                           #field applied along the z axis
                           Hext = 0.0, 
                           #maximum temperature for domain wall width computation
                           Tmax_dw = 0.0):
        
        meshrect = self.ns.meshrect()
        length = meshrect[3] - meshrect[0]
        
        self.stage = 0
        
        #starting state : Bloch domain wall along x axis
        self.ns.dwall('z', 'y', length, 0.0)
        self.ns.setfield(Hext, 0, 0)
        self.ns.pbc(self.ns.meshfocus(), 'x', -1)
        
        dT = self.ns.setdt()
        
        #relax Bloch wall at zero Kelvin      
        self.ns.temperature(0)
        self.ns.setstage('Relax')
        self.ns.setode('LLGStatic', 'SDesc')
        self.ns.editstagestop(0, 'mxh', 1e-5)
        #run with cuda 0: make sure this runs in DP mode, otherwise might not reach low enough mxh due to limited precision
        self.ns.cuda(0)
        self.ns.Run()
        
        self.ns.reset()
        self.ns.cuda(1)
        
        #first simulation stage: thermalize, no data save
        if self.use_sllg:
            self.ns.setstage('Relax')
            self.ns.editstagestop(0, 'iter', self.iter_thrm_sllg[self.stage])
            self.ns.setode('sLLG', 'TEuler')
            self.ns.setdt(dT)
        else:
            self.ns.setstage('MonteCarlo')
            self.ns.editstagestop(0, 'iter', self.iter_thrm_dw[self.stage])
        
        #second simulation stage: save enough points to average
        if self.use_sllg:
            self.ns.addstage('Relax')
            self.ns.editstagestop(1, 'iter', self.iter_avrg_sllg[self.stage])
        else:
            self.ns.addstage('MonteCarlo')
            self.ns.editstagestop(1, 'iter', self.iter_avrg_dw[self.stage])
            
        self.ns.editdatasave(1, 'iter', self.iter_avrg_spacing)
        
        self.ns.setdata('<M>')
        self.ns.adddata('e_anis')
        self.ns.adddata('e_exch')
        self.ns.adddata('dwpos_x')
        self.ns.savedatafile('temporary.txt')
        
        self.ns.mcconstrain(0)
        self.ns.mcserial(0, self.ns.meshfocus())
        self.ns.mccomputefields(1)
        self.ns.cuda(1)
        
        if self.use_sllg:
            return self.__sweep_temperature_up(Tsweep_range, use_default_step, self.iter_thrm_sllg, self.iter_avrg_sllg, outputFile, Tmax_dw, True)
        else:
            return self.__sweep_temperature_up(Tsweep_range, use_default_step, self.iter_thrm_dw, self.iter_avrg_dw, outputFile, Tmax_dw, True)
    
    ######################################################################################################
    # Calculate m scaling
    
    def compute_magnetization_scaling(self, Temperature_range = [0, 700], use_default_step = True, outputFile = '', Hext = 0.0, direction = [90, 0], get_energies = False):
        """m scaling for given temperature range, external field and starting orientation along given direction (polar, azimuthal in degrees), save to outputFile.
           if use_default_step = False then expecting Temperature_range to have temperature values to simulate.
           outputFile saves: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
           return: Trange, mrange, susrel_l, (Krange) calculated (Krange is not normalized and included in return only if get_energies == True)
           Trange: simulated temperature values
           mrange: magnetization length scaling
           susrel_l: relative longitudinal susceptibility (units 1/T)
           Krange: anisotropy energy density (J/m^3) as a function of temperature along given direction
        """

        mu_s = self.ns.setparam(self.ns.meshfocus(), 'mu_s')
        J = self.ns.setparam(self.ns.meshfocus(), 'J')
        K = self.ns.setparam(self.ns.meshfocus(), 'K1')
        meshrect = self.ns.meshrect()
        cellsize = self.ns.cellsize()

        if len(outputFile): 
            self.ns.dp_newfile(outputFile)
            self.ns.savecomment(outputFile, 'MC sweep. mu_s = %f (muB), J = %E (J), K = %E (J), size = %f, %f, %f (nm), cellsize = %f, %f, %f (nm), Hext = %f (A/m), theta = %f, phi = %f (deg.)' % 
                                (mu_s, J, K, meshrect[3] / 1e-9, meshrect[4] / 1e-9, meshrect[5] / 1e-9, cellsize[0] / 1e-9, cellsize[1] / 1e-9, cellsize[2] / 1e-9, Hext, direction[0], direction[1]))
            self.ns.savecomment(outputFile, 'T\tm\tm_std\tsusrel_l\tsusrel_x\tsusrel_y\tsusrel_z\te_anis (J/m^3)\te_exch (J/m^3)\tdw_width (m)')
        
        if get_energies:
            Trange, mrange, susrel, Krange = self.__run_MC_classical(Temperature_range, use_default_step, outputFile, Hext, direction, get_energies)
            return Trange, mrange, susrel, Krange
        else:
            Trange, mrange, susrel, Krange = self.__run_MC_classical(Temperature_range, use_default_step, outputFile, Hext, direction, get_energies)
            return Trange, mrange, susrel
    
    def compute_magnetization_scaling_constrained(self, Temperature_range = [0, 700], use_default_step = True, outputFile = '', Hext = 0.0, direction = [90, 0]):
        """m scaling for given temperature range using constraining direction, external field along constraining direction (polar, azimuthal in degrees), save to outputFile
           if use_default_step = False then expecting Temperature_range to have temperature values to simulate
           outputFile saves: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
           return: Trange, mrange, susrel_l, Krange calculated (Krange is not normalized)
           Trange: simulated temperature values
           mrange: magnetization length scaling
           susrel_l: relative longitudinal susceptibility (units 1/T)
           Krange: anisotropy energy density (J/m^3) as a function of temperature along given direction
        """
        
        mu_s = self.ns.setparam(self.ns.meshfocus(), 'mu_s')
        J = self.ns.setparam(self.ns.meshfocus(), 'J')
        K = self.ns.setparam(self.ns.meshfocus(), 'K1')
        meshrect = self.ns.meshrect()
        cellsize = self.ns.cellsize()

        if len(outputFile): 
            self.ns.dp_newfile(outputFile)
            self.ns.savecomment(outputFile, 'CMC sweep. mu_s = %f (muB), J = %E (J), K = %E (J), size = %f, %f, %f (nm), cellsize = %f, %f, %f (nm), Hext = %f (A/m), theta = %f, phi = %f (deg.)' % 
                                (mu_s, J, K, meshrect[3] / 1e-9, meshrect[4] / 1e-9, meshrect[5] / 1e-9, cellsize[0] / 1e-9, cellsize[1] / 1e-9, cellsize[2] / 1e-9, Hext, direction[0], direction[1]))
            self.ns.savecomment(outputFile, 'T\tm\tm_std\tsusrel_l\tsusrel_x\tsusrel_y\tsusrel_z\te_anis (J/m^3)\te_exch (J/m^3)\tdw_width (m)')
        
        return self.__run_MC_constrained(Temperature_range, use_default_step, outputFile, Hext, direction)
    
    ######################################################################################################
    # Calculate uniaxial k scaling
    
    def compute_uniaxial_anisotropy_scaling(self, Temperature_range = [0, 700], use_default_step = True, outputFile = '', hard_axis = [90, 90], easy_axis = [90, 0]):
        """k scaling for given temperature range using constrained MC, save to outputFile
           if use_default_step = False then expecting Temperature_range to have temperature values to simulate
           outputFile saves: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width, first pass along hard axis, second pass along easy axis
           return: Trange, mrange, Krange calculated (Krange is not normalized)
           Trange: simulated temperature values
           mrange: magnetization length scaling
           Krange: anisotropy energy density (J/m^3) as a function of temperature, obtained as the difference between the hard and easy axes directions
        """
        
        #this should already be added, but make sure
        self.ns.addmodule(self.ns.meshfocus(), 'aniuni')
        
        mu_s = self.ns.setparam(self.ns.meshfocus(), 'mu_s')
        J = self.ns.setparam(self.ns.meshfocus(), 'J')
        K = self.ns.setparam(self.ns.meshfocus(), 'K1')
        meshrect = self.ns.meshrect()
        cellsize = self.ns.cellsize()

        if len(outputFile): 
            self.ns.dp_newfile(outputFile)
            self.ns.savecomment(outputFile, 'CMC sweeps, 1) hard, 2) easy. mu_s = %f (muB), J = %E (J), K = %E (J), size = %f, %f, %f (nm), cellsize = %f, %f, %f (nm), hard_axis = %f, %f (deg.), easy_axis = %f, %f (deg.)' % 
                                (mu_s, J, K, meshrect[3] / 1e-9, meshrect[4] / 1e-9, meshrect[5] / 1e-9, cellsize[0] / 1e-9, cellsize[1] / 1e-9, cellsize[2] / 1e-9, hard_axis[0], hard_axis[1], easy_axis[0], easy_axis[1]))
            self.ns.savecomment(outputFile, 'T\tm\tm_std\tsusrel_l\tsusrel_x\tsusrel_y\tsusrel_z\te_anis (J/m^3)\te_exch (J/m^3)\tdw_width (m)')
        
        Trange, mrange, susrel, Krange_hard = self.__run_MC_constrained(Temperature_range, use_default_step, outputFile, 0.0, hard_axis)
        Trange, mrange, susrel, Krange_easy = self.__run_MC_constrained(Trange, False, outputFile, 0.0, easy_axis)
        
        Krange = [K_hard - K_easy for (K_hard, K_easy) in zip(Krange_hard, Krange_easy)]
        return Trange, mrange, Krange
    
    ######################################################################################################
    # Calculate cubic k scaling
    
    def compute_cubic_anisotropy_scaling(self, Temperature_range = [0, 700], use_default_step = True, outputFile = '', hard_axis = [np.degrees(np.arctan(np.sqrt(2))), 45], easy_axis = [90, 0]):
        """k scaling for given temperature range using constrained MC, save to outputFile
           if use_default_step = False then expecting Temperature_range to have temperature values to simulate
           outputFile saves: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width, first pass along hard axis, second pass along easy axis
           return: Trange, mrange, Krange calculated (Krange is not normalized)
           Trange: simulated temperature values
           mrange: magnetization length scaling
           Krange: anisotropy energy density (J/m^3) as a function of temperature, obtained as 3x the difference between the hard (cube diagonal) and easy axes directions
        """
        
        #this should already be added, but make sure
        self.ns.addmodule(self.ns.meshfocus(), 'anicubi')
        
        mu_s = self.ns.setparam(self.ns.meshfocus(), 'mu_s')
        J = self.ns.setparam(self.ns.meshfocus(), 'J')
        K = self.ns.setparam(self.ns.meshfocus(), 'K1')
        meshrect = self.ns.meshrect()
        cellsize = self.ns.cellsize()

        if len(outputFile): 
            self.ns.dp_newfile(outputFile)
            self.ns.savecomment(outputFile, 'CMC sweeps, 1) hard, 2) easy. mu_s = %f (muB), J = %E (J), K = %E (J), size = %f, %f, %f (nm), cellsize = %f, %f, %f (nm), hard_axis = %f, %f (deg.), easy_axis = %f, %f (deg.)' % 
                                (mu_s, J, K, meshrect[3] / 1e-9, meshrect[4] / 1e-9, meshrect[5] / 1e-9, cellsize[0] / 1e-9, cellsize[1] / 1e-9, cellsize[2] / 1e-9, hard_axis[0], hard_axis[1], easy_axis[0], easy_axis[1]))
            self.ns.savecomment(outputFile, 'T\tm\tm_std\tsusrel_l\tsusrel_x\tsusrel_y\tsusrel_z\te_anis (J/m^3)\te_exch (J/m^3)\tdw_width (m)')
        
        Trange, mrange, susrel, Krange_hard = self.__run_MC_constrained(Temperature_range, use_default_step, outputFile, 0.0, hard_axis)
        Trange, mrange, susrel, Krange_easy = self.__run_MC_constrained(Trange, False, outputFile, 0.0, easy_axis)
        
        Krange = [K_hard - K_easy for (K_hard, K_easy) in zip(Krange_hard, Krange_easy)]
        return Trange, mrange, Krange
    
    ######################################################################################################
    # Calculate exchange stiffness a scaling
    
    def compute_exchange_stiffness_scaling(self, Temperature_range = [0, 700], use_default_step = True, outputFile = '', Hext = 0.0, K = 5e-23):
        """a scaling for given temperature range using domain wall free energy method with MC, save to outputFile
           if use_default_step = False then expecting Temperature_range to have temperature values to simulate
           outputFile saves: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width, first pass for uniform easy axis state, second pass for Bloch domain wall state transverse to easy axis
           return: Trange, mrange, Arange calculated, Krange calculated (Arange and Krange not normalized)
           Trange: simulated temperature values
           mrange: magnetization length scaling
           Krange: micromagnetic anisotropy energy density (J/m^3) as a function of temperature
           Arange: micromagnetic exchange stiffness (J/m) as a function of temperature
        """
        
        #Find exchange stiffness scaling, A = A0 * a(T), with the domain wall free energy method (e.g. [Hinzke et al., PRB 77, 094407 (2008)])
        #
        # For a Bloch domain wall with uniaxial anisotropy axis perpendicular to it, we have:
        #
        # DW_width = PI * SQRT(A / K)
        # dF = 4 * C * SQRT(A * K)
        #
        # Here dF is the difference in free energy between the state with a domain wall and uniform state along easy axis
        #
        # dF is obtained from the total internal energy difference dE as 
        # F(1/T) = -T * Integral(1 / T to 0) { dE(1/T') d(1/T')}. NOTE : this can be truncated to just above Tc, since dE(T') ~= 0 for T' > Tc.
        # Evaluate as: F(T) ~= -T * Integral(T, Tmax) { dE(T') d(1/T') }, where Tmax > Tc. When numerically evaluating integral, better evaluate it using this form than with the substitution d(1/T') -> -dT' / T'^2.
        # dE is the internal energy difference between the state with a domain wall and uniform state along easy axis, where the internal energy is just the sum of the exchange and anisotropy energies (assuming others relatively negligible)
        #
        # C is the cross-section area, where the domain wall longitudinal direction is normal to it
        #
        # Eliminating K gives:
        # 
        # A = dF * DW_width / (4 * C * PI)
        #
        # This formula works until close to Tc, as long as we can fit for a DW width
        # When close to Tc and above, use the dF formula with K(T) computed separately
        
        #1.
        #Compute K(T) using the method in simulate_m_k_uniaxial_scaling : run algorithm along easy axis, and obtain K scaling as 1.0 - 3*K_easy / 2*K0_hard
        #At the same time save internal energy data (as energy density) for uniform state along easy axis: e_anis and e_exch
                
        meshname = self.ns.meshfocus()
        
        self.ns.pbc(meshname, 'x', 1)
        self.ns.addmodule(meshname, 'aniuni')
        self.ns.setparam(meshname, 'K1', K)
        
        #easy axis along z : Bloch domain wall along x axis
        self.ns.setparam(meshname, 'ea1', [0, 0, 1])
        
        self.ns.setangle(90, 0, meshname)
        self.ns.computefields()
        K0 = self.ns.showdata('e_anis')

        Trange, mrange, susrel, Krange = self.compute_magnetization_scaling(Temperature_range, use_default_step, outputFile, Hext, direction = [0, 0], get_energies = True)
        #data = self.ns.Get_Data_Columns(outputFile, [0, 1, 7])
        #Trange = data[0][:int(len(data[0])/2)]
        #mrange = data[1][:int(len(data[1])/2)]
        #Krange = data[2][:int(len(data[2])/2)]
        
        #K(T) = K0 * k(T) data
        Krange_mm = [K0 - 1.5 * Keasy for Keasy in Krange]
        
        #2.
        
        #Compute internal energy as a function of temperature for a Bloch domain wall, as well as domain wall width
        Tmax_dw = 0.0
        for (m, T) in zip(mrange, Trange):
            if m <= self.dw_width_m_setpoint:
                Tmax_dw = T
                break
            
        self.__run_MC_BlochWall(Trange, False, outputFile, Hext, Tmax_dw)
        
        #3. Compute free energy
        
        data = self.ns.Get_Data_Columns(outputFile, [7, 8, 9])
        e_anis_unif = data[0][:int(len(data[0])/2)]
        e_exch_unif = data[1][:int(len(data[1])/2)]
        e_anis_dw = data[0][int(len(data[0])/2):]
        e_exch_dw = data[1][int(len(data[1])/2):]
        dw_width = data[2][int(len(data[2])/2):]
        
        #compute internal energy
        dE = [ea_dw + ex_dw - ea_u - ex_u for (ea_u, ea_dw, ex_u, ex_dw) in zip(e_anis_unif, e_anis_dw, e_exch_unif, e_exch_dw)]
        
        trapezoidal_sum_terms = []
        for idx in range(len(Trange) - 1):
            
            dE_j1, dE_j2 = dE[idx], dE[idx + 1]
            T_j1, T_j2 = Trange[idx], Trange[idx + 1]
            
            if T_j1 > 0:
                trapezoidal_sum_terms.append((dE_j1 + dE_j2) * (1.0 / T_j2 - 1.0 / T_j1) / 2.0)
            else:
                trapezoidal_sum_terms.append(0.0)

        #last point should be zero : dE, dF zero above Tc
        trapezoidal_sum_terms.append(0.0)
            
        dF = []
        #dF at T = 0 needs special treatment (still slightly inaccurate)
        if Trange[0] == 0.0:
            dF.append((dE[0] + dE[1]) / 2)
            for idx in range(1, len(Trange)): 
                dF.append(-Trange[idx] * sum(trapezoidal_sum_terms[idx:]))
        else:
            for idx in range(len(Trange)): 
                dF.append(-Trange[idx] * sum(trapezoidal_sum_terms[idx:]))
        
        meshrect = self.ns.meshrect()
        L = meshrect[3] - meshrect[0]
        
        #4. Compute A(T)
        
        #A calculated 3 ways:
        #1. Arange_dw: dw width only, switching to free energy only below dw_width_m_setpoint where dw width not available
        #2. Arange_dF: free energy only
        #3. Arange_dFdw: combined method, except at below dw_width_m_setpoint where free energy only is used
        Arange_dw, Arange_dF, Arange_dFdw = [], [], []
        for (dFv, dw, m, Kv) in zip(dF, dw_width, mrange, Krange_mm):
         
            Adw = (dw / np.pi)**2 * Kv
            AdF = (dFv * L)**2 / (16.0 * Kv)
            AdFdw = dFv * L * dw / (4 * np.pi)
            
            if m > self.dw_width_m_setpoint:
                Arange_dw.append(Adw)
                Arange_dF.append(AdF)
                Arange_dFdw.append(AdFdw)
            else:
                #below m dw setpoint only use dF since the dw width is too inaccurate in general.
                #With larger cross-sections and more averages the setpoint can be decreased, but default 0.1 is very close to Tc (T/Tc ~ 0.99).
                Arange_dw.append(AdF)
                Arange_dF.append(AdF)
                Arange_dFdw.append(AdF)
        """
        #A(0) found directly
        if Trange[0] == 0.0: 
            Arange_dw[0] = self.ns.showa()
            Arange_dF[0] = Arange_dw[0]
            Arange_dFdw[0] = Arange_dw[0]
        """
            
        return Trange, mrange, Arange_dw, Arange_dF, Arange_dFdw, Krange_mm
    
    ######################################################################################################
    # Calculate interfacial surface exchange scaling
    
    def compute_RKKY_scaling(self, 
                             bMesh, tMesh, Temperature_range = [0, 700], use_default_step = True, outputFile = '', 
                             use_H_method = False, mrange = [], krange = []):
        """RKKY scaling for given temperature range, save to outputFile.
           if use_default_step = False then expecting Temperature_range to have temperature values to simulate.
           outputFile saves: T m e_anis e_surfexch
        """

        mu_s = self.ns.setparam(self.ns.meshfocus(), 'mu_s')
        J = self.ns.setparam(self.ns.meshfocus(), 'J')
        K = self.ns.setparam(self.ns.meshfocus(), 'K1')
        Js = self.ns.setparam(self.ns.meshfocus(), 'Js')
        meshrect = self.ns.meshrect()
        cellsize = self.ns.cellsize()

        if len(outputFile): 
            self.ns.dp_newfile(outputFile)
            self.ns.savecomment(outputFile, 'MC sweep, RKKY. mu_s = %f (muB), J = %E (J), K = %E (J), Js = %E (J), size = %f, %f, %f (nm), cellsize = %f, %f, %f (nm).' % 
                                (mu_s, J, K, Js, (meshrect[3] - meshrect[0]) / 1e-9, (meshrect[4] - meshrect[1]) / 1e-9, (meshrect[5] - meshrect[2]) / 1e-9, cellsize[0] / 1e-9, cellsize[1] / 1e-9, cellsize[2] / 1e-9))
            self.ns.savecomment(outputFile, 'T\tm\te_anis (J/m^3)\te_surfexch (J/m^3)')
        
        self.__run_MC_classical_RKKY(bMesh, tMesh, Temperature_range, use_default_step, outputFile, use_H_method, mrange, krange)
    
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
        
######################################################
# Calculate m and fit 

def simulate_m_scaling(ns, fileName, Tmax = 700, Hext = 0.0, averaging_level = 2):
    """
    This routine is intended for efficient simulation of m scaling only. It cannot be used to simulate k dependence also.
    Sweep temperature from 0 to Tmax along x axis, with Hext applied field (A/m).
    Save results in fileName as: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
    Return: Tc, beta, Tc_std, beta_std, obtained from Bloch fit.
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level)
    
    #Set Euler to save memory - LLG not used with Monte-Carlo
    ns.setode('LLG', 'Euler')
    
    mc.compute_magnetization_scaling([0, Tmax], True, fileName, Hext)
    data = ns.Get_Data_Columns(fileName, [0, 1])
    Trange, mrange = data[0], data[1]

    #m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
    
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    return Tc, beta, Tc_std, beta_std

def simulate_m_scaling_constrained(ns, fileName, Tmax = 700, Hext = 0.0, direction = [90, 0], averaging_level = 2):
    """
    This routine is intended for testing of m scaling simulation using constrained MC only. It will also calculate the anisotropy energy density along constraining direction.
    Sweep temperature from 0 to Tmax along given direction, with Hext applied field (A/m).
    Save results in fileName as: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
    Return: Tc, beta, Tc_std, beta_std, obtained from Bloch fit.
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level)
    
    #Set Euler to save memory - LLG not used with Monte-Carlo
    ns.setode('LLG', 'Euler')
    
    mc.compute_magnetization_scaling_constrained([0, Tmax], True, fileName, Hext, direction)
    data = ns.Get_Data_Columns(fileName, [0, 1])
    Trange, mrange = data[0], data[1]

    #m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
        
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    return Tc, beta, Tc_std, beta_std

######################################################
# Calculate susceptibilities

def simulate_susceptibilities(ns, fileName, Tmax = 700, Hext = 0.0, averaging_level = 2):
    """
    This routine is intended mainly for simulation of susceptibilities. It cannot be used to simulate k dependence also.
    Sweep temperature from 0 to Tmax along x axis, with Hext applied field (A/m).
    Save results in fileName as: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
    Return: Tc, beta, Tc_std, beta_std, obtained from Bloch fit.
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level, fixed_average = True)
    
    #Set Euler to save memory - LLG not used with Monte-Carlo
    ns.setode('LLG', 'Euler')
    
    mc.compute_magnetization_scaling([0, Tmax], True, fileName, Hext)
    data = ns.Get_Data_Columns(fileName, [0, 1, 3, 4, 5, 6])
    Trange, mrange, susrel_l, susrel_x, susrel_y, susrel_z = data[0], data[1], data[2], data[3], data[4], data[5]

    #1. m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
        
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    #2. susceptibilities
    plt.axes(xlabel = 'Temperature (K)', ylabel = r'$\tilde{\chi}$ (1/T)')
    plt.yscale('log')
    
    plt.plot(Trange, susrel_l, 's', label = 'l')
    plt.plot(Trange, susrel_x, '>', label = 'x')
    plt.plot(Trange, susrel_y, '^', label = 'y')
    plt.plot(Trange, susrel_z, 'o', label = 'z')
    
    plt.legend()
    plt.savefig(fileName_noextension + '_susrel.png')
    plt.show()
    
    return Tc, beta, Tc_std, beta_std

######################################################
# Calculate k and fit using constrained Monte Carlo

def simulate_k_uniaxial_scaling(ns, fileName, Tmax = 700, hard_axis = [90, 90], easy_axis = [90, 0], averaging_level = 2):
    """
    This routine is intended for simulation of uniaxial anisotropy scaling. It uses the constrained MC algorithm.
    Sweep temperature from 0 to Tmax along hard axis, then along easy axis.
    Save results in fileName as: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
    Return: Tc, beta, power, Tc_std, beta_std, power_std obtained from Bloch fit, with power being k(T) ~= m(T) ^ power obtained by fitting m to k.
    Additionally it saves in fileName + '_T_m_K.txt' the T, m, and K computed values.
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level)
    
    #Set Euler to save memory - LLG not used with Monte-Carlo
    ns.setode('LLG', 'Euler')
    
    ns.addmodule(ns.meshfocus(), 'aniuni')
    theta = np.radians(easy_axis[0]); phi = np.radians(easy_axis[1])
    ns.setparam(ns.meshfocus(), 'ea1', [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    
    ns.setangle(hard_axis[0], hard_axis[1])
    ns.computefields()
    K0 = ns.showdata('e_anis')
    
    Trange, mrange, Krange = mc.compute_uniaxial_anisotropy_scaling([0, Tmax], True, fileName, hard_axis, easy_axis)
    krange = [value / K0 for value in Krange]
    
    #1. m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
        
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    #2. k scaling plot
    power, power_std = bfit.fit_m_power(Trange, krange, mrange)
    print('Anisotropy m scaling power: %0.2f +/- %0.2f' % (power, power_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'k')
    plt.plot(Trange, krange, 'o', label = 'computed')
    plt.plot(Trange, bfit.m_power_fitted(mrange, power), '--', color = 'black', label = 'fitted m$^{%0.2f}$' % power)
    plt.legend()
    plt.savefig(fileName_noextension + '_k.png')
    plt.show()
    
    fileName_processed = fileName_noextension + '_T_m_K.txt'
    ns.Save_Data_Columns(fileName_processed, [Trange, mrange, Krange])
    
    return Tc, beta, power, Tc_std, beta_std, power_std  
    
def simulate_k_cubic_scaling(ns, fileName, Tmax = 700, hard_axis = [np.degrees(np.arctan(np.sqrt(2))), 45], easy_axis1 = [90, 0], easy_axis2 = [90, 90], averaging_level = 2):
    """
    This routine is intended for simulation of cubic anisotropy scaling. It uses the constrained MC algorithm.
    Sweep temperature from 0 to Tmax along hard axis, then along easy axis.
    Save results in fileName as: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
    Return: Tc, beta, power, Tc_std, beta_std, power_std obtained from Bloch fit, with power being k(T) ~= m(T) ^ power obtained by fitting m to k.
    Additionally it saves in fileName + '_T_m_K.txt' the T, m, and K computed values.
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level)
    
    #Set Euler to save memory - LLG not used with Monte-Carlo
    ns.setode('LLG', 'Euler')
    
    ns.addmodule(ns.meshfocus(), 'anicubi')
    theta = np.radians(easy_axis1[0]); phi = np.radians(easy_axis1[1])
    ns.setparam(ns.meshfocus(), 'ea1', [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    theta = np.radians(easy_axis2[0]); phi = np.radians(easy_axis2[1])
    ns.setparam(ns.meshfocus(), 'ea2', [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    
    ns.setangle(hard_axis[0], hard_axis[1])
    ns.computefields()
    K0 = ns.showdata('e_anis')
    
    Trange, mrange, Krange = mc.compute_cubic_anisotropy_scaling([0, Tmax], True, fileName, hard_axis, easy_axis1)
    krange = [value / K0 for value in Krange]
    
    #1. m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
        
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    #2. k scaling plot
    power, power_std = bfit.fit_m_power(Trange, krange, mrange)
    print('Anisotropy m scaling power: %0.2f +/- %0.2f' % (power, power_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'k')
    plt.plot(Trange, krange, 'o', label = 'computed')
    plt.plot(Trange, bfit.m_power_fitted(mrange, power), '--', color = 'black', label = 'fitted m$^{%0.2f}$' % power)
    plt.legend()
    plt.savefig(fileName_noextension + '_k.png')
    plt.show()
    
    fileName_processed = fileName_noextension + '_T_m_K.txt'
    ns.Save_Data_Columns(fileName_processed, [Trange, mrange, Krange])
    
    return Tc, beta, power, Tc_std, beta_std, power_std
    
######################################################
# Calculate both m and k, and fit 

def simulate_m_k_uniaxial_scaling(ns, fileName, Tmax = 700, hard_axis = [90, 90], easy_axis = [90, 0], averaging_level = 2):
    """
    This routine is intended for efficient simulation of both m and uniaxial k scaling. 
    It uses the standard MC algorithm to simulate along the easy axis only, with anisotropy energy density temperature scaling obtained as 1.0 - 3*K_easy / 2*K0_hard.
    Sweep temperature from 0 to Tmax along easy axis.
    Save results in fileName as: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
    Return: Tc, beta, power, Tc_std, beta_std, power_std obtained from Bloch fit, with power being k(T) ~= m(T) ^ power obtained by fitting m to k.
    Additionally it saves in fileName + '_T_m_K.txt' the T, m, and K computed values.
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level)
    
    #Set Euler to save memory - LLG not used with Monte-Carlo
    ns.setode('LLG', 'Euler')
    
    ns.addmodule(ns.meshfocus(), 'aniuni')
    theta = np.radians(easy_axis[0]); phi = np.radians(easy_axis[1])
    ns.setparam(ns.meshfocus(), 'ea1', [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    
    ns.setangle(hard_axis[0], hard_axis[1])
    ns.computefields()
    K0 = ns.showdata('e_anis')
    
    Trange, mrange, susrel, Krange = mc.compute_magnetization_scaling([0, Tmax], True, fileName, Hext = 0.0, direction = easy_axis, get_energies = True)
    krange = [1.0 - K * (3.0/2) / K0 for K in Krange]
    Krange = [K0 * k for k in krange]

    #1. m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
        
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    #2. k scaling plot
    power, power_std = bfit.fit_m_power(Trange, krange, mrange)
    print('Anisotropy m scaling power: %0.2f +/- %0.2f' % (power, power_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'k')
    plt.plot(Trange, krange, 'o', label = 'computed')
    plt.plot(Trange, bfit.m_power_fitted(mrange, power), '--', color = 'black', label = 'fitted m$^{%0.2f}$' % power)
    plt.legend()
    plt.savefig(fileName_noextension + '_k.png')
    plt.show()
    
    fileName_processed = fileName_noextension + '_T_m_K.txt'
    ns.Save_Data_Columns(fileName_processed, [Trange, mrange, Krange])
    
    return Tc, beta, power, Tc_std, beta_std, power_std

def simulate_m_k_cubic_scaling(ns, fileName, Tmax = 700, hard_axis = [np.degrees(np.arctan(np.sqrt(2))), 45], easy_axis1 = [90, 0], easy_axis2 = [90, 90], averaging_level = 2):
    """
    This routine is intended for efficient simulation of both m and uniaxial k scaling. 
    It uses the standard MC algorithm to simulate along the easy axis only, with anisotropy energy density temperature scaling obtained as 1.0 - 5*K_easy / 3*K0_hard.
    Sweep temperature from 0 to Tmax along easy axis.
    Save results in fileName as: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
    Return: Tc, beta, power, Tc_std, beta_std, power_std obtained from Bloch fit, with power being k(T) ~= m(T) ^ power obtained by fitting m to k.
    Additionally it saves in fileName + '_T_m_K.txt' the T, m, and K computed values.
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level)
    
    #Set Euler to save memory - LLG not used with Monte-Carlo
    ns.setode('LLG', 'Euler')
    
    ns.addmodule(ns.meshfocus(), 'anicubi')
    theta = np.radians(easy_axis1[0]); phi = np.radians(easy_axis1[1])
    ns.setparam(ns.meshfocus(), 'ea1', [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    theta = np.radians(easy_axis2[0]); phi = np.radians(easy_axis2[1])
    ns.setparam(ns.meshfocus(), 'ea2', [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    
    ns.setangle(hard_axis[0], hard_axis[1])
    ns.computefields()
    K0 = ns.showdata('e_anis')
    
    Trange, mrange, susrel, Krange = mc.compute_magnetization_scaling([0, Tmax], True, fileName, Hext = 0.0, direction = easy_axis1, get_energies = True)
    krange = [1.0 - K * (5.0/3) / K0 for K in Krange]
    Krange = [K0 * k for k in krange]

    #1. m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
        
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    #2. k scaling plot
    power, power_std = bfit.fit_m_power(Trange, krange, mrange)
    print('Anisotropy m scaling power: %0.2f +/- %0.2f' % (power, power_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'k')
    plt.plot(Trange, krange, 'o', label = 'computed')
    plt.plot(Trange, bfit.m_power_fitted(mrange, power), '--', color = 'black', label = 'fitted m$^{%0.2f}$' % power)
    plt.legend()
    plt.savefig(fileName_noextension + '_k.png')
    plt.show()
    
    fileName_processed = fileName_noextension + '_T_m_K.txt'
    ns.Save_Data_Columns(fileName_processed, [Trange, mrange, Krange])
    
    return Tc, beta, power, Tc_std, beta_std, power_std
        
######################################################
# Calculate exchange stiffness
    
def simulate_exchange_stiffness(ns, fileName, Tmax = 700, Hext = 0.0, K = 5e-23, sllg = False, averaging_level = 2):
    """
    This routine is intended for simulation of exchange stiffness using the domain wall method. It can also obtain the k (uniaxial) and m scalings at the same time.
    Sweep temperature from 0 to Tmax along easy axis with applied field Hext (A/m).
    Save results in fileName as: T m m_std susrel_l susrel_x susrel_y susrel_z e_anis e_exch dw_width
    Return: Tc, beta, power_k, power_a, Tc_std, beta_std, power_k_std, power_a_std obtained from Bloch fit, with power being k(T) ~= m(T) ^ power_k, resp. a(T) ~= m(T) ^ power_a, obtained by fitting m to k.
    Additionally it saves in fileName + '_T_m_A_K.txt' the T, m, A, and K computed values.
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level, use_sllg = sllg)
    
    if sllg:
        ns.setode('sLLG', 'TEuler')
        ns.setdt(0.5e-15)
    else:
        #Set Euler to save memory - LLG not used with Monte-Carlo
        ns.setode('LLG', 'Euler')
    
    Trange, mrange, Arange_dw, Arange_dF, Arange_dFdw, Krange = mc.compute_exchange_stiffness_scaling([0, Tmax], True, fileName, Hext, K)
    arange_dw = [value / Arange_dw[0] for value in Arange_dw]
    arange_dF = [value / Arange_dF[0] for value in Arange_dF]
    arange_dFdw = [value / Arange_dFdw[0] for value in Arange_dFdw]
    krange = [value / Krange[0] for value in Krange]
    
    #1. m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
        
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    #2. k scaling plot
    power_k, power_k_std = bfit.fit_m_power(Trange, krange, mrange)
    print('Anisotropy m scaling power: %0.2f +/- %0.2f' % (power_k, power_k_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'k')
    plt.plot(Trange, krange, 'o', label = 'computed')
    plt.plot(Trange, bfit.m_power_fitted(mrange, power_k), '--', color = 'black', label = 'fitted m$^{%0.2f}$' % power_k)
    plt.legend()
    plt.savefig(fileName_noextension + '_k.png')
    plt.show()
    
    #3. a scaling plot
    power_a, power_a_std = bfit.fit_m_power(Trange, arange_dw, mrange)
    print('Exchange stiffness m scaling power: %0.2f +/- %0.2f' % (power_a, power_a_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'a')
    #plt.plot(Trange, arange_dFdw, 's', label = 'computed (dFdw)')
    #plt.plot(Trange, arange_dF, '.', label = 'computed (dF)')
    plt.plot(Trange, arange_dw, 'o', label = 'computed (dw)')
    plt.plot(Trange, bfit.m_power_fitted(mrange, power_a), '--', color = 'black', label = 'fitted m$^{%0.2f}$' % power_a)
    plt.legend()
    plt.savefig(fileName_noextension + '_a.png')
    plt.show()
    
    fileName_processed = fileName_noextension + '_T_m_A_K.txt'
    ns.Save_Data_Columns(fileName_processed, [Trange, mrange, Arange_dw, Krange])
    
    return Tc, beta, power_k, power_a, Tc_std, beta_std, power_k_std, power_a_std

######################################################
# Calculate RKKY surface exchange coupling
    
def simulate_RKKY(ns, bMesh, tMesh, fileName, Tmax = 700, averaging_level = 2, use_H_method = False, skip_normal_method = False):
    """
    This routine is intended for simulation of surface exchange coupling temperature dependence. It can also obtain m scaling at the same time.
    Sweep temperature from 0 to Tmax along easy axis with applied field Hext (A/m).
    Save results in fileName as: T m e_anis e_surfexch
    """
    
    mc = MonteCarlo(ns, averaging_level = averaging_level)
    
    #Set Euler to save memory - LLG not used with Monte-Carlo
    ns.setode('LLG', 'Euler')
    
    if not skip_normal_method: mc.compute_RKKY_scaling(bMesh, tMesh, [0, Tmax], True, fileName)
    
    data = ns.Get_Data_Columns(fileName, [0, 1, 2, 3])
    Trange, mrange, Krange, J1range = data[0], data[1], data[2], data[3]
    j1range = [value / J1range[0] if value / J1range[0] > 0 else 0.0 for value in J1range]
    K0 = ns.showk()
    krange = [1.0 - K * (3.0/2) / K0 for K in Krange]
    
    #1. m scaling plot
    bfit = BlochFit()
    Tc, beta, Tc_std, beta_std = bfit.fit_Bloch(Trange, mrange)
    print('Bloch fit Curie temperature: %0.2f K +/- %0.2f, with beta = %0.3f +/- %0.3f' % (Tc, Tc_std, beta, beta_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
    plt.plot(Trange, mrange, 'o', label = 'computed', zorder = -1)
    plt.plot(Trange, bfit.Bloch_fitted(Tc, beta, Trange), '--', color = 'black', label = 'Bloch. Tc = %0.2f K' % Tc)
    plt.legend()
    
    fileName_noextension = fileName
    if len(fileName) > 4 and fileName[-4:] == '.txt':
        fileName_noextension = fileName[:-4]
        
    plt.savefig(fileName_noextension + '_m.png')
    plt.show()
    
    #2. J1 scaling plot
    data_red = np.asarray([[T, m ,j] for (T, m, j) in zip(Trange, mrange, j1range) if T < Tc]).T
    Tred, mred, j1red = data_red[0], data_red[1], data_red[2]
    
    popt, pcov = curve_fit(lambda T, x, y: 1 - x * (T/Tc)**y, np.asarray(Tred), j1red)
    x, y = popt[0], popt[1]
    power_J1, power_J1_std = bfit.fit_m_power(Trange, j1range, mrange)
    print('RKKY m scaling power: %0.2f +/- %0.2f' % (power_J1, power_J1_std))
    
    plt.axes(xlabel = 'Temperature (K)', ylabel = 'j1')
    plt.plot(Trange, j1range, 'o', label = 'computed')
    plt.plot(Trange, bfit.m_power_fitted(mrange, power_J1), '--', color = 'black', label = 'fitted m$^{%0.2f}$' % power_J1)
    plt.plot(Tred, 1 - x * (np.asarray(Tred) / Tc)**y, '--', color = 'black', label = 'fitted 1 - %0.2ft^%0.2f' % (x, y))
    plt.legend()
    plt.savefig(fileName_noextension + '_j1.png')
    plt.show()
    
    #3. j1 vs m plot    
    plt.axes(xlabel = 'm', ylabel = 'j1')
    plt.plot(mred, j1red, 'o', label = 'computed')
    
    plt.plot(mred, 1 - x * (np.asarray(Tred) / Tc)**y, '--', color = 'black', label = 'fitted 1 - %0.2ft^%0.2f' % (x, y))
    plt.plot(mred, bfit.m_power_fitted(mred, power_J1), '--', color = 'blue', label = 'fitted m$^{%0.2f}$' % power_J1)
    plt.legend()
    plt.yscale('log')
    plt.savefig(fileName_noextension + '_j1vsm.png')
    plt.show()
    
    if use_H_method:
        outputFile = fileName_noextension + '_Hmethod.txt'
        mc.compute_RKKY_scaling(bMesh, tMesh, Trange, True, outputFile, True, mrange, krange)
        data = ns.Get_Data_Columns(outputFile, [3])
        j1red_Hmethod = [value / data[0][0] for (value, T) in zip(data[0], Trange) if T < Tc]        
        #4. J1 scaling plot comparing the 2 methods
        plt.axes(xlabel = 'm', ylabel = 'j1')
        plt.plot(mred, j1red, 'o', label = 'computed (energy)')
        plt.plot(mred, j1red_Hmethod, 's', label = 'computed (Hsat)')
        plt.plot(mred, bfit.m_power_fitted(mred, power_J1), '--', color = 'black', label = 'fitted m$^{%0.2f}$' % power_J1)
        plt.legend()
        plt.yscale('log')
        plt.savefig(fileName_noextension + '_j1_2methods.png')
        plt.show()
        
    
    
    
    