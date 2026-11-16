import numpy as np

class conv:
    def __init__(self, pos, exp, mode, A3 = 0):
        self.pos = pos
        self.exp = exp
        self.mode = mode
        self.A3 = A3
        self.min = 0
        self.max = -1

        self.Lf = {'upper':3.7, 'cent':2.5+1.4115219637390113,'lower':3.7}

        self.He3 = 300

        self.t_scale = {'upper':np.linspace(0.137,0.144, 500),'cent':np.linspace(0.134,0.14, 500),'lower':np.linspace(0.137,0.144, 500)}

        self.I = np.array([])
        self.I_err = np.array([])
        self.dE = np.array([])

        self.t = np.array([])
        self.Ei = np.array([])
        self.Ef = np.array([])

        self.tSx =  {'upper':-6*np.pi/180, 'cent':0,'lower':6*np.pi/180}#{'upper':-8.75*np.pi/180, 'cent':0,'lower':8.75*np.pi/180} #{'upper':0, 'cent':0,'lower':0}#

        self.tDy = 0.009668565550687947
        self.tDy_array = np.array([i*self.tDy for i in range(self.He3)])

        self.ki = np.array([])
        self.kf = np.array([])

        self.qh = np.array([])
        self.qk = np.array([])

    def dataset(self):
        self.I = np.genfromtxt(f'Data_reduced/data_{self.exp}_{self.mode}_{self.pos}.txt')

        self.t = np.vstack([np.hstack([np.transpose(np.vstack([self.t_scale[self.pos]  for _ in range(80)]))]) for _ in range(self.He3)])#np.vstack([np.transpose((np.linspace(0.26, 0.29, 500))) for _ in range(self.He3)])
        self.Ef = np.vstack([np.hstack([np.vstack([np.genfromtxt(f'Storage/Ef_model_{self.pos}_{self.mode}.txt')  for _ in range(500)])]) for _ in range(self.He3)]) #[self.min:self.max]
        self.Ei = self.ToF_conversion()
        self.dE = self.Ei - self.Ef
            
        self.ki = 0.69469*np.sqrt(self.Ei)
        self.kf =  0.69469*np.sqrt(self.Ef)

        self.ql = - self.kf*np.sin(self.tSx[self.pos])
        self.qh = np.vstack([self.ki[500*x:500*(x+1),:] - self.kf[500*x:500*(x+1),:]*np.cos(self.tSx[self.pos])*np.cos(self.tDy_array[x]) for x in range(self.He3)])
        self.qk = np.vstack([ -self.kf[500*x:500*(x+1),:]*np.cos(self.tSx[self.pos])*np.sin(self.tDy_array[x]) for x in range(self.He3)])
        # self.qh = np.cos(self.A3)*self.qh - np.sin(self.A3)*self.qk
        # self.qk = np.sin(self.A3)*self.qh + np.cos(self.A3)*self.qk
        
        self.name = f'Data_reduced/{self.exp}_{self.mode}_dataredu_{self.pos}.txt'
        print('SAVING REDUCED FILE:', self.name)
        np.savetxt(self.name, np.array([self.I.flatten(), np.sqrt(self.I.flatten()), self.dE.flatten(), self.qh.flatten(), self.qk.flatten(),self.ql.flatten()]).T) 


    def ToF_conversion(self):
        J_to_meV = 6.24150907*1e21
        m_n = 1.674927471e-27 #kg        

        Li =  162.001-6.36941 #m 
        d_Li = {'upper': 0.1, 'cent':0.36,'lower':0.1}#{'upper': 0.45, 'cent':0.22,'lower':0.45}
        t_offset = 0.00791506 #s

        self.t_to_Ei=((Li+d_Li[self.pos])**2 * m_n) /(2*((self.t-t_offset)-np.sqrt((np.array(self.Lf[self.pos])**2*m_n)/(2*np.array(self.Ef)*meV_to_J)))**2) *J_to_meV 
        return self.t_to_Ei
    