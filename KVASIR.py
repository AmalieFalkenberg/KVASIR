import numpy as np
import mcstasscript as ms 

def KVASIR_instrument(exp, mode, start = 0, verbrose = 0):
    KVASIR = ms.McStas_instr("KVASIR")

    t  = KVASIR.add_parameter("t", value=0)     
    A3 = KVASIR.add_parameter("A3", value=0)

    l_analyzers = [np.float64(2.75039750397504), np.float64(2.7366273662736624), np.float64(2.722467224672247), np.float64(2.707917079170792), np.float64(2.6929469294692945), np.float64(2.677556775567756), np.float64(2.6617766177661775), np.float64(2.645576455764558), np.float64(2.6289562895628955), np.float64(2.6118861188611886), np.float64(2.5943959439594395), np.float64(2.5851858518585185), np.float64(2.5943959439594395), np.float64(2.6032760327603275), np.float64(2.6117961179611795), np.float64(2.6200162001620013), np.float64(2.6278762787627876), np.float64(2.6354363543635433), np.float64(2.6426364263642634), np.float64(2.649536495364954), np.float64(2.6561065610656107), np.float64(2.6118861188611886), np.float64(2.6289562895628955), np.float64(2.645576455764558), np.float64(2.6617766177661775), np.float64(2.677556775567756), np.float64(2.6929469294692945), np.float64(2.707917079170792), np.float64(2.722467224672247), np.float64(2.7366273662736624), np.float64(2.75039750397504)]
    cent_frac = 1/3
    d_d_SDy_cent = 0.4

    N_analyzer = 31
    width = 0.0127 #m
    dt = 2*np.arcsin(width/(2*2.5))*180/np.pi
    theta =  np.linspace(dt*N_analyzer+dt/2, -dt*N_analyzer-dt/2,N_analyzer)


    #### PARAMETERS ####
    print('SETTING PARAMETERS')

    width= KVASIR.add_declare_var("double", "width", value=0.0127)
    d_PG= KVASIR.add_declare_var("double", "d_PG", value=3.355) 
    u = KVASIR.add_declare_var("double","u", value=1e-5)
    h_D = KVASIR.add_declare_var("double", "h_D", value=0.4) #heigt of detector 

    KVASIR.add_declare_var("int", "wedge")

    Ef = KVASIR.add_declare_var("double","Ef", value=1.9)
    d_AD = KVASIR.add_declare_var("double","d_AD", value=1.2)
    d_SA = KVASIR.add_declare_var("double","d_SA", value=2.5)


    Lf = KVASIR.add_declare_var("double","Lf")
    KVASIR.append_initialize("Lf = d_AD + d_SA;")

    lamda_f = KVASIR.add_declare_var("double", "lambda_f") #final wavelenght
    KVASIR.append_initialize("lambda_f=1/(0.11056*sqrt(Ef));")


    #ANALYZER CALCULATION
    tA = KVASIR.add_declare_var("double", "tA") #Bragg angle of analyzer
    KVASIR.append_initialize("tA=asin(lambda_f/(2*d_PG))*RAD2DEG;")

    dtS = KVASIR.add_declare_var("double", "dtS", value=[*theta], array = N_analyzer) #

    len = KVASIR.add_declare_var("double", "len", value=l_analyzers, array = N_analyzer)

    #DETECTOR CALCULATION
    d_SDz = KVASIR.add_declare_var("double", "d_SDz") #z-component of distance between sample and detector
    KVASIR.append_initialize("d_SDz=d_SA+d_AD*cos(2*tA*DEG2RAD);")

    d_SDy = KVASIR.add_declare_var("double", "d_SDy") #heigt of detector relative arm Bragg 
    KVASIR.append_initialize("d_SDy=d_AD*sin(2*tA*DEG2RAD);")

    tDy= KVASIR.add_declare_var("double", "tDy")
    KVASIR.append_initialize("tDy=asin(width/(d_SDz))*RAD2DEG;")

    #### ADDING COMPONENTS ####
    print('GENERATING INSTRUMENT')
    ##SAMPLE##

    source = KVASIR.add_component(f"source", "MCPL_input", AT=[0,0,0], RELATIVE="ABSOLUTE")
    source.filename= f'"Storage/Virtual_output_{mode}.mcpl"'
    source.pos_smear = 0.005
    source.dir_smear = 0.05

    Arm_sample = KVASIR.add_component(f"Arm_sample", "Arm", AT=[0,0,0.58000], ROTATED = [0,f"tDy*t",0], RELATIVE="ABSOLUTE")

    if exp in ('incoherent', 'gen_prismatic_Ef'):
        sample = KVASIR.add_component(f"sample", "Incoherent", AT=[0,0,0], RELATIVE="Arm_sample")
        sample.radius=0.005
        sample.yheight=0.01
        sample.thickness=0.001
        sample.target_index = 1
        sample.focus_ah =20
        sample.focus_aw = 1

    if exp == 'powder':
        sample = KVASIR.add_component("sample", "PowderN", AT=[0,0,0], RELATIVE="Arm_sample")
        sample.reflections = '"Na2Ca3Al2F14.laz"'
        sample.radius = 0.005
        sample.yheight = 0.01 
        sample.thickness = 0.001
        sample.p_inc = 0
        sample.d_phi = 25
        sample.target_index = 1

    if exp == 'magnon':
        sample = KVASIR.add_component(f"sample", "Magnon_oI", AT=[0,0,0], ROTATED=[0,"90-A3",0], AT_RELATIVE="Arm_sample", ROTATED_RELATIVE = "ABSOLUTE")
        sample.radius=0.005
        sample.yheight=0.01
        sample.sigma_abs = 0
        sample.sigma_inc = 0
        sample.T = 0.05
        sample.a1 = 4.873
        sample.a2 = 4.873
        sample.a3 = 3.130
        sample.j = 0.00304
        sample.jc = -0.0056
        sample.ja = 0.0008
        sample.jb = 0.0008
        sample.S = 2.5
        sample.B = 0
        sample.D = -0.0023
        sample.FM = 0
        sample.mode_input = 2
        sample.e_steps_low = 50
        sample.e_steps_high = 50
        sample.target_index = 1
        sample.focus_yh ="1"
        sample.focus_xw = "0.03"#f"dtS[{int(N_analyzer/2)}]"
        sample.verbose = 0

    if exp == 'single_crystal':
        sample = KVASIR.add_component(f"sample","Single_crystal", AT = [0,0,0], ROTATED=[0,"90-A3",0], AT_RELATIVE=f"Arm_sample", ROTATED_RELATIVE = "ABSOLUTE")
        sample.xwidth=0.01
        sample.yheight=0.01
        sample.zdepth = 0.0001
        sample.mosaic = 30
        sample.p_transmit = 0.001
        sample.delta_d_d = "0.001/2.35"
        sample.reflections = '"ybgo.cif"'
        sample.sigma_inc = -1
        sample.sigma_abs = -1
        sample.order = 1
        sample.powder = 0

    
    target = KVASIR.add_component(f"target", "Arm", AT=[0,f"len[{int(N_analyzer/2)}]*sin(dtS[{int(N_analyzer/2)}]*DEG2RAD)",f"len[{int(N_analyzer/2)}]*cos(dtS[{int(N_analyzer/2)}]*DEG2RAD)"], RELATIVE="Arm_sample") # Arm at the center of analyser, target_index aims at it

    if verbrose > 0: #sanity checks
        Lmon = KVASIR.add_component(f"Lmon","L_monitor", AT = [0,0,"-u"], ROTATED=[0,0,0], RELATIVE="PREVIOUS")
        Lmon.xwidth = 0.25#f"2*len[{n}]*sin(tDy/2*DEG2RAD)"
        Lmon.yheight = 1.1
        Lmon.Lmin = 3
        Lmon.Lmax = 4
        Lmon.restore_neutron = 1

        fourpi = KVASIR.add_component(f"fourpi", "PSD_monitor_4PI", AT=[0,0,0], RELATIVE="Arm_sample")
        fourpi.radius = 0.5
        fourpi.restore_neutron =1

    ## Analyzer section

    for n in range(N_analyzer): # Placing and rotating each analyzer
        if theta[n] > (dt*N_analyzer+dt/2)*cent_frac:
            bz = KVASIR.add_declare_var("double", f"bz_{n}")
            KVASIR.append_initialize(f"bz_{n} = len[{n}]*cos(dtS[{n}]*DEG2RAD);")
            
            by = KVASIR.add_declare_var("double", f"by_{n}")
            KVASIR.append_initialize(f"by_{n} = len[{n}]*sin(dtS[{n}]*DEG2RAD);")

            az = KVASIR.add_declare_var("double", f"az_{n}")
            KVASIR.append_initialize(f"az_{n} = (bz_{n}-d_SDz);")
            
            ay = KVASIR.add_declare_var("double", f"ay_{n}")
            KVASIR.append_initialize(f"ay_{n} = (by_{n}+d_SDy);")
            
            dtA = KVASIR.add_declare_var("double", f"dtA_{n}")
            KVASIR.append_initialize(f"dtA_{n} = -acos( (az_{n}*bz_{n} + ay_{n}*by_{n}) / (len[{n}]* sqrt(az_{n}*az_{n} + ay_{n}*ay_{n}) ) )*RAD2DEG;")

        if theta[n] < (dt*N_analyzer+dt/2)*cent_frac and theta[n] > -(dt*N_analyzer+dt/2)*cent_frac:

            az = KVASIR.add_declare_var("double", f"az_{n}")
            KVASIR.append_initialize(f"az_{n} = len[{n}]*cos(dtS[{n}]*DEG2RAD)-(d_SDz);")
            
            ay = KVASIR.add_declare_var("double", f"ay_{n}")
            KVASIR.append_initialize(f"ay_{n} = len[{n}]*sin(dtS[{n}]*DEG2RAD)-d_SDy-{d_d_SDy_cent};")

            bz = KVASIR.add_declare_var("double", f"bz_{n}")
            KVASIR.append_initialize(f"bz_{n} = len[{n}]*cos(dtS[{n}]*DEG2RAD);")
            
            by = KVASIR.add_declare_var("double", f"by_{n}")
            KVASIR.append_initialize(f"by_{n} = len[{n}]*sin(dtS[{n}]*DEG2RAD);")
            
            dtA = KVASIR.add_declare_var("double", f"dtA_{n}")
            KVASIR.append_initialize(f"dtA_{n} = acos( (az_{n}*bz_{n} + ay_{n}*by_{n}) / (len[{n}]* sqrt(az_{n}*az_{n} + ay_{n}*ay_{n}) ) )*RAD2DEG;")
        
        if theta[n] < -(dt*N_analyzer+dt/2)*cent_frac:
            az = KVASIR.add_declare_var("double", f"az_{n}")
            KVASIR.append_initialize(f"az_{n} = (len[{n}]*cos(dtS[{n}]*DEG2RAD)-d_SDz);")
            
            ay = KVASIR.add_declare_var("double", f"ay_{n}")
            KVASIR.append_initialize(f"ay_{n} = (len[{n}]*sin(dtS[{n}]*DEG2RAD)-d_SDy);")

            bz = KVASIR.add_declare_var("double", f"bz_{n}")
            KVASIR.append_initialize(f"bz_{n} = len[{n}]*cos(dtS[{n}]*DEG2RAD);")
            
            by = KVASIR.add_declare_var("double", f"by_{n}")
            KVASIR.append_initialize(f"by_{n} = len[{n}]*sin(dtS[{n}]*DEG2RAD);")
            
            dtA = KVASIR.add_declare_var("double", f"dtA_{n}")
            KVASIR.append_initialize(f"dtA_{n} = acos( (az_{n}*bz_{n} + ay_{n}*by_{n}) / (len[{n}]* sqrt(az_{n}*az_{n} + ay_{n}*ay_{n}) ) )*RAD2DEG;")
        
        Arm = KVASIR.add_component(f"Arm_{n}", "Arm", AT=[0,0,0], ROTATED=[f"-dtS[{n}]",0,0], RELATIVE="Arm_sample") 
        Analyser = KVASIR.add_component(f"Analyzer_{n}","Single_crystal", AT = [0,0,f"len[{n}]"], ROTATED=[f"dtA_{n}/2",90,0], RELATIVE=f"Arm_{n}")
        Analyser.zdepth = f"2*len[{n}]*sin(tDy/2*DEG2RAD)"
        Analyser.yheight = width
        Analyser.xwidth = 0.0005
        Analyser.mosaic = 90
        Analyser.delta_d_d = "0.001/2.35"
        Analyser.reflections = '"C_graphite.lau"'
        Analyser.sigma_inc = -1
        Analyser.order = 1
        Analyser.append_EXTEND("if (SCATTERED) {wedge = 1;};")

        Analyser.ax = 0
        Analyser.ay = 2.14
        Analyser.az = -1.24

        Analyser.bx = 0
        Analyser.by = 0
        Analyser.bz = 2.47

        Analyser.cx = 6.71
        Analyser.cy = 0
        Analyser.cz = 0

        ## flagging neutrons - simulating potential collimators
        if theta[n] > (dt*N_analyzer+dt/2)*cent_frac:
            Analyser.append_EXTEND("if (SCATTERED) {wedge = 2;};")

        if theta[n] < (dt*N_analyzer+dt/2)*cent_frac and theta[n] > -(dt*N_analyzer+dt/2)*cent_frac:
            Analyser.append_EXTEND("if (SCATTERED) {wedge = 1;};")

        if theta[n] < -(dt*N_analyzer+dt/2)*cent_frac:
            Analyser.append_EXTEND("if (SCATTERED) {wedge = 0;};")

    ## Detector bank
    
    if exp == 'gen_prismatic_Ef':
        if mode == 'high_res':
            E_scale = {'upper':[1.84,1.98],'cent':[1.925,2.1],'lower':[1.84,1.98]}
        if mode == 'wide_dyn_range':
            E_scale = {'upper':[7.35,7.9],'cent':[7.7,8.4],'lower':[7.35,7.9]}

        KVASIR.add_component(f"arm_detector", "Arm", AT = [0,0,0], RELATIVE="Arm_sample")  
        det_He3 = KVASIR.add_component(f"det_He3_0_upper", "Monitor_nD", AT = [0,"d_SDy", d_SDz], AT_RELATIVE="PREVIOUS", ROTATED=[0,f"tDy*{n}",0], ROTATED_RELATIVE="Arm_sample")
        det_He3.options = f'"cylinder, y, limits=[-0.2,0.2], bins=80, energy, limits = {E_scale['upper']}, bins=500"' 
        det_He3.xwidth = width
        det_He3.yheight = "h_D"
        det_He3.filename = f'"Detector_upper.dat"'
        det_He3.set_WHEN("wedge == 0")

        det_He3 = KVASIR.add_component(f"det_He3_0_cent", "Monitor_nD", AT = [0,f"d_SDy+{d_d_SDy_cent}", f"d_SDz"], AT_RELATIVE=f"arm_detector", ROTATED=[0,f"tDy*{n}",0], ROTATED_RELATIVE="Arm_sample")
        det_He3.options = f'"cylinder, y, limits=[-0.2,0.2], bins=80, energy, limits = {E_scale['cent']}, bins=500"' 
        det_He3.xwidth = width
        det_He3.yheight = "h_D"
        det_He3.filename = f'"Detector_cent.dat"'
        det_He3.set_WHEN("wedge == 1")

        KVASIR.add_component(f"arm_detector_lower", "Arm", AT = [0,0,0], RELATIVE="Arm_sample")
        det_He3 = KVASIR.add_component(f"det_He3_0_lower", "Monitor_nD",  AT = [0,"-d_SDy", d_SDz], AT_RELATIVE="PREVIOUS", ROTATED_RELATIVE="Arm_sample")
        det_He3.options = f'"cylinder, y, limits=[-0.2,0.2], bins=80, energy, limits = {E_scale['lower']}, bins=500"'
        det_He3.xwidth = width
        det_He3.yheight = "h_D"
        det_He3.filename = f'"Detector_lower.dat"'
        det_He3.set_WHEN("wedge == 2")

    else:
        N_detector  =  24
        if mode == 'high_res':
            t_scale = {'upper':[0.2725,0.282],'cent':[0.264,0.277],'lower':[0.2725,0.282]}
        if mode == 'wide_dyn_range':
            t_scale = {'upper':[0.137,0.144],'cent':[0.134,0.14],'lower':[0.137,0.144]}

        for n in range(N_detector):
            KVASIR.add_component(f"arm_detector_{n}", "Arm", AT = [0,0,0], RELATIVE="Arm_sample", ROTATED=[0,f"tDy*{n}-tDy*{N_detector/2}",0])
            det_He3 = KVASIR.add_component(f"det_He3_{n}_ToF_upper", "Monitor_nD", AT = [0,"d_SDy", d_SDz], AT_RELATIVE="PREVIOUS", ROTATED=[0,f"tDy*{n}-tDy*{N_detector/2}",0], ROTATED_RELATIVE="Arm_sample")
            det_He3.options = f'"cylinder, y, limits=[-0.2,0.2], bins=80, t, limits =  {t_scale['upper']}, bins=500"' # t, limits = [0.2725,0.282]
            det_He3.xwidth = width
            det_He3.yheight = "h_D"
            det_He3.filename = f'"Detector_{n}_ToF_upper.dat"'
            det_He3.set_WHEN("wedge == 0")

            det_He3 = KVASIR.add_component(f"det_He3_{n}_ToF_cent", "Monitor_nD", AT = [0,f"d_SDy+{d_d_SDy_cent}", f"d_SDz"], AT_RELATIVE=f"arm_detector_{n}", ROTATED=[0,f"tDy*{n}-tDy*{N_detector/2}",0], ROTATED_RELATIVE="Arm_sample")
            det_He3.options = f'"cylinder, y, limits=[-0.2,0.2], bins=80, t, limits = {t_scale['cent']}, bins=500"' # t, limits = [0.264,0.277]
            det_He3.xwidth = width
            det_He3.yheight = "h_D"
            det_He3.filename = f'"Detector_{n}_ToF_cent.dat"'
            det_He3.set_WHEN("wedge == 1")

            KVASIR.add_component(f"arm_detector_{n}_lower", "Arm", AT = [0,0,0], RELATIVE="Arm_sample", ROTATED=[0,f"tDy*{n}-tDy*{N_detector/2}",0])
            det_He3 = KVASIR.add_component(f"det_He3_{n}_ToF_lower", "Monitor_nD", AT = [0,"-d_SDy", d_SDz], AT_RELATIVE="PREVIOUS", ROTATED=[0,f"tDy*{n}-tDy*{N_detector/2}",0], ROTATED_RELATIVE="Arm_sample")
            det_He3.options = f'"cylinder, y, limits=[-0.2,0.2], bins=80, t, limits =  {t_scale['lower']}, bins=500"' 
            det_He3.xwidth = width
            det_He3.yheight = "h_D"
            det_He3.set_WHEN("wedge == 2")
            det_He3.filename = f'"Detector_{n}_ToF_lower.dat"'


    #### SIMULATE EXPERIMENT ####
    print(f'RUNNING SIMULATION - {exp}')

    if exp == 'gen_prismatic_Ef':
        KVASIR.set_parameters(t=50)
        sample.SPLIT=100
        KVASIR.settings(ncount=1e7, mpi=6, output_path=f'Data_raw/signal_on_det_{mode}')
        data = KVASIR.backengine()

    elif exp in ('incoherent', 'powder'):
        for i in range(start,330):
            KVASIR.set_parameters(t=i)
            KVASIR.settings(ncount=1e5, mpi=6, output_path=f'Data_raw/{exp}_{mode}_n{i}')
            data = KVASIR.backengine()

    elif exp in ('magnon', 'single_crystal'):
        for j in range(90):
            for i in range(start,330):
                KVASIR.set_parameters(t=i)
                KVASIR.set_parameters(A3=j)
                KVASIR.settings(ncount=1e7, mpi=6, output_path=f'Data_raw/{exp}_{mode}_n{i}_A3{j}')
                data = KVASIR.backengine()




        



            



