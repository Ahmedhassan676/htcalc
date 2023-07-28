import streamlit as st 
from st_aggrid import AgGrid
import pandas as pd
import pandas as pd
import numpy as np
import streamlit as st
import ht
import openpyxl

from thermo import ChemicalConstantsPackage, PRMIX, CEOSLiquid, CEOSGas, FlashPureVLS,IAPWS95Gas,IAPWS95Liquid, GibbsExcessLiquid
from thermo.interaction_parameters import IPDB
from thermo.nrtl import NRTL
from physical_prop import *
@st.cache_data
def load_table():
    url ='http://raw.githubusercontent.com/Ahmedhassan676/htcalc/main/heat_table.csv'

    return pd.read_csv(url, index_col=[0])

@st.cache_data
def load_data_table():
    url ='http://raw.githubusercontent.com/Ahmedhassan676/htcalc/main/data_tables.csv'

    return pd.read_csv(url)
if "rating_table" not in st.session_state:
    st.session_state.rating_table = load_table().iloc[2:12,:]
@st.cache_data
def load_const_table():
    url ='http://raw.githubusercontent.com/Ahmedhassan676/htcalc/main/j_consts.csv'

    return pd.read_csv(url)  
calc_list = ['Duty','LMTD','Ft','Corrected LMTD','Surface Area','Tube Heat transfer Coef.','Shell Heat transfer Coef.','Uclean','Udirty','Uservice','Over Design','Over Surface','Shell Pressure Drop','Tube Pressure Drop','Shell Reynolds Number','Tube Reynolds Number','Tube Velocity','Shell Velocity']
geo_input_list = ['Shell D','Baffle Spacing','Number of baffles','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut']
para_input_list = ['Tube Flow rate','Tube inlet temperature','Tube outlet temperature','Tube Density','Tube Heat Capacity','Tube Viscosity', 'Tube Thermal conductivity','Tube Fouling factor','Shell flow rate','Shell inlet temperature','Shell outlet temperature','Shell Density','Shell Heat Capacity','Shell Viscosity', 'Shell Thermal conductivity','Shell Fouling factor'] 
       
j_const = load_const_table()
def Heat_balance(shell_side, Tube_list, Shell_list,s2,s3):
    m_t,t1_t,t2_t,rho_t,Cp_t,mu_t,k_t,fouling_t = Tube_list[0], Tube_list[1], Tube_list[2], Tube_list[3], Tube_list[4], Tube_list[5], Tube_list[6], Tube_list[7]
    m_s,t1_s,t2_s,rho_s,Cp_s,mu_s,k_s,fouling_s = Shell_list[0], Shell_list[1], Shell_list[2], Shell_list[3], Shell_list[4], Shell_list[5], Shell_list[6], Shell_list[7]
    if shell_side == 'Hot Side':
        T1 = t1_s  
        T2 =t2_s
        m_c = m_t
        m_h = m_s
        t1 = t1_t
        t2 = t2_t
        Cp_h = Cp_s
        Cp_c = Cp_t
    else:
        T1 = t1_t  
        T2 =t2_t
        m_c = m_s
        m_h = m_t
        t1 = t1_s
        t2 = t2_s
        Cp_h = Cp_t
        Cp_c = Cp_s
    if s2 == 'Hot side mass flow':
        Q = m_c * Cp_c * (t2-t1)
        m_h = Q/(Cp_h*(T1-T2))
    elif s2 == 'Hot side T1':
        Q = m_c * Cp_c * (t2-t1)
        T1 = T2 + (Q/m_h*Cp_h)
    elif s2 == 'Hot side T2':
        Q = m_c * Cp_c * (t2-t1)
        T2 = T1 - (Q/m_h*Cp_h)
    elif s2 == 'Cold side mass flow':
        Q = m_h * Cp_h * (T1-T2)
        m_c = Q/(Cp_c*(t2-t1)) 
    elif s2 == 'Cold side T1':
        Q = m_h * Cp_h * (T1-T2)
        t1 = t2 - (Q/m_c*Cp_c)
    else: #cold side T2
        Q = m_h * Cp_h * (T1-T2)
        t2 = t1 + (Q/m_c*Cp_c)
    dTlm = ht.LMTD(T1,T2,t1,t2)
    ft = ht.F_LMTD_Fakheri(t1,t2,T1,T2,s3)
    UA = Q/(dTlm*ft)
    Q = Q *1.163 # Kcal to W
    
    ntu_calc = ht.effectiveness_NTU_method(mh=m_h/3600, mc=m_c/3600, Cph=Cp_h*4184, Cpc=Cp_c*4184,subtype='S&T', Tci=t1, Thi=T1, UA=UA,n_shell_tube=s3)
    
    HB_data = [Q,dTlm,ft]  
    return HB_data,ntu_calc
def kern(Tube_list, Shell_list, geo_list,s3,HB_data,geo_input_df,calculations_df):
            m_t,t1_t,t2_t,rho_t,Cp_t,mu_t,k_t,fouling_t = Tube_list[0], Tube_list[1], Tube_list[2], Tube_list[3], Tube_list[4], Tube_list[5], Tube_list[6], Tube_list[7]
            m_s,t1_s,t2_s,rho_s,Cp_s,mu_s,k_s,fouling_s = Shell_list[0], Shell_list[1], Shell_list[2], Shell_list[3], Shell_list[4], Shell_list[5], Shell_list[6], Shell_list[7]
            Di,Do,tn,pn,L,tpitch,pitch,b_cut,shell_D,b_space = geo_list[3], geo_list[2], geo_list[0], geo_list[1], geo_list[6], geo_list[5], geo_list[4], geo_list[8], geo_list[-1], geo_list[7]
            Q, dTlm, ft = HB_data[0], HB_data[1], HB_data[2]
            if pitch == 'square' or 'rotated square 45':
              De = 4*(((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.25))/(3.14*Do*0.001)
            else:
              De = 8*(0.43301*((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.125))/(3.14*Do*0.001)
            
            C = tpitch-Do
            As = (shell_D*b_space*C)/(tpitch*1000)
            Gs = m_s/(As*3600)
            velocity_s = Gs/rho_s
            Res = (De*Gs)/mu_s
            f = np.exp(0.576-(0.19*np.log(Res)))
            Nb = (L/b_space)-1
            dp_s = ((f*(Gs**2)*(Nb+1)*shell_D)/(2*rho_s*De))*0.000010197
            L = L/1000
            A = np.pi*L*Do*0.001*s3*tn
            Cp_t = Cp_t*4184
            Cp_s = Cp_s*4184
            cross_A=(np.pi*0.25*(Di**2))*(tn/pn)
            velocity_t = m_t/(rho_t*3600*cross_A)
            Ret=(rho_t*velocity_t*Di)/mu_t
            f_t =1/(1.58*np.log(Ret)-3.28)**2 # valid for Re 2300 to 5,000,000 and Pr 0.5 to 2000
            port_1 = f_t*L*pn/Di
            port_2 = rho_t*(velocity_t**2)/2
            dp_t = (4*(port_1)+4*(pn))*port_2*0.000010197
            h_shell = (0.36*((De*Gs/mu_s)**0.55)*((Cp_s*mu_s/k_s)**(1/3)))*k_s/De #for Re between 2000 and 1,000,000
            Pr = Cp_t*mu_t/k_t
            Nu = ((0.5*f_t*(Ret-1000)*Pr))/(1+12.7*((0.5*f_t)**0.5)*((Pr**(2/3))-1)) # valid for Re 2300 to 5,000,000 (Gnielinski)
            h_t = Nu *k_t/Di
            d_ratio = Do/(Di*1000)
            Uc = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell))
            Ud = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell)+fouling_s+(d_ratio*fouling_t))
            U_calc = Q/(ft*dTlm*A)
            Rdesign = - (1/Uc) + (1/Ud)
            Rsevice = - (1/Uc) + (1/U_calc)
            OD_k=100*((Ud-U_calc)/Ud)
            OV = ((Uc/U_calc)-1)*100
            geo_input_df.loc['Number of baffles','Kern_summary'] = Nb
            calculations_df.loc[['Surface Area','Tube Heat transfer Coef.','Shell Heat transfer Coef.','Uclean','Udirty','Uservice','Shell Pressure Drop','Tube Pressure Drop','Shell Reynolds Number','Tube Reynolds Number','Tube Velocity','Shell Velocity','Over Design','Over Surface'],'Kern_summary']= A,h_t,h_shell,Uc,Ud,U_calc,dp_s,dp_t,Res,Ret,velocity_t,velocity_s,OD_k, OV
            return dp_s, dp_t, h_shell, h_t, Uc, Ud, U_calc, Rdesign, Rsevice

def bell_delaware(Tube_list, Shell_list ,h_t,h_shell,geo_list,s3,HB_data, geo_input_df, calculations_df):
    no_of_shells = s3
    m_t,t1_t,t2_t,rho_t,Cp_t,mu_t,k_t,fouling_t = Tube_list[0], Tube_list[1], Tube_list[2], Tube_list[3], Tube_list[4], Tube_list[5], Tube_list[6], Tube_list[7]
    m_s,t1_s,t2_s,rho_s,Cp_s,mu_s,k_s,fouling_s = Shell_list[0], Shell_list[1], Shell_list[2], Shell_list[3], Shell_list[4], Shell_list[5], Shell_list[6], Shell_list[7]
    Di,Do,tn,pn,L,tpitch,pitch,b_cut,shell_D,b_space = geo_list[3], geo_list[2], geo_list[0], geo_list[1], geo_list[6]/1000, geo_list[5], geo_list[4], geo_list[8], geo_list[-1], geo_list[7]
    # Tube pitch layout
    if pitch == 'square':
      t_p_angle = 90
    elif pitch == 'rotated square 45':
      t_p_angle = 45
    else:
      t_p_angle = 30
    p_ratio = [1.25,1.285,1.33,1.5]
    
    T_wall = (t1_t+t2_t)*0.5+h_shell*((t2_t+t2_s)*0.5-(t1_t+t2_t)*0.5)/(h_shell+h_t)
    Di = Di*1000
    print('dp for T_wall '+str(T_wall))
    mu_s_w = mu_s # assuming negligible T_wall effects for shell
    mu_t_w = mu_t # assuming negligible T_wall effects for tubes
    shell_D = shell_D*1000
  
    k_s = k_s/1.16
    k_t = k_t /1.16
    fouling_s = fouling_s*1.16
    fouling_t = fouling_t*1.16
    #Shell outer tube limit
    D_otl = shell_D - (12.5+(shell_D/200))
    print('dp for D_otl '+str(D_otl))
    #Height of baffle cut
    L_c = b_cut/(100*shell_D)
    # Diametral Shell baffle clearance
    D_sb = 3.1+0.004*shell_D
    # Tubesheet thickness
    L_s = 0.1*shell_D
    Lb_cut = b_space # central baffle spacing
    LB_in = b_space # inlet baffle spacing
    LB_out = b_space # Outlet baffle spacing
    L_tb = 0.4 # Diametral Tube-Baffle Clearance
    
    # tube arrangement
    if t_p_angle == 30:
      t_arrg = np.sqrt(3)/2
    elif t_p_angle == 45:
      t_arrg = 1/np.sqrt(2)
    else : t_arrg = 1


    #LMTD and Q
    LMTD = HB_data[1] 
    f_t = HB_data[2]
    corrected_LMTD = LMTD * f_t
    Q = HB_data[0]

    #Tube wall temperature
    T_wall = (t1_t+t2_t)*0.5+h_shell*((t1_s+t2_s)*0.5-(t1_t+t2_t)*0.5)/(h_shell+h_t)


    # Shell Side heat transfer coefficient
    t_p = tpitch #p_ratio *Do
    #t_p_effective 
    if t_p_angle == 45:
      t_p_effective = t_p/(2**0.5)
    else: t_p_effective = t_p

    #cross flow area 
    S_m = (Lb_cut/1000)*((shell_D-D_otl)+(D_otl-Do)*(t_p-Do)/t_p_effective)/1000
    Gs = (m_s/3600)/S_m
    Re_s = (Do/1000)*Gs/(mu_s*0.001)
    print('dp for Re_s '+str(Re_s))
    Pr_s = (mu_s/1000)*Cp_s/k_s*3600
    # j- ideal factor coefficients
    a1 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s < j_const['Reynolds_max']) & (j_const['Layout'] == t_p_angle) ,j_const['a_{1}'],0).sum()
    a2 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s < j_const['Reynolds_max']) & (j_const['Layout'] == t_p_angle) ,j_const['a_{2}'],0).sum()
    a3 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s < j_const['Reynolds_max']) & (j_const['Layout'] == t_p_angle) ,j_const['a_{3}'],0).sum()
    a4 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s < j_const['Reynolds_max']) & (j_const['Layout'] == t_p_angle) ,j_const['a_{4}'],0).sum()
    a=a3/(1+0.14*((Re_s)**a4))
    #print(a)
    #print(Re_s)
    j=a1*((1.33/(t_p/Do))**a)*(Re_s**a2)
    print('dp for j '+str(j))
    #print(j)
    h_s_i = Cp_s*(m_s/S_m)*j*(Pr_s**(-2/3))*((mu_s/mu_s_w)**0.14)
    print(h_s_i)
    # 1. correction factor for baffle window flow

    theta_CTL = 2*np.arccos(shell_D*(1-2*b_cut/100)/(D_otl-Do)) #b_cut baffle cut, shell_D isnide shell diameter mm
    F_w = (theta_CTL-np.sin(theta_CTL))/(2*np.pi)

    F_c = 1-2*F_w # Fraction of tubes in cross flow
    j_c = 0.55+0.72*F_c

    # 2. correction factor for baffle leakage
    theta_Ds = 2*np.arccos(1-(2*b_cut)/100) # Baffle window angle
    S_sb = (shell_D/1000)*(D_sb/1000)*(np.pi-0.5*theta_Ds)    # Shell to baffle leakage area
    S_tb = (np.pi/4)*((((Do/1000)+(L_tb/1000))**2)-((Do/1000)**2))*tn*(1-F_w)  # Tube to baffle leakage area

    r_L = (S_sb+S_tb)/S_m # ratio of leakage to cross flow
    r_S = S_sb/(S_sb+S_tb) # ratio of shell-baffle to total area

    j_l = 0.44*(1-r_S)+(1-0.44*(1-r_S))*np.exp(-2.2*r_L) # Baffle leakage correction factor

    # 3. correction factor for bundle bypass
     # number_of_sealing_strips   
    P_p = t_p * t_arrg # Tube row distance in flow direction
    N_TCC = (shell_D/P_p)*(1-(2*b_cut)/100) #N_TCC number of tube rows between baffle
    N_ss = int(N_TCC/6)
    r_ss = N_ss/N_TCC # Nss Number of sealing strips
    print('N_ss is '+str(N_ss)+' While N_Tcc ia '+str(N_TCC))
    S_b = (Lb_cut/1000)*(shell_D-D_otl-(Do/2))/1000 #bundle pybass area
    if Re_s < 100:
      C_j = 1.35
    else: C_j = 1.25 # Correlation constant

    if r_ss >= 0.5:
      j_b = 1
    else: j_b = np.exp(-1*C_j*(S_b/S_m)*(1-((2*r_ss)**(1/3))))

    # 4. Correction factor for adverse temperature gradient
    N_tcw = (0.8/P_p)*(shell_D*b_cut/100-(shell_D-(D_otl-Do))/2)
    N_b = 1 +int((L-(2*L_s*0.001)-(LB_in+LB_out)*0.001)/(Lb_cut*0.001)) # number of baffles
    print('value for N_b '+str(N_b))
    N_c = (N_tcw +N_TCC)*(1+N_b) # tube rows crossed in entire exchanger
    j_RL = (10/N_c)**0.18
    if Re_s <= 20:
      j_R = j_RL
    elif Re_s <100:
      j_R = j_RL+((20-Re_s)/80)*(j_RL-1)
    else: j_R = 1

    # 5. correction factor for unequal baffle spacing
    if Re_s < 100:
      n1 = 1/3
    else: n1 = 0.6
    j_s =((N_b-1)+(LB_in/Lb_cut)**(1-n1)+(LB_out/Lb_cut)**(1-n1))/((N_b-1)+(LB_in/Lb_cut)+(LB_out/Lb_cut))


    h_shell = j_s * j_R *j_b * j_l * j_c * h_s_i
    print(h_shell)

    ### Tube Side Heat transfer coeficient
    a_tube = (np.pi*((Di/1000)**2)*tn)/(4*pn) # Flow area
    v_t = (m_t/3600)/(rho_t*a_tube) #velocity through tube
    print('value for v_t '+str(v_t))
    Re_t = (Di/1000)*rho_t*v_t/(mu_t/1000)
    print('value for Re_t '+str(Re_t))
    Pr_t = Cp_t*(mu_t/1000)/k_t*3600
    print('value for Pr_t '+str(Pr_t))
    L_eff = L-2*L_s/1000 # Effective tube length

    # Nusselt Number Calculation
    Nu_laminar = 1.86*(Re_t*Pr_t*(Di/1000)/L_eff)**(1/3)
    # Turbulent flow Petukhov-Kirillov
    f_turbulent = (1.58*np.log(Re_t)-3.28)**-2
    print('value for f_turbulent '+str(f_turbulent))
    Nu_turb = (f_turbulent/2)*Re_t*Pr_t/(1.07+12.7*((f_turbulent/2)**0.5)*((Pr_t**(2/3))-1))
    # Transition flow Nu
    Nu_2300 = 1.86*(2300*Pr_t*(Di/1000)/L_eff)**(1/3)

    f_Re_10000 = (1.58*np.log(10000)-3.28)**-2
    Nu_10000 = (f_Re_10000/2)*10000*Pr_t/(1.07+12.7*((f_Re_10000/2)**0.5)*((Pr_t**(2/3))-1))

    Nu_trans = Nu_2300+(Nu_10000-Nu_2300)*(Re_t-2300)/(10000-2300)


    if Re_t <= 2300:
      Nu_tube = Nu_laminar
    elif Re_t < 10000:
      Nu_tube = Nu_trans
    else: Nu_tube = Nu_turb
    print('value for Nu_tube '+str(Nu_tube))
    h_t_i=Nu_tube*k_t/(Di/1000)*(mu_t/mu_t_w)**0.14
    print('value for h_t_i '+str(h_t_i))

    # Overall Heat tansfer coefficient
    dict_of_conductivity = {'Carbon Steel':38.69,'Copper':324.42,'Inconel':12.95,'Monel':21.28,'Nickel':52.09,'Stainless Steel':13.54}
    k_w_t = dict_of_conductivity['Carbon Steel']
    T_wall = (t1_t+t2_t)/2+h_s_i*((t1_s+t2_s)/2-(t1_t+t2_t)/2)/(h_s_i+h_t_i)
    wall_resistance = (Do/2000)*np.log(Do/Di)/k_w_t


    U_clean = 1/((1/h_shell)+(Do/(h_t_i*Di))+wall_resistance)
    print('dp for U_clean '+str(U_clean))
    U_dirty = 1/((1/U_clean)+fouling_t+fouling_s)
    print('dp for U_dirty '+str(U_dirty))
    ### Shell side pressure drop
    b1 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s < j_const['Reynolds_max']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{1}'],0).sum()
    b2 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s < j_const['Reynolds_max']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{2}'],0).sum()
    b3 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s < j_const['Reynolds_max']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{3}'],0).sum()
    b4 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s < j_const['Reynolds_max']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{4}'],0).sum()
    b = b3/(1+0.14*(Re_s**b4))
    f_s = b1*((1.33/(t_p/Do))**b)*(Re_s**b2)
    dp_shell_ideal = 2*f_s*(Gs**2)*N_TCC*(mu_s_w/mu_s)**0.14/(rho_s)/100000 #N_TCC number of tube rows between baffle
    #print(dp_shell_ideal)

    # 1. correction factor for baffle leakage
    f_p = 0.8-0.15*(1+r_S) # r_s ratio of shell-baffle to total area
    R_L = np.exp(-1.33*(1+r_S)*(r_L**f_p))# r_l ratio of  leakage area to cross flow

    # pressure drop for an ideal window section
    S_wg = (((shell_D/1000)**2)/8)*(theta_Ds-np.sin(theta_Ds)) #Gross window area
    S_wt = tn*F_w*np.pi*((Do/1000)**2)/4 #Window area occupied with tubes
    S_w = S_wg - S_wt # Net Cross flow area
    G_w = (m_s/3600)/np.sqrt(S_m*S_w)
    v_w = G_w/rho_s
    print('dp for G_w '+str(G_w))
    D_w = 4*S_w/(np.pi*(Do/1000)*tn*F_w+(shell_D/1000)*theta_Ds)

    #pressure drop for turbluent flow in ideal window section
    if Re_s >= 100:
      dp_window = N_b*R_L*(2+0.6*N_tcw)*(G_w**2)/(2*rho_s)/100000

      
    else:
      dp_window = N_b*R_L*((26*G_w*Cp_s/rho_s)*(N_tcw/(tpitch-Do)/1000 + (Lb_cut/1000)/(D_w**2)) + (G_w**2)/rho_s)/100000
        

    # 2. Correction factor for bundle bypass effect
    if Re_s < 100:
      C_r = 4.5
    else: C_r = 3.7

    if r_ss >= 0.5:
      R_b = 1
    else:
      R_b = np.exp(-1*C_r*(S_b/S_m)*(1-(2*r_ss)**(1/3)))

    # 3. Correction for unequal baffle spacing inlet/ outlet
    if Re_s < 100:
      n = 1 # Slope of friction factor curve
    else: n = 0.2
    R_s = (Lb_cut/LB_in)**(2-n) + (Lb_cut/LB_out)**(2-n) # correction factor

    # Shell side pressure drop (Excluding nozzles)
    dp_cent_baff = dp_shell_ideal*(N_b-1)*R_L*R_b # pressure drop in baffle windows
    dp_baff_window = dp_window # pressure drop in baffle windows
    dp_entrance_exit = dp_shell_ideal*R_s*R_b*(1+N_tcw/N_TCC) # pressure drop in entrance / exit baffles
    print('dp for cent baffle '+str(dp_cent_baff))
    print('dp for baffle wind '+str(dp_baff_window))
    print('dp for entrance '+str(dp_entrance_exit))

    total_dp_shell = dp_cent_baff + dp_baff_window + dp_entrance_exit
    print('overall shell dp is {}'.format(total_dp_shell))
    # Tube side pressure drop
    total_dp_tube =((4*f_turbulent*L*pn/(Di/1000))+4*pn)*rho_t*(v_t**2)/2/100000
    print('overall Di is {}'.format(Di))
    ### Area Calculations
    print('dp for corrected_LMTD '+str(corrected_LMTD))
    print('dp for LMTD '+str(LMTD))
    print(corrected_LMTD/LMTD)
    A_required = Q/(corrected_LMTD*U_dirty*1.163)
    A_available = L_eff * tn *np.pi*(Do/1000)*no_of_shells
    U_required = Q/(corrected_LMTD*A_available*1.163*no_of_shells)
    A_over = ((U_clean/U_required)-1)*100
    OD = ((U_dirty/U_required)-1)*100
    print('dp for A_required '+str(A_required))
    print('dp for A_available '+str(A_available))
    print('dp for U_required '+str(U_required))
    print('od is {} while OV is {}'.format(OD,A_over))
    geo_input_df.loc[['Number of baffles','Length'],'Bell_summary'] = N_b,L_eff*1000
    calculations_df.loc[['Surface Area','Tube Heat transfer Coef.','Shell Heat transfer Coef.','Uclean','Udirty','Uservice','Shell Pressure Drop','Tube Pressure Drop','Shell Reynolds Number','Tube Reynolds Number','Tube Velocity','Shell Velocity','Over Design','Over Surface'],'Bell_summary']= A_available,h_t_i,h_s_i,U_clean,U_dirty,U_required,total_dp_shell,total_dp_tube,Re_s,Re_t,v_t,v_w,OD,A_over
            
    return U_clean,U_dirty,U_required,OD,total_dp_shell,total_dp_tube


def main():
    html_temp="""
    <div style="background-color:lightblue;padding:16px">
    <h2 style="color:black"; text-align:center> Heat Exchangers Calculation </h2>
    </div>
    
        """
   
    st.markdown(html_temp, unsafe_allow_html=True)
    if 'calculations_df' not in st.session_state:
        st.session_state.calculations_df = pd.DataFrame(index=calc_list)
    if 'geo_input_df' not in st.session_state:
        st.session_state.geo_input_df = pd.DataFrame(index=geo_input_list)
    if 'para_input_df' not in st.session_state:
        st.session_state.para_input_df = pd.DataFrame(index=para_input_list)  
    s1 = st.selectbox('Select Calculations required',('Heat Exchanger Assessment','HEx Rating from a TEMA datasheet','Prelaminary Design','Perform Trials'), key = 'type')
    if s1 == 'Heat Exchanger Assessment':
        wizard_form_header()
        st.markdown('---')
        if st.session_state['current_step'] == 1:
         
            s_prop = st.selectbox('Estimate Shell & Tube Fluids properties?',('No','Yes'), key = 'prop') 
            if "rating_var" not in st.session_state:
                st.session_state.rating_var = pd.DataFrame() 
            if s_prop == 'No':
                st.session_state.rating_var = st.data_editor(st.session_state.rating_table)
            else:
                st.session_state.rating_var = st.data_editor(main_prop())
        if st.session_state['current_step'] == 2:
            st.session_state.rating_var = st.data_editor(st.session_state.rating_var)
        if st.session_state['current_step'] == 3:
            rating_df = st.session_state.rating_var
            s2 = st.selectbox('Select Heat Balance variable',('Hot side mass flow','Hot side T1','Hot side T2','Cold side mass flow','Cold side T1','Cold side T2'), key = 'HB') 
            s3 = st.selectbox('Number of Shells',(1,2,3,4,5,6,7,8), key='shells')
            if 's3' not in st.session_state:
                 st.session_state.s3 = s3
            if 'dp_calc_check' not in st.session_state:
                st.session_state.dp_calc_check = st.checkbox("Calculate pressure drop?")
            else: st.session_state.dp_calc_check = st.checkbox("Calculate pressure drop?", value=st.session_state.dp_calc_check)
            shell_side = st.selectbox('Shell Side is the..?',('Cold Side','Hot Side'), key = 'shell_side')
            
            try:
              t1_s = float(rating_df.iloc[1,1])
              t2_s = float(rating_df.iloc[2,1] )
              m_s = float(rating_df.iloc[0,1])
              Cp_s = float(rating_df.iloc[6,1]) 
              mu_s = float(rating_df.iloc[7,1])/1000
              rho_s =  float(rating_df.iloc[5,1])
              mu_s = float(rating_df.iloc[7,1])/1000
              k_s = float(rating_df.iloc[8,1])
              fouling_s = float(rating_df.iloc[9,1])
              mu_t = float(rating_df.iloc[7,2])/1000
              fouling_t = float(rating_df.iloc[9,2])
              rho_t =  float(rating_df.iloc[5,1])
              m_t = float(rating_df.iloc[0,2])
              t1_t = float(rating_df.iloc[1,2]) 
              t2_t = float(rating_df.iloc[2,2]) 
              Cp_t = float(rating_df.iloc[6,2])
              k_t = float(rating_df.iloc[8,2]) 
              
              Shell_list = [m_s, t1_s, t2_s, rho_s, Cp_s, mu_s, k_s, fouling_s]
              Tube_list = [m_t, t1_t, t2_t, rho_t, Cp_t, mu_t, k_t, fouling_t]
              
          
              
              
              HB_data,ntu_calc = Heat_balance(shell_side, Tube_list, Shell_list,s2,s3)
              if 'HB_data' not in st.session_state:
                 st.session_state.HB_data = HB_data
              if 'ntu_calc' not in st.session_state:
                st.session_state.ntu_calc = ntu_calc
              Q, dTlm, ft = HB_data[0], HB_data[1], HB_data[2]
              
                
            except (UnboundLocalError,IndexError,ZeroDivisionError): pass
            #except IndexError: pass
            if not st.session_state.dp_calc_check:
                A = st.number_input('Total Heat Exchanger(s) Area', key = 'a')
                U = st.number_input('Service U Kcal/hr.m2.C', key = 'U')
                try:
                    U_calc = Q/(ft*dTlm*A)
                    st.session_state.calculations_df.loc[['Surface Area','Udirty','Uservice','Over Design'],'summary'] = A, U, U_calc, 100*(U-U_calc)/U
                    st.session_state.para_input_df.loc[:,'summary'] = [m_t, t1_t, t2_t, rho_t, Cp_t, mu_t*1000, k_t, fouling_t,m_s, t1_s, t2_s, rho_s, Cp_s, mu_s*1000, k_s, fouling_s]
                    st.session_state.calculations_df.loc[['Duty','LMTD','Ft','Corrected LMTD'],'summary'] = Q,dTlm,ft,dTlm*ft
                    
                except (ZeroDivisionError,UnboundLocalError): pass
            def get_index(series, value):
                      n = 0 
                      series = series.reset_index()
                      del series['index']
                      try:
                        n = int(series[series.iloc[:,0]==str(value)].index[0])
                      except IndexError:
                        
                        n = int(series[series.iloc[:,0]==str(int(value*1000))].index[0]) 
                      return n    
            if st.session_state.dp_calc_check:
                st.session_state.para_input_df.loc[:,'Kern_summary'] = st.session_state.para_input_df.loc[:,'Bell_summary'] = [m_t, t1_t, t2_t, rho_t, Cp_t, mu_t*1000, k_t, fouling_t,m_s, t1_s, t2_s, rho_s, Cp_s, mu_s*1000, k_s, fouling_s]
                st.session_state.calculations_df.loc[['Duty','LMTD','Ft','Corrected LMTD'],'Kern_summary'] = st.session_state.calculations_df.loc[['Duty','LMTD','Ft','Corrected LMTD'],'Bell_summary'] = Q,dTlm,ft,dTlm*ft
                geo_table = load_table().iloc[26:,:2].rename(columns={'Shell Fluid':'Value'})
                tube_table = load_data_table().iloc[1:11,1:5]
                thickness_table = load_data_table().iloc[11:36,1:4]
                shell_table = load_data_table().iloc[37:67,1]
                pitch_options = ['square','rotated square 45','triangle 30','triangle 60']
                if 'Do_ind' not in st.session_state:
                   st.session_state.Do_ind,st.session_state.thick_ind,st.session_state.pitch_ind,st.session_state.shell_D_ind = 0,0,0,0
                 
                Do = float(st.selectbox('tube OD (mm)?',tube_table.iloc[1:10,0],index=st.session_state.Do_ind, key = 'Do_st'))
                st.session_state.Do_ind = get_index(tube_table.iloc[1:10,0],Do)
                
                
                thick = float(st.selectbox('tube gauge (thickness)?',thickness_table.iloc[1:25,2],index=st.session_state.thick_ind, key = 'thick_st'))
                st.session_state.thick_ind = get_index(thickness_table.iloc[1:25,2],thick)
                Di = (Do - 2*thick)*0.001
                
                pitch = st.selectbox('pitch type?',pitch_options,index=st.session_state.pitch_ind, key = 'pitch_st')
                st.session_state.pitch_ind = pitch_options.index(pitch)
                shell_D = float(st.selectbox('Shell diameter (mm)?',shell_table, index=st.session_state.shell_D_ind,key = 'Shell ID_st'))/1000
                st.session_state.shell_D_ind = get_index(shell_table,shell_D)
                
                try:
                    if 'geo_df' not in st.session_state:
                       st.session_state.geo_df = geo_table.iloc[[0,1,4,6,7,8],:]
                    
                    st.session_state.geo_df = st.data_editor(st.session_state.geo_df)
                    tn = float(st.session_state.geo_df.loc['Number of tubes','Value'])
                    pn = float(st.session_state.geo_df.loc['number of passes','Value'])
                    #Do = float(st.session_state.geo_df.iloc[2,1])
                    #Di = (Do - 2*float(st.session_state.geo_df.iloc[3,1]))*0.001
                    L = float(st.session_state.geo_df.loc['Tube length','Value'])
                    tpitch = float(st.session_state.geo_df.loc['pitch','Value'])
                    b_space = float(st.session_state.geo_df.loc['baffle spacing','Value'])
                    b_cut = float(st.session_state.geo_df.loc['baffle cut','Value'])
                    geo_list = [tn ,pn,Do, Di, pitch, tpitch,L, b_space, b_cut,shell_D]

                except ValueError: pass
                try:
                    st.session_state.geo_input_df.loc[['Shell D','Baffle Spacing','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut'],'Kern_summary']=st.session_state.geo_input_df.loc[['Shell D','Baffle Spacing','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut'],'Bell_summary']=shell_D,b_space,Do,Di*1000,L,tn,pn,tpitch,pitch,b_cut
                    dp_s, dp_t, h_shell, h_t, Uc, Ud, U_calc, Rdesign, Rsevice = kern(Tube_list, Shell_list, geo_list,s3,HB_data,st.session_state.geo_input_df,st.session_state.calculations_df)

                    #L = L*1000
                    geo_list = [tn ,pn,Do, Di, pitch, tpitch,L, b_space, b_cut,shell_D]
                    Shell_list = [m_s, t1_s, t2_s, rho_s, Cp_s, mu_s*1000, k_s, fouling_s]
                    Tube_list = [m_t, t1_t, t2_t, rho_t, Cp_t, mu_t*1000, k_t, fouling_t]
                    
                    U_clean,U_dirty,U_required,OD,total_dp_shell,total_dp_tube=bell_delaware(Tube_list, Shell_list ,h_t,h_shell,geo_list,s3,HB_data,st.session_state.geo_input_df,st.session_state.calculations_df)
                except UnboundLocalError: pass 
                except ValueError: pass
        if st.session_state['current_step'] == 4:
            if 'submitted' not in st.session_state:
                st.session_state.submitted = False
            submit(st.session_state.submitted,st.session_state.ntu_calc)
            
        if st.session_state['current_step'] == 5:
            shell_table = load_data_table().iloc[37:67,1]
            options_list = ['Shell D','Baffle Spacing','Do','Length','Number of tubes','Number of passes','Tube pitch','pitch type']
            opt_dict = {}
            #for i in options_list:
            #    opt_dict[i] = None
            trials_options = st.multiselect('select trials basis',options_list) 
            with st.form("my_form"):
              
              try:
                  if 'Shell D' in trials_options:
                    shell_D =  float(st.selectbox('Shell diameter (mm)?',shell_table, key = 'ShellD trial'))/1000  
                    opt_dict['Shell D'] = shell_D
                  #else: shell_D=st.session_state.geo_input_df.loc['Shell D','Kern_summary']
                  if 'Baffle Spacing' in trials_options:
                    b_space = st.number_input('Input baffle Spacing (mm)', key='B_space')
                    opt_dict['Baffle Spacing'] = b_space
                  #else: b_space=st.session_state.geo_input_df.loc['Shell D','Kern_summary']
                  if 'Length' in trials_options:
                    L = float(st.number_input('Input Length (mm)', key='tubelength'))  
                    opt_dict['Length'] = L 
                  if 'Number of tubes' in trials_options:
                    tn = float(st.number_input('Input number of tubes', key='tn')) 
                    opt_dict['Number of tubes'] = tn
                  if 'Number of passes' in trials_options:
                    pn = float(st.number_input('Input Number of passes', key='pn'))
                    opt_dict['Number of passes'] = pn
                  if 'Tube pitch' in trials_options:
                    tpitch = float(st.number_input('Input Tube pitch (mm)', key='pn'))  
                    opt_dict['Tube pitch'] = tpitch   
                  if 'pitch type' in trials_options:
                    pitch = st.selectbox('pitch type?',('square','rotated square 45','triangle 30','triangle 60'), key = 'pitch_trials') 
                    opt_dict['pitch type'] = pitch 
              except UnboundLocalError: pass
              trials_bttn = st.form_submit_button("Try!")    
    
                    
              if trials_bttn:
                
                  #st.session_state.para_input_df.loc[:,'Kern_summary'] = st.session_state.para_input_df.loc[:,'Bell_summary'] = [m_t, t1_t, t2_t, rho_t, Cp_t, mu_t*1000, k_t, fouling_t,m_s, t1_s, t2_s, rho_s, Cp_s, mu_s*1000, k_s, fouling_s]
                  #st.session_state.calculations_df.loc[['Duty','LMTD','Ft','Corrected LMTD'],'Kern_summary'] = st.session_state.calculations_df.loc[['Duty','LMTD','Ft','Corrected LMTD'],'Bell_summary'] = Q,dTlm,ft,dTlm*ft
                  for i in opt_dict.keys():
                    if i in st.session_state.geo_input_df.index:
                        st.session_state.geo_input_df.loc[i,'Kern_summary']=st.session_state.geo_input_df.loc[i,'Bell_summary']=opt_dict[i]
                  geo_list =list(st.session_state.geo_input_df.loc[['Number of tubes','Number of passes','Do','Di','pitch type','Tube pitch','Length','Baffle Spacing','baffle cut','Shell D'],'Kern_summary'].values) #[tn ,pn,Do, Di, pitch, tpitch,L, b_space, b_cut,shell_D]
                  Shell_list = list(st.session_state.para_input_df.iloc[8:16,0].values)
                  Tube_list = list(st.session_state.para_input_df.iloc[:8,0].values)
                  #inside tube diameter
                  geo_list[3]=geo_list[3]/1000
                  #shell diameter
                  #geo_list[-1]=geo_list[-1]*1000
                  Shell_list[5]=Shell_list[5]/1000
                  Tube_list[5]=Tube_list[5]/1000
                  st.write(Tube_list, Shell_list, geo_list)
                  #st.session_state.geo_input_df.loc[['Shell D','Baffle Spacing','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type'],'Kern_summary']=st.session_state.geo_input_df.loc[['Shell D','Baffle Spacing','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type'],'Bell_summary']=shell_D,b_space,Do,Di*1000,L,tn,pn,tpitch,pitch
                  dp_s, dp_t, h_shell, h_t, Uc, Ud, U_calc, Rdesign, Rsevice = kern(Tube_list, Shell_list, geo_list,st.session_state.s3,st.session_state.HB_data,st.session_state.geo_input_df,st.session_state.calculations_df)
                  Shell_list[5]=Shell_list[5]*1000
                  Tube_list[5]=Tube_list[5]*1000
                  st.write(Tube_list, Shell_list, geo_list)
                  U_clean,U_dirty,U_required,OD,total_dp_shell,total_dp_tube=bell_delaware(Tube_list, Shell_list ,h_t,h_shell,geo_list,st.session_state.s3,st.session_state.HB_data,st.session_state.geo_input_df,st.session_state.calculations_df)
                  st.session_state.summary = pd.concat([st.session_state.calculations_df,st.session_state.para_input_df,st.session_state.geo_input_df])
                  st.session_state.summary['Kern_summary'] = st.session_state.summary['Kern_summary'].apply(lambda x: convert_to_float_or_string(x))
                  st.session_state.summary['Bell_summary'] = st.session_state.summary['Bell_summary'].apply(lambda x: convert_to_float_or_string(x))
                  st.write(st.session_state.summary)
        st.markdown('---')
        st.session_state.submitted = wizard_form_footer()      
    elif s1 == 'HEx Rating from a TEMA datasheet':      
     
      try:
          st.session_state.calculations_df = pd.DataFrame(index=calc_list)
          st.session_state.geo_input_df = pd.DataFrame(index=geo_input_list)
          st.session_state.para_input_df = pd.DataFrame(index=para_input_list) 
          uploaded_file_power = st.file_uploader('Choose a file', key = 2)
          shell_side = st.selectbox('Shell Side is the..?',('Cold Side','Hot Side'), key = 'shell_side')
          s2 = st.selectbox('Select Heat Balance variable',('Hot side mass flow','Hot side T1','Hot side T2','Cold side mass flow','Cold side T1','Cold side T2'), key = 'HB') 
          s3 = st.selectbox('Number of Shells',(1,2,3,4,5,6,7,8), key='shells')
          
          if uploaded_file_power:
              workbook= openpyxl.load_workbook(uploaded_file_power, data_only=True)
              try:
                    thickness_table = load_data_table().iloc[11:36,1:4]
                    thickness_table.columns = thickness_table.iloc[0]
                    thickness_table = thickness_table[1:]
                    worksheet = workbook['Sheet1']     
                    s3 = worksheet['I10'].value
                    t1_s = worksheet['H20'].value
                    t2_s =worksheet['I20'].value
                    m_s = worksheet['H14'].value
                    Cp_s = (worksheet['H24'].value+worksheet['I24'].value)/2
                    mu_s = (worksheet['H22'].value+worksheet['I22'].value)*0.001/2
                    rho_s =  (worksheet['H21'].value+worksheet['I21'].value)/2
                    k_s = (worksheet['H25'].value+worksheet['I25'].value)/2
                    fouling_s = worksheet['H30'].value
                    mu_t = (worksheet['K22'].value+worksheet['L22'].value)*0.001/2
                    fouling_t =  worksheet['K30'].value
                    rho_t =  (worksheet['K21'].value+worksheet['L21'].value)/2
                    m_t = worksheet['K14'].value
                    t1_t =  worksheet['K20'].value
                    t2_t = worksheet['L20'].value
                    Cp_t = (worksheet['K24'].value+worksheet['L24'].value)/2
                    k_t =(worksheet['K25'].value+worksheet['L25'].value)/2
                    #shell_side ='Cold Side'
                    #s2 = 'Cold side T2'

                    Shell_list = [m_s, t1_s, t2_s, rho_s, Cp_s, mu_s, k_s, fouling_s]
                    Tube_list = [m_t, t1_t, t2_t, rho_t, Cp_t, mu_t, k_t, fouling_t]

                    HB_data,ntu_calc = Heat_balance(shell_side, Tube_list, Shell_list,s2,s3)
                    Q, dTlm, ft = HB_data[0], HB_data[1], HB_data[2]

                    Do = worksheet['F42'].value
                    thick = float(thickness_table[thickness_table['Gauge']==str(worksheet['H42'].value)]['mm']) #2.108
                    print(worksheet['H42'].value)
                    print(thick)
                    Di = (Do - 2*thick)*0.001

                    shell_D = worksheet['H44'].value/1000
                    tn = worksheet['D42'].value

                    #Do = float(st.session_state.geo_df.iloc[2,1])
                    #Di = (Do - 2*float(st.session_state.geo_df.iloc[3,1]))*0.001
                    L = worksheet['J42'].value
                    tpitch = Do*worksheet['M43'].value
                    b_space = worksheet['M48'].value
                    b_cut = worksheet['I48'].value
                    

                    print(tn,shell_D,b_cut,b_space,tpitch,b_cut)
                    
                    pn = worksheet['H37'].value

                    pitch =worksheet['M42'].value
                    if pitch == 'square':
                      t_p_angle = 90
                    elif pitch == 'rotated square 45':
                      t_p_angle = 45
                    elif pitch == 'triangle 30':
                      t_p_angle = 30
                    else: t_p_angle = 60
                    geo_list = [tn ,pn,Do, Di, pitch, tpitch,L, b_space, b_cut,shell_D]

                    try:
                        st.session_state.para_input_df.loc[:,'Kern_summary'] = st.session_state.para_input_df.loc[:,'Bell_summary'] = [m_t, t1_t, t2_t, rho_t, Cp_t, mu_t*1000, k_t, fouling_t,m_s, t1_s, t2_s, rho_s, Cp_s, mu_s*1000, k_s, fouling_s]
                        st.session_state.calculations_df.loc[['Duty','LMTD','Ft','Corrected LMTD'],'Kern_summary'] = st.session_state.calculations_df.loc[['Duty','LMTD','Ft','Corrected LMTD'],'Bell_summary'] = Q,dTlm,ft,dTlm*ft
                       
                        st.session_state.geo_input_df.loc[['Shell D','Baffle Spacing','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut'],'Kern_summary']=st.session_state.geo_input_df.loc[['Shell D','Baffle Spacing','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut'],'Bell_summary']=shell_D,b_space,Do,Di*1000,L,tn,pn,tpitch,pitch,b_cut
           
                        dp_s, dp_t, h_shell, h_t, Uc, Ud, U_calc, Rdesign, Rsevice = kern(Tube_list, Shell_list, geo_list,s3,HB_data,st.session_state.geo_input_df,st.session_state.calculations_df)
                        geo_list = [tn ,pn,Do, Di, pitch, tpitch,L, b_space, b_cut,shell_D]
                        Shell_list = [m_s, t1_s, t2_s, rho_s, Cp_s, mu_s*1000, k_s, fouling_s]
                        Tube_list = [m_t, t1_t, t2_t, rho_t, Cp_t, mu_t*1000, k_t, fouling_t]
                        U_clean,U_dirty,U_required,OD,total_dp_shell,total_dp_tube=bell_delaware(Tube_list, Shell_list ,h_t,h_shell,geo_list,s3,HB_data,st.session_state.geo_input_df,st.session_state.calculations_df)
                        print(U_clean,U_dirty,U_required,OD,total_dp_shell,total_dp_tube)
                    #except UnboundLocalError: pass 
                    except ValueError: pass
            
              except TypeError: st.write('Please Check your dataset')
      
      except ValueError:
        st.write('Error in file')
      st.session_state.dp_calc_check = st.checkbox("Caclulate pressure drop?")
      max_tubes = ht.hx.Ntubes(DBundle=shell_D,Do=Do/1000,pitch=tpitch/1000,Ntp=pn,angle=t_p_angle)
      tube_mask = (max_tubes < tn)
      if st.button("Reveal Calculations", key = 'calculations_table22'):
        if not st.session_state.dp_calc_check:
          st.session_state.calculations_df = st.session_state.calculations_df.dropna(how='any')
          st.session_state.summary = pd.concat([st.session_state.calculations_df, st.session_state.para_input_df])
          st.session_state.summary['Kern_summary'] = st.session_state.summary['Kern_summary'].apply(lambda x: convert_to_float_or_string(x))
          st.session_state.summary['Bell_summary'] = st.session_state.summary['Bell_summary'].apply(lambda x: convert_to_float_or_string(x))
          st.write(st.session_state.summary)
          st.write(pd.DataFrame([ntu_calc]).transpose().rename(columns={0:'NTU Calculations'}))
        else:
          st.session_state.summary = pd.concat([st.session_state.calculations_df,st.session_state.para_input_df,st.session_state.geo_input_df])
          st.session_state.summary['Kern_summary'] = st.session_state.summary['Kern_summary'].apply(lambda x: convert_to_float_or_string(x))
          st.session_state.summary['Bell_summary'] = st.session_state.summary['Bell_summary'].apply(lambda x: convert_to_float_or_string(x))
          st.write(st.session_state.summary)
          st.write(pd.DataFrame([ntu_calc]).transpose().rename(columns={0:'NTU Calculations'}))
          if tube_mask:
            st.warning('Max Tube count for the selected shell diameter is '+str(max_tubes))
         
     
      


if 'current_view' not in st.session_state:
    st.session_state['current_view'] = 'Grid'

if 'current_step' not in st.session_state:
    st.session_state['current_step'] = 1


def set_page_view(page):
    st.session_state['current_view'] = page
    st.session_state['current_step'] = 1         

def set_form_step(action,step=None):
    if action == 'Next':
        st.session_state['current_step'] = st.session_state['current_step'] + 1
    if action == 'Back':
        st.session_state['current_step'] = st.session_state['current_step'] - 1
    if action == 'Jump':
        st.session_state['current_step'] = step


##### wizard functions ####
def wizard_form_header():
    sf_header_cols = st.columns([1,1.75,1])
        
    with sf_header_cols[1]:            
        st.subheader('Heat Exchanger Rating')
            
    # determines button color which should be red when user is on that given step
    wh_type = 'primary' if st.session_state['current_step'] == 1 else 'secondary'
    ff_type = 'primary' if st.session_state['current_step'] == 2 else 'secondary'
    lo_type = 'primary' if st.session_state['current_step'] == 3 else 'secondary'
    sf_type = 'primary' if st.session_state['current_step'] == 4 else 'secondary'
    fb_type = 'primary' if st.session_state['current_step'] == 5 else 'secondary'

    step_cols = st.columns([.5,.85,.85,.85,.85,.85])    
    step_cols[1].button('Fluids Properties',on_click=set_form_step,args=['Jump',1],type=wh_type)
    step_cols[2].button('Complete table',on_click=set_form_step,args=['Jump',2],type=ff_type)        
    step_cols[3].button('Exchanger Geometry',on_click=set_form_step,args=['Jump',3],type=lo_type)      
    step_cols[4].button('Calculations summary',on_click=set_form_step,args=['Jump',4],type=sf_type)
    step_cols[5].button('Perform Trials',on_click=set_form_step,args=['Jump',5],type=fb_type)     
### Replace Wizard Form Body with this ###
                 

def wizard_form_footer():    
    form_footer_container = st.empty()
    with form_footer_container.container():
        
        disable_back_button = True if st.session_state['current_step'] == 1 else False
        disable_next_button = True if st.session_state['current_step'] == 5 else False
        
        form_footer_cols = st.columns([5,1,1,2.5])
        
        form_footer_cols[0].button('Cancel',on_click=set_page_view,args=['Grid'])
        form_footer_cols[1].button('Back',on_click=set_form_step,args=['Back'],disabled=disable_back_button)
        form_footer_cols[2].button('Next',on_click=set_form_step,args=['Next'],disabled=disable_next_button)
        file_ready = False if st.session_state['current_step'] == 4 else True
        submitted = form_footer_cols[3].button('📤 Double Click!',disabled=file_ready)  
    return submitted
        
def submit(button,ntu_calc):
   if button :
      try:
              if not st.session_state.dp_calc_check:
              
                  st.session_state.calculations_df = st.session_state.calculations_df.dropna(how='any')
                  st.session_state.summary = pd.concat([st.session_state.calculations_df, st.session_state.para_input_df])
                  st.session_state.summary['summary'] = st.session_state.summary['summary'].apply(lambda x: convert_to_float_or_string(x))
                  
                  st.write(st.session_state.summary)
                  
              else:
                  st.session_state.summary = pd.concat([st.session_state.calculations_df,st.session_state.para_input_df,st.session_state.geo_input_df])
                  st.session_state.summary['Kern_summary'] = st.session_state.summary['Kern_summary'].apply(lambda x: convert_to_float_or_string(x))
                  st.session_state.summary['Bell_summary'] = st.session_state.summary['Bell_summary'].apply(lambda x: convert_to_float_or_string(x))
                  st.write(st.session_state.summary)
              
              st.write(pd.DataFrame([ntu_calc]).transpose().rename(columns={0:'NTU Calculations'}))
          #except UnboundLocalError: pass 
      except IndexError: pass     

### Replace Render Wizard View With This ###
def render_wizard_view():
    with st.expander('',expanded=True):
        wizard_form_header()
        st.markdown('---')
        main()
        st.markdown('---')
        wizard_form_footer()


##### grid functions ####


#render_wizard_view()
if __name__ == '__main__':
    main()