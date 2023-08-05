import streamlit as st
import numpy as np
import pandas as pd
from scipy.optimize import fsolve,root
import ht

def main_kern(Tube_list, Shell_list,HB_data,j_const,Do,thick,L,geo_input_list,dp_sin,dp_tin,s3):
    m_t,t1_t,t2_t,rho_t,Cp_t,mu_t,k_t,fouling_t = Tube_list[0], Tube_list[1], Tube_list[2], Tube_list[3], Tube_list[4], Tube_list[5]*0.001, Tube_list[6], Tube_list[7]
    m_s,t1_s,t2_s,rho_s,Cp_s,mu_s,k_s,fouling_s = Shell_list[0], Shell_list[1], Shell_list[2], Shell_list[3], Shell_list[4], Shell_list[5]*0.001, Shell_list[6], Shell_list[7]
    Cp_t = Cp_t*4184
    Cp_s = Cp_s*4184
    L = L/1000
    Di = (Do - 2*thick)*0.001
    pn = 2 # assumed
    tpitch = 23.81 # assumed
    Q, dTlm, ft = HB_data[0], HB_data[1], HB_data[2]
    U_assumed = 350
    corrected_LMTD = dTlm*ft
    A_required = Q/(corrected_LMTD*U_assumed)
    velocity_t = 0
    while velocity_t < 2:
        tn = int (A_required/(np.pi*L*Do*0.001*s3))
        
        cross_A=(np.pi*0.25*(Di**2))*(tn/pn)
        velocity_t = m_t/(rho_t*3600*cross_A)
        bundle = ht.hx.DBundle_for_Ntubes_Phadkeb(tn, Do/1000, tpitch/1000, pn, angle=30)
        shell_D = 522 #bundle*1000+ht.shell_clearance(DBundle=bundle)*1000
        Ret=(rho_t*velocity_t*Di)/mu_t
        f_t =1/(1.58*np.log(Ret)-3.28)**2 # valid for Re 2300 to 5,000,000 and Pr 0.5 to 2000
        #port_1 = f_t*L*pn/Di
        #port_2 = rho_t*(velocity_t**2)/2
        #dp_t = (4*(port_1)+4*(pn))*port_2*0.000010197
        Pr = Cp_t*mu_t/k_t
        Nu = ((0.5*f_t*(Ret-1000)*Pr))/(1+12.7*((0.5*f_t)**0.5)*((Pr**(2/3))-1)) # valid for Re 2300 to 5,000,000 (Gnielinski)
        h_t = Nu *k_t/Di
        if velocity_t < 2:
            pn +=2
        st.write(pn,velocity_t,shell_D,h_t,Ret,Nu)
    b_space = shell_D/5 # assumed
    C = tpitch-Do
    As = (shell_D*b_space*C)/(tpitch*1000000)
    st.warning(As)
    st.warning(shell_D)
    Gs = m_s/(As*3600)
    velocity_s = Gs/rho_s 
    pitch = 'triangle 30'
    if pitch == 'square' or 'rotated square 45':
        De = 4*(((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.25))/(3.14*Do*0.001)
    else:
        De = 8*(0.43301*((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.125))/(3.14*Do*0.001)  
    
    h_shell = (0.36*((De*Gs/mu_s)**0.55)*((Cp_s*mu_s/k_s)**(1/3)))*k_s/De
    d_ratio = Do/(Di*1000)
    Uc = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell))
    Ud = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell)+fouling_s+(d_ratio*fouling_t))
    percentage_diff = ((Ud-U_assumed)/U_assumed)*100

    Res = (De*Gs)/mu_s
    f = np.exp(0.576-(0.19*np.log(Res)))
    Nb = int((L*1000/b_space)-1)
    st.warning(Nb)
    dp_s = ((f*(Gs**2)*(Nb+1)*shell_D)/(2*rho_s*De*1000))*0.000010197
    Ret=(rho_t*velocity_t*Di)/mu_t
    f_t =1/(1.58*np.log(Ret)-3.28)**2 # valid for Re 2300 to 5,000,000 and Pr 0.5 to 2000
    port_1 = f_t*L*pn/Di
    port_2 = rho_t*(velocity_t**2)/2
    dp_t = (4*(port_1)+4*(pn))*port_2*0.000010197
    error_dp_s = dp_s-(dp_sin/10**8)
    error_dp_t = dp_t-(dp_tin/10**8)
    if error_dp_s < 0:
        b_space -=(shell_D/5)*0.1
    st.write(dp_s,Res,Gs,velocity_s,De)
