import pandas as pd
import numpy as np
import streamlit as st
import ht
@st.cache_data
def load_table():
    url ='http://raw.githubusercontent.com/Ahmedhassan676/htcalc/main/heat_table.csv'

    return pd.read_csv(url, index_col=[0])
@st.cache_data
def load_summary_table():
    url ='http://raw.githubusercontent.com/Ahmedhassan676/htcalc/main/summary_table.csv'

    return pd.read_csv(url, index_col=[0])
@st.cache_data
def load_data_table():
    url ='http://raw.githubusercontent.com/Ahmedhassan676/htcalc/main/data_tables.csv'

    return pd.read_csv(url)
rating_table = load_table().iloc[2:12,:]

def main():
    html_temp="""
    <div style="background-color:lightblue;padding:16px">
    <h2 style="color:black"; text-align:center> Heat Exchangers Calculation </h2>
    </div>
    
        """
    st.markdown(html_temp, unsafe_allow_html=True)
    
    s1 = st.selectbox('Select Calculations required',('Heat Exchanger Assessment','Heat Exchanger Rating','Quick Heat Exchanger Sizing','Calculate fouling'), key = 'type')
    if s1 == 'Heat Exchanger Assessment':
      s2 = st.selectbox('Select Heat Balance variable',('Hot side mass flow','Hot side T1','Hot side T2','Cold side mass flow','Cold side T1','Cold side T2'), key = 'HB')  
      rating_df = st.experimental_data_editor(rating_table)
      dp_calc_check = st.checkbox("Caclulate pressure drop?")
      shell_side = st.selectbox('Shell Side is the..?',('Cold Side','Hot Side'), key = 'shell_side')
      
      s3 = st.selectbox('Number of Shells',(1,2,3,4,5,6,7,8), key='shells')
      
      
      try:
        t1_s = float(rating_df.iloc[1,1])
        t2_s = float(rating_df.iloc[2,1] )
        m_s = float(rating_df.iloc[0,1])
        Cp_s = float(rating_df.iloc[6,1]) 
        mu_s = float(rating_df.iloc[7,1])
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
          
      except UnboundLocalError: pass
      if not dp_calc_check:
        A = st.number_input('Total Heat Exchanger(s) Area', key = 'a')
        U = st.number_input('Service U Kcal/hr.m2.C', key = 'U')
        U_calc = Q/(ft*dTlm*A)
      if dp_calc_check:
        geo_table = load_table().iloc[26:,:2].rename(columns={'Shell Fluid':'Value'})
        tube_table = load_data_table().iloc[1:11,1:5]
        thickness_table = load_data_table().iloc[11:36,1:4]
        shell_table = load_data_table().iloc[37:67,1]
        #st.write(thickness_table.iloc[1:25,2])
        Do = float(st.selectbox('tube OD (mm)?',tube_table.iloc[1:10,0], key = 'Do'))
        thick = float(st.selectbox('tube gauge (thickness)?',thickness_table.iloc[1:25,2], key = 'thick'))
        Di = (Do - 2*thick)*0.001
        pitch = st.selectbox('pitch type?',('square','triangle'), key = 'pitch')
        shell_D = float(st.selectbox('Shell diameter (mm)?',shell_table, key = 'Shell ID'))/1000
        geo_df = st.experimental_data_editor(geo_table.iloc[[0,1,4,6,7,8],:])
        try:
            tn = float(geo_df.loc['Number of tubes','Value'])
            pn = float(geo_df.loc['number of passes','Value'])
            #Do = float(geo_df.iloc[2,1])
            #Di = (Do - 2*float(geo_df.iloc[3,1]))*0.001
            L = float(geo_df.loc['Tube length','Value'])
            tpitch = float(geo_df.loc['pitch','Value'])
            b_space = float(geo_df.loc['baffle spacing','Value'])
            b_cut = float(geo_df.loc['baffle cut','Value'])
            #shell_D = float(geo_df.iloc[-1,1])/1000
        
            if pitch == 'square':
              De = 4*(((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.25))/(3.14*Do*0.001)
            else:
              De = 8*(0.43301*((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.125))/(3.14*Do*0.001)
            
            C = tpitch-Do
            As = (shell_D*b_space*C)/(tpitch*1000)
            Gs = m_h/(As*3600)
            Res = (1000*De*Gs)/mu_s
            f = np.exp(0.576-(0.19*np.log(Res)))
            Nb = (L/b_space)-1
            dp_s = ((f*(Gs**2)*(Nb+1)*shell_D)/(2*rho_s*De))*0.000010197
            L = L/1000
            A = np.pi*L*Do*0.001*s3*tn
            Cp_t = Cp_t*4184
            Cp_s = Cp_s*4184
            cross_A=(np.pi*0.25*(Di**2))*(tn/pn)
            velocity_t = m_c/(rho_t*3600*cross_A)
            Ret=(rho_t*velocity_t*Di)/mu_t
            f_t =1/(1.58*np.log(Ret)-3.28)**2
            port_1 = f_t*L*pn/Di
            port_2 = rho_t*(velocity_t**2)/2
            dp_t = (4*(port_1)+4*(pn))*port_2*0.000010197
            h_shell = (0.36*((De*Gs/mu_s)**0.55)*((Cp_s*mu_s/k_s)**(1/3)))*k_s/De
            Pr = Cp_t*mu_t/k_t
            Nu = ((0.5*f_t*(Ret-1000)*Pr))/(1+12.7*((0.5*f_t)**0.5)*((Pr**(2/3))-1))
            h_t = Nu *k_t/Di
            d_ratio = Do/(Di*1000)
            Uc = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell))
            Ud = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell)+fouling_s+(d_ratio*fouling_t))
            U_calc = Q/(ft*dTlm*A)
            Rdesign = - (1/Uc) + (1/Ud)
            Rsevice = - (1/Uc) + (1/U_calc)
        except UnboundLocalError: pass 
      if st.button("Reveal Calculations", key = 'calculations_table22'):
        if not dp_calc_check:
          summary_df = load_summary_table()
          summary_df.iloc[0,1] = Q
          summary_df.iloc[2,1] = U
          summary_df.iloc[3,1] = U_calc
          summary_df.iloc[4,1] = 100*(U-U_calc)/U
          
          st.write(summary_df.iloc[[0,2,3,4],1])
        else:
          summary_df = load_summary_table()
          summary_df.iloc[0,1] = Q
          summary_df.iloc[1,1] = Uc
          summary_df.iloc[2,1] = Ud
          summary_df.iloc[3,1] = U_calc
          summary_df.iloc[4,1] = 100*(Ud-U_calc)/Ud
          summary_df.iloc[5,1] = Rdesign
          summary_df.iloc[6,1] = Rsevice
          summary_df.iloc[7,1] = dp_s
          summary_df.iloc[8,1] = dp_t
          st.write(summary_df)
          
             
          
if __name__ == '__main__':
    main()
