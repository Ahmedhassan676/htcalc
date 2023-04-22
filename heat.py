import pandas as pd
import numpy as np
import streamlit as st
import ht
@st.cache_data
def load_table():
    url ='heat_table.csv'

    return pd.read_csv(url, index_col=[0])
rating_table = load_table().iloc[:11,:]

def main():
    html_temp="""
    <div style="background-color:lightblue;padding:16px">
    <h2 style="color:black"; text-align:center> Heat Exchangers Calculation </h2>
    </div>
    
        """
    st.markdown(html_temp, unsafe_allow_html=True)
    
    s1 = st.selectbox('Select Calculations required',('Heat Exchanger Assessment','Heat Exchanger Rating','Quick Heat Exchanger Sizing','Calculate fouling'), key = 'type')
    if s1 == 'Heat Exchanger Assessment':
      s2 = st.selectbox('Select Heat Balance variable',('Hot side mass flow','Hot side T1','Hot side T2',' Cold side mass flow','Cold side T1','Cold side T2'), key = 'HB')  
      rating_df = st.experimental_data_editor(rating_table)
      dp_calc_check = st.checkbox("Caclulate pressure drop?")
      shell_side = st.selectbox('Shell Side is the..?',('Cold Side','Hot Side'), key = 'shell_side')
      if s2 == 'Hot side mass flow':
        s3 = st.selectbox('Number of Shells',(1,2,3,4,5,6,7,8), key='shells')
        A = st.number_input('Total Heat Exchanger(s) Area', key = 'a')
        U = st.number_input('Service U Kcal/hr.m2.C', key = 'U')
        if shell_side == 'Hot Side':
          try:
            T1 = float(rating_df.iloc[3,1])
            T2 = float(rating_df.iloc[4,1] )
            Cp_h = float(rating_df.iloc[8,1]) 
            mu_s = float(rating_df.iloc[9,1])
            mu_t = float(rating_df.iloc[9,2])/1000
            rho_s =  float(rating_df.iloc[7,1])
            rho_t =  float(rating_df.iloc[7,1])
            m_c = float(rating_df.iloc[2,2])
            t1 = float(rating_df.iloc[3,2]) 
            t2 = float(rating_df.iloc[4,2]) 
            Cp_c = float(rating_df.iloc[8,2]) 
            Q = m_c * Cp_c * (t2-t1)
            m_h = Q/(Cp_h*(T1-T2))
            dTlm = ht.LMTD(T1,T2,t1,t2)
            
            ft = ht.F_LMTD_Fakheri(t1,t2,T1,T2,s3)
            U_calc = Q/(ft*dTlm*A)
          except UnboundLocalError: pass
        if dp_calc_check:
          geo_table = load_table().iloc[26:,:2]
          pitch = st.selectbox('pitch type?',('square','triangle'), key = 'pitch')
          geo_df = st.experimental_data_editor(geo_table)
          try:
              tn = float(geo_df.iloc[0,1])
              pn = float(geo_df.iloc[1,1])
              Do = float(geo_df.iloc[2,1])
              Di = (Do - 2*float(geo_df.iloc[3,1]))*0.001
              L = float(geo_df.iloc[4,1])
              tpitch = float(geo_df.iloc[6,1])
              b_space = float(geo_df.iloc[7,1])
              b_cut = float(geo_df.iloc[8,1])
              shell_D = float(geo_df.iloc[-1,1])/1000
          
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
              cross_A=(np.pi*0.25*(Di**2))*(tn/pn)
              velocity_t = m_c/(rho_t*3600*cross_A)
              Ret=(rho_t*velocity_t*Di)/mu_t
              f_t =1/(1.58*np.log(Ret)-3.28)**2
              port_1 = f_t*L*pn/Di
              port_2 = rho_t*(velocity_t**2)/2
              dp_t = (4*(port_1)+4*(pn))*port_2*0.000010197
          except UnboundLocalError: pass 
        if st.button("Reveal Calculations", key = 'calculations_table22'):
          st.write(U_calc,ft,Q,ft*dTlm,dp_s,dp_t,port_2)

if __name__ == '__main__':
    main()
