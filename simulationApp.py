import streamlit as st
from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#st.set_page_config(layout="wide")
st.set_page_config(page_title='SimulationApp',page_icon=':100:')
st.title('Simulador de Pandemias :earth_americas:')
@st.cache
def deriv(y, t, N, beta, gamma, sigma, xi, mu, nu):
    #Ecuaciones del cambio para solucionar
    S, E, I, R, D = y
    dSdt = -beta(t) * S * I / N + xi * R -nu(t)*S #nu(t)*S vacunados en tiempo t, solo vacunación sobre susceptibles es efectiva
    dEdt = beta(t) * S * I / N - sigma * E
    dIdt = sigma * E - (1-mu)*gamma * I - mu*I
    dRdt = (1-mu)*gamma * I - xi*R +nu(t) *S
    dDdt = mu * I
    return dSdt, dEdt, dIdt, dRdt,dDdt

#La sidebar NO se actualiza cada vez que recarga el archivo!!!
N = st.sidebar.number_input('Población Total',min_value=1_000,max_value=100_000_000,value=1_000_000,step=1_000)
D = st.sidebar.number_input('Duración de la Infección',min_value=1,max_value=20,value=10,step=1) # infections lasts four days
gamma = 1.0 / D
Inc = st.sidebar.number_input('Periodo de Incubación',min_value=0,max_value=15,value=2,step=1)
sigma = 1.0 / Inc  # incubation period of five days

R_0 = st.sidebar.number_input('Tasa de Reproducción Base (R0)',min_value=0.0,max_value=10.0,value=2.5,step=0.5)

tap = st.sidebar.slider('Porcentaje Uso del tapabocas',0,100,0,1)
tap = tap*0.35/100

hig = st.sidebar.slider('Porcentaje de Higiene (lavado de manos)',0,100,0,1)
hig = hig*0.25/100


mu = st.sidebar.slider('Porcentaje de Mortalidad',0,10,0,1)
mu = mu/100

#beta = (1-tap)*(1-hig)*R_0 * gamma   # R_0 = beta / gamma, so beta = R_0 * gamma


if st.sidebar.checkbox('¿Pierde la Inmunidad?'):

    Inm = st.sidebar.number_input('Periodo de Perdida de Inmunidad',min_value=1,value=365,step=1)
    xi = 1.0/Inm
else: 
    xi=0.0

if st.sidebar.checkbox('¿Hay Encerramiento?'):

    L = st.sidebar.number_input('Comienzo del Encerramiento',min_value=1,value=200,step=1)
    dur = st.sidebar.slider('Duración',min_value=1,max_value=60,value=20,step=1)
else: 
    L = 1000
    dur=0
#Para el encerramiento:
def R__0(t):
    return R_0 if ((t < L)|(t>(L+dur))) else 0.9
def beta(t):
    return (1-tap)*(1-hig)*R__0(t) * gamma

if st.sidebar.checkbox('¿Hay vacunación?'):
    vac = st.sidebar.slider('Porcentaje de Vacunados',0,100,0,1)
    dia = st.sidebar.number_input('Día de Vacunación',1,None,10,1)
else:
    vac=0
    dia=0
def nu(t):
    return vac/100 if (t<dia+1+vac/100)&(t>=dia) else 0.0   #Tasa de vacunación

t_max = st.sidebar.number_input('Tiempo Máximo de Simulación',min_value=100,value=400,step=10)

st.write('Parametros: ','\u03C3:',str(sigma),'\u03B3:',str(gamma),'\u03B2:',str(beta(0)),'\u03BE:',str(xi),'\u03BC',str(mu),'\u03BD',str(nu(0)))
S0, E0, I0, R0, D0 = N-1 , 1, 0, 0, 0  # initial conditions: one exposed
#A los susceptibles les quitamos los vacunados, y a los recuperados se los sumamos (En las ecuaciones)
t = np.linspace(0, t_max, t_max)
#Es necesario que sean t_max pasos para que el df sea acorde

y0 = S0, E0, I0, R0, D0 # Initial conditions vector

# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, sigma, xi, mu, nu))
S, E, I, R, D = ret.T

population = pd.DataFrame(columns=['Susceptibles','Expuestos','Infectados','Recuperados','Fallecidos','Total'])
population['Susceptibles'] = S
population['Expuestos'] = E
population['Infectados'] = I
population['Recuperados'] = R
population['Fallecidos'] = D
population['Total'] = S+E+I+R+D

if st.selectbox('Tipo de Gráfico',('Area','Linea'))=='Area':
    st.area_chart(population[['Total','Recuperados','Fallecidos','Infectados','Expuestos']])

else:
    st.line_chart(population)

    
#to json
# result = population.to_json(orient="records")
# parsed = json.loads(result)
# st.write(json.dumps(parsed, indent=4))
