
# coding: utf-8

# # Equations to derive leaf energy balance components from wind tunnel measurements and compare against leaf model</h1>
# 
# 

# In[1]:

get_ipython().run_cell_magic(u'capture', u'storage', u"# The above redirects all output of the below commands to the variable 'storage' instead of displaying it.\n# It can be viewed by typing: 'storage()'\n# Setting up worksheet and importing equations from other worksheets\nload('temp/Worksheet_setup.sage')\nload_session('temp/leaf_enbalance_eqs.sobj')\nfun_loadvars(dict_vars) # any variables defined using var2() are saved with their attributes in dict_vars. Here we re-define them based on dict_vars from the other worksheet.")


# In[2]:

# %load -s fun_SS 'temp/leaf_enbalance_eqs.sage'
def fun_SS(vdict1):
    '''
    Steady-state T_l, R_ll, H_l and E_l under forced conditions.
    Parameters are given in a dictionary (vdict) with the following entries:
    a_s, a_sh, L_l, P_a, P_wa, R_s, Re_c, T_a, g_sw, v_w
    ''' 
    vdict = vdict1.copy()

    # Nusselt number
    vdict[nu_a] = eq_nua.rhs().subs(vdict)
    vdict[Re] = eq_Re.rhs().subs(vdict)
    vdict[Nu] = eq_Nu_forced_all.rhs().subs(vdict)
    
    # h_c
    vdict[k_a] = eq_ka.rhs().subs(vdict)
    vdict[h_c] = eq_hc.rhs().subs(vdict)
 
    # gbw
    vdict[D_va] = eq_Dva.rhs().subs(vdict)
    vdict[alpha_a] = eq_alphaa.rhs().subs(vdict)
    vdict[rho_a] =  eq_rhoa.rhs().subs(vdict)
    vdict[Le] =  eq_Le.rhs().subs(vdict)
    vdict[g_bw] = eq_gbw_hc.rhs().subs(vdict)   
    
    # Hl, Rll
    vdict[R_ll] = eq_Rll.rhs().subs(vdict)
    vdict[H_l] = eq_Hl.rhs().subs(vdict)   

    # El
    vdict[g_tw] =  eq_gtw.rhs().subs(vdict)
    vdict[C_wa] = eq_Cwl.rhs()(P_wl = P_wa, T_l = T_a).subs(vdict)
    vdict[P_wl] = eq_Pwl.rhs().subs(vdict)
    vdict[C_wl] = eq_Cwl.rhs().subs(vdict)
    vdict[E_lmol] = eq_Elmol.rhs().subs(vdict)
    vdict[E_l] = eq_El.rhs().subs(eq_Elmol).subs(vdict)    

    # Tl
    try:
        vdict[T_l] = find_root((eq_Rs_enbal - R_s).rhs().subs(vdict), 273, 373)
    except: 
        print 'too many unknowns for finding T_l: ' + str((eq_Rs_enbal - R_s).rhs().subs(vdict).args())
    
    # Re-inserting T_l
    Tlss = vdict[T_l]
    for name1 in [C_wl, P_wl, R_ll, H_l, E_l, E_lmol]:
        vdict[name1] = vdict[name1].subs(T_l = Tlss)
    
    # Test for steady state
    if n((E_l + H_l + R_ll - R_s).subs(vdict))>1.:
        return 'error in energy balance: El + Hl + Rll - R_s = ' + str(n((E_l + H_l + R_ll - R_s).subs(vdict))) 
    return vdict


# In[27]:

# Test


# In[3]:

# Test
vdict = cdict.copy()
vdict[a_s] = 1.0    # one sided stomata
vdict[g_sw] = 0.01    
vdict[T_a] = 273 + 25.5
vdict[T_w] = vdict[T_a] # Wall temperature equal to air temperature
vdict[P_a] = 101325
rha = 0.5
vdict[P_wa] = rha*eq_Pwl.rhs()(T_l = T_a).subs(vdict)
vdict[L_l] = 0.03
#vdict[L_A] = vdict[L_l]^2
vdict[Re_c] = 3000
vdict[R_s] = 0.
#vdict[Q_in] = 0
vdict[v_w] = 1

dict_print(fun_SS(vdict))


# ## Gas and energy exchange in a leaf chamber
# Calculations based on leaf_capacitance_steady_state1. However, following the LI-6400XT user manual (Eq. 17-3), we replace the air temperature by wall temperature in the calculation of the net longwave balance of the leaf, as wall temperature can be measured in the chamber. Following the same equation, we also add the leaf thermal emissivity of 0.95 (P. 17-3). 
# **Note that in order to measure sensible heat flux from the leaf, wall temperature must be equal to air temperature!**

# In[4]:

width = 0.05
height = 0.03
volume = 310/(100^3)
print 'Volume = ' + str((volume*1000).n()) + ' l'

print 'min flow rate for flushing = ' + str((volume*3/100^3/60).n()) + ' m3/s'
print 'min flow rate for flushing = ' + str((volume*3*1000).n()) + ' l/min'
print 'flow rate for 1 m/s direct wind = ' + str(width*height*1) + ' m3/s'
print 'flow rate for 1 m/s direct wind = ' + str(width*height*1*1000*60) + ' l/m'
print 'flow rate for 5 m/s direct wind = ' + str(width*height*5) + ' m3/s'
print 'flow rate for 5 m/s direct wind = ' + str(width*height*5*1000*60) + ' l/m'


# In[5]:

width = 0.05
height = 0.03
volume = 400/(100^3)
print 'Volume = ' + str((volume*1000).n()) + ' l'

print 'min flow rate for flushing = ' + str((volume*3/100^3/60).n()) + ' m3/s'
print 'min flow rate for flushing = ' + str((volume*3*1000).n()) + ' l/min'
print 'flow rate for 1 m/s direct wind = ' + str(width*height*1) + ' m3/s'
print 'flow rate for 1 m/s direct wind = ' + str(width*height*1*1000*60) + ' l/m'


# <h2>Chamber insulation material</h2>

# In[6]:

var2('c_pi', 'Heat capacity of insulation material', joule/kilogram/kelvin, latexname='c_{pi}')
var2('lambda_i', 'Heat conductivity of insulation material', joule/second/meter/kelvin)
var2('rho_i', 'Density of insulation material', kilogram/meter^3)
var2('L_i', 'Thickness of insulation material', meter)
var2('A_i', 'Conducting area of insulation material', meter^2)
var2('Q_i', 'Heat conduction through insulation material', joule/second)
var2('dT_i', 'Temperature increment of insulation material', kelvin)


# In[7]:

assumptions(c_pi)


# In[8]:

eq_Qi = Q_i == dT_i*lambda_i*A_i/L_i
units_check(eq_Qi)


# In[9]:

eq_Li = solve(eq_Qi, L_i)[0]
units_check(eq_Li)


# In[10]:

# The amount of heat absorbed by the insulation material per unit area to increase the wall temperature by the same amount as dT_i for given heat flux Q_i
units_check(c_pi*rho_i*dT_i*L_i)


# In[11]:

(c_pi*rho_i*dT_i*L_i).subs(eq_Li)


# In[12]:

# From http://www.sager.ch/default.aspx?navid=25, PIR
vdict = {}
vdict[lambda_i] = 0.022
vdict[c_pi] = 1400
vdict[rho_i] = 30
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[13]:

# From http://www.sager.ch/default.aspx?navid=25, Sagex 15
vdict = {}
vdict[lambda_i] = 0.038
vdict[c_pi] = 1400
vdict[rho_i] = 15
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[14]:

units_check(lambda_i*A_i*dT_i*L_i)


# In[15]:

# Assuming a 30x10x5 cm chamber, how thick would the insulation have to be in order to lose 
# less than 0.01 W heat for 0.1 K dT_i?
eq_Li(A_i = 0.3*0.1*2 + 0.3*0.05*2, dT_i = 0.1, Q_i = 0.01).subs(vdict)


# In[16]:

# From http://www.sager.ch/default.aspx?navid=25, Sagex 30
vdict = {}
vdict[lambda_i] = 0.033
vdict[c_pi] = 1400
vdict[rho_i] = 30
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# ## Leaf radiation balance

# In[17]:

# Additional variables
var2('alpha_l', 'Leaf albedo, fraction of shortwave radiation reflected by the leaf', watt/watt)
var2('R_d', 'Downwelling global radiation', joule/second/meter^2)
var2('R_la', 'Longwave radiation absorbed by leaf', joule/second/meter^2, latexname='R_{la}')
var2('R_ld', 'Downwards emitted/reflected global radiation from leaf', joule/second/meter^2, latexname='R_{ld}')
var2('R_lu', 'Upwards emitted/reflected global radiation from leaf', joule/second/meter^2, latexname='R_{lu}')
var2('R_u', 'Upwelling global radiation', joule/second/meter^2)
var2('S_a', 'Radiation sensor above leaf reading', joule/second/meter^2)
var2('S_b', 'Radiation sensor below leaf reading', joule/second/meter^2)
var2('S_s', 'Radiation sensor beside leaf reading', joule/second/meter^2)


# ### Measuring radiative exchange
# 
# <p>The leaf is exposed to downwelling radiation ($R_d$) originating from shortwave irradiance entering through the glass window plus the longwave irradiance transmitted througha and emitted by the glass window, plus the upwelling radiation ($R_u$) emitted by the bottom glass window.</p>
# <p>The leaf itself reflects some of the radiation in both direction and emits its own black body longwave radiation. The sum of refelcted and emitted radiation away from the leaf is denoted as $R_{lu}$ and $R_{ld}$ for upward and downwards respectively. We have three net radiation sensors in place, one above the leaf measuring $S_a$, one below the leaf measureing $S_b$ and one at the same level beside the leaf measureing $S_s$. These sensor measure:</p>
# <p><img style="float: right;" src="figures/Leaf_radbalance.png" alt="" width="400" height="300" /></p>
# <p>$S_a = R_d - R_{lu}$</p>
# <p>$S_b = R_{ld} - R_u$</p>
# <p>$S_s = R_d - R_u$</p>
# <p>This leaves us with 3 equations with 4 unknows, so we either have to assume that $R_{lu} = R_{ld}$, assuming that both sides of the leaf have the same temperature or $R_u = 0$ to solve the algebraic problem. In daylight, $R_d >> R_u$, so this assumption should not lead to a big bias, however this would imply that $S_b = R_{ld}$, which is certainly incorrect.</p>
# <p>Unfortunately, the assumption $R_{lu} = R_{ld}$ does not help solve the problem as it just implies that $S_s = S_a + S_b$:</p>

# In[18]:

eq_Rs_Rd = R_s == (1-alpha_l)*R_d
eq_Sa = S_a == R_d - R_lu
eq_Sb = S_b == R_ld - R_u
eq_Ss = S_s == R_d - R_u


# In[19]:

# Assuming R_lu = R_ld
eq_assRldRlu = R_ld == R_lu
solve([eq_Sa, eq_Sb.subs(eq_assRldRlu), eq_Ss], R_d, R_lu, R_u)


# In[20]:

# More specifically,
eq1 = solve(eq_Sa, R_d)[0]
eq2 = solve(eq_Sb, R_ld)[0].subs(eq_assRldRlu)
eq3 = solve(eq_Ss, R_u)[0] 
solve(eq1.subs(eq2).subs(eq3), S_s)


# In[21]:

# Assuming that R_u = 0
eq_assRu0 = R_u == 0
soln = solve([eq_Sa, eq_Sb, eq_Ss, eq_assRu0], R_d, R_lu, R_ld, R_u)
print soln


# In[22]:

#eq_Rd = soln[0][0]
#eq_Rlu = soln[0][1]
#eq_Rld = soln[0][2]
#eq_Rlu


# <p>However, what we can do in any case is to quantify the net radiative energy absorbed by the leaf as</p>
# <p>$\alpha_l R_s - R_{ll} = S_a - S_b$:</p>

# In[23]:

# Leaf radiation balance
eq_Rs_Rll = R_s - R_ll == R_d + R_u - R_lu - R_ld
eq_Rbalance = R_s - R_ll == S_a - S_b


# In[24]:

solve([eq_Sa, eq_Sb, eq_Ss, R_d + R_u - R_lu - R_ld == S_a - S_b], R_d, R_lu, R_ld, R_u)


# ## Leaf water vapour exchange and energy balace
# The leaf water vapour exchange and energy balance equations were imported from [leaf_enbalance_eqs](Leaf_enbalance_eqs.ipynb). Here we will use an additional equation to estimate the thickness of the leaf boundary layer and get a feeling for the minimum distance between leaf and sensors to avoid interference with the boundary layer conditions.

# In[25]:

# Blasius solution for BL thickness (http://en.wikipedia.org/wiki/Boundary-layer_thickness)
var2('B_l', 'Boundary layer thickness', meter)
vdict = cdict.copy()
Ta = 300
vdict[T_a] = Ta
vdict[T_l] = Ta
vdict[L_l] = 0.15
vdict[v_w] = 0.5
vdict[Re_c] = 3000
vdict[a_s] = 1
vdict[P_a] = 101325
vdict[P_wa] = 0
eq_Bl = B_l == 4.91*sqrt(nu_a*L_l/v_w)
print eq_Bl.subs(eq_nua).subs(vdict)
eq_gbw_hc.subs(eq_hc).subs(eq_Nu_forced_all).subs(eq_Re).subs(eq_rhoa).subs(eq_Le).subs(eq_Dva).subs(eq_alphaa, eq_nua, eq_ka).subs(vdict)


# In[26]:

vdict[L_l] = 0.03
vdict[v_w] = 8
print eq_Bl.subs(eq_nua).subs(vdict)
eq_gbw_hc.subs(eq_hc).subs(eq_Nu_forced_all).subs(eq_Re).subs(eq_rhoa).subs(eq_Le).subs(eq_Dva).subs(eq_alphaa, eq_nua, eq_ka).subs(vdict)


# In[27]:

# How does B_l scale with wind speed?
vdict[v_w] = v_w
P = plot(eq_Bl.rhs().subs(eq_nua).subs(vdict), (v_w, 0.5,5))
P.axes_labels(['Wind speed (m s$^{-1}$)', 'Boundary layer thickness (m)'])
P


# In[28]:

# Maximum sensible heat flux of a 3x3 cm leaf irradiated by 600 W/m2
600*0.03^2


# <h1>Chamber mass and energy balance</h1>
# <p>Usually, we know the volumetric inflow into the chamber, so to convert to molar inflow (mol s$^{-1}$), we will use the ideal gas law: $P_a V_c = n R_{mol} T_{in}$, where $n$ is the amount of matter in the chamber (mol). To convert from a volume to a flow rate, we replace $V_c$ by $F_{in,v}$. Note that partial pressures of dry air and vapour are additive, such that</p>
# <p>$P_a = P_w + P_d$</p>
# <p>However, the volumes are not additive, meaning that:</p>
# <p>$P_d V_a = n_d R_{mol} T_{a}$</p>
# <p>$(P_a - P_d) V_a = n_a R_{mol} T_{a}$</p>
# <p>i.e. we use the same volume ($V_a$) for both the vapour and the dry air. This is because both the vapour and the dry air are well mixed and occupy the same total volume. Their different amounts are expressed in their partial pressures. If we wanted to calculate the partial volumes they would take up in isolation from each other, we would need to specify at which pressure this volume is taken up and if we used the same pressure for both (e.g. $P_a$), we would obtain a volume fraction for water vapour equivalent to its partial pressure fraction in the former case.</p>
# <p>Therefore, we will distinguish the molar flow rates of water vapour ($F_{in,mol,v}$) and dry air ($F_{in,mol,a}$) but they share a common volumetric flow rate ($F_{in,v}$).</p>

# In[29]:

var2('W_c', 'Chamber width', meter)
var2('L_c', 'Chamber length', meter)
var2('H_c', 'Chamber height', meter)
var2('V_c', 'Chamber volume', meter^3)
var2('n_c', 'molar mass of gas in chamber', mole)
var2('F_in_v', 'Volumetric flow rate into chamber', meter^3/second, latexname='F_{in,v}')
var2('F_in_mola', 'Molar flow rate of dry air into chamber', mole/second, latexname='F_{in,mol,a}')
var2('F_in_molw', 'Molar flow rate of water vapour into chamber', mole/second, latexname='F_{in,mol,w}')
var2('F_out_mola', 'Molar flow rate of dry air out of chamber', mole/second, latexname='F_{out,mol,a}')
var2('F_out_molw', 'Molar flow rate of water vapour out of chamber', mole/second, latexname='F_{out,mol,w}')
var2('F_out_v', 'Volumetric flow rate out of chamber', meter^3/second, latexname='F_{out,v}')
var2('T_d', 'Dew point temperature of incoming air', kelvin)
var2('T_in', 'Temperature of incoming air', kelvin, latexname='T_{in}')
var2('T_out', 'Temperature of outgoing air (= chamber T_a)', kelvin, latexname='T_{out}')
var2('T_room', 'Lab air temperature', kelvin, latexname='T_{room}')
var2('P_w_in', 'Vapour pressure of incoming air', pascal, latexname='P_{w,in}')
var2('P_w_out', 'Vapour pressure of outgoing air', pascal, latexname='P_{w,out}')
var2('R_H_in', 'Relative humidity of incoming air', latexname='R_{H,in}')
var2('L_A', 'Leaf area', meter^2)


# In[30]:

eq_V_c = fun_eq(V_c == W_c*L_c*H_c)
eq_F_in_mola = fun_eq(F_in_mola == (P_a - P_w_in)*F_in_v/(R_mol*T_in))
eq_F_in_molw = fun_eq(F_in_molw == (P_w_in)*F_in_v/(R_mol*T_in))
eq_F_out_mola = fun_eq(F_out_mola == (P_a - P_w_out)*F_out_v/(R_mol*T_out))
eq_F_out_molw = fun_eq(F_out_molw == (P_w_out)*F_out_v/(R_mol*T_out))
eq_F_out_v = fun_eq(F_out_v == (F_out_mola + F_out_molw)*R_mol*T_out/P_a)


# <p>At steady state, $F_{out,mola} = F_{in,mola}$ and $F_{out,molw} = F_{in,molw} + E_{l,mol} L_A$. In the presence of evaporation, we can simply add Elmol to get F_out_v as a function of F_in_v</p>
# <p>Assuming that the pressure inside the chamber is constant and equal to the pressure outside, we compute the change in volumetric outflow due to a change in temperature and due to the input of water vapour by transpiration as:</p>

# In[31]:

eq_Foutv_Finv_Tout = eq_F_out_v.subs(F_out_mola = F_in_mola, F_out_molw = F_in_molw + E_lmol*L_A).subs(eq_F_in_mola, eq_F_in_molw).simplify_full()
units_check(eq_Foutv_Finv_Tout)


# In[32]:

eq_Foutv_Finv_Tout


# In[33]:

# Other way, using molar in and outflow
eq_Foutmolw_Finmolw_Elmol = F_out_molw == (F_in_molw + E_lmol*L_A)
units_check(eq_Foutmolw_Finmolw_Elmol)


# In[ ]:




# <h2>Change in air temperature</h2>
# <p>See also <a href="http://www.engineeringtoolbox.com/mixing-humid-air-d_694.html">http://www.engineeringtoolbox.com/mixing-humid-air-d_694.html</a> and <a href="http://www.engineeringtoolbox.com/enthalpy-moist-air-d_683.html">http://www.engineeringtoolbox.com/enthalpy-moist-air-d_683.html</a> for reference.</p>
# <p>We will assume that the air entering the chamber mixes with the air inside the chamber at constant pressure, i.e. the volume of the mixed air becomes the chamber volume plus the volume of the air that entered. The temperature of the mixed air is then the sum of their enthalpies plus the heat added by the fan and by sensible heaflux, divided by the sum of their heat capacities. The addition of water vapour through evaporation by itself should not affect the air temperature, but the volume of the air.</p>
# <p> </p>
# <p>Alternatively, we could assume that a given amount of air is added to a constant volume, leading to an increase in pressure. Addition of water vapour would lead to an additional increase in pressure. In addition, addition/removal of heat by sensible heat flux and the chamber fan would affect both temperature and pressure.To calculate both temperature and pressure, we need to track the internal energy in addition to the mole number. According to Eq. 6.1.3 in Kondepudi & Prigogine (2006), the internal energy of an ideal gas is given by (see also Eq. 2.2.15 in Kondepuid & Prigogine):</p>
# <p>$U = N(U_0 + C_v T)$</p>
# <p>where</p>
# <p>$U_0 = M c^2$</p>
# <p>The relation between molar heat capacities at constant pressure and volume is given as :</p>
# <p>$C_v = C_p - R_{mol}$</p>
# <p>Any heat exchanged by sensible heat flux, across the walls and the fan can be added to total $U$, and then knowledge about total $C_v$ will let us calculate air temperature inside the chamber. After that, we can use the ideal gas law to calculate volume or pressure, depending in which of those we fixed:</p>
# <p>$P V = n R T$</p>
# <p> </p>
# <p>The difference in water vapour pressure and temperature between the incoming and outgoing air is a function of the latent and sensible heat flux, as well as the flow rate. The heat fluxes associated with the incoming and outgoing air are $T_{in} (c_{pa} F_{in,mola} M_{air} + c_{pv} F_{in,molw} M_{w})$ and $T_{out} (c_{pa} F_{out,mola} M_{air} + c_{pv} F_{out,molw} M_{w})$ respectively. The difference between the two plus any additional heat sources/sinks ($Q_{in}$) equals the sensible heat flux at constant air temperature (steady state).</p>

# In[34]:

units_check(eq_F_out_v)


# In[35]:

var2('M_air', 'Molar mass of air (kg mol-1)', kilogram/mole, value = 0.02897, latexname='M_{air}')  # http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
var2('c_pv', 'Specific heat of water vapour at 300 K', joule/kilogram/kelvin, latexname = 'c_{pv}', value = 1864) # source: http://www.engineeringtoolbox.com/water-vapor-d_979.html
var2('Q_in', 'Internal heat sources, such as fan', joule/second, latexname = 'Q_{in}')
eq_chamber_energy_balance = 0 == H_l*L_A + Q_in + T_in*(c_pa*M_air*F_in_mola + c_pv*M_w*F_in_molw) - (T_out*(c_pa*M_air*F_out_mola + c_pv*M_w*F_out_molw)); show(eq_chamber_energy_balance)
eq_Hl_enbal = solve(eq_chamber_energy_balance.subs(F_out_mola == F_in_mola, F_out_molw == F_in_molw + L_A*E_lmol), H_l)[0].expand()
print units_check(eq_Hl_enbal)


# In[36]:

soln = solve(eq_chamber_energy_balance.subs(eq_Foutmolw_Finmolw_Elmol, F_out_mola == F_in_mola), T_out)
print soln
eq_Tout_Finmol_Tin = soln[0]
units_check(eq_Tout_Finmol_Tin).simplify_full().convert().simplify_full()


# In[37]:

eq_Tout_Finv_Tin = eq_Tout_Finmol_Tin.subs(eq_F_in_mola, eq_F_in_molw).simplify_full()
print eq_Tout_Finv_Tin
show(eq_Tout_Finv_Tin)
units_check(eq_Tout_Finv_Tin).simplify_full().convert()


# In[ ]:




# <p>The molar outflux of dry air equals the molar influx of dry air, while the molar outflux of water vapour equals the molar influx plus the evaporation rate. The sum of both can be used to obtain the volumetric outflow:</p>

# In[38]:

# F_out_v as function of F_inv and T_in
eq1 = (eq_F_out_molw + eq_F_out_mola).simplify_full(); show(eq1)
eq2 = eq1.subs(F_out_mola == F_in_mola, eq_Foutmolw_Finmolw_Elmol).subs(eq_F_in_mola, eq_F_in_molw); show(eq2)
soln = solve(eq2,F_out_v); print soln
eq_Foutv_Finv = soln[0]; show(eq_Foutv_Finv)
units_check(eq_Foutv_Finv)


# In[39]:

eq_Foutv_Finv


# In[40]:

eq_Foutv_Finv_Tout


# In[41]:

eq_Tout_Finmol_Tin


# In[42]:

# Finding the T_in that would balance sensible heat release by the plate for given F_inv
assume(F_in_v > 0)
assume(F_out_v > 0)
assume(E_lmol >=0)
show(eq_chamber_energy_balance)
soln = solve(eq_chamber_energy_balance.subs(F_out_mola == F_in_mola).subs(eq_Foutmolw_Finmolw_Elmol).subs(eq_F_in_mola, eq_F_in_molw), T_in)
print soln
eq_T_in_ss = soln[0]
show(eq_T_in_ss)


# In[43]:

# Calculating Q_in from T_in, F_in_v and T_out
soln = solve(eq_T_in_ss,Q_in)
print soln
eq_Qin_Tin_Tout = soln[0]


# In[44]:

# Calculating H_l from T_in, F_in_v and T_out
soln = solve(eq_T_in_ss,H_l)
print soln
eq_Hl_Tin_Tout = soln[0]
eq_Hl_Tin_Tout.show()


# In[45]:

# Calculating H_l from T_in, T_out and Fmol
soln = solve(eq_chamber_energy_balance.subs(F_out_mola == F_in_mola), H_l)
print soln
eq_Hl_Tin_Tout_Fmol = soln[0].simplify_full()
eq_Hl_Tin_Tout_Fmol.show()


# <h2>Calculation of volumetric flow rate based on Cellkraft measurements</h2>
# <p>Cellcraft uses Arden-Buck equation to convert between vapour pressure and dew point (<a href="http://en.wikipedia.org/wiki/Arden_Buck_Equation">http://en.wikipedia.org/wiki/Arden_Buck_Equation</a>).<br />The air flow rate is given by the Cellkraft humidifier in l/min, but it refers to dry air at 0 $^o$C and 101300 Pa.</p>
# <p>We will use the reported dew point temperature to obtain the vapour pressure of the air coming out from the Cellkraft humidifier, then the ideal gas law to obtain the molar flow of dry air, leading to three equations with three unknowns:</p>
# <p>$F_{\mathit{in}_{\mathit{mola}}} = \frac{F_{\mathit{in}_{\mathit{va}_{n}}} P_{r}}{{R_{mol}} T_{r}}$</p>
# <p>$F_{\mathit{in}_{v}} = \frac{{\left(F_{\mathit{in}_{\mathit{mola}}} + F_{\mathit{in}_{\mathit{molw}}}\right)} {R_{mol}} T_{\mathit{in}}}{P_{a}}$</p>
# <p>$F_{\mathit{in}_{\mathit{molw}}} = \frac{F_{\mathit{in}_{v}} P_{v_{\mathit{in}}}}{{R_{mol}} T_{\mathit{in}}}$</p>

# In[46]:

# Calculating vapour pressure in the incoming air from reported dew point by the Cellkraft humidifier
var2('T0', 'Freezing point in kelvin',kelvin, value = 273.15)
eq_Pwin_Tdew = P_w_in == 611.21*exp((18.678 - (T_d - T0)/234.5)*((T_d - T0)/(257.14-T0+T_d))) # c
list_Tdew = np.array(srange(-40,50,6))+273.25
list_Pwin = [eq_Pwin_Tdew.rhs()(T_d = dummy).subs(cdict) for dummy in list_Tdew]
print list_Pwin
P = plot(eq_Pwin_Tdew.rhs().subs(cdict), (T_d, 253,303), frame = True, axes = False, legend_label = 'Arden Buck')
P += plot(eq_Pwl.rhs().subs(cdict), (T_l, 253,303), color = 'red', linestyle = '--', legend_label = 'Clausius-Clapeyron')
P += list_plot(zip(list_Tdew,list_Pwin), legend_label = 'Cellcraft')
P.axes_labels(['Dew point (K)', 'Saturation vapour pressure (Pa)'])
P


# In[47]:

eq_F_in_molw.show()


# In[ ]:

eq_Finv_Finmol = F_in_v == (F_in_mola + F_in_molw)*R_mol*T_in/P_a
units_check(eq_Finv_Finmol)


# In[ ]:

var2('F_in_va_n', 'Volumetric inflow of dry air at 0oC and 101325 Pa', meter^3/second, latexname='F_{in,v,a,n}') 
var2('P_r', 'Reference pressure',pascal, value = 101325)
var2('T_r', 'Reference temperature', kelvin)

eq_Finmola_Finva_ref = fun_eq(F_in_mola == F_in_va_n * P_r/(R_mol*T_r))


# <p>To get $F_{in,v}$ and $F_{in,mol,w}$, we will consider that:</p>
# <p>$P_d = P_a - P_w$</p>
# <p>$P_w F_{in,v} = F_{in,mol,w} R_{mol} T_{in}$</p>
# <p>$(P_a - P_w) F_{in,v} = F_{in,mol,a} R_{mol} T_{a}$</p>

# In[ ]:

eq_Finmolw_Finv = F_in_molw == (P_w_in*F_in_v)/(R_mol*T_in)
print units_check(eq_Finmolw_Finv)
eq_Finv_Finmola = F_in_v == F_in_mola*R_mol*T_in/(P_a - P_w_in)
print units_check(eq_Finv_Finmola)
eq_Finmolw_Finmola_Pwa = fun_eq(eq_Finmolw_Finv.subs(eq_Finv_Finmola))
eq_Finv_Finva_ref = fun_eq(eq_Finv_Finmola.subs(eq_Finmola_Finva_ref))


# In[ ]:

vdict = cdict.copy()
vdict[F_in_va_n] = 10e-3/60   # 10 l/min reported by Cellkraft
vdict[T_d] = 273.15 + 10    # 10oC dew point
vdict[P_a] = 101325.
vdict[T_r] = 273.15
vdict[P_w_in] = eq_Pwin_Tdew.rhs().subs(vdict)
print vdict[P_w_in]
print vdict[F_in_va_n]

vdict[T_in] = 273.15+0 
inflow0 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 0 oC: ' + str(inflow0) + ' m3/s'

vdict[T_in] = 273.15+25  
inflow25 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 25 oC: ' + str(inflow25) + ' m3/s'


print '25oC/0oC: ' + str(inflow25/inflow0)

vdict[T_in] = 273.15+25 
vdict[P_w_in] = 0. 
inflow25 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 25 oC without added vapour: ' + str(inflow25) + ' m3/s'


# <h2>Vapour pressure</h2>
# <p>The water fluxes associated with the incoming and the outgoing air according to the ideal gas law are $P_{v,in} F_{in,v}/(R_{mol} T_{in})$ and $P_{v,out}  F_{out,v}/(R_{mol} T_{out})$ respectively.</p>

# In[ ]:

eq_Foutv_Finv_Tout.show()


# In[ ]:

eq_F_out_v.show()


# In[ ]:

eq_F_in_mola.subs(eq_Finv_Finva_ref).show()


# In[ ]:

eq_Finv_Finva_ref.show()


# In[ ]:

eq_F_in_molw.show()
eq_F_out_molw.show()
eq_Foutmolw_Finmolw_Elmol.show()


# In[ ]:

eq1 = eq_F_out_molw.rhs() == eq_Foutmolw_Finmolw_Elmol.rhs().subs(eq_F_in_molw)
eq1.show()
soln = solve(eq1, P_w_out)
eq_Pwout_Elmol = fun_eq(soln[0].subs(eq_Foutv_Finv_Tout.subs(eq_F_in_molw)).simplify_full())
units_check(eq_Pwout_Elmol)


# <p><span style="color: #ff0000;">It is a bit surprising that steady-state $P_{w_{out}}$ does not depend on $T_{out}$.</span></p>

# In[ ]:

soln[0].subs(eq_Foutv_Finv_Tout.subs(eq_F_in_molw)).simplify_full().show()
soln[0].subs(eq_F_out_v).subs(F_out_mola = F_in_mola, F_out_molw = F_in_molw + E_lmol*L_A).simplify_full().show()


# <p><span style="color: #ff0000;">The above are equivalent, because $F_{in,v} P_a = (F_{in,mol,a} + F_{in,mol,w}) R_{mol} T_{in}$<br /></span></p>

# In[ ]:

eq_Pwout_Elmol.subs(eq_Finv_Finva_ref).simplify_full().show()


# In[ ]:

show(soln[0])
show(eq_Foutv_Finv_Tout)
show(eq_F_in_molw)


# In[ ]:

show(soln[0].subs(eq_Foutv_Finv_Tout.subs(eq_F_in_molw)).simplify())


# In[ ]:

# T_out cancels out when the above is expanded
show(soln[0].subs(eq_Foutv_Finv_Tout.subs(eq_F_in_molw)).expand())


# <p>To convert from energetic to molar units, we need to divide $E_l$ by $\lambda_E M_w$:</p>

# In[ ]:

eq_Elmol_El = E_lmol == E_l/(lambda_E*M_w)
print units_check(eq_Elmol_El)
eq_El_Elmol = E_l == E_lmol*lambda_E*M_w


# In[ ]:

# In order to keep P_w_out = P_wa = const., we need to adjust P_w_in accordingly.
soln = solve(eq_Pwout_Elmol, P_w_in)
eq_Pwin_Elmol = soln[0].simplify_full()
show(eq_Pwin_Elmol)


# In[ ]:

vdict = cdict.copy()
vdict[F_in_va_n] = 10e-3/60   # 10 l/min reported by Cellkraft
vdict[T_d] = 273.15 + 10    # 10oC dew point
vdict[P_a] = 101325.
vdict[T_r] = 273.15
vdict[P_w_in] = eq_Pwin_Tdew.rhs().subs(vdict)
print vdict[P_w_in]
print vdict[F_in_va_n]

vdict[T_in] = 273.15+0 
inflow0 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 0 oC: ' + str(inflow0) + ' m3/s'

vdict[T_in] = 273.15+25  
inflow25 = eq_Finv_Finva_ref.rhs().subs(vdict)
vdict[F_in_v] = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 25 oC: ' + str(inflow25) + ' m3/s'


print '25oC/0oC: ' + str(inflow25/inflow0)

vdict[T_in] = 273.15+25 
vdict[P_w_in] = 0. 
inflow25 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 25 oC without added vapour: ' + str(inflow25) + ' m3/s'

vdict[E_l] = 100. # assuming 100 W/m2 El
vdict[L_A] = 0.03^2
vdict[E_lmol] = eq_Elmol_El.rhs().subs(vdict)
vdict[T_out] = 273+20. # Assuming 20oC T in chamber

eq_Pwout_Elmol.subs(vdict)


# <h2>Net radiation measurement</h2>
# <p>According to Incropera_fundamentals, Table 13.1, the view factor (absorbed fraction of radiation emitted by another plate) of a small plate of width $w_i$ at a distance $L$ from a parallel larger plate of width $w_j$ is calculated as:</p>

# In[ ]:

var2('L_s', 'Width of net radiation sensor', meter)
var2('L_ls', 'Distance between leaf and net radiation sensor', meter)
var2('F_s', 'Fraction of radiation emitted by leaf, absorbed by sensor', 1)
Wi = L_s/L_ls
Wj = L_l/L_ls
eq_Fs = F_s == (sqrt((Wi + Wj)^2 + 4) - sqrt((Wj - Wi)^2 + 4))/(2*Wi)
units_check(eq_Fs)


# In[ ]:

vdict = cdict.copy()
vdict[L_l] = 0.03
vdict[L_s] = 0.01
print eq_Fs.rhs().subs(vdict)(L_ls = 0.01)
P = plot(eq_Fs.rhs().subs(vdict), (L_ls, 0.001, 0.02))
P.axes_labels(['Distance to leaf (m)', 'Fraction of radiation captured'])
P


# # Saving definitions to separate file
# In the below, we save the definitions and variables to separate files in the /temp directory, one with the extension .sage, from which we can selectively load functions using
# `%load fun_name filenam.sage`
# and one with the extension .sobj, to be loaded elsewhere using 
# `load_session()`

# In[ ]:

fun_export_ipynb('leaf_chamber_eqs', 'temp/')
save_session('temp/leaf_chamber_eqs.sobj')


# # Table of symbols

# In[ ]:

# Creating dictionary to substitute names of units with shorter forms
var('m s J Pa K kg mol')
subsdict = {meter: m, second: s, joule: J, pascal: Pa, kelvin: K, kilogram: kg, mole: mol}
var('N_Re_L N_Re_c N_Le N_Nu_L N_Gr_L N_Sh_L')
dict_varnew = {Re: N_Re_L, Re_c: N_Re_c, Le: N_Le, Nu: N_Nu_L, Gr: N_Gr_L, Sh: N_Sh_L}
dict_varold = {v: k for k, v in dict_varnew.iteritems()}
variables = sorted([str(variable.subs(dict_varnew)) for variable in udict.keys()],key=str.lower)
tableheader = [('Variable', 'Description (value)', 'Units')]
tabledata = [('Variable', 'Description (value)', 'Units')]
for variable1 in variables:
    variable2 = eval(variable1).subs(dict_varold)
    variable = str(variable2)
    tabledata.append((eval(variable),docdict[eval(variable)],fun_units_formatted(variable)))

table(tabledata, header_row=True)


# In[ ]:



