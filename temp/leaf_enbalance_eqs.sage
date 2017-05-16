
# coding: utf-8

# # Leaf energy balance equations
# Based on the following paper:
# Schymanski, S.J. and Or, D. (2016): [Leaf-scale experiments reveal important omission in the Penman-Monteith equation.](http://www.hydrol-earth-syst-sci-discuss.net/hess-2016-363/) Hydrology and Earth System Sciences Discussions, p.1–33. doi: 10.5194/hess-2016-363.
# 
# Author: Stan Schymanski (stan.schymanski@env.ethz.ch)
# 
# Note: This worksheet is prepared for the open source software [sage](http://www.sagemath.org). It relies on definitions provided in Worksheet [Worksheet_setup](Worksheet_setup.ipynb).

# In[1]:

get_ipython().run_cell_magic(u'capture', u'storage', u"# The above redirects all output of the below commands to the variable 'storage' instead of displaying it.\n# It can be viewed by typing: 'storage()'\n\n# Loading modules, settings and custom functions\nload('temp/Worksheet_setup.sage')")


# ## Definition of variables
# All variables are also listed in a [Table](#Table-of-symbols) at the end of this document.

# In[2]:

var2('alpha_a', 'Thermal diffusivity of dry air',meter^2/second)
var2('a_s', 'Fraction of one-sided leaf area covered by stomata (1 if stomata are on one side only, 2 if they are on both sides)', meter^2/meter^2)
var2('a_sh', 'Fraction of projected area exchanging sensible heat with the air (2)', meter^2/meter^2, value = 2, latexname = 'a_{sh}')
var2('c_pa', 'Specific heat of dry air (1010) ', joule/kilogram/kelvin, latexname = 'c_{pa}', value = 1010)
#var2('c_pamol', 'Molar specific heat of dry air (29.19) ', joule/mole/kelvin, latexname = 'c_{pa,mol}', value = 29.19) # https://en.wikipedia.org/wiki/Heat_capacity#Specific_heat_capacity
var2('C_wa', 'Concentration of water in the free air ',mole/meter^3, latexname = 'C_{wa}')
var2('C_wl', 'Concentration of water in the leaf air space ',mole/meter^3, latexname = 'C_{wl}')
var2('D_va', 'Binary diffusion coefficient of water vapour in air',meter^2/second, latexname = 'D_{va}')
var2('E_lmol', 'Transpiration rate in molar units',mole/second/meter^2,latexname='E_{l,mol}')
var2('E_l', 'Latent heat flux from leaf',joule/second/meter^2)
var2('epsilon_l', 'Longwave emmissivity of the leaf surface (1.0)', value = 1, units = 1/1)
var2('g', 'Gravitational acceleration (9.81)', meter/second^2, value= 9.81)
var2('g_bw', 'Boundary layer conductance to water vapour ',meter/second,latexname='g_{bw}')
var2('g_bwmol', 'Boundary layer conductance to water vapour ',mole/meter^2/second,latexname='g_{bw,mol}')
var2('Gr', 'Grashof number', latexname = 'N_{Gr_L}')
var2('g_sw', 'Stomatal conductance to water vapour',meter/second,latexname='g_{sw}')
var2('g_swmol', 'Stomatal conductance to water vapour',mole/meter^2/second,latexname='g_{sw,mol}')
var2('g_tw', 'Total leaf conductance to water vapour',meter/second,latexname='g_{tw}')
var2('g_twmol', 'Total leaf layer conductance to water vapour',mole/meter^2/second,latexname='g_{tw,mol}')
var2('h_c', 'Average 1-sided convective transfer coefficient', joule/meter^2/second/kelvin)
var2('H_l', 'Sensible heat flux from leaf',joule/second/meter^2)
var2('k_a', 'Thermal conductivity of dry air',joule/second/meter/kelvin)
var2('lambda_E', 'Latent heat of evaporation (2.45e6)',joule/kilogram,value = 2.45e6)
var2('Le', 'Lewis number', latexname = 'N_{Le}')
var2('L_l', 'Characteristic length scale for convection (size of leaf)',meter)
var2('M_N2', 'Molar mass of nitrogen (0.028)',kilogram/mole,value = 0.028)
var2('M_O2', 'Molar mass of oxygen (0.032)',kilogram/mole,value = 0.032)
var2('M_w', 'Molar mass of water (0.018)',kilogram/mole,value = 0.018)
var2('nu_a', 'Kinematic viscosity of dry air',meter^2/second)
var2('Nu', 'Nusselt number',latexname = 'N_{Nu_L}')
var2('P_a', 'Air pressure',pascal)
var2('Pr', 'Prandtl number (0.71)',latexname = 'N_{Pr}',value = 0.71)
var2('P_N2', 'Partial pressure of nitrogen in the atmosphere', pascal, latexname = 'P_{N2}')
var2('P_O2', 'Partial pressure of oxygen in the atmosphere', pascal, latexname = 'P_{O2}')
var2('P_wa', 'Vapour pressure in the atmosphere', pascal, latexname = 'P_{wa}')
var2('P_was', 'Saturation vapour pressure at air temperature', pascal, latexname = 'P_{was}')
var2('P_wl', 'Vapour pressure inside the leaf', pascal, latexname = 'P_{wl}')
var2('r_bw', 'Boundary layer resistance to water vapour, inverse of $g_{bw}$', second/meter, latexname = 'r_{bw}')
var2('r_sw', 'Stomatal resistance to water vapour, inverse of $g_{sw}$', second/meter, latexname = 'r_{sw}')
var2('r_tw', 'Total leaf resistance to water vapour, $r_{bv} + r_{sv}$', second/meter, latexname = 'r_{tw}')
var2('Re_c', 'Critical Reynolds number for the onset of turbulence',latexname='N_{Re_c}')
var2('Re', 'Reynolds number',latexname='N_{Re_L}')
var2('rho_a', 'Density of dry air', kilogram/meter^3)
var2('rho_al', 'Density of air at the leaf surface', kilogram/meter^3)
var2('R_ll', 'Longwave radiation away from leaf',joule/second/meter^2, latexname = 'R_{ll}')
var2('R_mol', 'Molar gas constant (8.314472)',joule/mole/kelvin,latexname='R_{mol}', value = 8.314472)
var2('R_s', 'Solar shortwave flux',joule/second/meter^2)
var2('Sh', 'Sherwood number',latexname = 'N_{Sh_L}')
var2('sigm', 'Stefan-Boltzmann constant (5.67e-8)',joule/second/meter^2/kelvin^4,'real',latexname='\\sigma',value=5.67e-8)
var2('T_a', 'Air temperature',kelvin)
var2('T_l', 'Leaf temperature',kelvin)
var2('T_w', 'Radiative temperature of objects surrounding the leaf', kelvin)
var2('v_w', 'Wind velocity',meter/second)


# # Mathematical derivations
# ## Leaf energy balance
# The material below follows closely derivations published previously \citep{schymanski_stomatal_2013}, with re-organisation of equations for greater consistence with the present paper.
# 
# The leaf energy balance is determined by the dominant energy fluxes between the leaf and its surroundings, including radiative, sensible, and latent energy exchange (linked to mass exchange).
# 
# 
# In this study we focus on steady-state conditions, in which the energy balance can be written as:
# ##### {eq_Rs_enbal}
# $$ R_s = R_{ll} + H_l + E_l $$ 
# where $R_s$ absorbed short wave radiation, $R_{ll}$ is the net longwave balance, i.e. the emitted minus the absorbed, $H_l$ is the sensible heat flux away from the leaf and $E_l$ is the latent heat flux away from the leaf. In the above, extensive variables are defined per unit leaf area. Following our previous work \citep{schymanski_stomatal_2013}, this study considers spatially homogeneous planar leaves, i.e. homogenous illumination and a negligible temperature gradient between the two sides of the leaf. 
# The net longwave emission is represented by the difference between blackbody radiation at leaf temperature ($T_l$) and that at the temperature of the surrounding objects ($T_w$) \citep{Monteith_principles_2007}:
# ##### {eq_Rll}
# $$R_{ll} = a_{sH} \epsilon_l \sigma (T_{l}^{4} - T_{w}^{4})$$
# 
# where $a_{sH}$ is the fraction of projected leaf area exchanging radiative and sensible heat (2 for a planar leaf, 1 for a soil surface), $\epsilon_l$ is the leaf's longwave emmissivity ($\approx 1$) and $\sigma$ is the Stefan-Boltzmann constant.
# Total convective heat transport away from the leaf is represented as:
# ##### {eq_Hl}
# $$H_l = a_{sH} h_c (T_l - T_a)$$
# 
# where $h_c$ is the average one-sided convective heat transfer coefficient, determined by properties of the leaf boundary layer.
# 
# Latent heat flux ($E_l$, W~m$^{-2}$) is directly related to the transpiration rate ($E_{l,mol}$) by:
# ##### {eq_El}
# $$E_l = E_{l,mol} M_w \lambda_E$$
# 
# where $M_w$ is the molar mass of water and $\lambda_E$ the latent heat of vaporisation. $E_{l,mol}$ (mol~m$^{-2}$~s$^{-1}$) was computed in molar units as a function of the concentration of water vapour within the leaf ($C_{wl}$, mol m$^{-3}$) and in the free air ($C_{wa}$, mol m$^{-3}$) \citep[Eq. 6.8]{incropera_fundamentals_2006}:
# ##### {eq_Elmol}
# $$E_{l,mol} = g_{tw}(C_{wl} - C_{wa}) $$
# 
# where $g_{tw}$ (m~s$^{-1}$) is the total leaf conductance for water vapour, dependent on stomatal ($g_{sw}$) and boundary layer conductance ($g_{bw}$) in the following way:
# ##### {eq_gtw}
# $$g_{tw} = \frac{1}{\frac{1}{g_{sw}} + \frac{1}{g_{bw}}}$$
# 

# In[3]:

eq_Rs_enbal = R_s == R_ll + H_l + E_l
units_check(eq_Rs_enbal).simplify_full()


# In[4]:

eq_Rll = R_ll == a_sh*sigm*epsilon_l*(T_l^4 - T_w^4)
units_check(eq_Rll).simplify_full()


# In[5]:

eq_Hl = H_l == a_sh*h_c*(T_l - T_a)
units_check(eq_Hl).simplify_full()


# In[6]:

eq_El = E_l == E_lmol*M_w*lambda_E
units_check(eq_El).simplify_full()


# In[7]:

ustr = str(units_check(M_w*lambda_E))
print str((M_w*lambda_E).subs(cdict)) + ustr


# In[8]:

eq_Elmol = E_lmol == g_tw*(C_wl - C_wa)
units_check(eq_Elmol).simplify_full()


# In[9]:

eq_gtw = g_tw == 1/(1/g_sw + 1/g_bw)
eq_gtw.simplify_full().show()
units_check(eq_gtw).simplify_full()


# ## Boundary layer conductance to water vapour
# The total leaf conductance to water vapour is determined by the boundary layer and stomatal conductances and equal to 1 over the sum of their respectife resistances ($g_{tw} = 1/(r_{sw} + r_{bw})$.
# The boundary layer conductance for water vapour is equivalent to the mass transfer coefficient for a wet surface \cite[Eq. 7.41]{incropera_fundamentals_2006}:
# ##### {eq_gbw}
# $$g_{bw} = N_{Sh_L} D_{va}/L_l$$
# 
# where $N_{Sh_L}$ is the dimensionless Sherwood number and $D_{va}$ is the diffusivity of water vapour in air.
# If the convection coefficient for heat is known, the one for mass ($g_{bw}$) can readily be calculated from the relation \cite[Eq. 6.60]{incropera_fundamentals_2006}:
# ##### {eq_gbw_hc}
# $$g_{bw} = \frac{a_s h_c}{\rho_a c_{pa} N_{Le}^{1-n}}$$
# 
# where $a_s$ is the fraction of one-sided transpiring surface area in relation to the surface area for sensible heat exchange, $c_{pa}$ is the constant-pressure heat capacity of air, $n$ is an empirical constant ($n=1/3$ for general purposes) and $N_{Le}$ is the dimensionless Lewis number, defined as \cite[Eq. 6.57]{incropera_fundamentals_2006}:
# ##### {eq_Le}
# $$N_{Le} = \alpha_a/D_{va}$$
# 
# where $\alpha_a$ is the thermal diffusivity of air. The value of $a_s$ was set to 1 for leaves with stomata on one side only, and to 2 for stomata on both sides. Other values could be used for leaves only partly covered by stomata.

# In[10]:

eq_gbw = g_bw == D_va*Sh/L_l
units_check(eq_gbw).simplify_full()
eq_gbw_hc = g_bw == a_s*h_c/(rho_a*c_pa*Le^(2/3))
print units_check(eq_gbw_hc).simplify_full()
eq_Le = Le == alpha_a/D_va
units_check(eq_Le).simplify_full()


# ## Effect of leaf temperature on the leaf-air vapour concentration gradient
# \label{sec_Tleaf}
# The concentration difference in Eq. [eq_Elmol](#{eq_Elmol}) is a function of the temperature and the vapour pressure differences between the leaf and the free air. Assuming that water vapour behaves like an ideal gas, we can express its concentration as:
# ##### {eq_Cwl}
# $$C_{wl} = \frac{P_{wl}}{R_{mol} T_l}
# $$
# 
# where $P_{wl}$ is the vapour pressure inside the leaf, $R_{mol}$ is the universal gas constant and $T_l$ is leaf temperature. A similar relation holds for the vapour concentration in free air, $C_{wa} = P_{wa}/(R_{mol} T_l)$. In this study the vapour pressure inside the leaf is assumed to be the saturation vapour pressure at leaf temperature, which is computed using the Clausius-Clapeyron relation \citep[Eq. B.3]{hartmann_global_1994}:
# ##### {eq_Pwl}
# $$P_{wl} = 611 \exp \left( \frac{\lambda_E M_w}{R_{mol}} \left( \frac{1}{273} - \frac{1}{T_l} \right) \right)
# $$
# where $\lambda_E$ is the latent heat of vaporisation and $M_w$ is the molar mass of water. 
# 
# Note that the dependence of the leaf-air water concentration difference ($C_{wl} - C_{wa}$) in Eq. [eq_Cwl](#{eq_Cwl}) is very sensitive to leaf temperature. For example if the leaf temperature increases by 5K relative to air temperature, $C_{wl} - C_{wa}$ would double, while if leaf temperature decreased by 6K, $C_{wl} - C_{wa}$ would go to 0 at 70\% relative humidity.

# <p> </p>
# <p> </p>

# In[11]:

eq_Cwl = C_wl == P_wl/(R_mol*T_l)
print units_check(eq_Cwl).simplify_full()
eq_Pwl = P_wl == 611*exp(lambda_E*M_w/R_mol*(1/273 - 1/T_l))
show(eq_Pwl)


# ## Concentration or vapour pressure gradient driving transpiration?
# 
# Note that $E_{l,mol}$ is commonly expressed as a function of the vapour pressure difference between the free air ($P_{wa}$) and the leaf ($P_{wl}$), in which the conductance ($g_{tw,mol}$) is expressed in molar units (mol~m$^{-2}$~s$^{-1}$):
# ##### {eq_Elmol_conv}
# $$ E_{l,mol} =  g_{tw,mol} \frac{P_{wl} - P_{wa}}{P_a}$$
# 
# For $P_{wl} = P_{wa}$, Eq. [eq_Elmol](#{eq_Elmol}) can still give a flux, whereas Eq. [eq_Elmol_conv](#{eq_Elmol_conv}) gives zero flux. This is because the concentrations of vapour in air (mol~m$^{-3}$) can differ due to differences in tempertaure, even if the partial vapour pressures are the same (see Eq. [eq_Cwl](#{eq_Cwl})). Therefore, the relation between $g_{tw}$ and $g_{v,mol}$ has an asymptote at the equivalent temperature. It can be obtained by combining Eqs. [eq_Elmol](#{eq_Elmol}) and [eq_Elmol_conv](#{eq_Elmol_conv}) and solving for $g_{tw,mol}$:
# ##### {eq_gtwmol_gtw}
# $$g_{tw,mol} = g_{tw}\frac{P_a(P_{wa} T_l - P_{wl} T_a)}{(P_{wa} - P_{wl}) R_{mol} T_a T_l)}$$
# 
# For $T_l = T_a$, the relation simplifies to:
# ##### {eq_gtwmol_gtw_iso}
# $$g_{tw,mol} = g_{tw}\frac{P_a }{R_{mol} T_a}$$
# 
# which, for typical values of $P_a$ and $T_a$ amounts to $g_{tw,mol}\approx40$~mol~m$^{-3}g_{tw}$.  For all practical purposes, we found that Eqs. [eq_Elmol](#{eq_Elmol}) and [eq_Elmol_conv](#{eq_Elmol_conv}) with $g_{tw,mol} = g_{tw}\frac{P_a }{R_{mol} T_a}$ give similar results when plotted as functions of leaf temperature.

# In[12]:

eq_Elmol_conv = E_lmol == g_twmol*(P_wl-P_wa)/P_a
print units_check(eq_Elmol_conv).simplify_full()
soln = solve([eq_Elmol.subs(eq_Cwl, eq_Cwl(C_wl = C_wa, P_wl = P_wa, T_l = T_a)), eq_Elmol_conv], g_twmol, E_lmol)
eq_gtwmol_gtw = soln[0][0]
print units_check(eq_gtwmol_gtw).simplify_full()
eq_gtwmol_gtw_iso = eq_gtwmol_gtw(T_l = T_a).simplify_full()
print units_check(eq_gtwmol_gtw_iso).simplify_full()


# ## Model closure
# 
# Given climatic forcing as $P_a$, $T_a$, $R_s$, $P_{wa}$ and $v_w$, and leaf-specific parameters $a_s$, $a_{sH}$, $L_l$ and $g_{sw}$, we need to compute $C_{wa}$, $h_c$, $g_{bw}$ and a series of other derived variables, as described below. 
# 
# The vapour concentration in the free air can be computed from vapour pressure analogously to Eq. [eq_Cwl](#{eq_Cwl}):
# ##### {eq_Cwa}
# $$C_{wa} = \frac{P_{wa}}{R_{mol} T_a}$$
# 
# The heat transfer coefficient ($h_c$) for a flat plate can be determined using the non-dimensional Nusselt number ($N_{Nu_L}$):
# ##### {eq_hc}
# $$h_{c} = k_a\frac{N_{Nu_L}}{L_l}$$
# 
# where $k_a$ is the thermal conductivity of the air in the boundary layer and $L_l$ is a characteristic length scale of the leaf. 
# 
# For sufficiently high wind speeds, inertial forces drive the convective heat transport (forced convection) and the relevant dimensionless number is the Reynolds number ($N_{Re_L}$), which defines the balance between inertial and viscous forces \cite[Eq. 6.41]{incropera_fundamentals_2006}:
# ##### {eq_Re}
# $$N_{Re_L} = \frac{v_w L_l}{\nu_a}$$
# 
# where $v_w$ is the wind velocity (m s$^{-1}$), $\nu_a$ is the kinematic viscosity of the air and $L_l$ is taken as the length of the leaf in wind direction.
# 
# In the absence of wind, buoyancy forces, driven by the density gradient between the air at the surface of the leaf and the free air dominate convective heat exchange (``free'' or ``natural convection''). The relevant dimensionless number here is the Grashof number ($N_{Gr_L}$), which defines the balance between buoyancy and viscous forces \cite[Eqs. 9.3 and 9.65]{incropera_fundamentals_2006}:
# ##### {eq_Gr}
# $$N_{Gr_L} =\frac{g(\frac{\rho_a - \rho_{al}}{\rho_{al}}) L_l^3}{\nu_a^2}$$
# 
# where $g$ is gravity, while $\rho_a$ and $\rho_{al}$ are the densities of the gas in the atmosphere and at the leaf surface respectively. 
# 
# For $N_{Gr_L} \ll N_{Re_L}^2$, forced convection is dominant and free convection can be neglected, whereas for $N_{Gr_L} \gg N_{Re_L}^2$ free convection is dominant and forced convection can be neglected \cite[P. 565]{incropera_fundamentals_2006}. For simplicity, the analysis is limited to forced conditions, which is satisfied by considering wind speeds greater than 0.5~m~s$^{-1}$ for $5 \times 5$ cm leaves.  %See leaf_chamber4.sws or leaf_capacitance_steady_state3.sws.
# 
# The average Nusselt number under forced convection was calculated as a function of the average Reynolds number and a critical Reynolds number ($N_{Re_c}$) that determines the onset of turbulence and depends on the level of turbulence in the free air stream or leaf surface properties \cite[P. 412]{incropera_fundamentals_2006}
# 
# ##### {eq_Nu_forced_all}
# $$N_{Nu_L} = (0.037 N_{Re_L}^{4/5} - C_1)N_{Pr}^{1/3}$$
# 
# with
# ##### {eq_C1}
# $$C_1 = 0.037 C_2^{4/5} - 0.664 C_2^{1/2}$$
# 
# and
# ##### {eq_C2}
# $$C_2 = \frac{N_{Re_L} + N_{Re_c} - |N_{Re_c} - N_{Re_L}|}{2}$$
# 
# Eq. [eq_C2](#{eq_C2}) was introduced to make Eq.  [eq_Nu_forced_all](#{eq_Nu_forced_all}) valid for all Reynolds numbers, and following considerations explained in our previous work \citep{schymanski_stomatal_2013}, we chose $N_{Re_c}=3000$ in the present simulations.

# In[13]:

eq_Cwa = C_wa == P_wa/(R_mol*T_a)
print units_check(eq_Cwa)
eq_hc = h_c == k_a*Nu/L_l
print units_check(eq_hc)
eq_Re = Re == v_w*L_l/nu_a
print units_check(eq_Re)
eq_Gr = Gr == g*(rho_a - rho_al)/rho_al * L_l^3/nu_a^2
print units_check(eq_Gr)
C2 = (Re + Re_c - (abs(Re - Re_c))/2)
C1 = 0.037*C2^(4/5) - 0.664*C2^(1/2)
eq_Nu_forced_all = Nu == (0.037*Re^(4/5) - C1)*Pr^(1/3)
eq_Nu_forced_all.show()


# ### Thermodynamic variables
# In order to simulate steady state leaf temperatures and the leaf energy balance terms using the above equations, it is necessary to calculate $\rho_a$, $D_{wa}$, $\alpha_a$, $k_a$, and $\nu_a$, while $L_l$, $Re_c$ and $g_{sv}$ are input parameters, and $P_{wa}$ and $v_w$ (vapour pressure and wind speed) are part of the environmental forcing.
# $D_{wa}$, $\alpha_a$, $k_a$ and $\nu_a$ were parameterised as functions of air temperature ($T_a$) only, by fitting linear curves to published data \cite[Table A.3]{Monteith_principles_2007}:
# ##### {eq_Dwa}
# $$D_{wa} = (1.49\times 10^{-7})T_a - 1.96\times 10^{-5}$$
# 
# ##### {eq_alphaa}
# $$\alpha_a = (1.32\times 10^{-7})T_a - 1.73\times 10^{-5}$$
# 
# ##### {eq_ka}
# $$k_a = (6.84\times 10^{-5})T_a + 5.62\times 10^{-3}$$
#  
# ##### {eq_nua}
# $$\nu_a = (9\times 10^{-8})T_a - 1.13\times 10^{-5}$$
# 
# 
# Assuming that air and water vapour behave like an ideal gas, and that dry air is composed of 79\% N$_2$ and 21\% O$_2$, we calculated air density as a function of temperature, vapour pressure and the partial pressures of the other two components using the ideal gas law:
# ##### {eq_rhoa_Pa_Ta}
# $$\rho_a = \frac{n_a M_a}{V_a} = M_a \frac{P_a}{R_{mol} T_a}$$
# 
# where $n_a$ is the amount of matter (mol), $M_a$ is the molar mass (kg~mol$^{-1}$), $P_a$ the pressure, $T_a$ the temperature and $R_{mol}$ the molar universal gas constant.
# This equation was used for each component, i.e. water vapour, N$_2$ and O$_2$, where the partial pressures of N$_2$ and O$_2$ are calculated from atmospheric pressure minus vapour pressure, yielding:
# ##### {eq_rhoa_Pwa_Ta}
# $$\rho_a = \frac{M_w P_{wa} + M_{N_2} P_{N_2} + M_{O_2} P_{O_2}}{R_{mol} T_a}$$
# 
# where $M_{N_2}$ and $M_{O_2}$ are the molar masses of nitrogen and oxygen respectively, while $P_{N_2}$ and $P_{O_2}$ are their partial pressures, calculated as:
# 
# ##### {eq_PN2}
# $$P_{N_2} = 0.79(P_a - P_{wa})$$
# 
# and
# ##### {eq_PO2}
# $$P_{O_2} = 0.21(P_a - P_{wa})$$
# 

# In[14]:

eq_Dva = D_va == 1.49e-07*T_a - 1.96e-05
eq_alphaa = alpha_a == (1.32e-07)*T_a - 1.73e-05
eq_ka = k_a == 6.84e-05*T_a + 5.63e-03
eq_nua = nu_a == 9.e-08*T_a - 1.13e-05
eq_rhoa_Pwa_Ta = rho_a == (M_w*P_wa + M_N2*P_N2 + M_O2*P_O2)/(R_mol*T_a)
print units_check(eq_rhoa_Pwa_Ta)


# In[15]:

eq_Pa = P_a == P_N2 + P_O2 + P_wa
print units_check(eq_Pa)
eq_PN2_PO2 = P_N2 == 0.79/0.21*P_O2
soln = solve([eq_Pa, eq_PN2_PO2], P_O2, P_N2)
eq_PO2 = soln[0][0]
print units_check(eq_PO2)
eq_PN2 = soln[0][1]
units_check(eq_PN2)


# In[16]:

eq_rhoa = eq_rhoa_Pwa_Ta.subs(eq_PN2, eq_PO2)
units_check(eq_rhoa).simplify_full()


# ## Example calculations

# In[17]:

# Energy balance as a function of variables that are independent of leaf temperature
eqenbalTl = (eq_Rs_enbal - R_s).rhs().subs(eq_El, eq_Hl, eq_Rll).subs(eq_Elmol).subs(eq_Cwl).subs(eq_Pwl)
eqenbalTl.subs(cdict).args()


# Below, we create a copy of cdict, which is a dictionary with general constants, generated during the definition of variables at the top of the worksheet, and some example values for external forcing and leaf properties needed to compute steady-state fluxes and leaf temperature:

# In[18]:

# Input values
vdict = cdict.copy()
vdict[a_s] = 1.0    # one sided stomata
vdict[g_sw] = 0.01    
vdict[T_a] = 273 + 25.5
vdict[T_w] = vdict[T_a] # Wall temperature equal to air temperature
vdict[P_a] = 101325
rha = 1
vdict[P_wa] = rha*eq_Pwl.rhs()(T_l = T_a).subs(vdict)
vdict[L_l] = 0.03
vdict[Re_c] = 3000
vdict[R_s] = 600
vdict[v_w] = 1


# Now, we compute all derived variables, using the equations described above.

# In[19]:

# Nusselt number
vdict[nu_a] = eq_nua.rhs().subs(vdict)
vdict[Re] = eq_Re.rhs().subs(vdict)
vdict[Nu] = eq_Nu_forced_all.rhs().subs(vdict)
vdict[Nu]


# In[20]:

# h_c
vdict[k_a] = eq_ka.rhs().subs(vdict)
vdict[h_c] = eq_hc.rhs().subs(vdict)
vdict[h_c]


# In[21]:

# gbw
vdict[D_va] = eq_Dva.rhs().subs(vdict)
vdict[alpha_a] = eq_alphaa.rhs().subs(vdict)
vdict[rho_a] =  eq_rhoa.rhs().subs(vdict)
vdict[Le] =  eq_Le.rhs().subs(vdict)
vdict[g_bw] = eq_gbw_hc.rhs().subs(vdict)
vdict[g_bw]


# In[22]:

# Hl, Rll
vdict[R_ll] = eq_Rll.rhs().subs(vdict)
print vdict[R_ll]
vdict[H_l] = eq_Hl.rhs().subs(vdict)
vdict[H_l]


# In[23]:

# El
vdict[g_tw] =  eq_gtw.rhs().subs(vdict)
vdict[C_wa] = eq_Cwl.rhs()(P_wl = P_wa, T_l = T_a).subs(vdict)
vdict[P_wl] = eq_Pwl.rhs().subs(vdict)
vdict[C_wl] = eq_Cwl.rhs().subs(vdict)
vdict[E_lmol] = eq_Elmol.rhs().subs(vdict)
vdict[E_l] = eq_El.rhs().subs(eq_Elmol).subs(vdict)
vdict[E_l]


# Now, we have a dictionary with $E_l$, $H_l$ and $R_{ll}$ as functions of leaf temperature ($T_l$) only. At steady state, [eq_Rs_enbal](#{eq_Rs_enbal}) must be satisfied, so we will substitute the above dictionary into eq_Rs_enbal, subtract $R_s$ and do a numerical search for the root of the equation to obtain steady-state leaf temperature:

# In[27]:

# Numerical search for steady-state Tl
vdict[T_l] = find_root((eq_Rs_enbal - R_s).rhs().subs(vdict), 273, 373)
vdict[T_l]


# Below, we define a function that performs all the above calculations given a dictionary with the instantaneous forcing:

# In[25]:

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


# In[26]:

# Test
vdict = cdict.copy()
#vdict= {}
vdict[a_s] = 1.0    # one sided stomata
vdict[g_sw] = 0.01    
vdict[T_a] = 273 + 25.5
vdict[T_w] = vdict[T_a] # Wall temperature equal to air temperature
vdict[P_a] = 101325
rha = 1
vdict[P_wa] = rha*eq_Pwl.rhs()(T_l = T_a).subs(vdict)
vdict[L_l] = 0.03
vdict[Re_c] = 3000
vdict[R_s] = 600
vdict[v_w] = 1
resdict = fun_SS(vdict)

fun_dict_print(resdict)


# ### Calculation using known $T_l$
# If leaf temperature ($T_l$) is known (e.g. measured), but stomatal conductance $g_{sw}$ is not known, we can still calculate $H_l$ and $R_{ll}$, and then obtain $E_l$ from the energy balance equation:

# In[27]:

eq_El_enbal = solve(eq_Rs_enbal, E_l)[0]
units_check(eq_El_enbal)


# In[28]:

# Test with known Tl but no g_sw
vdict = cdict.copy()
vdict[a_s] = 1.0    # one sided stomata
   
vdict[T_a] = 273 + 25.5
vdict[T_w] = vdict[T_a] # Wall temperature equal to air temperature
vdict[P_a] = 101325
rha = 1
vdict[P_wa] = rha*eq_Pwl.rhs()(T_l = T_a).subs(vdict)
vdict[L_l] = 0.03
vdict[Re_c] = 3000
vdict[R_s] = 600
vdict[v_w] = 1

vdict[C_wa] = eq_Cwl.rhs()(P_wl = P_wa, T_l = T_a).subs(vdict)
vdict[nu_a] = eq_nua.rhs().subs(vdict)
vdict[Re] = eq_Re.rhs().subs(vdict)
vdict[Nu] = eq_Nu_forced_all.rhs().subs(vdict)
vdict[k_a] = eq_ka.rhs().subs(vdict)
vdict[h_c] = eq_hc.rhs().subs(vdict)


vdict[T_l] = 305.65
vdict[H_l] = eq_Hl.rhs().subs(vdict)
vdict[R_ll] = eq_Rll.rhs().subs(vdict)

eq_El_enbal.subs(vdict)


# # Saving definitions to separate file
# In the below, we save the definitions and variables to separate files in the /temp directory, one with the extension .sage, from which we can selectively load functions using
# `%load fun_name filenam.sage`
# and one with the extension .sobj, to be loaded elsewhere using 
# `load_session()`

# In[29]:

fun_export_ipynb('leaf_enbalance_eqs', 'temp/')
save_session('temp/leaf_enbalance_eqs')


# # Table of symbols

# In[30]:

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



