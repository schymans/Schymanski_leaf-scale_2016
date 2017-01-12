
# coding: utf-8

# In[1]:

get_ipython().run_cell_magic(u'capture', u'storage', u"# The above redirects all output of the below commands to the variable 'storage' instead of displaying it.\n# It can be viewed by typing: 'storage()'\n# Setting up worksheet and importing equations from other worksheets\nload('temp/Worksheet_setup.sage')\nload_session('temp/leaf_enbalance_eqs.sobj')\nfun_loadvars(dict_vars) # any variables defined using var2() are saved with their attributes in dict_vars. Here we re-define them based on dict_vars from the other worksheet.")


# # Equations to compute stomatal conductance based on sizes and densities of stomata
# Equation numbers in comments refer to Lehmann & Or, 2013

# <h2>Definitions of additional variables</h2>

# In[2]:

# Following Lehmann_2015_Effects
var2('B_l', 'Boundary layer thickness', meter, domain1='positive')
var2('g_sp', 'Diffusive conductance of a stomatal pore', mole/second/meter^2, domain1='positive')
var2('r_sp', 'Diffusive resistance of a stomatal pore', meter^2*second/mole, domain1='positive')
var2('r_vs', 'Diffusive resistance of a stomatal vapour shell', meter^2*second/mole, domain1='positive')
var2('r_bwmol', 'Leaf BL resistance in molar units',  meter^2*second/mole, domain1='positive', latexname = 'r_{bw,mol}')
var2('r_end', 'End correction, representing resistance between evaporating sites and pores', meter^2*second/mole, domain1='positive')
var2('A_p', 'Cross-sectional pore area', meter^2, domain1='positive')
var2('d_p', 'Pore depth', meter, domain1='positive')
var2('V_m', 'Molar volume of air', meter^3/mole, domain1='positive')
var2('n_p', 'Pore density', 1/meter^2, domain1='positive')
var2('r_p', 'Pore radius (for ellipsoidal pores, half the pore width)', meter, domain1='positive')
var2('l_p', 'Pore length', meter, domain1='positive')
var2('k_dv', 'Ratio $D_{va}/V_m$', mole/meter/second, domain1='positive', latexname='k_{dv}')
var2('s_p', 'Spacing between stomata', meter, domain1='positive')
var2('F_p', 'Fractional pore area (pore area per unit leaf area)', meter^2/meter^2, domain1='positive')


# In[3]:

docdict[V_m]


# ## Equations from Lehmann & Or, 2013

# In[4]:

# Eqn1
eq_kdv = k_dv == D_va/V_m
print units_check(eq_kdv)
eq_gsp_A = g_sp == A_p/d_p*k_dv*n_p
print units_check(eq_gsp_A)
eq_Ap_rp = A_p == pi*r_p^2
print units_check(eq_Ap_rp)
eq_rsp_A = r_sp == 1/eq_gsp_A.rhs()
print units_check(eq_rsp_A)
eq_rsp_rp = r_sp == 1/eq_gsp_A.rhs().subs(eq_Ap_rp)
print units_check(eq_rsp_rp)


# In[5]:

# Eq. 2(c)
eq_rvs_B = r_vs == (1/(4*r_p) - 1/(pi*s_p))*1/(k_dv*n_p)
units_check(eq_rvs_B)


# In[6]:

# Eq. 2(d)
eq_rvs_S = r_vs == (1/(2*pi*r_p) - 1/(pi*s_p))*1/(k_dv*n_p)
units_check(eq_rvs_B)


# In[7]:

# Below Eq. 3(b)
assume(n_p>0)
eq_sp_np = s_p == 1/sqrt(n_p)
units_check(eq_sp_np)


# In[8]:

# Eq. 2(a)
eq_rend_rp = r_end == 1/(4*r_p*k_dv*n_p)
units_check(eq_rend_rp)


# In[9]:

# Eq. 2(b)
eq_rend_PW = r_end == log(4*(l_p/2)/r_p)/(2*pi*l_p/2*k_dv*n_p)
units_check(eq_rend_PW)


# In[10]:

eq_gswmol = g_swmol == 1/(r_end + r_sp + r_vs)
units_check(eq_gswmol)


# <p>It actually does not seem to make much sense to add r_end to the resistances, as it does not reflect the geometry of the inter-cellular air space! Also eq_rvs_B is the correct one to use for stomata, as eq_rvs_S is valid for droplets, not holes!</p>

# In[11]:

assume(s_p>0)
eq_np_sp = solve(eq_sp_np, n_p)[0]
units_check(eq_np_sp)


# In[12]:

# Eq. 5
eq_rbwmol = r_bwmol == B_l/k_dv
units_check(eq_rbwmol)


# In[13]:

docdict[k_dv]


# Below, we compare our system of equations to Eq. 7a in Lehmann & Or (2015), i.e. we check if $\frac{r_{bwmol}}{r_{tot}} = 1/(1 + s_p^2/B_l*(1/(2*r_p) + d_p/(\pi*r_p^2)))$:

# In[14]:

# Verification of Eq. 7a
eq1 = (r_bwmol/(2*r_end + r_sp + r_bwmol)).subs(eq_rend_rp, eq_rsp_rp, eq_rvs_B, eq_rbwmol).subs(eq_np_sp)
eq1.simplify_full().show()
eq7a = 1/(1 + s_p^2/B_l*(1/(2*r_p) + d_p/(pi*r_p^2)))
eq7a.show()
(eq1 - eq7a).simplify_full()


# ## Additional equations for calculating conductances of laser perforated foils

# In[15]:

# Pore area for circular pores
eq_Ap = A_p == pi*r_p^2
units_check(eq_Ap)


# In[16]:

# To deduce pore radius from pore area, assuming circular pore:
eq_rp_Ap = solve(eq_Ap_rp, r_p)[0]
units_check(eq_rp_Ap)


# In[17]:

# For regular pore distribution we can derive the pore density from porosity and pore area
eq_np = n_p == F_p/A_p


# In[18]:

# Conversion from mol/m2/s to m/s, i.e. from g_swmol to g_sw
eq_gsw_gswmol = solve(eq_gtwmol_gtw(g_tw = g_sw, g_twmol = g_swmol), g_sw)[0]
units_check(eq_gsw_gswmol)


# In[19]:

# Simplification, neglecting effect of leaf-air temperature difference
eq_gsw_gswmol_iso = eq_gsw_gswmol(T_l = T_a).simplify_full()
units_check(eq_gsw_gswmol_iso)


# In[20]:

# Saving session so that it can be loaded using load_session().
save_session('temp/stomatal_cond_eqs.sobj')


# <h2>Table of symbols</h2>

# In[21]:

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

