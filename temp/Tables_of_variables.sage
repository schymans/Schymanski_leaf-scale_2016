
# coding: utf-8

# # Table of symbolic variables used in other worksheets
# Below, we import variables and their definitions and units from other worksheets and display them in a sorted table. We also generate latex code for inclusion in manuscript.

# In[1]:

get_ipython().run_cell_magic(u'capture', u'storage', u"# The above redirects all output of the below commands to the variable 'storage' instead of displaying it.\n# It can be viewed by typing: 'storage()'\n# Setting up worksheet and importing equations for explicit leaf energy balance\nload('temp/Worksheet_setup.sage')")


# In[2]:

# List all .ipynb files
list_files = os.listdir('.')
for fname in list_files:
    if fname[-5:] == 'ipynb':
        print fname


# In[3]:

# From leaf_enbalance_eqs, E_PM_eqs and stomatal_cond_eqs

load_session('temp/leaf_enbalance_eqs.sobj')
dict_vars1 = dict_vars.copy()

load_session('temp/stomatal_cond_eqs.sobj')
dict_vars1.update(dict_vars)

load_session('temp/E_PM_eqs.sobj')
dict_vars1.update(dict_vars)

dict_vars = dict_vars1.copy()
fun_loadvars(vardict=dict_vars)   # re-loading variable definitions
udict = {}
for key1 in dict_vars.keys():
    udict[key1] = dict_vars[key1]['units']     # exporting units information from dict_vars to udict, which will be used below.


# In[4]:

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


# In[5]:

latex(table(tabledata))


# In[ ]:



