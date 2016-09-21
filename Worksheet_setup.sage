
# coding: utf-8

# # Importing modules, setting options and defining custom functions for use in other worksheets
# To import into another worksheet:
# 
# ```
# fun_include_ipynb('Worksheet_setup')
# ```
# or
# ```
# load('Worksheet_setup.sage')
# ```
# Note: This worksheet is prepared for the open source software [sage](http://www.sagemath.org). In sage 7.3, nbconvert is missing dependencies. To rectify, enter the following commands in a terminal:
# 
# ```
# sage --sh
# pip install entrypoints
# pip install configparser
# ```
# 

# In[1]:

get_ipython().magic(u'matplotlib inline')
import numpy as np  # for data array operations
import pylab   # for plotting
import json   # for reading notebook metadata

list_plot.options["frame"]=True
list_plot.options["axes"]=False
plot.options["frame"]=True
plot.options["axes"]=False

try:
    for s in units.amount_of_substance.trait_names():
        globals()[s] = getattr(units.amount_of_substance,s) 
    for s in units.energy.trait_names():
        globals()[s] = getattr(units.energy,s)
    for s in units.length.trait_names():
        globals()[s] = getattr(units.length,s)    
    for s in units.mass.trait_names():
        globals()[s] = getattr(units.mass,s)
    for s in units.pressure.trait_names():
        globals()[s] = getattr(units.pressure,s)
    for s in units.temperature.trait_names():
        globals()[s] = getattr(units.temperature,s)
    for s in units.time.trait_names():
        globals()[s] = getattr(units.time,s)
    for s in units.power.trait_names():
        globals()[s] = getattr(units.power,s)
except:     # in newer versions of sage, need to use _tab_completion() instead of trait_names()
    for s in units.amount_of_substance._tab_completion():
        globals()[s] = getattr(units.amount_of_substance,s) 
    for s in units.energy._tab_completion():
        globals()[s] = getattr(units.energy,s)
    for s in units.length._tab_completion():
        globals()[s] = getattr(units.length,s)    
    for s in units.mass._tab_completion():
        globals()[s] = getattr(units.mass,s)
    for s in units.pressure._tab_completion():
        globals()[s] = getattr(units.pressure,s)
    for s in units.temperature._tab_completion():
        globals()[s] = getattr(units.temperature,s)
    for s in units.time._tab_completion():
        globals()[s] = getattr(units.time,s)
    for s in units.power._tab_completion():
        globals()[s] = getattr(units.power,s)        


udict = {}
cdict = {}
docdict = {}

def var2(name,doc='',units=1,domain1='real',latexname='',value = 1234):
    '''
    Creates a symbolic variable in the given domain (standard:'real') and with the given 
    latexname. Further, it adds the string doc to docdict, the units to udict and the value
    to cdict. All entries are optional except for name.
    Usage example: 
        sage: var2('name','doc',watt/meter^2,'positive','N_{ame}',1.0)
    '''
    if len(latexname)>0:
        z = var(name,domain = domain1,latex_name = latexname)   
    else:
        z = var(name,domain = domain1) 
    if len(doc)>0:
        docdict[var(name)] = doc
    udict[var(name)] = var(name)*units
    if value != 1234:
        cdict[var(name)] = value
    return z
    
# Workaround from http://ask.sagemath.org/question/10260/convert-derived-si-units-to-base-units/
def convert(expr):
     op = expr.operator()
     ops = expr.operands()
     if op:
        if len(ops) == 2:
            return op(*map(convert, ops))
        else:
            return op(convert(ops[0]), reduce(op, map(convert, ops[1:])))
     else:
        return expr.convert() if hasattr(expr, 'convert') else expr
    
def units_check(eq, sfull = True):
    '''
    Checks whether all arguments are keys in udict and returns simplified units
    '''
    for blah in eq.arguments():
        udict[blah]
    eq.show()
    if sfull:
        return convert(eq.subs(udict)/eq).simplify_full()
    else:
        return convert(eq.subs(udict)/eq)

    
def plot_fit(xdata, ydata, model, initial_guess = None, parameters = None, variables = None):
    '''Fits the parameters of model to lists of xdata and ydata and plots the results
       Example:
          #linear regression through origin
        var('x m')
        model(x) = m*x
        plot_fit(xdata,ydata,model)
    '''
    data = zip(xdata,ydata)
    xdata1 = vector(xdata)
    ydata1 = vector(ydata)
    p1 = find_fit(data, model,solution_dict=True)
    pkeys = p1.keys()
    p2 = {}
    for i in pkeys: p2[i] = n(p1[i],digits=3)
    fitvals = vector([n(model(x1).subs(p1)) for x1 in xdata])
    # root mean square deviation
    rmsd = sqrt(mean(ydata1 - fitvals)^2)
    nrmsd = rmsd/(max(ydata) - min(ydata))
    blah = text("y="+str(model(x).subs(p2)),(0,0),fontsize=10,rgbcolor=(1,0,0))
    print "y="+str(model(x).subs(p2))
    print "NRMSD = "+str(n(nrmsd,digits=5))
    lentext = blah.xmax() - blah.xmin()
    # plot
    P = list_plot(data, frame=True, axes = False, faceted=True, color = "white")
    P += plot(model.subs(p1),(x,min(xdata),max(xdata)),color = 'black', legend_label = "y="+str(model(x).subs(p2)) + '\n ' + "NRMSD = "+str(n(nrmsd,digits=5)))   
    return P

def dict_print(vdict1, list_names = None):
    dtype1 = [('name', 'S10'), ('val', object)]
    vdict2 = vdict1.copy()
    if list_names:
        vdict2 = {}
        for name1 in list_names:
            vdict2[name1] = vdict1[name1]
    list1 = np.array(vdict2.items(), dtype = dtype1)
    list2 = np.sort(list1,order = 'name')
    for k, v in list2:
            print "{:<15} {:<15}".format(k,v)    
            
def fun_dict_print(vdict1):
    dtype1 = [('name', 'S10'), ('val', object)]
    list1 = np.array(vdict1.items(), dtype = dtype1)
    list2 = np.sort(list1,order = 'name')
    for k, v in list2:
            print "{:<15} {:<15}".format(k,v)
            
def fun_dict(dict1,dict2):
    '''
    Updates all entries of dict1 with entries co-occurring in dict2
    and returns the udated dict1.
    '''
    dict1.update(dict2)
    return dict1

def fun_comp_array(expr1, nparray):
    """
    Returns list of values for expr1 
    corresponding to each row of nparray
    where all arguments in expr1 are substituted
    by the corresponding values in nparray.
    """
    list_args = expr1.args()
    return [expr1.subs(dict(zip(list_args, [nparray[str(dummy)][i] for dummy in list_args]))) for i in srange(len(nparray))]
    
def fun_printvar(varname):
    pstring = str(varname)+': ' + docdict[varname]
    try: 
        value = cdict[varname]
        pstring1 = pstring + ', ' + str(n(value, prec = 16)) + ' ' + str(udict[varname]/varname)
    except: pstring1 = pstring + ', ' + str(udict[varname]/varname)
    print pstring1
    return

def fun_eq(eq, sfull = True):
    print units_check(eq)
    return eq

def fun_include_ipynb(worksheet, del1 = True, output = True):
    '''
    Exports worksheet.ipynb to a sage file and loads it into current worksheet.
    Example:
    fun_include_ipynb('Worksheet_setup')
    Need to import json before using this function. Unless option output=True is given,
    each line of text has a semicolon added in the .py file, preventing any output.
    '''
    str1 = 'jupyter nbconvert  --to=python \'' + worksheet+'.ipynb\''
    print str1
    print 'Exporting specified worksheet to .py file...'
    
    try:
        retcode = os.system(str1)
        if retcode < 0:
            print >>sys.stderr, "nbconvert was terminated by signal", -retcode
        else:
            print >>sys.stderr, "nbconvert returned", retcode
    except OSError as e:
        print >>sys.stderr, "Execution failed:", e
        print >>sys.stderr, "Trying ipython nbconvert instead..."
        str1 = 'ipython nbconvert  --to script \'' + worksheet+'.ipynb\''    # for new version of python
        try:      
            retcode = os.system(str1)
            if retcode < 0:
                print >>sys.stderr, "nbconvert was terminated by signal", -retcode
            else:
                print >>sys.stderr, "nbconvert returned", retcode
        except OSError as e:
            print >>sys.stderr, "Execution failed:", e              
        
    str1 = worksheet + '.py'
    
    print 'Checking if specified ipynb file was run with sage kernel...'
    with open(worksheet+'.ipynb') as data_file:    
        data = json.load(data_file)
    
    if data['metadata']['kernelspec']['name'][0:4] == 'sage':
        print 'Renaming .py file to .sage if notebook kernel was sage (to avoid exponent error)'
        str2 = worksheet + '.sage'
        os.rename(str1, str2)
        str1 = str2
    if output == False:
        input_file = open(str1, "r")
        str2 = '0_'+ str1
        output_file = open(str2, "w")
        for line in input_file:
            line = line.strip() # remove line end marks
            final = line+'\n'
            if len(line)>0:
                if line[-1] != ';':
                    final= line+';\n'
            output_file.write(final)
        input_file.close()
        os.remove(str1)
        output_file.close()
        str1 = str2
        
        
    print 'Loading file into current worksheet...'
    load(str1)
    if del1:
        print 'Removing file generated above...'
        os.remove(str1)


# Below, we export the above commands to a file called Worksheet_setup.sage, which can be loaded in other worksheets by typing:
# 
# `load('Worksheet_setup.sage')`

# In[2]:

str1 = 'ipython nbconvert  --to=python \'./Worksheet_setup.ipynb\''
print str1
print 'Exporting Worksheet_setup.ipynb to .py file...'

try:
    retcode = os.system(str1)
    if retcode < 0:
        print >>sys.stderr, "nbconvert was terminated by signal", -retcode
    else:
        print >>sys.stderr, "nbconvert returned", retcode
except OSError as e:
    print >>sys.stderr, "Execution failed:", e
    print >>sys.stderr, "Trying jupyter nbconvert instead..."
    str1 = 'jupyter nbconvert  --to script \'./Worksheet_setup.ipynb\''    # for new version of python
    try:      
        retcode = os.system(str1)
        if retcode < 0:
            print >>sys.stderr, "nbconvert was terminated by signal", -retcode
        else:
            print >>sys.stderr, "nbconvert returned", retcode
    except OSError as e:
        print >>sys.stderr, "Execution failed:", e        


os.system(str1)

str1 = 'Worksheet_setup.py'
str2 = 'Worksheet_setup.sage'
os.rename(str1, str2)


# In[ ]:



