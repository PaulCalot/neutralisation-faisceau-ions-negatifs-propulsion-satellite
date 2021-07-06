import matplotlib.pyplot as plt
import numpy as np

# bokeh plotting
from bokeh.io import push_notebook, show, output_notebook
from bokeh.plotting import figure
from bokeh.layouts import row, column
from bokeh.models import PrintfTickFormatter
from bokeh.models import CustomJS, DateRangeSlider
from bokeh.models import Legend, ColumnDataSource, Label, LabelSet, Range1d
from bokeh.palettes import magma, viridis
output_notebook(hide_banner=True)

from ipywidgets import interact, IntSlider, Dropdown, FloatSlider

# ---------------------------- Q2) MC ---------------------------- #

def plot_hist(array, bins, name, law = None, name_law = '', savefig = None):
    mean, std = np.mean(array), np.std(array)
    m, M = np.min(array), np.max(array)
    q10, q90 = np.quantile(array, q=0.1),np.quantile(array, q=0.9)
    
    dico = {
        'Mean' : mean,
        'Standard deviation' : std,
        'Minimum' : m,
        'Maximum' : M,
        'Quantile 0.1':q10,
        'Quantile 0.9':q90
    }
    
    N = array.shape[0]
    
    fig, ax = plt.subplots(figsize = (15,10))

    if(law!=None):
        x = np.arange(m, M)
        y = law(x, mean, std) # scale = 1/lambda = Esp(U)
        ax.plot(x,y, label = name_law)
        
    ax.hist(array, bins = bins , density = True)
    
    ax.set_title('Law of {} ; {} draws.'.format(name, N), fontsize=12)
    ax.set_xlabel(name, fontsize=10)
    ax.set_ylabel('quantity', fontsize=10)
    
    for key, value in dico.items():
        ax.plot([],[], '.', label = '{} : {}'.format(key, value), color = 'k')
        
    plt.legend(loc='upper left', fontsize='x-large', title_fontsize='x-small', ncol=1, bbox_to_anchor=(1.05, 1))

    if(savefig!=None):
        plt.savefig(savefig+'.png', dpi = 40)
    else:
        plt.show()
        
# --------------------------- Q3) Morris ------------------------- # 
def plot_morris(mu, sigma, names = 'default', verbose = False):
    if(names == 'default'):
        names = ['pa', 'pIH', 'pIU', 'pHD', 'pHU', 'pUD', 'NI', 'NH', 'NU', 'R0', 'mu', 'N', 't0', 'Im0', 'lambda1']
    
    if(verbose):
        dtype = [('names', 'S10'), ('mu', float), ('sigma', float)]
        values = [(n, m, s) for n,m,s in zip(names, mu, sigma)]
        arr = np.array(values, dtype=dtype)
        arr = np.sort(arr, order = 'mu')
        print('Parameters sorted by increasing mu:')
        pprint(arr)
        
    colors = magma(len(mu)) # https://www.geeksforgeeks.org/python-bokeh-colors-class/

    source = ColumnDataSource(data=dict(sigma=sigma,
                                    mu=mu,
                                    names=names,
                                    colors=colors))
    
    p = figure(plot_height=450, plot_width=900, title="Morris - sensibility analysis")

        
    p.scatter(x='mu', y='sigma', size=5, source=source, color ='colors')
    
    p.xaxis[0].axis_label = 'Mean'
    p.yaxis[0].axis_label = 'Standard deviation'

    labels = LabelSet(x='mu', y='sigma', x_offset=2, y_offset=2, text='names', source=source, render_mode='canvas')

    
    p.add_layout(labels)

    show(p)
    
# --------------------------- Q4,5) Sobol  ------------------------- # 
def plot_sobol(S, St, names= 'default', title = 'Sobol analysis', verbose = False):
    if(names == 'default'):
        names = ['Ei', 'Ti']
    
    if(verbose):
        dtype = [('names', 'S10'), ('S', float), ('St', float)]
        values = [(n, m, s) for n,m,s in zip(names, S, St)]
        arr = np.array(values, dtype=dtype)
        arr = np.sort(arr, order = 'St')
        print('Parameters sorted by increasing $S_t$:')
        pprint(arr)
    source = ColumnDataSource(data=dict(S=S,
                                    St=St,
                                    names=names))
    
    p = figure(x_range = names, plot_height=450, plot_width=900, title=title)

    p.circle(x='names', y='S', size=5, source=source, color='blue', legend_label = 'S')
    p.circle(x='names', y='St', size=5, source=source, color='red', legend_label = 'St')

    p.xaxis[0].axis_label = 'parameters'
    p.yaxis[0].axis_label = 'Index value'

    show(p)
    
def plot_sobol_time(S, St, names= 'default'):
    if(names == 'default'):
        names = ['pa', 'pIH', 'pIU', 'pHD', 'pHU', 'pUD', 'NI', 'NH', 'NU', 'R0', 'mu', 'N', 't0', 'Im0', 'lambda1']

    times = np.arange(0,St.shape[-1],1)
    
    p = figure(plot_height=450, plot_width=900, title="Sobol analysis")
    
    colors = magma(St.shape[0]) # https://www.geeksforgeeks.org/python-bokeh-colors-class/

    legend_it = []

    for k in range(St.shape[0]):
        source = ColumnDataSource(data=dict(S=S[k],
                                        St=St[k],
                                           times = times))
        c1 = p.x(x='times', y='S', source=source, color=colors[k], line_width=1.5, alpha=1.0, muted_color=colors[k], muted_alpha=0.01)
        c2 = p.line(x='times', y='St', source=source, color=colors[k], line_width=1.5, alpha=1.0, muted_color=colors[k], muted_alpha=0.01)

        p.xaxis[0].axis_label = 'time (days)'
        p.yaxis[0].axis_label = 'Index value'
        
        legend_it.append((names[k], [c1,c2]))
        
    c0 = p.line([],[],color = 'black')
    c1 = p.x([],[],color='blue')
    c2 = p.line([],[],color='blue')
    legend_it.append(('',[c0]))
    legend_it.append(('S',[c1]))
    legend_it.append(('St',[c2]))
    
    legend = Legend(items=legend_it)

    legend.click_policy="mute"
    p.add_layout(legend, 'right')
    
    show(p)