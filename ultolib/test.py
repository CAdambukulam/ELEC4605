import qcodes as qc
import qcodes.utils.validators as vals
from qcodes.plots.qcmatplotlib import MatPlot
from qcodes.loops import Loop
from qcodes.data.data_set import load_data
from qcodes.actions import Task, Wait
from qcodes.measure import Measure


fit_values = qc.ManualParameter(name='fit_values', unit='V')
def set_fit_curve(ax_val):
    fit_values.set(ax_val ** 2)
    
fit_axis = qc.Parameter(name='fit_axis', label='tau', 
                        unit='s', set_cmd=set_fit_curve)

loop = Loop(fit_axis.sweep(0, 1, num=10)).each(fit_values)
data_T1_fit = loop.get_data_set(name='T1_fit')