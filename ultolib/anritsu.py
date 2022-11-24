import pyvisa as visa
import numpy as np
import logging

import qcodes as qc
from qcodes import (Instrument, VisaInstrument, Parameter,
                    ManualParameter, MultiParameter,
                    validators as vals)
from qcodes.instrument.channel import InstrumentChannel

class MG3681A_base():
    instr_type='microwave source'
    
    def __init__(self) -> None:
        self.rm = visa.ResourceManager()
        self.rl = self.rm.list_resources()
        
        self._frequency=2.87e9
        self._power=-10
        self._output=False
        
        return
        
    def open(self, gpib_addr : str=r'GPIB0::3::INSTR') -> int:
        self.handle=self.rm.open_resource(gpib_addr)
        self.idn=self.handle.query('*IDN?')
        conn_str='Connected to '+self.instr_type+': '+self.idn
        print(conn_str)
        
        self.handle.write('OLDBM')
        self.handle.write('LVL ON')
        self.handle.write('PULM:STAT ON')
        
        self.output=False
        self.frequency=2.85e9
        self.power=-10
        
        return 0
    
    @property
    def frequency(self) -> float:
        f=self.handle.query('FREQ?')
        F=float(f)
        self._frequency=F
        return F
    
    @frequency.setter
    def frequency(self, f : float) -> None:
        F=str(f)
        self.handle.write('FREQ '+F+'HZ')
        self._frequency=f
        return
    
    @property
    def power(self) -> float:
        p=self.handle.query('OLVL?')
        P=float(p)
        self._power=P
        return P
        
    @power.setter
    def power(self, p : float) -> None:
        P=str(p)
        self.handle.write('OLVL '+P+'dBm')
        self._power=p
        return
    
    @property
    def output(self) -> bool:
        o=self.handle.query('LVL?')
        if o=='ON\n':
            self._output=True
            return True
        elif o=='OFF\n':
            self._output=False
            return False
        else:
            print('Communication error: Cannot read MG3681A status.')
            
    @output.setter 
    def output(self, o : bool) -> None:
        if o:
            self.handle.write('LVL ON')
            self._output=True
        else:
            self.handle.write('LVL OFF')
            self._output=False
            
            
class MG3681A(VisaInstrument):
    log=logging.getLogger(__name__)
    
    def __init__(self,name, address, **kwargs):
        super().__init__(name, address, terminator='', **kwargs)
        
        self.frequency = Parameter(
                name='frequency',
                unit='Hz',
                get_parser=float,
                get_cmd='FREQ?',
                set_cmd='FREQ {:.2f}HZ',
                instrument=self
            )
        
        self.power = Parameter(
                name='power',
                unit='dBm',
                get_parser=float,
                get_cmd='OLVL?',
                set_cmd='OLVL {:.1f}dBm',
                instrument=self
            )
        
        self.output = Parameter(
                name='output',
                unit=None,
                val_mapping={'ON': 'ON\n', 'OFF': 'OFF\n'},
                get_cmd='LVL?',
                set_cmd='LVL {}',
                instrument=self
            )
        
        self.output_level_unit = Parameter(
                name='output_level_unit',
                unit=None,
                val_mapping={'dBm': 'DBM', 'dBu': 'DBU', 
                             'Volt': 'V', 'Watt': 'W'},
                set_cmd='OL{}',
                get_cmd=None,
                instrument=self
            )
        
        self.pulse_modulation = Parameter(
                name='pulse_modulation',
                unit=None,
                val_mapping={'ON': 'ON\n', 'OFF': 'OFF\n',
                             'INT': 'INT\n', 'EXT': 'EXT\n'},
                set_cmd='PMO {}',
                get_cmd='PMO?',
                instrument=self
            )
        
        self.connect_message()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

