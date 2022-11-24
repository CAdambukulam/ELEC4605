# -*- coding: utf-8 -*-
import pyvisa as visa
import numpy as np
import logging

import qcodes as qc
from qcodes import (Instrument, VisaInstrument, Parameter,
                    ManualParameter, MultiParameter, InstrumentChannel,
                    validators as vals)

class KD3305P(VisaInstrument):
    log=logging.getLogger(__name__)
    
    def __init__(self, name, address, **kwargs):
        super().__init__(name, address, terminator='\n', **kwargs)
        
        channels=[1,2]
        for ch_num in channels:
            ch_name='ch{:d}'.format(ch_num)
            ch=KD3305P_channel(parent=self, name=ch_name,
                                channel=ch_num)
            self.add_submodule(name=ch_name, submodule=ch)
        
        self.connect_message()
            
        
        
class KD3305P_channel(InstrumentChannel):
    
    def __init__(self, parent : Instrument, name : str, channel : int):
        
        if channel not in [1, 2]:
            raise ValueError('channel must be either "(int) 1" or "(int) 2"')
            
        super().__init__(parent, name)
    
        self.voltage_setpoint = Parameter(
            name='voltage_setpoint',
            unit='V',
            get_cmd='VSET{:d}?'.format(channel),
            get_parser=float,
            set_cmd='VSET{:d}:{}'.format(channel, '{:.2f}'),
            set_parser=float,
            instrument=self
            )
        
        self.voltage_output = Parameter(
            name='voltage_output',
            label='Voltage',
            unit='V',
            get_cmd='VOUT{:d}?'.format(channel),
            get_parser=float,
            instrument=self
            )
        
        self.current_setpoint = Parameter(
            name='current_setpoint',
            unit='A',
            get_cmd='ISET{:d}?'.format(channel),
            get_parser=float,
            set_cmd='ISET{:d}:{}'.format(channel, '{:.2f}'),
            set_parser=float,
            instrument=self
            )
        
        self.current_output = Parameter(
            name='current_output',
            label='Current',
            unit='V',
            get_cmd='IOUT{:d}?'.format(channel),
            get_parser=float,
            instrument=self
            )
        
        
        