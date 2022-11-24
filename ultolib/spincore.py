import logging

import qcodes as qc
import ctypes
import numpy as np
import matplotlib.pyplot as plt


from qcodes import Instrument, Parameter
from qcodes.instrument.channel import InstrumentChannel
import qcodes.utils.validators as vals

from ultolib import spincore_base as api
#import spincore_base as api

class PulseBlasterESRPRO(Instrument):
    """
    This is the qcodes driver for the PulseBlaster ESR-PRO.
    The driver communicates with the underlying python wrapper spinapi.py,
    which in turn communicates with spinapi.dll.
    Both can be obtained from the manufacturer's website.
    Args:
        name (str): name of the instrument
        api_path(str): Path of the spinapi.py file
    Note that only a single instance can communicate with the PulseBlaster.
    To avoid a locked instrument, always close connection with the Pulseblaster.
    """
    log=logging.getLogger(__name__)
    program_instructions_map = {
        'CONTINUE': 0,  #inst_data=Not used
        'STOP': 1,      #inst_data=Not used
        'LOOP': 2,      #inst_data=Number of desired loops
        'END_LOOP': 3,  #inst_data=Address of instruction originating loop
        'JSR': 4,       #inst_data=Address of first instruction in subroutine
        'RTS': 5,       #inst_data=Not Used
        'BRANCH': 6,    #inst_data=Address of instruction to branch to
        'LONG_DELAY': 7,#inst_data=Number of desired repetitions
        'WAIT': 8}      #inst_data=Not used
    
    FLAG_MASK = 0b111 << 21

    def __init__(self, name, board_number=0, **kwargs):
        super().__init__(name, **kwargs)
    
        
        # It seems that the core_clock is not the same as the sampling rate.
        # At core_clock(500), the PulseBlaster uses 1 ns per wait duration.
        # The wait duration is inversely proportional to the core clock, in contrast to the sampling rate
        self.add_parameter('core_clock',
                           label='Core clock',
                           unit='MHz',
                           set_cmd=self.set_core_clock,
                           vals=vals.Numbers(0, 500))

        self.add_parameter('board_number',
                           set_cmd=None,
                           initial_value=board_number)

        self.add_function('initialize',
                          call_cmd=self.initialize)

        self.add_function('detect_boards',
                          call_cmd=self.detect_boards)

        self.add_function('select_board',
                          call_cmd=api.pb_select_board,
                          args=[vals.Enum(0, 1, 2, 3, 4)])

        self.add_function('start_programming',
                          call_cmd=self.start_programming)

        self.add_function('send_instruction',
                          call_cmd=self.send_instruction,
                          args=[vals.Ints(), vals.Strings(),
                                vals.Ints(), vals.Ints()])

        self.add_function('stop_programming',
                          call_cmd=self.stop_programming)

        self.add_function('start',
                          call_cmd=self.start)

        self.add_function('stop',
                          call_cmd=self.stop)

        self.add_function('close',
                          call_cmd=self.close)

        self.add_function('get_error',
                          call_cmd=api.pb_get_error)
        
        self.add_function('reset',
                          call_cmd=api.pb_reset)
        
        self.add_function('set_debug',
                          call_cmd=api.pb_set_debug,
                          args=[vals.Ints(0,1)])

        self.add_parameter('instruction_sequence',
                           set_cmd=None,
                           initial_value=[],
                           vals=vals.Anything(),
                           snapshot_value=False
                           )
        
        self.add_parameter('n_channels',
                           set_cmd=None,
                           initial_value=4,
                           snapshot_value=None)
        
        self._channel_list=[]
        for ch_num in range(0, self.n_channels()):
            ch_name='ch{:d}'.format(ch_num)
            ch=PulseBlasterESRPRO_channel(parent=self, name=ch_name,
                                channel=ch_num)
            self.add_submodule(ch_name, ch)
            self._channel_list.append(ch)

        self.setup(initialize=True)

    def initialize(self, debug : int=1):
        '''
        Initializes board. This needs to be performed before any communication with the board is possible
        Raises error if return message indicates an error.
        Returns:
            return_msg
        '''
        self.select_board(self.board_number())
        self.set_debug(debug)
        return_msg = api.pb_init()
        assert return_msg == 0, \
            'Error initializing board: {}'.format(api.pb_get_error())
        return return_msg

    def detect_boards(self):
        '''
        Detects the number of boards.
        Raises an error if the number of boards is zero, or an error has occurred.
        Returns:
            return_msg
        '''
        return_msg = api.pb_count_boards()
        assert return_msg > 0, 'No boards detected'
        return return_msg

    def setup(self, initialize=False):
        """
        Sets up the board, must be called before programming it
        Args:
            initialize: Whether to initialize (should only be done once at
            the start). False by default
        Returns:
        """
        self.detect_boards()
        if initialize:
            self.initialize()


    def set_core_clock(self, core_clock):
        self.select_board(self.board_number())
        # Does not return value
        api.pb_core_clock(core_clock)

    def start_programming(self):
        '''
        Indicate the start of device programming, after which instruction commands can be sent using PB.send_instruction
        Raises error if return message indicates an error.
        Returns:
            return_msg
        '''

        # Reset instruction sequence
        self.instruction_sequence([])

        # Needs constant PULSE_PROGRAM, which is set to equal 0 in the api)

        self.select_board(self.board_number())
        return_msg = api.pb_start_programming(api.PULSE_PROGRAM)
        assert return_msg == 0, \
            'Error starting programming: {}'.format(api.pb_get_error())
        return return_msg

    def send_instruction(self, flags, instruction, inst_args, length, log=True):
        '''
        Send programming instruction to Pulseblaster.
        Programming instructions can only be sent after the initial command pb.start_programming.
        Raises error if return message indicates an error.
        The different types of instructions are:
            ['CONTINUE', 'LOOP', 'END_LOOP', 'JSR', "RTS', 'BRANCH', 'LONG_DELAY', "WAIT']
        See manual for detailed description of each of the commands
        Args:
            flags: Bit representation of state of each output 0=low, 1=high
            instruction: Instruction to be sent, case-insensitive (see above for possible instructions)
            inst_args: Accompanying instruction argument, dependent on instruction type
            length: Number of clock cycles to perform instruction
            log: Whether to log to instruction_sequence (True by default)
        Returns:
            return_msg, which contains instruction address
        '''
        flags |= self.FLAG_MASK
        # Add instruction to log
        if log:
            self.instruction_sequence(self.instruction_sequence() +
                              [(bin(flags), instruction, inst_args, length)])

        instruction_int = self.program_instructions_map[instruction.upper()]
        #Need to call underlying spinapi because function does not exist in wrapper
        self.select_board(self.board_number())
        return_msg = api.spinapi.pb_inst_pbonly(ctypes.c_uint64(flags),
                                                ctypes.c_int(instruction_int),
                                                ctypes.c_int(inst_args),
                                                ctypes.c_double(length))
        assert return_msg >= 0, \
            'Error sending instruction: {}'.format(api.pb_get_error())
        return return_msg

    def send_instructions(self, *instructions):
        for instruction in instructions:
            self.send_instruction(*instruction, log=False)

        # Add instructions to log
        self.instruction_sequence(list(self.instruction_sequence()) + list(instructions))

    def stop_programming(self):
        '''
        Stop programming. After this function, instructions may not be sent to PulseBlaster
        Raises error if return message indicates an error.
        Returns:
            return_msg
        '''
        self.select_board(self.board_number())
        return_msg = api.pb_stop_programming()
        if return_msg != 0:
            print('Error starting programming: {}'.format(api.pb_get_error()))
        return return_msg

    def start(self):
        '''
        Start PulseBlaster sequence.
        Raises error if return message indicates an error.
        Returns:
            return_msg
        '''
        self.select_board(self.board_number())
        return_msg = api.pb_start()
        try:
            assert return_msg == 0, 'Error starting: {}'.format(api.pb_get_error())
        except AssertionError:
            api.pb_stop()
            api.pb_start()
            assert return_msg == 0, 'Error starting: {}'.format(
                api.pb_get_error())
        return return_msg

    def stop(self):
        '''
        Stop PulseBlaster sequence
        Raises error if return message indicates an error.
        Returns:
            return_msg
        '''
        self.select_board(self.board_number())
        return_msg = api.pb_stop()
        assert return_msg == 0, 'Error Stopping: {}'.format(api.pb_get_error())
        return return_msg
    
    def flush_channel_buffer(self, repeat : bool=True, contiuous : bool=True):
        
        N=self.n_channels()
        #Extract bit levels and durations from the PulseBlasterESRPRO_channel class.
        l_table=[[p.level for p in ch.pulse_sequence_buffer()] \
                 for ch in self._channel_list]
        d_table=[[]] * N
        for i in range(0, N):
            pulse_sequence=self._channel_list[i].pulse_sequence_buffer()
            M=len(pulse_sequence)
            d_table[i]=[0] * M
            for j in range(0, M):
                pulse_duration=pulse_sequence[j].duration
                if isinstance(pulse_duration, float):
                    d_table[i][j]=int(pulse_duration * 1e8) * 10
                elif isinstance(pulse_duration, Parameter):
                    d_table[i][j]=round(pulse_duration() * 1e8) * 10
                    if pulse_duration.unit != 's':
                        raise ValueError('PulseBlasterESRPRO: flush: Pulse duration'+
                                         ' unit must be "s" (seconds).')
                else:
                    raise TypeError('PulseBlasterESRPRO: flush: Pulse duration'+
                                    ' must be either of type "float" or "parameter".')
                    
        #Convert the pulse durations into time from start instruction.
        t_table=[[]] * N
        for i in range(0,N):
            M=len(d_table[i])
            t_line=[float] * (M+1)
            t_line[0]=0
            for j in range(1,M+1):
                t_line[j]=d_table[i][j-1]+t_line[j-1]
            t_table[i]=t_line
        
        end_time=t_table[0][-1]
        for i in range(0, N):
            if t_table[i][-1] > end_time:
                end_time=t_table[i][-1]
        
        #Get the bit flags and global time from the level and time tables.
        time=[0]
        flags=[int(0)]
        for i in range(0,N):
            M=len(d_table[i])
            k=0
            for j in range(0,M):
                K=len(time)
                while time[k] < t_table[i][j]:
                    k+=1
                    if k==K:
                        break
                if k==K:
                    time.append(t_table[i][j])
                    flags.append(flags[k-1] & ~(1 << i))
                elif time[k] > t_table[i][j]:
                    time.insert(k, t_table[i][j])
                    flags.insert(k, flags[k-1] & ~(1 << i))
                else:
                    flags[k]=flags[k] & ~(1 << i)
                flags[k]=flags[k] | (l_table[i][j] << i)
        time.append(end_time)
        duration=np.diff(time)  
        #Convert duration into nanoseconds
        #duration=np.floor(duration * 1e8) * 10
        
        inst=['CONTINUE'] * len(duration)
        inst_args=[0] * len(duration)
        
        #Check for whether a long delay is required.
        for i in range(0, len(duration)-1):
            delay_bits=int(np.ceil(duration[i] / 4))
            if delay_bits > 2**32-1:
                inst[i]='LONG_DELAY'
                inst_args[i]=delay_bits >> 32
                duration[i]/=float(delay_bits)
                
        self.start_programming()
        if len(duration) > 1:
            start=self.send_instruction(flags[0], inst[0], inst_args[0],
                                        length=float(duration[0]), log=True)
        else:
            start=0
        if start < 0: 
            print('Error starting programming: {}'.format(api.pb_get_error()))
        
        delay_bits=int(np.ceil(duration[-1] / 4))
        if delay_bits > 2**32-1:
            duration[-1]-=20
            delay_bits=int(np.ceil(duration[-1] / 4))
            inst[-1]='LONG_DELAY'
            inst_args[-1]=delay_bits >> 32
            duration[-1]/=float(delay_bits)
            flags.append(flags[-1])
            inst.append('BRANCH')
            inst_args.append(start)
            duration.append(20)
        else:
            inst[-1]='BRANCH'
            inst_args[-1]=start
        
        if len(duration)==1:
            sval=0
        else:
            sval=1
        for i in range(sval,len(duration)):
            self.send_instruction(flags[i], inst[i], inst_args[i], 
                                             float(duration[i]), log=True)
        
        self.stop_programming()
        
        if contiuous==True:
            self.start()
        return 
    
    def reset_channel_buffer(self):
        for i in range(0, self.n_channels()):
            self._channel_list[i].pulse_sequence_buffer([])
        return
    
    def plot_channel_buffer(self):
        #Get non-empty channel buffers
        valid_channel_list=[]
        for i in range(0, self.n_channels()):
            if self._channel_list[i].pulse_sequence_buffer() != []:
                valid_channel_list.append(self._channel_list[i])
        
        (fig,ax)=plt.subplots(nrows=len(valid_channel_list), ncols=1)
        fig.tight_layout()
        for i in range(len(valid_channel_list)):
            time_array=[]
            level_array=[]
            current_level=0
            next_level=0
            current_time=0
            for pulse in valid_channel_list[i].pulse_sequence_buffer():
                next_level=pulse.level
                time_array+=[current_time, current_time+10e-9]
                level_array+=[current_level, next_level]
                current_level=next_level
                if isinstance(pulse.duration, Parameter):
                    current_time+=pulse.duration()
                else:
                    current_time+=pulse.duration
            time_array+=[current_time]
            level_array+=[current_level]
            if len(valid_channel_list)==1:
                ax.plot(time_array, level_array)
                ax.set_xlabel('Time (s)')
                ax.set_ylabel('Level')
                ax.set_title('Ch{:d}'.format(valid_channel_list[i].channel))
            else:
                ax[i].plot(time_array, level_array)
                plt.subplots_adjust(hspace=0.5 + \
                                    (len(valid_channel_list)**2) * 0.15)
                ax[i].set_xlabel('Time (s)')
                ax[i].set_ylabel('Level')
                ax[i].set_title('Ch{:d}'.format(valid_channel_list[i].channel))
            
    def close(self):
        '''
        Terminate communication with PulseBlaster.
        After this command, communication is no longer possible with PulseBlaster
        Returns:
            None
        '''
        self.select_board(self.board_number())
        api.pb_close()
        super().close()
        
        
        
       
class PulseBlasterESRPRO_channel(InstrumentChannel):
    
    def __init__(self, parent : Instrument, name : str, channel : int):
        super().__init__(parent, name)
        
        self.channel=channel
        
        self._pulse_sequence=[]
        
        self.pulse_sequence_buffer=Parameter(
                name='pulse_sequence_buffer',
                initial_value=[],
                set_cmd=self._set_pulse_sequence,
                get_cmd=self._get_pulse_sequence,
                instrument=self
            )
        
    def _set_pulse_sequence(self, pulse_sequence):
        if len(pulse_sequence)==0:
            self._pulse_sequence=[]
            return
    
        p=pulse_sequence
        self._flatten_pulse_sequence(p)
        self._pulse_sequence=p
        return
        
    def _get_pulse_sequence(self):
        return self._pulse_sequence
        
    
    def _flatten_pulse_sequence(self, pulse_sequence, seqtypes=(list, tuple)):
        try:
            for i, x in enumerate(pulse_sequence):
                while isinstance(pulse_sequence[i], seqtypes):    
                    pulse_sequence[i:i+1]=pulse_sequence[i]
        except IndexError:
            pass
        return pulse_sequence
    
    
class pulse():
    
    def __init__(self, level : int, duration : (float, Parameter)):
        if level not in [0,1]:
            raise ValueError("pulse: 'level' must be 0 or 1")
        self.level=level
        self.duration=duration
        return
    
    
    
    
    
    
    
    
    
    
    
    
    