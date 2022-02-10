DEBUGGING_R = False
import sys
import pkg_resources
import re

# rpy2 for Python and R integration
import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)   # will display errors, but not warnings

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()




class Error_Handler:
    
    """
    Handles how errors are printed in Jupyter notebooks. By default, errors that
    are handled by AqEquil are printed with an error message, but no traceback.
    Errors that are not handled by AqEquil, such as those thrown if the user
    encounters a bug, will display a full traceback.
    
    If the error handler prints an error message without traceback, all future
    errors regardless of origin will be shown without traceback until the
    notebook kernel is restarted.
    
    Attributes
    ----------
    clean : bool
        Report exceptions without traceback? If True, only the error message is
        shown. If False, the entire error message, including traceback, is
        shown. Ignored if AqEquil is not being run in a Jupyter notebook.
    
    """
    def __init__(self, clean=True):
        self.clean = clean # bool: hide traceback?
        pass
    
    
    @staticmethod
    def hide_traceback(exc_tuple=None, filename=None, tb_offset=None,
                       exception_only=False, running_compiled_code=False):
        
        """
        Return a modified ipython showtraceback function that does not display
        traceback when encountering an error.
        """
        
        ipython = get_ipython()
        etype, value, tb = sys.exc_info()
        value.__cause__ = None  # suppress chained exceptions
        return ipython._showtraceback(etype, value, ipython.InteractiveTB.get_exception_only(etype, value))
        

    def raise_exception(self, msg):
        
        """
        Raise an exception that displays the error message without traceback. This
        happens only when the exception is predicted by the AqEquil package
        (e.g., for common user errors).
        """
        
        if self.clean and isnotebook():
            ipython = get_ipython()
            ipython.showtraceback = self.hide_traceback
        raise Exception(msg)


class Printer:
    def __init__(self, hide_traceback=True):
        self.hide_traceback = hide_traceback
        self.err_handler = Error_Handler(clean=self.hide_traceback)
        
    def __capture_r_output(self):
        """
        Capture and create a list of R console messages
        """
        
        # Record output #
        self.stdout = []
        self.stderr = []
        
        # print lines from R nicely if not debugging. Will ugly print otherwise.
        if not DEBUGGING_R:
        
            # Dummy functions #
            def add_to_stdout(line): self.stdout.append(line)
            def add_to_stderr(line): self.stderr.append(line)

            # Keep the old functions #
            self.stdout_orig = rpy2.rinterface_lib.callbacks.consolewrite_print
            self.stderr_orig = rpy2.rinterface_lib.callbacks.consolewrite_warnerror

            # Set the call backs #
            rpy2.rinterface_lib.callbacks.consolewrite_print     = add_to_stdout
            rpy2.rinterface_lib.callbacks.consolewrite_warnerror = add_to_stderr
    
    def __print_captured_r_output(self):
        printable_lines = [line for line in self.stdout if line not in ['[1]', '\n']]
        printable_lines = [line for line in printable_lines if re.search("^ \[[0-9]+\]$", line) is None]
        printable_lines = [re.sub(r" \\n", "", line) for line in printable_lines]
    
        [print(line[2:-1]) for line in self.stdout if line not in ['[1]', '\n'] and re.search("^ \[[0-9]+\]$", line) is None]
            
    def print_test_py(self):
        
        self.__capture_r_output()
        
        print_test_r_contents = pkg_resources.resource_string(
            __name__, 'print_test_r.r').decode("utf-8")

        ro.r(print_test_r_contents)

        ro.r.print_test_r()
        
        self.__print_captured_r_output()