class Error_Handler:
    
    """
    Handles how errors are printed in Jupyter notebooks. By default, errors that
    are handled by AqEquil are printed with an error message, but no traceback.
    Errors that are not handled by AqEquil, such as those thrown if the user
    encounters a bug, will display a full traceback.
    
    If the error handler prints an error message without traceback, all future
    errors regardless of origin will be shown without traceback until the
    notebook kernel is restarted.
    
    Parameters
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
        
        if self.clean and _isnotebook():
            ipython = get_ipython()
            ipython.showtraceback = self.hide_traceback
        raise Exception(msg)