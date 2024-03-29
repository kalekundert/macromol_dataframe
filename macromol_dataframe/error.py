class TidyError(Exception):
    # This class mimics `tidyexc`.  I don't want to use `tidyexc` directly, 
    # because it has weird incompatibilities with multiprocessing, and this 
    # code is meant to be used in pytorch data loaders.  Although I don't know 
    # the exact cause of the incompatibilities, I assume they have to do either 
    # with class-level state or with maintaining a data structure of arbitrary 
    # objects.  This code gets rid of those aspects of `tidyexc` and just keeps 
    # track of a few strings that will be used to make an error message.

    def __init__(self, brief=None, *, info=None, blame=None):
        self.brief = brief
        self.info = info or []
        self.blame = blame or []

    def __str__(self):
        info_strs = ['• ' + x for x in self.info]
        blame_strs = ['✖ ' + x for x in self.blame]
        msg_strs = [self.brief, *info_strs, *blame_strs]
        return '\n'.join(msg_strs)
