# -*- coding: utf-8 -*-
"""
mwtab.duplicates_dict
~~~~~~~~~~~

This module provides the :class:`~mwtab.duplicates_dict.DuplicatesDict` class
that is a dictionary that can handle duplicate keys.
"""

from collections import OrderedDict

import json_duplicate_keys as jdks



class DuplicatesDict(jdks.JSON_DUPLICATE_KEYS):
    def __init__(self, Jobj=None):
        if Jobj is None:
            Jobj = OrderedDict()
        super(DuplicatesDict, self).__init__(Jobj)
    
    def __iter__(self):
        return self.__Jobj.__iter__()

    def __getitem__(self, key):
        return self.__Jobj[key]
    
    def __setitem__(self, key, value):
        self.set(key, value, ordered_dict=True)
    
    def keys(self):
        return [key for key in self.__Jobj.keys() if not key.endswith('}}}')]
    
    def raw_keys(self):
        return self.__Jobj.keys()

