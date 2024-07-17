# -*- coding: utf-8 -*-
"""
mwtab.duplicates_dict
~~~~~~~~~~~

This module provides the :class:`~mwtab.duplicates_dict.DuplicatesDict` class
that is a dictionary that can handle duplicate keys.
"""

from collections import OrderedDict
import re

import json_duplicate_keys as jdks


DUPLICATE_KEY_REGEX = r'(.*)\{\{\{.*\}\}\}'

class DuplicatesDict(jdks.JSON_DUPLICATE_KEYS, OrderedDict):
    def __init__(self, Jobj=None):
        if Jobj is None:
            Jobj = OrderedDict()
        super(DuplicatesDict, self).__init__(Jobj)
    
    def __iter__(self):
        return self._JSON_DUPLICATE_KEYS__Jobj.__iter__()

    def __getitem__(self, key):
        return self._JSON_DUPLICATE_KEYS__Jobj[key]
    
    def __setitem__(self, key, value):
        self.set(key, value, ordered_dict=True)
    
    def __contains__(self, key):
        return key in self._JSON_DUPLICATE_KEYS__Jobj
    
    def __eq__(self, compare):
        return self._JSON_DUPLICATE_KEYS__Jobj == compare
    
    def __repr__(self):
        return self._JSON_DUPLICATE_KEYS__Jobj.__repr__().replace('OrderedDict', 'DuplicatesDict')
    
    def keys(self):
        return [key if not key.endswith('}}}') else re.match(DUPLICATE_KEY_REGEX, key).group(1) for key in self._JSON_DUPLICATE_KEYS__Jobj.keys()]
    
    def raw_keys(self):
        return self._JSON_DUPLICATE_KEYS__Jobj.keys()
    
    def items(self):
        items = []
        for key, value in self._JSON_DUPLICATE_KEYS__Jobj.items():
            key = key if not key.endswith('}}}') else re.match(DUPLICATE_KEY_REGEX, key).group(1)
            items.append((key, value))
        return items
    
    def raw_items(self):
        return self._JSON_DUPLICATE_KEYS__Jobj.items()
    
    def values(self):
        return self._JSON_DUPLICATE_KEYS__Jobj.values()

