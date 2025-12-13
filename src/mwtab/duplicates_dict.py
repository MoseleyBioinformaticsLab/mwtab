# -*- coding: utf-8 -*-
"""
mwtab.duplicates_dict
~~~~~~~~~~~

This module provides the :class:`~mwtab.duplicates_dict.DuplicatesDict` class
that is a dictionary that can handle duplicate keys.
"""

import re
import copy

DUPLICATE_KEY_REGEX = r'(.*)\{\{\{.*\}\}\}'

class DuplicatesDict(dict):
    """
    This needs to inherit from a dict type class in order for it to be printed out 
    correctly by the json module. i.e. json.dumps().
    
    Note:
        The C encoder implementation in the json module checks whether a dictionary 
        is empty, and if so, short cuts to just printing {} instead of iterating 
        over items. A dummy value must be set so that we don't trip that short cut. 
        This is set in the init under the key 'dummy' if set_dummy is True. It is 
        possible this could cause other issues trying to use this class with other 
        modules, if so, either set set_dummy to False, or delete the 'dummy' key.
    """
    def __init__(self, obj: dict|None = None, set_dummy: bool = True):
        obj = {} if obj is None else obj
        if isinstance(obj, DuplicatesDict):
            self.data = obj.data
            self._key_count = obj._key_count
        else:
            self.data = obj
            self._key_count = {key:0 for key in obj.keys()}
        
        self._set_dummy = set_dummy
        if set_dummy:
            self.setdefault('dummy', None)
    
    def __len__(self):
        return len(self.data)
    
    def __bool__(self):
        return len(self.data) > 0
    
    def __iter__(self):
        return self.data.__iter__()

    def __getitem__(self, key):
        return self.data[key]
        
    def __setitem__(self, key, value):
        if key in self.data:
            self._key_count[key] += 1
            self.data[key + '{{{_' + str(self._key_count[key]) + '_}}}'] = value
        else:
            self._key_count[key] = 0
            self.data[key] = value    
    
    def __contains__(self, key):
        return key in self.data
    
    def __eq__(self, compare):
        return self.data == compare
    
    def __ne__(self, compare):
        return self.data != compare
    
    def __repr__(self):
        return self.data.__repr__().replace('OrderedDict', 'DuplicatesDict')
    
    def __delitem__(self, key):
        del self.data[key]
        
    def keys(self):
        return [key if not key.endswith('}}}') else re.match(DUPLICATE_KEY_REGEX, key).group(1) for key in self.data]
    
    def raw_keys(self):
        return self.data.keys()
    
    def items(self):
        items = []
        for key, value in self.data.items():
            key = key if not key.endswith('}}}') else re.match(DUPLICATE_KEY_REGEX, key).group(1)
            items.append((key, value))
        return items
    
    def raw_items(self):
        return self.data.items()
    
    def values(self):
        return self.data.values()
    
    def update(self, iterable):
        self.data.update(iterable)
        
    def __deepcopy__(self, memo):
        new_dict = DuplicatesDict({}, self._set_dummy)
        memo[id(new_dict)] = new_dict
        for key, value in self.raw_items():
            new_dict[key] = copy.deepcopy(value, memo)
        return new_dict
    
    def __copy__(self):
        new_dict = DuplicatesDict({}, self._set_dummy)
        for key, value in self.raw_items():
            new_dict[key] = copy.copy(value)
        return new_dict



