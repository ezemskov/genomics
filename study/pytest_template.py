import pytest
import time

def setup():
    pass

def teardown():
    pass

def setup_module(module):
    pass
 
def teardown_module(module):
    pass
 
def setup_function(function):
    pass
 
def teardown_function(function):
    pass

def func_under_test(duration=0.000001):
    time.sleep(duration)
    return 123

def test_some_func(benchmark):
    result = benchmark(func_under_test)
    assert result == 123
  
class TestUM:
    def setup(self):
        pass
 
    def teardown(self):
        pass
 
    def setup_class(cls):
        pass
 
    def teardown_class(cls):
        pass
 
    def setup_method(self, method):
        pass
 
    def teardown_method(self, method):
        pass
 
    def test_some_method(self):
        pass
