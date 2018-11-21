"""
Unit and regression test for the dynamic_qm_properties package.
"""
import sys
sys.path.append("/home/david/Dropbox/projects17/DynQMProp/code/dynamic_qm_properties/dynamic_qm_properties")

# Import package, test suite, and other packages as needed
import dynamic_qm_properties
import pytest

def test_dynamic_qm_properties_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "dynamic_qm_properties" in sys.modules

a=amber.run_antechamber()
