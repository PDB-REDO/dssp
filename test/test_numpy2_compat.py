"""
DSSP numpy 2.x compatibility test script
This script tests whether the DSSP module works correctly with numpy 2.x
"""

import sys
import os
import gzip
import traceback

def find_test_file():
    """Find the test CIF file in the repository"""
    # Possible paths relative to script location
    possible_paths = [
        "test/1cbs.cif.gz",  # From repository root
        "../test/1cbs.cif.gz",  # From python-module directory
        "1cbs.cif.gz",  # From test directory
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    return None

def test_numpy_version():
    """Check numpy version"""
    try:
        import numpy as np
        print(f"✓ numpy version: {np.__version__}")
        major_version = int(np.__version__.split('.')[0])
        if major_version >= 2:
            print("✓ numpy 2.x detected")
        else:
            print(f"⚠ numpy {np.__version__} detected (not 2.x)")
        return True
    except ImportError as e:
        print(f"✗ Failed to import numpy: {e}")
        return False

def test_mkdssp_import():
    """Test mkdssp module import"""
    try:
        import mkdssp
        print("✓ Successfully imported mkdssp module")
        return True
    except ImportError as e:
        print(f"✗ Failed to import mkdssp: {e}")
        traceback.print_exc()
        return False

def test_basic_functionality():
    """Test basic functionality"""
    try:
        from mkdssp import dssp, TestBond, helix_type, chain_break_type, helix_position_type, structure_type
        print("✓ Successfully imported all mkdssp classes/functions")
        
        # Find test file
        test_file = find_test_file()
        if test_file is None:
            print("✗ Could not find test file (1cbs.cif.gz)")
            print("  Please run this script from the repository root or test directory")
            return False
        
        print(f"✓ Found test file: {test_file}")
        
        # Read test file
        try:
            with gzip.open(test_file, "rt") as f:
                file_content = f.read()
            print("✓ Successfully read test file")
        except Exception as e:
            print(f"✗ Failed to read test file: {e}")
            return False
        
        # Test DSSP object creation
        try:
            dssp_obj = dssp(file_content)
            print("✓ Successfully created DSSP object")
        except Exception as e:
            print(f"✗ Failed to create DSSP object: {e}")
            traceback.print_exc()
            return False
        
        # Test statistics
        try:
            stats = dssp_obj.statistics
            print(f"✓ Successfully retrieved statistics")
            print(f"  - Residues: {stats.residues}")
            print(f"  - Chains: {stats.chains}")
            print(f"  - H-bonds: {stats.H_bonds}")
        except Exception as e:
            print(f"✗ Failed to retrieve statistics: {e}")
            traceback.print_exc()
            return False
        
        # Test iteration
        try:
            count = 0
            for res in dssp_obj:
                count += 1
                if count == 1:
                    # Check first residue
                    print(f"  - First residue: {res.asym_id} {res.seq_id} {res.compound_id}")
            print(f"✓ Successfully iterated through {count} residues")
        except Exception as e:
            print(f"✗ Failed to iterate through residues: {e}")
            traceback.print_exc()
            return False
        
        # Test TestBond function
        try:
            a = dssp_obj.get('A', 137)
            b = dssp_obj.get('A', 6)
            result = TestBond(a, b)
            print(f"✓ Successfully tested TestBond function: {result}")
        except Exception as e:
            print(f"✗ Failed to test TestBond function: {e}")
            traceback.print_exc()
            return False
        
        return True
        
    except Exception as e:
        print(f"✗ Basic functionality test failed: {e}")
        traceback.print_exc()
        return False

def main():
    """Main test function"""
    print("=" * 60)
    print("DSSP numpy 2.x Compatibility Test")
    print("=" * 60)
    print()
    
    all_passed = True
    
    # Test 1: Check numpy version
    print("Test 1: Check numpy version")
    if not test_numpy_version():
        all_passed = False
    print()
    
    # Test 2: Import mkdssp
    print("Test 2: Import mkdssp module")
    if not test_mkdssp_import():
        all_passed = False
        print("\n⚠ Cannot import mkdssp, skipping remaining tests")
        print_summary(False)
        return 1
    print()
    
    # Test 3: Basic functionality
    print("Test 3: Test basic functionality")
    if not test_basic_functionality():
        all_passed = False
    print()
    
    # Summary
    print_summary(all_passed)
    
    return 0 if all_passed else 1

def print_summary(all_passed):
    """Print test summary"""
    print("=" * 60)
    if all_passed:
        print("✓ All tests passed!")
        print("✓ DSSP appears to be compatible with numpy 2.x")
    else:
        print("✗ Some tests failed")
        print("✗ There may be compatibility issues with numpy 2.x")
    print("=" * 60)

if __name__ == "__main__":
    sys.exit(main())
