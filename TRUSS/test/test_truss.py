#!/usr/bin/env python3
"""
Test script for TRUSS FEA program.
Creates test input, runs the program, and validates output.
Works with any language implementation that follows the same I/O format.
"""

import subprocess
import re
import os
import sys
from typing import Dict, List, Tuple


def create_test_input(filename: str = "test_input.txt"):
    """Create the test input file from Section 3.7."""
    content = """SPACE TRUSS EXAMPLE OF SECTION 3.7
3,4
1,0,1,0,72.0,0.,0.,0.,0.,-1000.0
2,1,1,1,0.0,36.0,0.,0.,0.,0.
3,1,1,1,0.0,36.0,72.0,0.,0.,0.
4,1,1,1,0.0,0.0,-48.0,0.,0.,0.
1,1,4,1.2E+6,0.187
2,1,2,1.2E+6,0.302
3,1,3,1.2E+6,0.729
"""
    with open(filename, 'w') as f:
        f.write(content)
    print(f"Created test input file: {filename}")


def run_program(executable: str, input_file: str) -> str:
    """Run the TRUSS program and capture output."""
    try:
        # Try different execution methods
        if executable.endswith('.py'):
            result = subprocess.run(['python3', executable, input_file],
                                  capture_output=True, text=True, timeout=10)
        elif executable.endswith('.exe') or os.access(executable, os.X_OK):
            result = subprocess.run([executable, input_file],
                                  capture_output=True, text=True, timeout=10)
        else:
            result = subprocess.run(['./' + executable, input_file],
                                  capture_output=True, text=True, timeout=10)
        
        if result.returncode != 0:
            print(f"Program exited with error code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return None
        
        return result.stdout
    except subprocess.TimeoutExpired:
        print("Program execution timed out")
        return None
    except Exception as e:
        print(f"Error running program: {e}")
        return None


def parse_displacements(output: str) -> Dict[int, Tuple[float, float, float]]:
    """Parse displacement values from output."""
    displacements = {}
    lines = output.split('\n')
    
    in_displacement_section = False
    for line in lines:
        if 'DISPLACEMENTS' in line and 'X' in line and 'Y' in line:
            in_displacement_section = True
            continue
        
        if in_displacement_section and 'NODE NUMBER' in line:
            # Parse line like: "NODE NUMBER 1   -0.7111E-01   0.0000E+00  -0.2662E+00"
            match = re.search(r'NODE NUMBER\s+(\d+)\s+([-+]?\d+\.\d+E[+-]\d+)\s+([-+]?\d+\.\d+E[+-]\d+)\s+([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                node = int(match.group(1))
                dx = float(match.group(2))
                dy = float(match.group(3))
                dz = float(match.group(4))
                displacements[node] = (dx, dy, dz)
        elif in_displacement_section and line.strip() == '':
            break
    
    return displacements


def parse_stresses(output: str) -> Dict[int, float]:
    """Parse stress values from output."""
    stresses = {}
    lines = output.split('\n')
    
    in_stress_section = False
    for line in lines:
        if 'STRESSES IN ELEMENTS' in line:
            in_stress_section = True
            continue
        
        if in_stress_section and '=' in line:
            # Parse line like: "          1 =    -0.28685E+04"
            match = re.search(r'(\d+)\s*=\s*([-+]?\d+\.\d+E[+-]\d+)', line)
            if match:
                element = int(match.group(1))
                stress = float(match.group(2))
                stresses[element] = stress
    
    return stresses


def parse_mud(output: str) -> int:
    """Parse MUD value from output."""
    match = re.search(r'NUMBER OF NONZERO UPPER CODIAGONALS\(MUD\)\s*=\s*(\d+)', output)
    if match:
        return int(match.group(1))
    return None


def validate_results(output: str, tolerance: float = 1e-3) -> bool:
    """Validate output against expected values."""
    print("\n" + "="*70)
    print("VALIDATION RESULTS")
    print("="*70)
    
    all_passed = True
    
    # Expected values from the documentation
    expected_displacements = {
        1: (-0.7111E-01, 0.0000E+00, -0.2662E+00),
        2: (0.0000E+00, 0.0000E+00, 0.0000E+00),
        3: (0.0000E+00, 0.0000E+00, 0.0000E+00),
        4: (0.0000E+00, 0.0000E+00, 0.0000E+00),
    }
    
    expected_stresses = {
        1: -0.28685E+04,
        2: -0.94819E+03,
        3: 0.14454E+04,
    }
    
    expected_mud = 11
    
    # Parse actual values
    actual_displacements = parse_displacements(output)
    actual_stresses = parse_stresses(output)
    actual_mud = parse_mud(output)
    
    # Validate MUD
    print("\nMUD (Maximum Upper Codiagonals):")
    if actual_mud == expected_mud:
        print(f"  ✓ PASS: MUD = {actual_mud}")
    else:
        print(f"  ✗ FAIL: Expected {expected_mud}, got {actual_mud}")
        all_passed = False
    
    # Validate displacements
    print("\nDisplacements:")
    for node in expected_displacements:
        if node not in actual_displacements:
            print(f"  ✗ FAIL: Node {node} displacement not found in output")
            all_passed = False
            continue
        
        exp_dx, exp_dy, exp_dz = expected_displacements[node]
        act_dx, act_dy, act_dz = actual_displacements[node]
        
        dx_diff = abs(exp_dx - act_dx)
        dy_diff = abs(exp_dy - act_dy)
        dz_diff = abs(exp_dz - act_dz)
        
        max_expected = max(abs(exp_dx), abs(exp_dy), abs(exp_dz), 1e-10)
        rel_tol = tolerance * max_expected
        
        if dx_diff <= rel_tol and dy_diff <= rel_tol and dz_diff <= rel_tol:
            print(f"  ✓ PASS: Node {node} displacements within tolerance")
        else:
            print(f"  ✗ FAIL: Node {node} displacements")
            print(f"    Expected: ({exp_dx:.4E}, {exp_dy:.4E}, {exp_dz:.4E})")
            print(f"    Actual:   ({act_dx:.4E}, {act_dy:.4E}, {act_dz:.4E})")
            print(f"    Diff:     ({dx_diff:.4E}, {dy_diff:.4E}, {dz_diff:.4E})")
            all_passed = False
    
    # Validate stresses
    print("\nStresses:")
    for element in expected_stresses:
        if element not in actual_stresses:
            print(f"  ✗ FAIL: Element {element} stress not found in output")
            all_passed = False
            continue
        
        expected = expected_stresses[element]
        actual = actual_stresses[element]
        diff = abs(expected - actual)
        rel_tol = tolerance * abs(expected)
        
        if diff <= rel_tol:
            print(f"  ✓ PASS: Element {element} stress within tolerance")
            print(f"    Expected: {expected:.5E}, Actual: {actual:.5E}")
        else:
            print(f"  ✗ FAIL: Element {element} stress")
            print(f"    Expected: {expected:.5E}")
            print(f"    Actual:   {actual:.5E}")
            print(f"    Diff:     {diff:.5E} (tolerance: {rel_tol:.5E})")
            all_passed = False
    
    print("\n" + "="*70)
    if all_passed:
        print("ALL TESTS PASSED ✓")
    else:
        print("SOME TESTS FAILED ✗")
    print("="*70 + "\n")
    
    return all_passed


def main():
    """Main test execution."""
    if len(sys.argv) < 2:
        print("Usage: test_truss.py <executable> [input_file] [tolerance]")
        print("Example: test_truss.py truss.py")
        print("Example: test_truss.py ./truss test_input.txt 0.01")
        sys.exit(1)
    
    executable = sys.argv[1]
    input_file = sys.argv[2] if len(sys.argv) > 2 else "test_input.txt"
    tolerance = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
    
    # Create test input if it doesn't exist
    if not os.path.exists(input_file):
        create_test_input(input_file)
    
    # Run the program
    print(f"\nRunning: {executable} {input_file}")
    print("-" * 70)
    output = run_program(executable, input_file)
    
    if output is None:
        print("Failed to run program or capture output")
        sys.exit(1)
    
    # Print the output
    print(output)
    print("-" * 70)
    
    # Validate results
    passed = validate_results(output, tolerance)
    
    # Exit with appropriate code
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
