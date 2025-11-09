#!/usr/bin/env python3
"""
TRUSS Finite Element Analysis Program
Solves 3D truss structures using the direct stiffness method.
"""

import sys
import numpy as np
from typing import Tuple, List


def read_input_file(filename: str) -> dict:
    """Read and parse the input file."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Line 1: Title
    title = lines[0].strip()
    
    # Line 2: NELE, NNODE
    nele, nnode = map(int, lines[1].strip().split(','))
    
    # Initialize arrays
    ifix = np.zeros((3, nnode + 1), dtype=int)  # 1-indexed
    xc = np.zeros(nnode + 1)
    yc = np.zeros(nnode + 1)
    zc = np.zeros(nnode + 1)
    force = np.zeros((3, nnode + 1))
    node_conn = np.zeros((2, nele + 1), dtype=int)
    E = np.zeros(nele + 1)
    A = np.zeros(nele + 1)
    
    # Lines 3 to 3+nnode-1: Node information
    for i in range(nnode):
        parts = lines[2 + i].strip().split(',')
        j = int(parts[0])
        ifix[0, j] = int(parts[1])
        ifix[1, j] = int(parts[2])
        ifix[2, j] = int(parts[3])
        xc[j] = float(parts[4])
        yc[j] = float(parts[5])
        zc[j] = float(parts[6])
        force[0, j] = float(parts[7])
        force[1, j] = float(parts[8])
        force[2, j] = float(parts[9])
    
    # Remaining lines: Element information
    for i in range(nele):
        parts = lines[2 + nnode + i].strip().split(',')
        k = int(parts[0])
        node_conn[0, k] = int(parts[1])
        node_conn[1, k] = int(parts[2])
        E[k] = float(parts[3])
        A[k] = float(parts[4])
    
    return {
        'title': title,
        'nele': nele,
        'nnode': nnode,
        'ifix': ifix,
        'xc': xc,
        'yc': yc,
        'zc': zc,
        'force': force,
        'node_conn': node_conn,
        'E': E,
        'A': A
    }


def compute_element_stiffness(x1: float, y1: float, z1: float,
                              x2: float, y2: float, z2: float,
                              E: float, A: float) -> Tuple[np.ndarray, float]:
    """Compute the 6x6 element stiffness matrix in global coordinates."""
    # Calculate element length
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    L = np.sqrt(dx**2 + dy**2 + dz**2)
    
    # Direction cosines
    cx = dx / L
    cy = dy / L
    cz = dz / L
    
    # Element stiffness in global coordinates
    k_elem = np.zeros((6, 6))
    factor = E * A / L
    
    # Create the direction cosine vector
    c = np.array([cx, cy, cz])
    
    # Build the 6x6 stiffness matrix
    for i in range(3):
        for j in range(3):
            k_elem[i, j] = factor * c[i] * c[j]
            k_elem[i, j+3] = -factor * c[i] * c[j]
            k_elem[i+3, j] = -factor * c[i] * c[j]
            k_elem[i+3, j+3] = factor * c[i] * c[j]
    
    return k_elem, L


def assemble_global_stiffness(data: dict) -> Tuple[np.ndarray, np.ndarray, int]:
    """Assemble the global stiffness matrix and force vector."""
    nnode = data['nnode']
    nele = data['nele']
    ndof = 3 * nnode
    
    # Initialize global stiffness matrix and force vector
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)
    
    # Assemble element stiffness matrices
    for k in range(1, nele + 1):
        node1 = data['node_conn'][0, k]
        node2 = data['node_conn'][1, k]
        
        x1, y1, z1 = data['xc'][node1], data['yc'][node1], data['zc'][node1]
        x2, y2, z2 = data['xc'][node2], data['yc'][node2], data['zc'][node2]
        
        k_elem, L = compute_element_stiffness(x1, y1, z1, x2, y2, z2,
                                              data['E'][k], data['A'][k])
        
        # Global DOF indices for this element
        dof1 = [(node1 - 1) * 3, (node1 - 1) * 3 + 1, (node1 - 1) * 3 + 2]
        dof2 = [(node2 - 1) * 3, (node2 - 1) * 3 + 1, (node2 - 1) * 3 + 2]
        dofs = dof1 + dof2
        
        # Add element stiffness to global stiffness
        for i in range(6):
            for j in range(6):
                K[dofs[i], dofs[j]] += k_elem[i, j]
    
    # Assemble force vector
    for j in range(1, nnode + 1):
        F[(j-1)*3] = data['force'][0, j]
        F[(j-1)*3 + 1] = data['force'][1, j]
        F[(j-1)*3 + 2] = data['force'][2, j]
    
    # Calculate MUD (maximum upper codiagonal)
    mud = 0
    for i in range(ndof):
        for j in range(i+1, ndof):
            if K[i, j] != 0:
                mud = max(mud, j - i)
    
    return K, F, mud


def apply_boundary_conditions(K: np.ndarray, F: np.ndarray, 
                             data: dict) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    """Apply fixed displacement boundary conditions."""
    nnode = data['nnode']
    fixed_dofs = []
    
    # Identify fixed DOFs
    for j in range(1, nnode + 1):
        for i in range(3):
            if data['ifix'][i, j] == 1:
                dof = (j - 1) * 3 + i
                fixed_dofs.append(dof)
    
    # Create reduced system by removing fixed DOFs
    free_dofs = [i for i in range(3 * nnode) if i not in fixed_dofs]
    
    K_reduced = K[np.ix_(free_dofs, free_dofs)]
    F_reduced = F[free_dofs]
    
    return K_reduced, F_reduced, free_dofs


def solve_displacements(K: np.ndarray, F: np.ndarray, 
                       free_dofs: List[int], nnode: int) -> np.ndarray:
    """Solve for nodal displacements."""
    # Solve reduced system
    d_reduced = np.linalg.solve(K, F)
    
    # Expand to full displacement vector
    d = np.zeros(3 * nnode)
    for i, dof in enumerate(free_dofs):
        d[dof] = d_reduced[i]
    
    return d


def compute_stresses(d: np.ndarray, data: dict) -> np.ndarray:
    """Compute element stresses."""
    nele = data['nele']
    stresses = np.zeros(nele + 1)
    
    for k in range(1, nele + 1):
        node1 = data['node_conn'][0, k]
        node2 = data['node_conn'][1, k]
        
        x1, y1, z1 = data['xc'][node1], data['yc'][node1], data['zc'][node1]
        x2, y2, z2 = data['xc'][node2], data['yc'][node2], data['zc'][node2]
        
        # Element length and direction cosines
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        L = np.sqrt(dx**2 + dy**2 + dz**2)
        
        cx = dx / L
        cy = dy / L
        cz = dz / L
        
        # Displacements at nodes
        d1 = d[(node1-1)*3:(node1-1)*3+3]
        d2 = d[(node2-1)*3:(node2-1)*3+3]
        
        # Axial strain
        delta_L = cx * (d2[0] - d1[0]) + cy * (d2[1] - d1[1]) + cz * (d2[2] - d1[2])
        strain = delta_L / L
        
        # Stress = E * strain
        stresses[k] = data['E'][k] * strain
    
    return stresses


def print_output(data: dict, d: np.ndarray, stresses: np.ndarray, mud: int):
    """Print the output in the specified format."""
    print(data['title'])
    print()
    print(f"NUMBER OF ELEMENTS(NELE) = {data['nele']}")
    print(f"NUMBER OF NODES(NNODE)   = {data['nnode']}")
    print()
    print("NODE POINTS")
    print("K     IFIX        XC(K)           YC(K)           ZC(K)")
    
    for j in range(1, data['nnode'] + 1):
        print(f"{j}    {data['ifix'][0,j]} {data['ifix'][1,j]} {data['ifix'][2,j]}    "
              f"{data['xc'][j]:12.6E}    {data['yc'][j]:12.6E}    {data['zc'][j]:12.6E}")
    
    print()
    print("               FORCE(1,K)      FORCE(2,K)      FORCE(3,K)")
    for j in range(1, data['nnode'] + 1):
        print(f"              {data['force'][0,j]:12.6E}    "
              f"{data['force'][1,j]:12.6E}    {data['force'][2,j]:12.6E}")
    
    print()
    print("   ELEMENTS")
    print("K      NODE(1,K)          E(K)               A(K)")
    for k in range(1, data['nele'] + 1):
        print(f"{k}        {data['node_conn'][0,k]}   {data['node_conn'][1,k]}         "
              f"{data['E'][k]:12.4E}        {data['A'][k]:12.4E}")
    
    print()
    print(f"NUMBER OF NONZERO UPPER CODIAGONALS(MUD) = {mud}")
    print()
    print(" DISPLACEMENTS       X            Y            Z")
    
    for j in range(1, data['nnode'] + 1):
        dx = d[(j-1)*3]
        dy = d[(j-1)*3 + 1]
        dz = d[(j-1)*3 + 2]
        print(f"NODE NUMBER {j}   {dx:12.4E}   {dy:12.4E}  {dz:12.4E}")
    
    print()
    print("STRESSES IN ELEMENTS (IN CURRENT UNITS)")
    print()
    print("ELEMENT NUMBER     STRESS")
    for k in range(1, data['nele'] + 1):
        print(f"          {k} =    {stresses[k]:12.5E}")


def main():
    """Main program execution."""
    if len(sys.argv) != 2:
        print("Usage: truss.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Read input data
    data = read_input_file(input_file)
    
    # Assemble global stiffness matrix
    K, F, mud = assemble_global_stiffness(data)
    
    # Apply boundary conditions
    K_reduced, F_reduced, free_dofs = apply_boundary_conditions(K, F, data)
    
    # Solve for displacements
    d = solve_displacements(K_reduced, F_reduced, free_dofs, data['nnode'])
    
    # Compute stresses
    stresses = compute_stresses(d, data)
    
    # Print output
    print_output(data, d, stresses, mud)


if __name__ == "__main__":
    main()
