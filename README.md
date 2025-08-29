# Mapped Tent Pitching Scheme for the 1D Wave Equation

This repository is part of the **NAPDE project**, where the **Mapped Tent Pitching (MTP) scheme** has been implemented and applied to the **1D wave equation**.

---

## Project Description

The Mapped Tent Pitching (MTP) scheme is a **locally implicit space-time method**.  
The computational domain is decomposed into *tents*, which are space-time elements allowing local progression in time.  
Each tent is mapped into a reference domain where the problem is solved, enabling **stable time integration** even on highly varying meshes.  

The MTP approach provides:
- Local adaptivity in both space and time.  
- Stability properties similar to implicit schemes.  
- Computational efficiency close to explicit schemes.  

In this project, the MTP method has been implemented and tested on the **1D wave equation**.  
The method handles Dirichlet boundary conditions and allows convergence tests by varying the spatial discretization.

---

## Installation

To install and use the repository, navigate to the folder where you want to clone it and use **SSH**:

```bash
git clone git@github.com:GaetaEmanuele/MTP_project_NAPDE.git

## Code Organization

The repository is structured as follows:

- `src/` – Main source code:
  - `FEM_1D/` – Routines for spatial discretization using the Finite Element Method.  
  - `ODE/` – Time integration schemes and boundary condition handling.  
  - `post_processing/` – Tools for analyzing and visualizing results.  
  - `Tents/` – Construction and management of tents for the mapped scheme.  

- `tests/` – Test scripts to validate the implementation.

---

## Running Tests

Two main tests are included:

1. **`prova1.m`**  
   - Solves the **1D wave equation** on \( I \times [0,T] \) with **Dirichlet boundary conditions**.  
   - Demonstrates the full implementation of the MTP scheme.  

2. **`prova2.m`**  
   - Convergence test.  
   - Modify the spatial discretization parameter in `C_dati.m` by setting:  
     ```matlab
     h = 0.1;
     ```  
   - Run `prova2.m` to observe the convergence behavior of the scheme.

## Authors

- **Emanuele Gaeta** ([GitHub](https://github.com/GaetaEmanuele))
