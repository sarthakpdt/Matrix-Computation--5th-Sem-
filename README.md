# 🔌 Solving Electrical Circuits Using SVD (Singular Value Decomposition)

This project implements a **manual SVD-based solver** to compute **currents in electrical circuits** using the matrix equation:

\[
A.x = b
\]

Where:  
- **A** → Resistance / Conductance matrix  
- **x** → Unknown current vector  
- **b** → Voltage vector  

The solver can handle **consistent**, **inconsistent**, **overdetermined**, and **underdetermined** systems using **SVD and the Moore–Penrose pseudoinverse**.

---

## 🚀 Project Overview

Electrical circuits often produce large systems of equations using:

- Kirchhoff’s Current Law (KCL)  
- Kirchhoff’s Voltage Law (KVL)  
- Ohm’s Law (V = I × R)

Complex or irregular circuits generate matrices that may be:

- Singular  
- Rank-deficient  
- Overdetermined (more equations than unknowns)  
- Underdetermined (more unknowns than equations)

Traditional inverse-based methods fail in these cases.  
**SVD always works**, even when the matrix is not invertible.

---

## 🔍 Why SVD?

SVD factorizes any matrix into:

\[
A = UΣV^T
\]

From this, the **pseudoinverse** is computed as:

\[
A^+ = VΣ^+U^T
\]

This guarantees:

- **Exact solution** (if system is consistent)  
- **Least-squares solution** (if inconsistent)  
- **Minimum-norm solution** (if underdetermined)

This makes SVD the most stable and general technique to solve circuit equations.

---

## 📁 UF Sparse Matrix (University of Florida Collection)

The project also solves a real circuit matrix downloaded from the  
**UF Sparse Matrix Collection (Bomhof / circuit_1)**.

The uploaded files include:

- `circuit_1_x.mtx` → Matrix **A**
- `circuit_1_b.mtx` → Vector **b**

These files come from a real electrical circuit simulation based on  
Modified Nodal Analysis (MNA).

### 🔹 What does this matrix represent?

In simple words:

- It is a **large collection of equations** generated during circuit simulation.
- Each row represents one electrical constraint (KCL/KVL).
- The matrix relates **node voltages**, **currents**, and **resistances**.
- The vector **b** contains known inputs (voltage sources, injections).
- Solving **A x = b** gives the current that satisfies all circuit equations.

---

## 🧠 Code Features

✔ Manual construction of **U**, **Σ**, and **V**  
✔ Custom pseudoinverse calculation  
✔ Consistency check using residual norm  
✔ Supports all matrix shapes:
- Square matrices (m = n)
- Overdetermined (m > n)
- Underdetermined (m < n)
✔ Integration with UF Sparse Matrix files  
✔ Residual computation:  
\[
r = b - Ax
\]

---

## 📌 How the Code Works

### **1️⃣ Loading the Matrix**
The function `load_uf_matrix()` loads `.mtx` files and converts them to dense format.

### **2️⃣ Applying SVD**
The solver manually computes:

- \(A^T.A\)
- Eigenvalues → Singular values
- Eigenvectors → U and V matrices
- Pseudoinverse → \(A^+ = VΣ^+U^T\)

### **3️⃣ Solving the Circuit**
\[
x = A^+.b
\]
---

## 🧪 Test Matrices Included

The code tests:

- ✔ Consistent square matrices  
- ✔ Inconsistent square matrices  
- ✔ Overdetermined consistent systems  
- ✔ Overdetermined inconsistent systems  
- ✔ Underdetermined consistent systems  
- ✔ Underdetermined inconsistent systems  
- ✔ Large 50×50 random systems  
- ✔ UF Sparse Matrix (real circuit)

---

## 📄 Full Code (svd_solver.py)

> The full code is available in `svd_solver.py`  
> It includes:
> - SVD solver  
> - UF matrix loader  
> - Multiple test cases  
> - Residual analysis  

---

## 📈 Output Summary

The program prints:

- Matrix A and vector b  
- Singular values  
- U, Σ, V matrices  
- Pseudoinverse  
- Final solution vector  
- Residual  
- Consistency check  

It also displays results for:

- UF circuit matrix  
- Hardcoded example matrices  
- 50x50 random matrices  

---

## 📘 Verification

The report demonstrates both:

- Manual step-by-step calculations  
- Program output comparison  

It verifies that the SVD solution matches the expected theoretical result.

---

## 🛠 Requirements

Install SciPy if not installed:
pip install scipy

Or for specific Python version:
<path_to_python.exe> -m pip install scipy
## 🏁 Conclusion

This project demonstrates how **SVD is a powerful and reliable tool** for solving electrical circuits, especially when the system is inconsistent, singular, or poorly conditioned.

It successfully solves both:

- Real-world UF circuit matrices  
- Custom test matrices  

using a complete manual implementation of SVD and the pseudoinverse.

---
