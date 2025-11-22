#importing module numpy as we need matrix operations
import numpy as np
from scipy.io import mmread
#defining function to solve linear system Ax=b using manual SVD
def solve_system_manual_svd(A,b,tol=1e-10):
    #Converting A and b to numpy arrays
    A=np.array(A,dtype=float)
    b=np.array(b,dtype=float)
    #getting dimensions of A
    m,n=A.shape
    print("\n==============================")
    #printing A and b
    print("A=\n",A)
    print("b=\n",b)
    #1)Compute A^T.A to caluclate V matrix 
    AtA=A.T@A
    print("\nA^T.A=\n",AtA)
    #2)Eigen-decomposition of ATA => V
    #eigenvalues and eigenvectors
    evalsA,evecsA=np.linalg.eig(AtA)
    #sorting the eigenvalues and eigenvectors in descending order
    idx=np.argsort(evalsA)[::-1]
    evalsA=evalsA[idx]
    evecsA=evecsA[:,idx]
    #we have V from eigenvectors of A^T A
    V=evecsA
    print("\nV from eigenvectors of A^T.A:\n",V)
    print("Eigenvalues λ_i=\n",evalsA)
    #3)Singular values
    sig=np.sqrt(np.maximum(evalsA,0))
    print("\nSingular values σ_i=\n",sig)
    #Constructing Sigma matrix
    Sigma=np.zeros((m,n))
    #Filling diagonal of Sigma with singular values
    for i in range(min(m,n)):
        Sigma[i,i]=sig[i]
    #printing Sigma
    print("\nΣ=\n",Sigma)
    #Sigma pseudoinverse
    Sigma_pinv=np.zeros((n,m))
    #Filling diagonal of Sigma_pinv with reciprocals of non-zero singular values
    for i in range(min(m,n)):
        if sig[i]>tol:
            Sigma_pinv[i,i]=1/sig[i]
    #printing Sigma_pinv
    print("\nΣ⁺=\n",Sigma_pinv)
    #4)Compute U using AA^T eigenvectors
    AAt=A@A.T
    print("\nA.A^T=\n",AAt)
    #eigen-decomposition of AAt => U
    evalsU,evecsU=np.linalg.eig(AAt)
    #sorting the eigenvalues and eigenvectors in descending order
    idx2=np.argsort(evalsU)[::-1]
    evalsU=evalsU[idx2]
    evecsU=evecsU[:,idx2]
    #we have U from eigenvectors of A A^T
    U=evecsU
    print("\nU from eigenvectors of A.A^T:\n",U)
    print("Eigenvalues λ_i=\n",evalsU)
    #A⁺=V Σ⁺ U^T
    A_pinv=V@Sigma_pinv@U.T
    print("\nA⁺=V.Σ⁺.U^T=\n",A_pinv)
    #x=A⁺b
    x=A_pinv@b
    print("\nx=A⁺.b=\n",x)
    #7)Residual
    r=b-A@x
    normr=np.linalg.norm(r)
    print("\nResidual (r=b-Ax)=\n",r)
    print("||r||₂=",normr)
    #Checking consistency
    if normr<1e-8:
        print("\nSYSTEM CONSISTENT")
    else:
        print("\nSYSTEM INCONSISTENT")
    print("==============================\n")
    return x,normr
#for loading UF sparse matrix files
def load_uf_matrix(A_path, b_path):
    print("\nLoading UF Sparse Matrix Files...")
    print("A file:", A_path)
    print("b file:", b_path)
    # Reading the matrix files using mmread
    A = mmread(A_path)
    b = mmread(b_path)
    # Convert to dense only if matrix is sparse
    if hasattr(A, "toarray"):
        A = A.toarray()
    if hasattr(b, "toarray"):
        b = b.toarray()
    b = np.array(b).flatten()
    #print shapes
    print("Loaded UF matrix shapes -> A:", A.shape, ", b:", b.shape)
    return A, b
# ================================
# LOAD CIRCUIT_1 FROM UF DATABASE
# ================================
UF_A_PATH = r"C:\All Projects\Matrix Computation (5th Sem)\circuit_1_x.mtx"
UF_B_PATH = r"C:\All Projects\Matrix Computation (5th Sem)\circuit_1_b.mtx"
A_uf, b_uf = load_uf_matrix(UF_A_PATH, UF_B_PATH)
print("\nSOLVING UF CIRCUIT_1 USING MANUAL SVD...\n")
solve_system_manual_svd(A_uf, b_uf)
# m = n CONSISTENT
A1 = [
    [5,-1,0,0,0],
    [-1,4,-1,0,0],
    [0,-1,3,-1,0],
    [0,0,-1,2,-1],
    [0,0,0,-1,1]
]
b1 = [10,5,0,0,0]

print("\nCONSISTENT (m=n)")
solve_system_manual_svd(A1, b1)
# m = n INCONSISTENT
A2 = [
    [1,2,3],
    [2,4,6],
    [1,1,1]
]
b2 = [6,12,5]  # inconsistent because eq1 & eq2 proportional but b not proportional

print("\nINCONSISTENT (m=n)")
solve_system_manual_svd(A2, b2)
# m > n CONSISTENT (overdetermined but consistent)
A3 = [
    [1,2],
    [2,4],
    [3,6]
]
b3 = [5,10,15]  # perfectly proportional → consistent

print("\nCONSISTENT (m>n)")
solve_system_manual_svd(A3, b3)
# m > n INCONSISTENT
A4 = [
    [1,2],
    [2,4],
    [3,6]
]
b4 = [5,10,20]   # breaks proportionality → inconsistent

print("\nINCONSISTENT (m>n)")
solve_system_manual_svd(A4, b4)
# m < n CONSISTENT (infinite solutions)
A5 = [
    [1,2,0,0],
    [0,1,1,1]
]
b5 = [5,3]  # always solvable → consistent

print("\nCONSISTENT (m<n)")
solve_system_manual_svd(A5, b5)
# m < n INCONSISTENT
A6 = [
    [1,2,0,0],
    [0,1,1,1]
]
b6 = [5,300]   # second equation cannot satisfy → inconsistent

print("\nINCONSISTENT (m<n)")
solve_system_manual_svd(A6, b6)
# 50x50 CONSISTENT MATRIX
np.random.seed(0)
R = np.random.randint(1,5,(50,50))
A7 = (R + R.T) + 50*np.eye(50)     # SPD → invertible → consistent
b7 = np.random.randint(1,10,50)

print("\nCONSISTENT (50x50)")
solve_system_manual_svd(A7, b7)
# 50x50 INCONSISTENT MATRIX
np.random.seed(1)
A8 = np.random.randint(1, 10, (50, 50))
A8[0] = A8[1].copy()
b8 = np.random.randint(1, 10, 50)
b8[0] = 1
b8[1] = 9     
print("\nINCONSISTENT (50x50)")
solve_system_manual_svd(A8, b8)

