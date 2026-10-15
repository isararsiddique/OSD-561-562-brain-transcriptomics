"""Gain-constrained sensor regression with a convex optimality certificate.

For raw-unit sensor matrix S and outside response Y, minimize
  ||Y-SC||_F^2/(2n) + alpha ||D C||_F^2/2,
subject to ||C||_2 <= sqrt(L^2-1), where D is training sensor RMS scaling.
The operator [0,0;-C.T,I] has gain at most L and annihilates sensor coordinates.
This is a standard convex spectral constraint, not a new projection theorem.
"""
import numpy as np

def solve_constrained(S,Y,scale,alpha=.1,max_gain=2.,tolerance=1e-7,max_iter=2000):
    n,p=S.shape;cap=np.sqrt(max_gain**2-1)
    # Optimize in the response row space. This is exact because orthogonal
    # response directions only increase loss, ridge penalty and spectral norm.
    eig,U=np.linalg.eigh(Y@Y.T);keep=eig>max(float(eig.max()),1.)*1e-12
    U=U[:,keep];sing=np.sqrt(eig[keep]);E=U*sing
    H=S.T@S/n+alpha*np.diag(scale**2);F=S.T@E/n
    lipschitz=float(np.linalg.eigvalsh(H)[-1]);step=1/lipschitz
    C=np.zeros((p,len(sing)));Z=C.copy();acceleration=1.
    constant=float(np.sum(E*E)/(2*n));gap=np.inf
    def objective(A):return constant+.5*float(np.sum(A*(H@A)))-float(np.sum(A*F))
    def project(A):
        left,s,right=np.linalg.svd(A,full_matrices=False)
        return (left*np.minimum(s,cap))@right
    previous=objective(C)
    for iteration in range(1,max_iter+1):
        candidate=project(Z-step*(H@Z-F));value=objective(candidate)
        # Monotone restart prevents accelerated oscillation near the boundary.
        if value>previous+1e-12*max(1,constant):
            Z=C.copy();acceleration=1.;candidate=project(Z-step*(H@Z-F));value=objective(candidate)
        grad=H@candidate-F
        gap=float(np.sum(grad*candidate)+cap*np.linalg.svd(grad,compute_uv=False).sum())
        new_acceleration=(1+np.sqrt(1+4*acceleration**2))/2
        Z=candidate+(acceleration-1)/new_acceleration*(candidate-C)
        C=candidate;acceleration=new_acceleration;previous=value
        if gap<=tolerance*max(1,constant):break
    Q=(C/sing)@U.T
    result={'iterations':iteration,'objective':value,'objective_baseline':constant,
            'Frank_Wolfe_duality_gap':max(0.,gap),'relative_gap':max(0.,gap)/max(1,constant),
            'converged':gap<=tolerance*max(1,constant),
            'coefficient_spectral_norm':float(np.linalg.svd(C,compute_uv=False)[0]),
            'response_rank':len(sing),'alpha':alpha,'gain_cap':max_gain}
    return Q,result
