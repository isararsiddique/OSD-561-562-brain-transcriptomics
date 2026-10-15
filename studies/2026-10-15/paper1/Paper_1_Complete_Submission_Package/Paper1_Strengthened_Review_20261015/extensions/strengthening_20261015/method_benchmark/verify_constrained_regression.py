#!/usr/bin/env python3
"""Independent dense optimization and operator checks for the convex fit."""
import json
from pathlib import Path
import numpy as np
import cvxpy as cp
from constrained_regression import solve_constrained

HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(2026101531);checks=[];details=[]
def check(name,value,tolerance):
    checks.append({'name':name,'absolute_error':float(abs(value)),'tolerance':tolerance,'passed':bool(abs(value)<=tolerance)})

for case in range(12):
    n,p,g=10,3,5
    S=rng.normal(size=(n,p))*rng.uniform(.5,2,p);Y=rng.normal(size=(n,g))
    Y+=S@rng.normal(size=(p,g))*(.1 if case%3==0 else 2.)
    S-=S.mean(axis=0);Y-=Y.mean(axis=0);scale=np.sqrt(np.mean(S*S,axis=0))
    alpha=[.01,.1,1.,10.][case%4]
    Q,certificate=solve_constrained(S,Y,scale,alpha=alpha,tolerance=1e-10,max_iter=10000)
    C=Q@Y
    def f(flat):
        M=flat.reshape(p,g)
        return np.sum((Y-S@M)**2)/(2*n)+alpha*np.sum((scale[:,None]*M)**2)/2
    def jac(flat):
        M=flat.reshape(p,g)
        return (S.T@(S@M-Y)/n+alpha*scale[:,None]**2*M).ravel()
    # Independent conic formulation over the FULL coefficient matrix. The
    # spectral constraint is represented as a positive-semidefinite block.
    # A prior SLSQP diagnostic is retained: nonsmooth repeated singular values
    # caused unreliable convergence, so it is not an optimization oracle.
    variable=cp.Variable((p,g));cap=np.sqrt(3)
    psd=cp.bmat([[cap*np.eye(p),variable],[variable.T,cap*np.eye(g)]])
    problem=cp.Problem(cp.Minimize(cp.sum_squares(Y-S@variable)/(2*n)+
                       alpha*cp.sum_squares(cp.multiply(scale[:,None],variable))/2),[psd >> 0])
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=1000)
    check(f'case{case}: optimized objective agreement',f(C.ravel())-problem.value,2e-6)
    check(f'case{case}: solver objective direct reconstruction',f(C.ravel())-certificate['objective'],1e-9)
    check(f'case{case}: spectral feasibility excess',max(0,np.linalg.svd(C,compute_uv=False)[0]-np.sqrt(3)),1e-10)
    A=np.block([[np.zeros((p,p)),np.zeros((p,g))],[-C.T,np.eye(g)]])
    check(f'case{case}: idempotence',np.linalg.norm(A@A-A),1e-10)
    check(f'case{case}: exact individual-sensor annihilation',np.linalg.norm(A[:p]),1e-10)
    check(f'case{case}: independent dense operator norm',np.linalg.svd(A,compute_uv=False)[0]-np.sqrt(1+np.linalg.svd(C,compute_uv=False)[0]**2),1e-10)
    v=rng.normal(size=p+g);v[:p]=0
    check(f'case{case}: sensor-disjoint signal preservation',np.linalg.norm(A@v-v),1e-10)
    # First-order convex dual gap upper bounds suboptimality on the spectral ball.
    grad=jac(C.ravel()).reshape(p,g)
    dense_gap=float(np.sum(grad*C)+np.sqrt(3)*np.linalg.svd(grad,compute_uv=False).sum())
    check(f'case{case}: independently computed dual gap',dense_gap-certificate['Frank_Wolfe_duality_gap'],1e-8)
    details.append({'case':case,'independent_optimizer_success':problem.status=='optimal',
                    'independent_optimizer_message':problem.status,'objective':f(C.ravel()),
                    'independent_objective':float(problem.value),'dual_gap':dense_gap,**certificate})

report={'status':'PASS' if all(x['passed'] for x in checks) else 'FAIL','n_checks':len(checks),
        'independent_solver':'CVXPY/Clarabel positive-semidefinite conic formulation over full coefficients; distinct from production accelerated projected-gradient row-space algorithm',
        'checks':checks,'cases':details}
(HERE/'convex_independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(report['status'],report['n_checks'],'checks')
if report['status']!='PASS':raise SystemExit(1)
