#!/usr/bin/env python3
"""Compute distance of an assembled surface force to the two-phase Q1 pressure range."""
import argparse, json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr

p = argparse.ArgumentParser()
p.add_argument('prefix', help='output prefix from cutfem_surface_range')
p.add_argument('--maxiter', type=int, default=6000)
p.add_argument('--max-relative', type=float, help='fail if the range distance exceeds this value')
a = p.parse_args()
n = len(np.loadtxt(a.prefix+'_nodes.csv', delimiter=',', skiprows=1))
t = np.loadtxt(a.prefix+'_bt.csv', delimiter=',', skiprows=1)
B = coo_matrix((t[:,2], (t[:,0].astype(int), t[:,1].astype(int))), shape=(3*n,2*n)).tocsr()
B.sum_duplicates(); B.eliminate_zeros()
d = np.genfromtxt(a.prefix+'_forces.csv', delimiter=',', names=True)
idx = 3*d['node'].astype(int)+d['component'].astype(int)
f = np.zeros(3*n); f[idx] = d['surface']
pressure = np.zeros(3*n); pressure[idx] = d['pressure']
colnorm = np.sqrt(np.asarray(B.power(2).sum(axis=0)).ravel())
active = np.flatnonzero(colnorm > 0)
A = B[:,active].multiply(1/colnorm[active]).tocsr()
B128 = B.astype(np.longdouble); f128 = f.astype(np.longdouble)
def norm(x): return float(np.sqrt(np.sum(np.asarray(x,dtype=np.longdouble)**2)))
pstar = np.zeros(2*n); pstar[:n] = 10
r0 = B128 @ pstar.astype(np.longdouble) + f128
mismatch = norm(B128 @ pstar.astype(np.longdouble) - pressure.astype(np.longdouble))
gauge = norm(B128 @ np.ones(2*n,dtype=np.longdouble))
if mismatch > 1e-12*norm(pressure) or gauge > 1e-13*norm(pressure):
    raise RuntimeError('Pressure matrix consistency or gauge check failed')
if norm(f) == 0: raise ValueError('Surface force has zero norm; relative distance is undefined')
print(json.dumps({'stage':'assembly','free_dofs':len(idx),'active_pressure_dofs':len(active),'surface_norm':norm(f),'pressure_mismatch':mismatch,'gauge_norm':gauge,'prescribed_relative':norm(r0)/norm(f)}),flush=True)
result = lsmr(A,-f,atol=0,btol=0,conlim=1e15,maxiter=a.maxiter)
popt = np.zeros(2*n); popt[active] = result[0]/colnorm[active]
r = B128 @ popt.astype(np.longdouble)+f128
print(json.dumps({'stage':'range','absolute':norm(r),'relative':norm(r)/norm(f),'iterations':result[2],'stop':result[1],'normal_residual':result[4]}),flush=True)
if result[1] in (3, 6, 7): raise RuntimeError('Least-squares iteration limit reached')
if a.max_relative is not None and norm(r)/norm(f) > a.max_relative: raise RuntimeError('Range residual exceeded requested threshold')
random_force = np.zeros(3*n); random_force[idx] = np.random.default_rng(12345).standard_normal(len(idx))
negative = lsmr(A,-random_force,atol=0,btol=0,conlim=1e15,maxiter=a.maxiter)
negative_relative = norm(random_force+A@negative[0])/norm(random_force)
print(json.dumps({'stage':'negative_control','relative':negative_relative,'iterations':negative[2]}),flush=True)
if negative[1] in (3, 6, 7) or negative_relative < 0.1: raise RuntimeError('Negative control failed')
