"""Scaled convex disk QP with certified projected-KKT Newton corrections.

Force splits may be nonunique. No physical regularization is added. Explicit
rank-revealing least squares acts ONLY on the Newton correction, not spatial M_i.
"""
import numpy as np


class ConeCertificateError(RuntimeError):
    pass


def project_disks(y,radii):
    out=y.copy();D=np.zeros((y.size,y.size))
    for k,(v,r) in enumerate(zip(y,radii)):
        length=np.linalg.norm(v);block=slice(2*k,2*k+2)
        if length<=r and r>0:
            D[block,block]=np.eye(2)
        elif length>0 and r>0:
            direction=v/length;out[k]=r*direction
            D[block,block]=(r/length)*(np.eye(2)-np.outer(direction,direction))
        else:out[k]=0.
    return out,D


def solve(base,bases,mobility,coefficient,limits,tolerance=2e-12):
    b=np.asarray(base,float);T=np.asarray(bases,float)
    M=np.asarray(mobility,float);caps=np.asarray(limits,float);count=len(caps)
    if b.shape!=(count,2) or T.shape!=(count,3,2) or M.shape!=(3,3):
        raise ValueError('cone shape mismatch')
    if not all(np.isfinite(a).all() for a in (b,T,M,caps)) or not np.isfinite(coefficient):
        raise ValueError('nonfinite cone input')
    if coefficient<0 or np.any(caps<0) or tolerance<=0:
        raise ValueError('invalid cone coefficient/capacity/tolerance')
    if count==0:return np.empty((0,3))
    ms=float(np.linalg.norm(M,2))
    if np.linalg.norm(M-M.T)>1e-12*ms or np.linalg.eigvalsh(M).min() < -1e-12*ms:
        raise ValueError('symmetric PSD mobility required')
    if np.linalg.norm(T.swapaxes(1,2)@T-np.eye(2))>1e-12*count:
        raise ValueError('orthonormal tangent columns required')
    fs=float(np.max(caps))
    if fs==0:return np.zeros((count,3))
    vs=max(float(np.linalg.norm(b)),coefficient*ms*fs)
    if vs==0:return np.zeros((count,3))
    A=T.transpose(1,0,2).reshape(3,2*count)
    H=(coefficient*fs/vs)*(A.T@M@A)
    q=b/vs;r=caps/fs
    if coefficient*ms==0:
        f=np.zeros_like(q)
        for k,u in enumerate(q):
            speed=np.linalg.norm(u)
            if speed:f[k]=-r[k]*u/speed
        return fs*np.einsum('nki,ni->nk',T,f)
    # A common step loses the slow face in z-Proj(z-step*g) when a different
    # face sets the velocity scale. Positive scalar steps PER DISK have exactly
    # the same KKT zeros: no change to H, capacities, mobility or friction law.
    # Keep a GLOBAL curvature/physical-mobility floor, including all faces.
    # A tangent projection may cancel to zero in a constrained direction; its
    # roundoff must not set a huge step. These are numerical scales only: H
    # is unchanged. All x components use the same force unit fs (not caps[k]).
    curvature=max(float(np.linalg.norm(H,2)),coefficient*ms*fs/vs)
    rates=np.maximum(curvature,np.linalg.norm(q,axis=1))
    step=(1./rates)[:,None]
    row_step=np.repeat(step[:,0],2)
    target=min(tolerance,2e-14)
    x=np.zeros_like(q);last=[]

    def evaluate(z,with_jac=False):
        g=q+(H@z.ravel()).reshape(count,2)
        p,D=project_disks(z-step*g,r);F=z-p
        gap=float(np.sum(g*z)+np.sum(r*np.linalg.norm(g,axis=1)))
        error=max(float(np.linalg.norm(F)),abs(gap))
        J=np.eye(z.size)-D@(np.eye(z.size)-row_step[:,None]*H) if with_jac else None
        return error,F,g,gap,J

    for iteration in range(200):
        error,F,g,gap,J=evaluate(x,True)
        if error<=target:return fs*np.einsum('nki,ni->nk',T,x)
        # Explicit minimum-norm Newton correction, not a preferred force law.
        delta,_,rank,_=np.linalg.lstsq(J,-F.ravel(),rcond=64*np.finfo(float).eps)
        direction=delta.reshape(count,2);merit=float(np.sum(F*F));accepted=False
        for backtrack in range(32):
            alpha=2.**(-backtrack);trial,_=project_disks(x+alpha*direction,r)
            next_error,next_F,_,_,_=evaluate(trial)
            if next_error<=target or np.sum(next_F*next_F)<=(1.-1e-4*alpha)*merit:
                x=trial;accepted=True;break
        if not accepted:x,_=project_disks(x-step*g,r)
        if iteration>=196:
            last.append(dict(iteration=iteration+1,residual=error,gap=gap,
                             newton_rank=int(rank),line_search_accepted=accepted))
    error,F,g,gap,_=evaluate(x)
    failure=ConeCertificateError('coupled cone KKT certificate not attained')
    failure.context=dict(base_tangent=b,tangent_bases=T,mobility=M,
        coefficient=coefficient,limits_N=caps,scaled_residual=error,
        signed_gap=gap,target=target,last_updates=last,force_tangent_N=fs*x,
        disk_projected_steps=step[:,0])
    raise failure
