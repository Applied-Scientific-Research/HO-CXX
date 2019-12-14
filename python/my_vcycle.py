"""Generic AMG solver."""


from warnings import warn

import sys
import scipy as sp
import numpy as np

try:
    import cupy
    #import cupyx as cupy
    import cupyx.scipy.sparse as cupy_sp
except ImportError:
    print("cuPy isn't available")

def _spmv(A,x):
    return A * x
SpMV = _spmv

try:
    import sp_solver
except ImportError:
    print("local sp_solver isn't available")
else:
    SpMV = sp_solver.SpMV_viaMKL
    print("Using MKL SpMV")

class amg_solver:

    def __init__(self, levels):

        self.levels = levels

        for level in levels[:-1]:
            if not hasattr(level, 'R'):
                level.R = level.P.H

        self.coarse_PLU  = sp.linalg.lu_factor( self.levels[-1].A.todense() )
        self.coarse_Ainv = sp.linalg.inv( self.levels[-1].A.todense() )

    def __repr__(self):
        """Print basic statistics about the multigrid hierarchy."""
        output = 'multilevel_solver\n'
        output += 'Number of Levels:     %d\n' % len(self.levels)
        output += 'Operator Complexity: %6.3f\n' % self.operator_complexity()
        output += 'Grid Complexity:     %6.3f\n' % self.grid_complexity()

        total_nnz = sum([level.A.nnz for level in self.levels])

        output += '  level   unknowns     nonzeros\n'
        for n, level in enumerate(self.levels):
            A = level.A
            output += '   %2d   %10d   %10d [%5.2f%%]\n' %\
                (n, A.shape[1], A.nnz,
                 (100 * float(A.nnz) / float(total_nnz)))

        return output

    def cycle_complexity(self):

        nnz = [level.A.nnz for level in self.levels]

        def V(level):
            if len(self.levels) == 1:
                return nnz[0]
            elif level == len(self.levels) - 2:
                return 2 * nnz[level] + nnz[level + 1]
            else:
                return 2 * nnz[level] + V(level + 1)

        flops = V(0)

        return float(flops) / float(nnz[0])

    def operator_complexity(self):
        """Operator complexity of this multigrid hierarchy.

        Defined as:
            Number of nonzeros in the matrix on all levels /
            Number of nonzeros in the matrix on the finest level

        """
        return sum([level.A.nnz for level in self.levels]) /\
            float(self.levels[0].A.nnz)

    def grid_complexity(self):
        """Grid complexity of this multigrid hierarchy.

        Defined as:
            Number of unknowns on all levels /
            Number of unknowns on the finest level

        """
        return sum([level.A.shape[0] for level in self.levels]) /\
            float(self.levels[0].A.shape[0])

    #def psolve(self, b):
    #    """Lagacy solve interface."""
    #    return self.solve(b, maxiter=1)

    def aspreconditioner(self):
        from scipy.sparse.linalg import LinearOperator

        shape = self.levels[0].A.shape
        dtype = self.levels[0].A.dtype

        def matvec(b):
            return self.solve(b, maxiter=1)

        return LinearOperator(shape, matvec, dtype=dtype)

    def solve(self, b, maxiter=1):

        x = np.zeros_like(b)

        # Create uniform types for A, x and b
        # Clearly, this logic doesn't handle the case of real A and complex b
        from scipy.sparse.sputils import upcast
        from pyamg.util.utils import to_type

        tp = upcast(b.dtype, x.dtype, self.levels[0].A.dtype)
        [b, x] = to_type(tp, [b, x])
        b = np.ravel(b)
        x = np.ravel(x)

        A = self.levels[0].A

        for iter in range(maxiter):
            if len(self.levels) == 1:
                # hierarchy has only 1 level
                x = self.coarse_solve(b)
            else:
                self.__solve(0, x, b)

        return x

    def __solve(self, lvl, x, b):

        #print('starting: ', lvl)
        A = self.levels[lvl].A

        self.levels[lvl].presmoother(A, x, b)

        residual = b - A * x

        coarse_b = self.levels[lvl].R * residual
        coarse_x = np.zeros_like(coarse_b)

        if lvl == len(self.levels) - 2:
            coarse_x[:] = self.coarse_solve(coarse_b)
        else:
            self.__solve(lvl + 1, coarse_x, coarse_b)

        x += self.levels[lvl].P * coarse_x   # coarse grid correction

        self.levels[lvl].postsmoother(A, x, b)

        #print('finished: ', lvl)

    def coarse_solve(self, b):
        #print('coarse')
        #return sp.linalg.lu_solve(self.coarse_PLU, b)
        return self.coarse_Ainv.dot( b )

def get_exec_space(A):
    if 'cupy' in sys.modules:
        return cupy.get_array_module(A)
    else:
        return np
    
class ml_level:
    def __init__(self, A, D_inv, R = None, P = None, n_pre = 1, n_post = 1, omega = 1.0, on_device = False):
        self.A = A
        self.D_inv = D_inv
        self.R = R
        self.P = P
        self.n_pre = n_pre
        self.n_post = n_post
        self.omega = omega
        self.on_device = on_device

        self.exec_space = get_exec_space(A)
        #self.exec_space = np
        #if 'cupyx.scipy.sparse' in sys.modules:
        #    self.exec_space = cupy.get_array_module(A)

        #print("exec_space: ", self.exec_space)

class coarse_solver:
    def __init__(self, A_inv, on_device = False):
        self.A_inv = A_inv

        if on_device:
            self.A_inv = cupy.asarray(A_inv)

            _b = cupy.ones( A_inv.shape[0], dtype=A_inv.dtype )
            _x = cupy.matmul( self.A_inv, _b )

        #self.exec_space = cupy.get_array_module(self.A_inv)
        self.exec_space = get_exec_space(self.A_inv)

        #print(self.A_inv.shape)

    def __call__(self, b, x = None):
        if x is None:
            x = self.exec_space.empty_like(b)
        #print(self.A_inv.shape, x.shape, b.shape)
        x[:] = self.exec_space.matmul( self.A_inv, b )
        return x

class amg_solver_gpu:

    def __init__(self, pyamg_levels, use_device = True):

        self.levels = []

        from pyamg.util.utils import scale_rows, get_block_diag, get_diagonal
        from pyamg.util.linalg import approximate_spectral_radius

        for lvl,level in enumerate(pyamg_levels):
            print('host: ', lvl)
            A = level.A
            D_inv = get_diagonal(level.A, inv=True)
            R = None
            P = None
            if hasattr(level, 'R'):
                R = level.R
            if hasattr(level, 'P'):
                P = level.P

            omega = 1.

            if lvl < len(pyamg_levels)-1:
                rho_D_inv = None
                if hasattr(level.A, 'rho_D_inv'):
                    rho_D_inv = level.A.rho_D_inv
                    #print(lvl,'reused')
                else:
                    D_inv_A = scale_rows(level.A, D_inv, copy=True)
                    rho_D_inv = approximate_spectral_radius(D_inv_A)

                # limit rho to [4/3,2)
                omega = 1.0 / rho_D_inv
                omega = 1.0 / max(4./3., min(2., rho_D_inv))
                #omega = (4./3.)/2.5

                print(lvl, omega, rho_D_inv)

            n_post = 2
            n_pre = 1
            new_level = ml_level( A, D_inv, R=R, P=P, omega=omega, on_device=False, n_post=2, n_pre=1 )
            self.levels.append( new_level )

        self.coarse_solver = coarse_solver( sp.linalg.inv( self.levels[-1].A.todense() ), on_device=False )

        self.gpu_levels = None
        self.gpu_coarse_solver = None

        if 'cupyx.scipy.sparse' in sys.modules and use_device:
            self.gpu_levels = []
            for lvl,level in enumerate(self.levels):
                R = None
                P = None
                A = cupy_sp.csr_matrix( level.A )
                D_inv = cupy.asarray( level.D_inv)
                if level.R is not None:
                    R = cupy_sp.csr_matrix(level.R)
                if level.P is not None:
                    P = cupy_sp.csr_matrix(level.P)

                self.gpu_levels.append( ml_level( A, D_inv, R=R, P=P, n_pre=level.n_pre, n_post=level.n_post, omega=level.omega, on_device = True) )
                #print("gpu level: ", lvl)

            self.gpu_coarse_solver = coarse_solver( self.coarse_solver.A_inv, on_device=True )
            #print("gpu coarse: ")

    def __repr__(self):
        """Print basic statistics about the multigrid hierarchy."""
        output = 'multilevel_solver\n'
        output += 'Number of Levels:     %d\n' % len(self.levels)
        output += 'Operator Complexity: %6.3f\n' % self.operator_complexity()
        output += 'Grid Complexity:     %6.3f\n' % self.grid_complexity()

        total_nnz = sum([level.A.nnz for level in self.levels])

        output += '  level   unknowns     nonzeros\n'
        for n, level in enumerate(self.levels):
            A = level.A
            output += '   %2d   %10d   %10d [%5.2f%%]\n' %\
                (n, A.shape[1], A.nnz,
                 (100 * float(A.nnz) / float(total_nnz)))

        if self.gpu_levels is not None:
            output += "\nEnabled GPU offloading\n"

        return output

    def cycle_complexity(self):

        nnz = [level.A.nnz for level in self.levels]

        def V(level):
            if len(self.levels) == 1:
                return nnz[0]
            elif level == len(self.levels) - 2:
                return 2 * nnz[level] + nnz[level + 1]
            else:
                return 2 * nnz[level] + V(level + 1)

        flops = V(0)

        return float(flops) / float(nnz[0])

    def operator_complexity(self):
        """Operator complexity of this multigrid hierarchy.

        Defined as:
            Number of nonzeros in the matrix on all levels /
            Number of nonzeros in the matrix on the finest level

        """
        return sum([level.A.nnz for level in self.levels]) /\
            float(self.levels[0].A.nnz)

    def grid_complexity(self):
        """Grid complexity of this multigrid hierarchy.

        Defined as:
            Number of unknowns on all levels /
            Number of unknowns on the finest level

        """
        return sum([level.A.shape[0] for level in self.levels]) /\
            float(self.levels[0].A.shape[0])

    def psolve(self, b):
        """Legacy solve interface."""
        return self.solve(b, maxiter=1)

    def aspreconditioner(self):
        from scipy.sparse.linalg import LinearOperator

        shape = self.levels[0].A.shape
        dtype = self.levels[0].A.dtype

        def matvec(b):
            return self.solve(b, maxiter=1)

        return LinearOperator(shape, matvec, dtype=dtype)

    def solve(self, b, maxiter=1):

        x = None

        if self.gpu_levels is not None:
            #print('gpu vcycle')

            A = self.gpu_levels[0].A

            exec_space = get_exec_space(A)
            #print(exec_space)
            #print(get_exec_space(b))

            assert A.dtype == b.dtype

            b_dev = None
            if exec_space == get_exec_space(b):
                #print('b is on device')
                b_dev = b
            else:
                b_dev = cupy.asarray(b)
            x_dev = exec_space.zeros_like(b_dev)
            #b_dev = cupy.asarray(b)
            #x_dev = cupy.zeros_like(b_dev)

            for iter in range(maxiter):
                if len(self.gpu_levels) == 1:
                    # hierarchy has only 1 level
                    self.gpu_coarse_solver( b_dev, x_dev )
                else:
                    self.__solve(self.gpu_levels, 0, x_dev, b_dev)

            if exec_space == get_exec_space(b):
                #print('leaving x on device')
                x = x_dev
            else:
                x = cupy.asnumpy(x_dev)
            #x = cupy.asnumpy(x_dev)

        else:

            #print('cpu vcycle')

            A = self.levels[0].A

            assert A.dtype == b.dtype

            b = np.asarray(b)
            x = np.zeros_like(b)

            for iter in range(maxiter):
                if len(self.levels) == 1:
                    # hierarchy has only 1 level
                    self.coarse_solver( b, x )
                else:
                    self.__solve(self.levels, 0, x, b)

        return x

    def residual(self, A, x, b):
        #return b - A.dot( x )
        return b - SpMV( A, x )

    def jacobi(self, level, x, b, zero_guess = False):
        if zero_guess: # Ax = 0
            x = level.omega * level.D_inv * b
        else:
            x += level.omega * level.D_inv * self.residual( level.A, x, b )

    def __solve(self, levels, lvl, x, b):

        level = levels[lvl]
        #print('starting: ', lvl, len(levels), level.R is None, level.P is None)
        A = level.A

        # presmoother (1 iteration)
        for i in range(level.n_pre):
            self.jacobi( level, x, b, (i==0) )

        #coarse_b = level.R * self.residual( A, x, b )
        coarse_b = SpMV( level.R, self.residual( A, x, b ) )
        coarse_x = level.exec_space.zeros_like(coarse_b)

        if lvl == len(levels) - 2:
            if levels[-1].on_device:
                self.gpu_coarse_solver( coarse_b, coarse_x )
            else:
                self.coarse_solver( coarse_b, coarse_x )
        else:
            self.__solve(levels, lvl + 1, coarse_x, coarse_b)

        #x += level.P * coarse_x   # coarse grid correction
        x += SpMV( level.P, coarse_x )   # coarse grid correction

        # postsmoother
        for i in range(level.n_post):
            self.jacobi( level, x, b )

        #print('finished: ', lvl)

def my_bicgstab(A, b, M=None, x0=None, tol=1e-5, maxiter=None, callback=None, residuals=None):
    """Biconjugate Gradient Algorithm with Stabilization.

    Solves the linear system Ax = b. Left preconditioning is supported.

    Parameters
    ----------
    A : array, matrix, sparse matrix, LinearOperator
        n x n, linear system to solve
    b : array, matrix
        right hand side, shape is (n,) or (n,1)
    M : array, matrix, sparse matrix, LinearOperator
        n x n, inverted preconditioner, i.e. solve M A A.H x = M b.
    x0 : array, matrix
        initial guess, default is a vector of zeros
    tol : float
        relative convergence tolerance, i.e. tol is scaled by ||r_0||_2
    maxiter : int
        maximum number of allowed iterations
    callback : function
        User-supplied function is called after each iteration as
        callback(xk), where xk is the current solution vector
    residuals : list
        residuals has the residual norm history,
        including the initial residual, appended to it

    Returns
    -------
    (xNew, info)
    xNew : an updated guess to the solution of Ax = b
    info : halting status of bicgstab

            ==  ======================================
            0   successful exit
            >0  convergence to tolerance not achieved,
                return iteration count instead.
            <0  numerical breakdown, or illegal input
            ==  ======================================

    Notes
    -----
    The LinearOperator class is in scipy.sparse.linalg.interface.
    Use this class if you prefer to define A or M as a mat-vec routine
    as opposed to explicitly constructing the matrix.  A.psolve(..) is
    still supported as a legacy.

    Examples
    --------
    >>> from pyamg.krylov.bicgstab import bicgstab
    >>> from pyamg.util.linalg import norm
    >>> import numpy as np
    >>> from pyamg.gallery import poisson
    >>> A = poisson((10,10))
    >>> b = np.ones((A.shape[0],))
    >>> (x,flag) = bicgstab(A,b, maxiter=2, tol=1e-8)
    >>> print norm(b - A*x)
    4.68163045309

    References
    ----------
    .. [1] Yousef Saad, "Iterative Methods for Sparse Linear Systems,
       Second Edition", SIAM, pp. 231-234, 2003
       http://www-users.cs.umn.edu/~saad/books.html

    """
    exec_space = get_exec_space(A)

    if b.dtype != A.dtype:
        raise TypeError('b and A types do not match: ' + b.dtype.char + ' ' + A.dtype.char)

    dtype = b.dtype

    n = A.shape[0]

    #print("bicgstab: ", n, dtype, exec_space, tol, M)

    x = None
    if x0 is None:
        x = exec_space.zeros( (n), dtype=dtype )
    else:
        x = exec_space.asarray( x0, dtype=dtype )

    if M is None:
        #print('Creating dummy idenity precon')
        from scipy.sparse.linalg import LinearOperator

        def matvec(b):
            return b

        M = LinearOperator( A.shape, matvec, dtype=dtype )

    # Check iteration numbers
    if maxiter is None:
        maxiter = n + 5
    elif maxiter < 1:
        raise ValueError('Number of iterations must be positive')

    def norm(v):
        return exec_space.linalg.norm(v)

    def dot(x,y):
        return exec_space.dot(x,y)

    def axpby(alpha, x, beta, y, z = None):
        # DAXPY-like function with optional output array.
        # z = alpha*x + beta*y ... or ...
        # y = alpha*x + beta*y
        if z is None:
            z = y

        z = alpha * x + beta * y

        return z

    psolve = None
    if hasattr( M, '__call__' ):
        def _psolve(M, b):
            return M(b)

        psolve = _psolve
    elif hasattr( M, 'psolve' ):
        def _psolve(M, b):
            return M.psolve(b)

        psolve = _psolve
    else:
        raise "Unknow psolve method"

    # Prep for method
    r = None
    if x0 is None:
        r = b.copy()
    else:
        #r = b - A*x
        r = b - SpMV(A,x)
    normr = norm(r)

    if residuals is not None:
        residuals[:] = [normr]

    # Check initial guess ( scaling by b, if b != 0,
    #   must account for case when norm(b) is very small)
    normb = norm(b)
    #print(normr, normb)
    if normb == 0.0:
        normb = 1.0
    if normr < tol*normb:
        return (x, 0)

    # Scale tol by ||r_0||_2
    if normr != 0.0:
        tol = tol*normr

    rhat = r.copy()
    p = r.copy()

    rho_old = 1.0
    alpha = 1.0
    omega = 1.0
    beta  = 1.0

    niters = 0

    # Begin BiCGStab
    while True:
        rho = dot(rhat, r)

        if niters > 0:
            beta = ( rho / rho_old ) * (alpha / omega)
            p = r + beta * (p - omega * v)

        #phat = M.psolve(p)
        phat = psolve( M, p )
        #v = A * phat
        v = SpMV(A, phat)

        alpha = rho / dot( rhat, v)
        s = r - alpha * v

        norms = norm(s)
        if norms < tol:
            x += alpha * phat
            return (x, 0)

        #shat = M.psolve(s)
        shat = psolve( M, s )
        #t = A * shat
        t = SpMV( A, shat )

        omega = dot(t, s) / dot( t, t )

        x += alpha*phat + omega*shat

        r = s - omega * t

        rho_old = rho

        niters += 1

        normr = norm(r)
        #if niters % 10 == 0:
        #    print(niters, normr, normr/normb)

        if residuals is not None:
            residuals.append(normr)

        if callback is not None:
            if hasattr(normr, 'get'): # If this is a CuPy object, pull it back.
                normr = normr.get().item()
            callback( None, normr=normr )

        if normr < tol:
            return (x, 0)

        if niters >= maxiter:
            return (x, niters)

    return (x, niters)

class bicgstab:

    def __init__(self, use_device = False):
        self.use_device = use_device

    def __call__(self, A, b, M=None, x0=None, tol=1e-5, maxiter=None, callback=None, residuals=None):
        x = None
        info = None

        if 'cupy' in sys.modules and self.use_device:
            print('bicgstab on device')
            _A = cupy_sp.csr_matrix( A )
            _b = cupy.asarray( b )

            _x, info = my_bicgstab(_A, _b, M, x0=x0, tol=tol, maxiter=maxiter, callback=callback, residuals=residuals)
            x = _x.get()
        else:
            print('bicgstab on host')
            x, info = my_bicgstab(A, b, M, x0=x0, tol=tol, maxiter=maxiter, callback=callback, residuals=residuals)

        return (x,info)
