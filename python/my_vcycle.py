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

class amg_vcycle_solver:

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
            return self.vcycle(b, ncycles=1)

        return LinearOperator(shape, matvec, dtype=dtype)

    def vcycle(self, b, ncycles=1):

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

        for iter in range(ncycles):
            if len(self.levels) == 1:
                # hierarchy has only 1 level
                x = self.coarse_solve(b)
            else:
                self.__vcycle(0, x, b)

        return x

    def __vcycle(self, lvl, x, b):

        #print('starting: ', lvl)
        A = self.levels[lvl].A

        self.levels[lvl].presmoother(A, x, b)

        residual = b - A * x

        coarse_b = self.levels[lvl].R * residual
        coarse_x = np.zeros_like(coarse_b)

        if lvl == len(self.levels) - 2:
            coarse_x[:] = self.coarse_solve(coarse_b)
        else:
            self.__vcycle(lvl + 1, coarse_x, coarse_b)

        x += self.levels[lvl].P * coarse_x   # coarse grid correction

        self.levels[lvl].postsmoother(A, x, b)

        #print('finished: ', lvl)

    def coarse_solve(self, b):
        #print('coarse')
        #return sp.linalg.lu_solve(self.coarse_PLU, b)
        return self.coarse_Ainv.dot( b )

class gpu_level:
    def __init__(self, A, D_inv, R = None, P = None):
        self.A = A
        self.D_inv = D_inv
        self.R = R
        self.P = P
        self.n_pre = 1
        self.n_post = 3
        self.omega = 1.

class amg_vcycle_solver_gpu:

    def __init__(self, levels):

        self.levels = levels

        self.omega = []

        from pyamg.util.utils import scale_rows, get_block_diag, get_diagonal
        from pyamg.util.linalg import approximate_spectral_radius

        if True:
            for lvl,level in enumerate(levels[:-1]):
                rho_D_inv = None
                if hasattr(level.A, 'rho_D_inv'):
                    rho_D_inv = level.A.rho_D_inv
                    print(lvl,'reused')
                else:
                    D_inv = get_diagonal(level.A, inv=True)
                    D_inv_A = scale_rows(level.A, D_inv, copy=True)
                    rho_D_inv = approximate_spectral_radius(D_inv_A)

                # limit rho to [4/3,2)
                omega = 1.0 / rho_D_inv
                omega = 1.0 / max(4./3., min(2., rho_D_inv))
                self.omega.append( omega )
                print(lvl, omega, rho_D_inv)

        self.coarse_PLU  = sp.linalg.lu_factor( self.levels[-1].A.todense() )
        self.coarse_Ainv = sp.linalg.inv( self.levels[-1].A.todense() )

        self.gpu_levels = None
        self.gpu_coarse_Ainv = None

        if 'cupyx.scipy.sparse' in sys.modules:
            self.gpu_levels = []
            for lvl,level in enumerate(levels):
                R = None
                P = None
                assert hasattr(level,'A')
                A = cupy_sp.csr_matrix(level.A)
                D_inv_host = get_diagonal(level.A, inv=True)
                D_inv = cupy.asarray( D_inv_host)
                if hasattr(level,'R'):
                    R = cupy_sp.csr_matrix(level.R)
                if hasattr(level,'P'):
                    P = cupy_sp.csr_matrix(level.P)

                self.gpu_levels.append( gpu_level( A, D_inv, R, P ) )
                print("gpu level: ", lvl)

            self.gpu_coarse_Ainv = cupy.asarray( self.coarse_Ainv )
            print("gpu coarse: ", lvl)

            _b = cupy.ones( self.coarse_Ainv.shape[0] )
            _x = cupy.matmul( self.gpu_coarse_Ainv, _b )

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
            return self.vcycle(b, ncycles=1)

        return LinearOperator(shape, matvec, dtype=dtype)

    def vcycle(self, b, ncycles=1):

        #print('gpu vcycle')

        A = self.gpu_levels[0].A

        assert A.dtype == b.dtype

        b_dev = cupy.asarray(b)
        x_dev = cupy.zeros_like(b_dev)

        for iter in range(ncycles):
            if len(self.gpu_levels) == 1:
                # hierarchy has only 1 level
                x_dev[:] = cupy.matmul( self.gpu_coarse_Ainv, b_dev )
            else:
                self.__vcycle(0, x_dev, b_dev)

        x = cupy.asnumpy(x_dev)

        return x

    def residual(self, A, x, b):
        return b - A.dot( x )

    def jacobi(self, lvl, x, b, zero_guess = False):
        if zero_guess: # Ax = 0
            x = self.omega[lvl] * self.gpu_levels[lvl].D_inv * b
        else:
            x += self.omega[lvl] * self.gpu_levels[lvl].D_inv * self.residual( self.gpu_levels[lvl].A, x, b )

    def __vcycle(self, lvl, x, b):

        #print('starting: ', lvl)
        A = self.gpu_levels[lvl].A

        # presmoother (1 iteration)
        for i in range(self.gpu_levels[lvl].n_pre):
            self.jacobi( lvl, x, b, (i==0) )

        coarse_b = self.gpu_levels[lvl].R * self.residual( A, x, b )
        coarse_x = cupy.zeros_like(coarse_b)

        if lvl == len(self.gpu_levels) - 2:
            coarse_x[:] = cupy.matmul( self.gpu_coarse_Ainv, coarse_b )
        else:
            self.__vcycle(lvl + 1, coarse_x, coarse_b)

        x += self.gpu_levels[lvl].P * coarse_x   # coarse grid correction

        # postsmoother
        for i in range(self.gpu_levels[lvl].n_post):
            self.jacobi( lvl, x, b )

        #print('finished: ', lvl)
