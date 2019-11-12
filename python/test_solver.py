import sys
import os
import numpy as np
import numpy.linalg as linalg
import scipy.sparse as sp
import scipy.sparse.linalg
import time
import getopt
#import block_sgs
import cblock_sgs
import sp_solver
import bilum
import itsol

has_pyamg = False
try:
    import pyamg
except ImportError:
    print("PyAMG isn't available")
else:
    print("PyAMG successfully loaded")
    has_pyamg = True

has_pyamgcl = False
try:
    import pyamgcl
except ImportError:
    print("PyAMGCL isn't available")
else:
    print("PyAMGCL successfully loaded")
    has_pyamgcl = True

packages = ['scipy','itsol']
if has_pyamg:
    packages.append('pyamg')
if has_pyamgcl:
    packages.append('pyamgcl')

def timestamp():
    return time.time()

def load_matrix ( fname ):

    if not os.path.exists( fname ):
        print("Error: file does not exist: " + fname)
        sys.exit(2)

    with open( fname, 'rb') as f:
        try:
           encoding = 'utf-8'
           #num_rows    = int( f.readline().decode(encoding).strip('\n').split()[-1] )
           #num_columns = int( f.readline().decode(encoding).strip('\n').split()[-1] )
           #num_values  = int( f.readline().decode(encoding).strip('\n').split()[-1] )
           #index_type  = f.readline().decode(encoding).strip('\n').split()[-1]
           #value_type  = f.readline().decode(encoding).strip('\n').split()[-1]
           #print("num_rows: {}".format(num_rows))
           #print("num_columns: {}".format(num_columns))
           #print("num_values: {}".format(num_values))
           #print("index_type: {}".format(index_type))
           #print("value_type: {}".format(value_type))

           num_rows = np.fromfile( f, dtype=np.int32, count=1 )[0]
           num_columns = num_rows
           print( num_rows )

           #print(f.readline().decode(encoding))
           rowptr = np.fromfile( f, dtype=np.int32, count=(num_rows+1) )
           #line = f.readline().decode(encoding) # Read trailing \n
           #print(rowptr[:15])
           #print(rowptr[-15:-1])
           #print(rowptr.size)

           num_values = rowptr[num_rows]
           print( num_values )

           #print(f.readline().decode(encoding))
           colidx = np.fromfile( f, dtype=np.int32, count=num_values )
           #line = f.readline().decode(encoding) # Read trailing \n
           #print(colidx[:15])
           #print(colidx[-15:-1])
           #print(colidx.size)

           #line = f.readline()
           #print(line)
           #print(line.decode(encoding))
           values = np.fromfile( f, dtype=np.float64, count=num_values )

           #print(values[:15])
           #print(values.size)

           #print(values.shape)
           #print(colidx.shape)
           #print(rowptr.shape)

           return sp.csr_matrix((values, colidx, rowptr), shape=(num_rows, num_columns))
    
        except IOError:
           print("Could not read file: " + fname)
           sys.exit(2)

def load_vector ( fname, num_rows ):

    if not os.path.exists( fname ):
        print("Error: file does not exist: " + fname)
        sys.exit(2)

    with open( fname, 'rb') as f:
        try:
           nr = np.fromfile( f, dtype=np.int32, count=1 )[0]
           if nr != num_rows:
               print("Error: vector length does not match matrix {} {}".format( nr, num_rows ))
               sys.exit(3)

           v = np.fromfile( f, dtype=np.float64, count=num_rows )
           if v.size != num_rows:
               print("Error: loaded vector length does not match matrix {} {}".format( v.size, num_rows ))
               sys.exit(4)

           return v

        except IOError:
           print("Could not read file: " + fname)
           sys.exit(2)

def test_bsr( A, b, bs ):

    shape = A.shape

    if not (shape[0] % bs == 0 and shape[1] % bs == 0):
        print("shape does not conform to blocksize {} {} {}".format( shape[0], shape[1], bs ))

    Ablk = sp.bsr_matrix( A, blocksize=(bs,bs) )

    print("Ablk: ", Ablk.shape, Ablk.nnz, Ablk.blocksize)

    #for row in range(9):
    #    print("row: ", row)
    #    n0 = A.indptr[row]
    #    nz = A.indptr[row+1] - n0
    #    for k in range(nz):
    #        print("  {}, {}".format( A.indices[n0+k], A.data[n0+k] ))

    #print("block: ", Ablk.indptr.shape, Ablk.indices.shape)
    #for row in range(1):
    #    n0 = Ablk.indptr[row]
    #    nz = Ablk.indptr[row+1] - n0
    #    print("row: ", row, n0, nz, type(Ablk.data), Ablk.data[0].shape)
    #    for k in range(nz):
    #        print("  {}, {}".format( Ablk.indices[n0+k], Ablk.data[n0+k] ))

    r = None
    rb = None

    for i in range(2):
        if i % 2 == 0:
            t_start = timestamp()
            sumr = 0
            for i in range(100):
                r  = A    * b
                sumr += r[0]
            t_stop = timestamp()
            print("scalar: ", t_stop - t_start)
        else:
            t_start = timestamp()
            for i in range(100):
                rb = Ablk * b
                sumr += r[0]
            t_stop = timestamp()
            print("block: ", t_stop - t_start)

    print("error: ", linalg.norm( r - rb ) )

    return Ablk

class block_diagonal_precond:

    def __init__ ( self, A, bs ):

        self.bs = bs

        shape = A.shape

        if not (shape[0] % bs == 0 and shape[1] % bs == 0):
            print("shape does not conform to blocksize {} {} {}".format( shape[0], shape[1], bs ))
            sys.exit(1)

        Ablk = sp.bsr_matrix( A, blocksize=(bs,bs) )

        print("Ablk: ", Ablk.shape, Ablk.nnz, Ablk.blocksize)

        t_start = timestamp()

        nbrows = Ablk.shape[0] // bs
        _Dinv = []
        for brow in range( nbrows ):
            k0 = Ablk.indptr[brow]
            nbcols = Ablk.indptr[brow+1] - k0
            #if brow % 100 == 0:
            #    print("brow: ", brow, k0, nbcols)
            for k in range(nbcols):
                #print("  {}".format( Ablk.indices[k0+k] ))
                if Ablk.indices[k0+k] == brow:
                    binv = linalg.inv( Ablk.data[k0+k] )
                    _Dinv.append(binv)
                    #print("  {}".format( binv ))

        Dinv = np.asarray(_Dinv)

        nrows = A.shape[0]
        indptr = np.arange( nbrows+1 )
        indices = np.arange( nbrows )
        print(indptr)
        print(indices)

        self.A = Ablk
        self.Dinv = sp.bsr_matrix( (Dinv, indices, indptr), blocksize=(bs,bs) )#shape=Ablk.shape, blocksize=bs )

        t_stop = timestamp()
        print("block-diagonal build: ", t_stop - t_start)
        #print(Dinv.shape)
        #print(self.Dinv)

    def M(self):

        def M_op(x):
            return self.Dinv * x

        M = sp.linalg.LinearOperator( self.Dinv.shape, M_op )

        return M

class block_Jacobi:

    def __init__ ( self, A, bs, doJacobi = True ):

        t_start = timestamp()

        # Split A in (D + LU) and then save D^(-1)

        shape = A.shape
        if not (shape[0] % bs == 0 and shape[1] % bs == 0):
            print("shape does not conform to blocksize {} {} {}".format( shape[0], shape[1], bs ))
            sys.exit(1)

        Ab = sp.bsr_matrix( A, blocksize=(bs,bs) )

        self.bs = bs

        self.nrows = A.shape[0]
        self.nbrows = self.nrows // self.bs

        print("block_Jacobi.__init__: ", Ab.shape, Ab.nnz, Ab.blocksize,
                                     self.nbrows, Ab.nnz // (self.bs**2), Ab.indptr.shape[0])

        value_type = A.data.dtype
        index_type = A.indices.dtype
        pointer_type = A.indptr.dtype

        Dinv = np.empty( (self.nbrows, bs, bs), dtype=value_type )

        nnzb = Ab.indices.shape[0] - self.nbrows

        LU_data = np.empty( (nnzb,bs,bs), dtype=value_type )
        LU_indices = np.empty( (nnzb), dtype=index_type )
        LU_indptr = np.empty( (self.nbrows+1), dtype=pointer_type )

        L_data = []
        L_indices = []
        L_indptr = np.empty( (self.nbrows+1), dtype=pointer_type )

        U_data = []
        U_indices = []
        U_indptr = np.empty( (self.nbrows+1), dtype=pointer_type )

        LU_nnzb = 0
        L_nnzb = 0
        U_nnzb = 0

        for brow in range( self.nbrows ):
            k0 = Ab.indptr[brow]
            nbcols = Ab.indptr[brow+1] - k0

            LU_indptr[brow] = LU_nnzb
            L_indptr[brow] = L_nnzb
            U_indptr[brow] = U_nnzb

            for k in range(nbcols):
                bcol = Ab.indices[k0+k]
                if bcol == brow:
                    Dinv[brow] = linalg.inv( Ab.data[k0+k] )
                else:
                    LU_data[LU_nnzb] = np.copy( Ab.data[k0+k] )
                    LU_indices[LU_nnzb] = bcol
                    LU_nnzb += 1

                if bcol < brow:
                    L_data.append( Ab.data[k0+k] )
                    L_indices.append( bcol )
                    L_nnzb += 1
                elif bcol > brow:
                    U_data.append( Ab.data[k0+k] )
                    U_indices.append( bcol )
                    U_nnzb += 1

        LU_indptr[self.nbrows] = LU_nnzb
        if LU_nnzb != (Ab.indices.shape[0] - self.nbrows):
            print("Error: LU_nnzb mismatch {} {}".format( LU_nnzb, Ab.indices.shape[0] - self.nbrows))
            sys.exit(-1)

        D_indptr = np.arange( self.nbrows+1 )
        D_indices = np.arange( self.nbrows )

        self.Dinv = sp.bsr_matrix( (Dinv, D_indices, D_indptr), blocksize=(bs,bs) )
        self.LU   = sp.bsr_matrix( (LU_data, LU_indices, LU_indptr), blocksize=(bs,bs) )

        L_indptr[self.nbrows] = L_nnzb
        L_data = np.asarray(L_data, dtype=value_type).reshape( (L_nnzb,bs,bs) )
        L_indices = np.asarray(L_indices, dtype=index_type )
        self.L    = sp.bsr_matrix( (L_data, L_indices, L_indptr), blocksize=(bs,bs), shape=A.shape )

        U_indptr[self.nbrows] = U_nnzb
        U_data = np.asarray(U_data, dtype=value_type).reshape( (U_nnzb,bs,bs) )
        U_indices = np.asarray(U_indices, dtype=index_type )
        self.U    = sp.bsr_matrix( (U_data, U_indices, U_indptr), blocksize=(bs,bs), shape=A.shape )

        print(self.Dinv.shape)
        print(self.LU.shape)
        print(self.L.shape)
        print(self.U.shape)

        self.Jacobi = doJacobi
        if self.Jacobi:
            self.niters = 2
        else:
            self.niters = 4
        self.time = 0.0
        self.nevals = 0

        self.cblock_sgs = cblock_sgs.cblock_sgs( niters = max(1, (self.niters // 2)) )

        t_stop = timestamp()
        print("block-Jacobi build: ", t_stop - t_start, self.Jacobi)

    def __del__(self):
        ms = self.time * 1000.0
        ave = 0.0
        if self.nevals > 0:
            ave = ms / self.nevals
        print("Jacobi precond time: {} (ms) {}".format( ms, ave ))

    def M(self):

        def M_op(x):

            t_start = timestamp()

            y = None

            if self.Jacobi == False:
                # Symmetric Gauss-Seidel iteration ...

                if True:

                    # # Forward sweep
                    # # Normally this is y = x - self.U * y but our initial guess is 0.
                    # # so we just copy the rhs.
                    # y = np.copy(x)

                    # # convert to block-vectors for easier matmuls
                    # y.shape = ( self.nbrows, self.bs )

                    # for i in range( self.nbrows ):
                    #     k0 = self.L.indptr[i]
                    #     nbcols = self.L.indptr[i+1] - k0
                    #     yi = np.copy( y[i] )
                    #     for k in range( nbcols ):
                    #         col = self.L.indices[k0+k]
                    #         yi -= np.matmul( self.L.data[k0+k], y[col] )

                    #     y[i] = np.matmul( self.Dinv.data[i], yi )

                    # # Backwards sweep
                    # y.shape = ( self.nrows )

                    # y = x - self.L * y

                    # # convert to block-vectors for easier matmuls
                    # y.shape = ( self.nbrows, self.bs )

                    # for i in range( self.nbrows-1, -1, -1 ):
                    #     k0 = self.U.indptr[i]
                    #     nbcols = self.U.indptr[i+1] - k0
                    #     yi = np.copy( y[i] )
                    #     for k in range( nbcols ):
                    #         col = self.U.indices[k0+k]
                    #         yi -= np.matmul( self.U.data[k0+k], y[col] )

                    #     y[i] = np.matmul( self.Dinv.data[i], yi )

                    # y.shape = ( self.nrows )

                    #y = np.copy(x)
                    #y = block_sgs.block_sgs( self.L, self.Dinv, self.U, x, y )

                    y = self.cblock_sgs( self.L, self.Dinv, self.U, x )

                else:

                    xb = np.reshape( x, (self.nbrows, self.bs) )
                    yb = np.zeros_like(xb)

                    def op( brow ):

                        k0 = self.LU.indptr[brow]
                        nbcols = self.LU.indptr[brow+1] - k0
                        ybb = np.copy( xb[brow] )
                        for k in range(nbcols):
                            col = self.LU.indices[k0+k]
                            ybb -= np.matmul( self.LU.data[k0+k], yb[col] )

                        yb[brow] = np.matmul( self.Dinv.data[brow], ybb )

                    for n in range( max(1, (self.niters // 2)) ):

                        if n % 2 == 0:
                            # Forward
                            for brow in range( self.nbrows ):
                                op(brow)
                        else:
                            # Backward
                            for brow in range( self.nbrows-1, -1, -1 ):
                                op(brow)

                    y = np.reshape( yb, (self.nrows) )

            else:

                # Jacobi iteration ...
                ytmp = [ np.zeros_like(x), np.empty_like(x) ]

                old = 0
                new = 1
                for n in range(self.niters):
                    new = 1 - old

                    if n == 0:
                       ytmp[new] = self.Dinv * x
                    else:
                       ytmp[new] = self.Dinv * ( x - self.LU * ytmp[old] )

                    old = (old + 1) % 2

                t_stop = timestamp()
                self.time += (t_stop - t_start)

                y = ytmp[new]

            t_stop = timestamp()
            self.time += (t_stop - t_start)
            self.nevals += 1

            return y

        M = sp.linalg.LinearOperator( self.Dinv.shape, M_op )

        return M

class Parameters:

    def __init__(self):

        self.tolerance = 1e-9
        self.absolute  = False
        self.maxiters  = 2000
        self.restart   = 16
        self.verbosity = 0
        self.numtrials = 1
        self.solver    = 'gmres'
        self.precon    = 'ilu'

    def __str__(self):

        s = "Parameters:\n"
        s += "  tolerance: {}\n".format(self.tolerance)
        s += "   absolute: {}\n".format(self.absolute)
        s += "   maxiters: {}\n".format(self.maxiters)
        s += "    restart: {}\n".format(self.restart)
        s += "  verbosity: {}\n".format(self.verbosity)
        s += "  numtrials: {}\n".format(self.numtrials)
        s += "     solver: {}\n".format(self.solver)
        s += "     precon: {}\n".format(self.precon)
        return s;

def help():
    print('{} -h (--help) -A (--A=) <CSR matrix> -b (--rhs) <rhs vector> -t (--tol) <tolerance>'.format(sys.argv[0]))

def get_options( argv ):

    A_file = None
    b_file = None
    x_file = None

    params = Parameters()

    try:
        opts, args = getopt.getopt( argv,"hA:b:x:t:vs:m:",["help","A=","b=","rhs=","x=","tol=","maxiters=","restart=","verbosity=","numtrials=","ntrials=",'solver=','precon=','M=','abstol='])
    except getopt.GetoptError:
        help()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            help()
            sys.exit()
        elif opt in ("-A", "--A"):
            A_file = arg
        elif opt in ("-b", "--b", "--rhs"):
            b_file = arg
        elif opt in ("-x", "--x"):
            x_file = arg
        elif opt in ("-t", "--tol"):
            params.tolerance = float(arg)
        elif opt in ("--abstol"):
            params.tolerance = float(arg)
            params.absolute  = True
        elif opt in ("-v"):
            params.verbosity += 1
        elif opt in ("--verbosity"):
            params.verbosity = int(arg)
        elif opt in ("--maxiters"):
            params.maxiters = int(arg)
        elif opt in ("--restart"):
            params.restart = int(arg)
        elif opt in ("--numtrials", "--ntrials"):
            params.numtrials = int(arg)
        elif opt in ("-s", "--solver"):
            params.solver = arg
        elif opt in ("-m", "--precon", "--M"):
            params.precon = arg

    print('A file is ', A_file)
    print('b file is ', b_file)
    print('x file is ', x_file)
    print(params)

    return A_file, b_file, x_file, params

class Monitor:

    def __init__(self, A=None, b=None, takes_residual=False, keep_residuals=False, frequency=0, normb=None):

        self.A = A
        self.b = b
        self.takes_residual = takes_residual
        #self.keep_residuals = keep_residuals
        self.frequency = frequency
        self.normb = normb

        self.residuals = []
        #if self.keep_residuals:
        #    if (self.takes_residual is False) and (A is None or b is None): # Must have A and b.
        #        print("Must provide both A and b to evaluate the residual norm if passing xk.")
        #        self.keep_residuals = False

        if self.frequency > 0 and (self.normb is None) and (self.b is not None):
            self.normb = linalg.norm( self.b )

        self.niters = 0

        #print(self)

    def op(self, xk_or_rk):

        #if self.keep_residuals or (self.frequency > 0 and self.niters % self.frequency == 0):
        normr = None
        if self.takes_residual:
            normr = linalg.norm( xk_or_rk )
        else:
            if self.b is None:
                print("Must specify b in Monitor")
                sys.exit(-1)
            normr = linalg.norm( self.b - self.A * xk_or_rk )

        #if self.keep_residuals:
        #    self.residuals.append( normr )
        self.residuals.append( normr )

        if self.frequency > 0 and self.niters % self.frequency == 0:
            if self.normb is not None:
                print("  iter: {}, {:10.5e}, {:10.5e}".format(self.niters, normr, normr / self.normb))
            else:
                print("  iter: {}, {:10.5e}".format(self.niters, normr))

        self.niters += 1

    def reset(self, b=None):
        self.niters = 0

        self.residuals = []

        self.b = b

        if self.b is None:
            print("Must specify b in Monitor")
            sys.exit(-1)
        else:
            self.normb = linalg.norm(b)

    def geometric_rate(self):
        l = len(self.residuals)
        if l > 1:
            if self.residuals[0] == 0:
                return -0.0
            return ( self.residuals[-1] /  self.residuals[0] ) ** (1.0 / float(l))
        else:
            return 0.0

    def average_rate(self):
        l = len(self.residuals)
        if l > 1:
            ave = 0.0
            for i in range(1,l):
                if self.residuals[i-1] == 0:
                    return -0.0
                ave += self.residuals[i] / self.residuals[i-1]
            return ave / float(l)
        else:
            return 0.0

    def __call__(self, xk_or_rk):
        self.op( xk_or_rk )

    def __str__(self):
        s  = "Monitor:\n"
        s += "  A: {}\n".format(self.A is not None)
        s += "  b: {}\n".format(self.b is not None)
        s += "  takes_residual: {}\n".format( self.takes_residual )
        s += "  keep_residuals: {}\n".format( self.keep_residuals )
        s += "  frequency: {}\n".format( self.frequency )
        s += "  normb: {}\n".format(self.normb is not None)
        return s

#def create_preconditioner( A = None, name = None, params = Parameters(), opts = {} ):
def create_preconditioner( A = None, params = Parameters(), opts = {} ):

    name = params.precon
    if name is None:
        return None

    package = None
    precon = name

    # Split the package name out of the request if present.
    names = name.split(':')

    if len(names) > 1:
        precon = names[1].lower()
        package = names[0].lower()
    else:
        precon = names[0].lower()

    if package is not None and not (package in( packages )):
        print('package {} is not available.'.format(package))
        sys.exit(1)

    if precon == 'identity':

        def M_op (xk):
            #print("Inside identity precon")
            return xk

        print("Instantiated identity preconditioner")

        return sp.linalg.LinearOperator( A.shape, M_op)

    elif precon == 'block-diag':

        M = block_diagonal_precond( A, bs=9 )

        print("Instantiated block diagonal preconditioner")
        return M.M()

    elif precon == 'block-jacobi' or precon == 'block-sgs':

        M = block_Jacobi( A, 9, precon == 'block-jacobi' )

        print("Instantiated block Jacobi preconditioner")
        return M.M()

    elif package == 'itsol':

        time_start = timestamp()

        print("Building ITSOL ({}.{}) Preconditioner".format(package,precon))

        M_itsol = itsol.itsol( A, precon )

        def M_op (xk):
            #print("Inside itsol precon")
            return M_itsol.solve( xk )

        M = sp.linalg.LinearOperator( A.shape, M_op)

        time_stop = timestamp()
        print("Instantiated ITSOL preconditioner in {:.4f} s".format( time_stop-time_start))

        return M

    elif precon == 'bilum':

        time_start = timestamp()

        print("Building BILUM")

        M_bilum = bilum.bilum( A )

        def M_op (xk):
            #print("Inside bilum precon")
            return M_bilum.solve( xk )

        M = sp.linalg.LinearOperator( A.shape, M_op)

        time_stop = timestamp()
        print("Instantiated BILUM preconditioner in {:.4f} s".format( time_stop-time_start))

        return M

    elif 'ilu' in precon:

        time_start = timestamp()
        A_csc = None
        if isinstance(A, sp.csc_matrix ):
            A_csc = A
        else:
            A_csc = sp.csc_matrix(A)

        drop_tol = None
        fill_factor = None

        if precon == 'ilu0':
           fill_factor = 1
        elif 'ilut' in precon:
           #drop_tol = 1e-3
           #fill_factor = 4.5

           # Look for ilut(#.#,#.#)
           first = precon.find('[')
           last = precon.find(']')
           if first != -1 and last != -1:
               tmp = precon[first+1:last]
               _opts = tmp.split(',')
               for o in _opts:
                   # split about the '='
                   kv = o.split('=')
                   if len(kv) != 2:
                       print("key/value pair not value {}".format(kv))
                   else:
                       key = kv[0]
                       val = kv[1]

                       if key == 'drop':
                           drop_tol = float(val)
                       elif key == 'fill':
                           fill_factor = float(val)
                       else:
                           print("Unknown ilut option: ", key, val)

        print("Building ILU ({}) with drop={} and fill={}".format(precon, drop_tol, fill_factor))

        if False:
            Ad = sp.csc_matrix([[1., 0., 0., 1.],
                                [5., 0., 2., 2.],
                                [1., 2., 3., 1.],
                                [0.,-1., 0., 3.]], dtype=np.double)
            print("Ad:", Ad.todense())
            print("Ad.data:", Ad.data)
            print("Ad.colptr:", Ad.indptr)
            print("Ad.rowind:", Ad.indices)
            LU = sp.linalg.splu( Ad )
            print(LU)
            print("Pc:", LU.perm_c)
            print("Pr:", LU.perm_r)
            print("L:", type(LU.L), LU.L)
            print("U:", type(LU.U), LU.U)
            b = np.array([1.,2.,3.,4.], dtype=np.double).reshape(4,1)
            x0 = LU.solve(b)
            print("b:", b)
            print("x:", x0)
            print("Ax:", Ad * x0)

            print("Ld:", LU.L.todense())
            print("Ud:", LU.U.todense())

            px = np.zeros_like(b)
            pb = np.zeros_like(b)
            y  = np.zeros_like(b)

            m = Ad.shape[0]

            Pr = np.zeros( (m,m), dtype=np.double)
            Pc = np.zeros( (m,m), dtype=np.double)
            for j in range(m):
                i = LU.perm_r[j]
                Pr[i,j] = 1.0

                i = LU.perm_c[j]
                Pc[j,i] = 1.0

            print("Pr:\n", Pr)
            print("inv(Pr):\n", linalg.inv(Pr))
            inv_perm_r = np.empty_like(LU.perm_r)
            for i in range(LU.perm_r.shape[0]):
                j = LU.perm_r[i]
                inv_perm_r[j] = i
            print("inv(r):\n", inv_perm_r)
            print("Pc:\n", Pc)
            print("Pr * b:\n", Pr.dot(b))

            # Forward solve (Pr)Ly = (Pr)b
            pb[ LU.perm_r[:] ] = b[:]
            print(pb)

            #pb = Pr.dot(b)
            #print(pb)
            pb.shape = (4,1)
            #print(pb.shape)

            Ld = LU.L.todense()
            Ud = LU.U.todense()

            y = np.matmul( linalg.inv(Ld), pb )
            #print(y.shape)
            #print(Ld.shape)
            #print(Ud.shape)
            #print(pb.shape)

            y[0,0] = pb[0,0]
            for i in range(1,m):
                xi = pb[i,0]
                for j in range(0,i):
                   xi -= Ld[i,j] * y[j,0]
                y[i,0] = xi # Lii == 1

            px = np.matmul( linalg.inv(Ud), y )

            #for i in range(m-1,-1,-1):
            #    xi = y[i,0]
            #    for j in range(i+1,m):
            #        xi -= Ud[i,j] * y[j,0]
            #    px[i,0] = xi / Ud[i,i]

            x = np.zeros_like(b)
            x = px[ LU.perm_c[:] ]
            #x = Pc.dot(px)
            print(x)
            print("Ax:", Ad * x)
            print("Ad:", Ad.todense())

            ipiv = np.arange(m)
            y, info = scipy.linalg.lapack.dgetrs( Ld, inv_perm_r, b )
            #y, info = scipy.linalg.lapack.dgetrs( Ld, LU.perm_r, b )
            pb[ LU.perm_r[:] ] = b[:]
            y, info = scipy.linalg.lapack.dgetrs( Ld, ipiv, pb )
            #ipiv = np.empty_like(inv_perm_r)
            #for i in range(m):
            #    ipiv[i] = i
            px, info = scipy.linalg.lapack.dgetrs( Ud, ipiv, y )
            x = px[ LU.perm_c[:] ]
            x = np.matmul( Pc, px )
            print("Lapack")
            print(x)
            print("Ax:", Ad * x)

            print("spsolve_triangular")
            # Apply the row permutation to the rhs vector
            pb[ LU.perm_r[:] ] = b[:]

            Ltmp = LU.L.tocsr()
            Utmp = LU.U.tocsr()
            y  = sp.linalg.spsolve_triangular( Ltmp, pb, lower=True, overwrite_A=True, overwrite_b=True)#, unit_diagonal=True)
            px = sp.linalg.spsolve_triangular( Utmp,  y, lower=False, overwrite_A=True, overwrite_b=True)#, unit_diagonal=False)
            x = px[ LU.perm_c[:] ]
            print(x)
            print("Ax:", Ad * x)

            #sys.exit(1)

        M_ilu = None

        use_superlu = False

        mixed_precision = False
        if mixed_precision:
            A_data = A_csc.data.astype(np.float32)
            A_mixed = sp.csc_matrix( (A_data, A_csc.indices, A_csc.indptr), shape=A_csc.shape)
            M_ilu = sp.linalg.spilu( A_mixed, fill_factor=fill_factor, drop_tol=drop_tol )
        else:
            M_ilu = sp.linalg.spilu( A_csc, fill_factor=fill_factor, drop_tol=drop_tol )

        time_end = timestamp()
        print("Finished ILU ( fill: {:.2f}% ) in {:.2f} ms".format( 100.*float(M_ilu.nnz) / A.nnz, 1000.*(time_end-time_start) ) )

        M_sp_solver = None
        if not use_superlu:
            M_L = M_ilu.L.tocsr()
            M_U = M_ilu.U.tocsr()
            M_pr = M_ilu.perm_r
            M_pc = M_ilu.perm_c

            print("M_L: {0:.2f}%".format(100.*(M_L.nnz / (M_L.shape[0]**2))))
            #print("M_L:\n", M_L[:100])
            #print("M_L.data:\n", M_L.data[:100])
            #print("M_L.indptr:\n", M_L.indptr[:100])
            #print("M_L.indices:\n", M_L.indices[:100])

            print("M_U: {0:.2f}%".format(100.*(M_U.nnz / (M_U.shape[0]**2))))
            #print("M_U:\n", M_U[:100])
            #print("M_U.data:\n", M_U.data[:100])
            #print("M_U.indptr:\n", M_U.indptr[:100])
            #print("M_U.indices:\n", M_U.indices[:100])

            #if True:
            #   bs = 4
            #   print(bs)
            #   M_Lb = sp.bsr_matrix( M_L, blocksize=(bs,bs))
            #   M_Ub = sp.bsr_matrix( M_U, blocksize=(bs,bs))
            #   print("M_Lb: {} {} {} {:.2f}".format(M_Lb.nnz, M_Lb.blocksize, M_Lb.indptr[-1], 100.*float(M_Lb.nnz)/(M_L.shape[0]**2)))
            #   print("M_Ub: {} {} {} {:.2f}".format(M_Ub.nnz, M_Ub.blocksize, M_Ub.indptr[-1], 100.*float(M_Ub.nnz)/(M_U.shape[0]**2)))
            #   print(M_Lb.data[0])
            #   print(linalg.inv(M_Lb.data[0]))
            #   print(M_Ub.data[0])
            #   print(linalg.inv(M_Ub.data[0]))

            #   M_L = M_Lb
            #   M_U = M_Ub

            #print(M_ilu)
            #print(M_ilu.perm_c.shape, M_ilu.perm_c[:10])
            #print(M_ilu.perm_r.shape, M_ilu.perm_r[:10])
            #print(M_ilu.L.shape, M_ilu.L.nnz, isinstance(M_ilu.L, sp.csc_matrix))
            #print(M_ilu.U.shape, M_ilu.U.nnz, isinstance(M_ilu.U, sp.csc_matrix))

            M_sp_solver = sp_solver.sp_solver( M_L, M_U, M_pr, M_pc )

        def M_op (xk):
            #print("Inside ilu precon")
            if use_superlu:
                return M_ilu.solve( xk )
            else:
                return M_sp_solver.solve(xk)

        M = sp.linalg.LinearOperator( A.shape, M_op)
        print("Instantiated ILU preconditioner")

        return M

    elif precon in ('krylov', 'gmres'):

        print("Inside krylov precon")

        solver = sp.linalg.gmres
        if package is not None:
            if package in ('scipy'):
                solver = sp.linalg.gmres
            elif package in ('pyamg'):
                solver = pyamg.krylov.gmres

        #solver = sp.linalg.bicgstab

        precon_opts = {}
        precon_opts['tol'] = 1e-3

        monitor = Monitor( takes_residual=True, keep_residuals=True, frequency=1 )

        precon_opts['maxiter'] = params.restart
        #precon_opts['restrt'] = params.restart
        #precon_opts['orthog'] = 'mgs'
        precon_opts['callback'] = monitor

        if solver == pyamg.krylov.gmres:
            precon_opts['maxiter'] = params.restart #// 2
            #precon_opts['restrt'] = None
            #precon_opts['callback'] = monitor
            #precon_opts['orthog'] = 'mgs'

        #for key in opts:
        #    precon_opts[k] = opts[key]

        #def I(x):
        #   return x

        #_M = sp.linalg.LinearOperator( A.shape, I)

        #_D = block_diagonal_precond(A,9)
        #_M = _D.M()
        _D = block_Jacobi(A,9,False)
        _M = _D.M()

        def M_op (xk):

            #zk = None
            #info = None

            #if solver == pyamg.krylov.gmres:
            #    zk, info = solver( A, xk, M=None, **precon_opts )
            #else:
            #    zk, info = solver( A, xk, M=None, **precon_opts )

            print("calling M: ", solver);
            #zk, info = solver( A, xk, M=None, **precon_opts )
            zk, info = solver( A, xk, M=_M, **precon_opts )
            print("info: ", info);

            return zk

        print("Instantiated krylov(gmres) preconditioner: ", M_op)
        return sp.linalg.LinearOperator( A.shape, M_op)

    elif 'amg' in precon:

        use_pyamg = True
        if package is not None and package == 'pyamgcl':
            use_pyamg = False

        amg_type = 'ruge_stuben'
        #amg_type = 'smoothed_aggregation'

        print("Matrix format: ", type(A))

        if use_pyamg:

            time_start = timestamp()

            # Set default options
            precon_opts = {}
            precon_opts['max_levels'] = 16
            precon_opts['max_coarse'] = 500
            #precon_opts['max_coarse'] = 50
            precon_opts['coarse_solver'] = 'lu'
            omega=0.533333333
            precon_opts['presmoother' ] = ('jacobi', {'omega':omega,'iterations':1})
            precon_opts['postsmoother'] = ('jacobi', {'omega':omega,'iterations':2})
            #precon_opts['presmoother' ] = ('jacobi', {'iterations':1})
            #precon_opts['postsmoother'] = ('jacobi', {'iterations':2})

            if amg_type == 'ruge_stuben':
                precon_opts['strength'] = ('classical', {'theta':0.5} )
                #precon_opts['strength'] = ('algebraic_distance')
                #precon_opts['CF'] = 'PMIS'
            elif amg_type == 'smoothed_aggregation':
                #precon_opts['strength'] = None
                #precon_opts['strength'] = 'classical'
                #precon_opts['strength'] = 'affinity'
                #precon_opts['strength'] = ('algebraic_distance')
                #precon_opts['strength'] = ('evolution')
                #precon_opts['smooth'] = ('energy')
                pass

            #for key in opts:
            #    precon_opts[k] = opts[key]

            print("Building {} AMG preconditioner".format(amg_type))
            print(precon_opts)

            if amg_type == 'ruge_stuben':
                ml = pyamg.ruge_stuben_solver( A, **precon_opts )
            else:
                ml = pyamg.smoothed_aggregation_solver( A, **precon_opts )
            print(ml)

            #M = ml.aspreconditioner()
            M = ml.aspreconditioner(cycle='V')
            #M = ml.aspreconditioner(cycle='W')
            #M = ml.aspreconditioner(cycle='F')

            time_end = timestamp()

            print("Instantiated PyAMG ({}) preconditioner in {} (ms)".format(amg_type, 1000*(time_end-time_start)))

#           for i in range(0, len(ml.levels)-1):
#               print(ml.levels[i])
#               print(ml.levels[i].presmoother)

            return M

        elif has_pyamgcl:

            time_start = timestamp()

            precon_opts = {}
            precon_opts['coarse_enough'] = 50
            precon_opts['direct_coarse'] = True

            if amg_type == 'ruge_stuben':
                precon_opts['coarsening.type'] = 'ruge_stuben'
                precon_opts['coarsening.eps_strong'] = 0.5
            else:
                precon_opts['coarsening.type'] = 'smoothed_aggregation'
                #precon_opts['coarsening.estimate_spectral_radius'] = 'true'
                precon_opts['coarsening.type'] = 'smoothed_aggr_emin'

            #precon_opts['relax.type'] = 'damped_jacobi'
            #precon_opts['relax.damping'] = 0.53333333
            precon_opts['npre'] = 1
            precon_opts['npost'] = 2

            #precon_opts['relax.type'] = 'gauss_seidel'

            #for key in opts:
            #    precon_opts[k] = opts[key]

            M = pyamgcl.amgcl(A, precon_opts)

            time_end = timestamp()

            print(M)

            print("Instantiated PyAMGCL ({}) preconditioner in {} (ms)".format(amg_type, 1000.*(time_end-time_start)))

            return M

        else:

            print('AMG preconditioner requested but not available.')
            sys.exit(1)

    else:

        print("Unknown precon: " + name)
        sys.exit(2)

class Solver:

    def __init__( self, A, params = Parameters(), opts = {} ):

        name = params.solver

        self.opts    = {}
        self.func    = None
        self.monitor = None
        self.name    = None
        self.package = None
        self.A       = A
        self.params  = params

        package = 'scipy'
        solver = name

        # Split the package name out of the request if present.
        names = name.split(':')

        if len(names) > 1:
            solver = names[1].lower()
            package = names[0].lower()
        else:
            solver = names[0].lower()

        if package and not ( package in(packages) ):
            print('package {} is not available.'.format(package))
            sys.exit(1)

        self.name    = solver
        self.package = package
        print('solver: ', solver, package)

        if 'gmres' in solver:

            if solver == 'fgmres':
                if has_pyamg is False:
                    print("fgmres requires pyamg")
                    sys.exit(1)
                else:
                    self.func = pyamg.krylov.fgmres
            else:
                if package == 'pyamg' and has_pyamg:
                    self.func = pyamg.krylov.gmres
                elif package == 'scipy':
                    self.func = sp.linalg.gmres

            self.opts = {}
            self.opts['tol']     = params.tolerance
            self.opts['maxiter'] = params.maxiters

            self.opts['restrt'] = params.restart
            #if solver == 'fgmres':
            #    self.opts['restrt'] = params.restart
            #else:
            #    self.opts['restart'] = params.restart

            for key in opts:
                self.opts[k] = opts[key]

            freq = 0
            if params.verbosity > 0:
                freq = params.restart
                if params.verbosity > 1 and solver == 'fgmres':
                    freq = 1

            self.monitor = Monitor( takes_residual=True, keep_residuals=True, frequency=freq )

        elif solver == 'bicgstab':

            if package == 'pyamg':
                self.func = pyamg.krylov.bicgstab
            elif package == 'scipy':
                self.func = sp.linalg.bicgstab

            self.opts = {}
            self.opts['tol']     = params.tolerance
            self.opts['maxiter'] = params.maxiters

            for key in opts:
                self.opts[k] = opts[key]

            freq = 0
            if params.verbosity > 0:
                freq = 1

            self.monitor = Monitor( A=A, takes_residual=False, keep_residuals=True, frequency=freq )

        elif solver == 'cg':

            if package == 'pyamg':
                self.func = pyamg.krylov.cg
            elif package == 'scipy':
                self.func = sp.linalg.cg

            self.opts = {}
            self.opts['tol']     = params.tolerance
            self.opts['maxiter'] = params.maxiters

            for key in opts:
                self.opts[k] = opts[key]

            freq = 0
            if params.verbosity > 0:
                freq = 1

            self.monitor = Monitor( A=A, takes_residual=False, keep_residuals=True, frequency=freq )

        if self.func is None:
            print("Solver={} in package={} is not available".format(solver,package))
            sys.exit(1)

        print("Instantiated solver={} from package={}".format(solver,package))

    def __str__(self):
        s  = "Solver:\n"
        s += "      name = {}\n".format(self.name)
        s += "   package = {}\n".format(self.package)
        s += "      opts = {}\n".format(self.opts)
        return s

    def __call__(self, b, M=None):

        if self.params.absolute:
            self.opts['tol'] = self.params.tolerance / linalg.norm(b)
            print("Reset tolerance to absolution: {}".format( self.opts['tol'] ))

        if self.monitor is not None:
            self.monitor.reset(b=b)

        return self.func( A=self.A, M=M, b=b, callback=self.monitor, **self.opts )

def main( argv ):

    # Get CLI options
    A_file, b_file, x_file, params = get_options( argv )
    
    if (A_file is None) or (b_file is None):
        print("Missing A or b file")
        help()
        sys.exit()
    
    A = load_matrix( A_file )
    #print('shape(A): ', A.shape)
    
    b = load_vector( b_file, A.shape[0] )
    #print('shape(b): ', b.shape)
    
    xref = None
    if x_file is not None:
        xref = load_vector( x_file, A.shape[0] )
        print('shape(xref): ', xref.shape)
    
    
    normb = linalg.norm( b )
    print("normb = {}".format(normb))
    
    verbose = (params.verbosity > 0)
    
    bs = 9
    Ab = test_bsr( A, b, bs )
    #A = Ab
    #BlkM = block_diagonal_precond( A, bs )
    #sys.exit(0)
    
    #if True:
    #    _D = block_diagonal_precond (A, bs)
    #    b = _D.Dinv * b
    #    #A = (_D.Dinv * _D.A).tocsr()
    #    A = (_D.Dinv * _D.A)
    
    def residual( x ):
        return b - A * x
    
    if True:
    
        solver = params.solver
        precon = params.precon
    
        M = create_preconditioner( A, params=params )#name=precon )
        S = Solver( A, params=params )#name=solver )
    
        print(S)
    
        print("Starting solver test: {}, {}".format(solver, precon) )
    
        run_time = 0.0
        info = None
    
        num_tests = 0
        max_tests = params.numtrials
    
    #   S.opts['tol'] = params.tolerance * normb
    
        while num_tests < max_tests:
    
            time_start = timestamp()
    
            S.monitor.reset(b=b)
    
            x, info = S( b, M=M )
    
            time_end = timestamp()
            run_time += (time_end - time_start)
            num_tests += 1
            if info != 0:
                break
    
            niters = S.monitor.niters
    
            if max_tests > 1:
                print("test {} finished in {:.2f} ms".format( num_tests, 1000.*(time_end-time_start) ))
    
        if info > 0:
            print("Failed to solve the problem in {} iterations {} and {} (ms)".format( info, linalg.norm( residual(x) ), 1000.*run_time) )
        elif info < 0:
            print("Error: invalid input parameter")
            sys.exit(10)
        else:
            print("Solved linear system {} times in {} iterations and {:.2f} ms".format( num_tests, niters, 1000.*run_time/num_tests ))
    
            if S.monitor.residuals is not None and params.verbosity > 1:
                print("Residuals:")
                print("iteration, norm(r), norm(r)/norm(b)")
                print("------------------------------------")
                rm1 = None
                for i in range(0,len(S.monitor.residuals)):
                    r = S.monitor.residuals[i]
                    print_this = True
                    if rm1 is not None:
                        if r == 0.0:
                            print_this = False
                        elif ( abs(r - rm1) / r) < 1e-6:
                            print_this = False
                    if print_this:
                        print("{}, {:10.5e}, {:10.5e}".format(i, r, r/normb))
                    rm1 = r
    
            print("Geometric convergence rate: {}".format( S.monitor.geometric_rate() ) )
            print("Average   convergence rate: {}".format( S.monitor.average_rate() ) )
    
            normr = linalg.norm( residual(x) )
            print("abserr= {}, relerr= {}".format( normr, normr / normb ))
            if xref is not None:
                norm_xref = linalg.norm( residual(xref) )
                print("xref: abserr= {}, relerr= {}".format( norm_xref, norm_xref / normb ) )

if __name__ == "__main__":
    main( sys.argv[1:] )

## Test out PyAMGCL
#if has_pyamgcl and False:
##try:
##    import pyamgcl
##except ImportError:
##    print("PyAMGCL isn't available")
##else:
#
#    time_start = timestamp()
#
#    opts = {}
#    opts['coarse_enough'] = 500
#    #opts['direct_coarse'] = True
#    #opts['coarsening.type'] = 'ruge_stuben'
#    #opts['coarsening.eps_strong'] = 0.5
#    opts['relax.type'] = 'damped_jacobi'
#    opts['relax.damping'] = 0.5333333333
#    #opts['relax.type'] = 'gauss_seidel'
#    opts['npre'] = 1
#    opts['npost'] = 2
#
#    M = pyamgcl.amgcl(A, opts)
#
#    time_end = timestamp()
#    print("PyAMGCL build took {} (ms)".format(1000.*(time_end-time_start)))
#
#    print(M)
#
#    time_start = timestamp()
#
#    #monitor = Monitor( takes_residual=True, keep_residuals=True, frequency=int(verbose), normb=normb )
#    monitor = Monitor( takes_residual=True, keep_residuals=True, frequency=1, normb=normb )
#
#    def monitor_callback(vec):
#        monitor.op(vec)
#
#    args = {}
#
#    #solver_name = ['scipy','gmres']
#    solver_name = ['pyamgcl','fgmres']
#    solver_name = ['pyamgcl','bicgstab']
#
#    x = None
#    info = None
#    niters = 0
#
#    if solver_name[0] == 'scipy':
#
#        args['restart'] = params.restart
#        args['callback'] = monitor_callback
#        args['tol'] = params.tolerance
#        args['maxiter'] = params.maxiters
#
#        solver = sp.linalg.gmres
#        x, info = solver( A, b, M=M, **args )
#
#        niters = monitor.niters
#
#    elif solver_name[0] == 'pyamgcl':
#
#        args['type']  = solver_name[1]
#        if 'gmres' in solver_name[1]:
#            args['M']     = params.restart
#        args['tol']   = params.tolerance
#        args['maxiter'] = params.maxiters
#
#        solver = pyamgcl.solver( P=M, prm=args )
#
#        x = solver( b )
#
#        print("err: {} iters: {}".format(solver.error, solver.iters ))
#        info = solver.error < params.tolerance
#        niters = solver.iters
#
#    time_end = timestamp()
#
#    print('info: {}, {}'.format( info, niters ))
#    print("{}.{} with PyAMGCL solve took {} (ms)".format(solver_name[0], solver_name[1], 1000.*(time_end-time_start)))
#
#    normr = linalg.norm( residual(x) )
#    print("abserr: {} relerr: {}".format(normr, normr/normb ))
#
#    if info == 0 and monitor.residuals is not None:
#        for i in range(0,len(monitor.residuals),params.restart):
#            r = monitor.residuals[i]
#            print("iter {}: {:10.5e}, {:10.5e}".format(i,r,r/normb))
