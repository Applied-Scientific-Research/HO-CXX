import sys
import os
import numpy as np
import numpy.linalg as linalg
import scipy.sparse as sp
import scipy.sparse.linalg
import time
import getopt

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

packages = ['scipy']
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
           print(rowptr[:15])
           print(rowptr[-15:-1])
           print(rowptr.size)

           num_values = rowptr[num_rows]
           print( num_values )

           #print(f.readline().decode(encoding))
           colidx = np.fromfile( f, dtype=np.int32, count=num_values )
           #line = f.readline().decode(encoding) # Read trailing \n
           print(colidx[:15])
           print(colidx[-15:-1])
           print(colidx.size)

           #line = f.readline()
           #print(line)
           #print(line.decode(encoding))
           values = np.fromfile( f, dtype=np.float64, count=num_values )

           print(values[:15])
           print(values.size)

           print(values.shape)
           print(colidx.shape)
           print(rowptr.shape)

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

class Parameters:

    def __init__(self):

        self.tolerance = 1e-9
        self.absolute  = False
        self.maxiters  = 2000
        self.restart   = 16
        self.verbosity = 1
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

def get_options():

    A_file = None
    b_file = None
    x_file = None

    params = Parameters()

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hA:b:x:t:vs:m:",["help","A=","b=","rhs=","x=","tol=","maxiters=","restart=","verbosity=","numtrials=","ntrials=",'solver=','precon=','M=','abstol='])
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

# Get CLI options
A_file, b_file, x_file, params = get_options()

if (A_file is None) or (b_file is None):
    print("Missing A or b file")
    help()
    sys.exit()

A = load_matrix( A_file )
print('shape(A): ', A.shape)

b = load_vector( b_file, A.shape[0] )
print('shape(b): ', b.shape)

xref = None
if x_file is not None:
    xref = load_vector( sys.argv[3], A.shape[0] )
    print('shape(xref): ', xref.shape)

class Monitor:

    def __init__(self, A=None, b=None, takes_residual=False, keep_residuals=False, frequency=0, normb=None):

        self.A = A
        self.b = b
        self.takes_residual = takes_residual
        self.keep_residuals = keep_residuals
        self.frequency = frequency
        self.normb = normb

        self.residuals = None
        if self.keep_residuals:
            if (self.takes_residual is False) and (A is None or b is None): # Must have A and b.
                print("Must provide both A and b to evaluate the residual norm if passing xk.")
                self.keep_residuals = False
        if self.keep_residuals:
            self.residuals = []

        if self.frequency > 0 and (self.normb is None) and (self.b is not None):
            self.normb = linalg.norm( self.b )

        self.niters = 0

        #print(self)

    def op(self, xk_or_rk):

        if self.keep_residuals or (self.frequency > 0 and self.niters % self.frequency == 0):
            normr = 0.0
            if self.takes_residual:
                normr = linalg.norm( xk_or_rk )
            else:
                normr = linalg.norm( self.b - self.A * xk_or_rk )
            if self.keep_residuals:
                self.residuals.append( normr )
            if self.frequency > 0 and self.niters % self.frequency == 0:
                if self.normb is not None:
                    print("  iter: {}, {:10.5e}, {:10.5e}".format(self.niters, normr, normr / normb))
                else:
                    print("  iter: {}, {:10.5e}".format(self.niters, normr))

        self.niters += 1

    def reset(self, b=None):
        self.niters = 0

        if self.residuals is not None:
           self.residuals = []

        if b is not None:
           self.b = b
           self.normb = linalg.norm(b)

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

def create_preconditioner( A = None, name = None, opts = {} ):

    if name is None:
        return None

    package = None
    precon = name

    # Split the package name out of the request if present.
    names = name.split('.')

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

    elif precon == 'block':

        M = block_diagonal_precond( A, bs=9 )

        print("Instantiated block diagonal preconditioner")
        return M.M()

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

        M_ilu = sp.linalg.spilu( A_csc, fill_factor=fill_factor, drop_tol=drop_tol )

        def M_op (xk):
            #print("Inside ilu precon")
            return M_ilu.solve( xk )

        M = sp.linalg.LinearOperator( A.shape, M_op)
        time_end = timestamp()
        print("Instantiated ILU ( fill: {:.2f}% ) in {:.2f} ms".format( 100.*float(M_ilu.nnz) / A.nnz, 1000.*(time_end-time_start) ) )

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

        _D = block_diagonal_precond(A,9)
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

    elif precon == 'amg':

        use_pyamg = True
        if package is not None and package == 'pyamgcl':
            use_pyamg = False

        amg_type = 'ruge_stuben'
        #amg_type = 'smoothed_aggregation'

        if use_pyamg:

            time_start = timestamp()

            # Set default options
            precon_opts = {}
            precon_opts['max_levels'] = 16
            precon_opts['max_coarse'] = 500
            precon_opts['coarse_solver'] = 'lu'
            omega=0.533333333
            precon_opts['presmoother' ] = ('jacobi', {'omega':omega,'iterations':1})
            precon_opts['postsmoother'] = ('jacobi', {'omega':omega,'iterations':2})

            if amg_type == 'ruge_stuben':
                precon_opts['strength'] = ('classical', {'theta':0.5} )

            for key in opts:
                precon_opts[k] = opts[key]

            print("Building {} AMG preconditioner".format(amg_type))
            print(precon_opts)

            if amg_type == 'ruge_stuben':
                ml = pyamg.ruge_stuben_solver( A, **precon_opts )
            else:
                ml = pyamg.smoothed_aggregation_solver( A, **precon_opts )
            print(ml)

            M = ml.aspreconditioner()

            time_end = timestamp()

            print("Instantiated PyAMG ({}) preconditioner in {} (ms)".format(amg_type, 1000*(time_end-time_start)))

#           for i in range(0, len(ml.levels)-1):
#               print(ml.levels[i])
#               print(ml.levels[i].presmoother)

            return M

        elif has_pyamgcl:

            time_start = timestamp()

            precon_opts = {}
            precon_opts['coarse_enough'] = 500
            precon_opts['direct_coarse'] = True

            if amg_type == 'ruge_stuben':
                precon_opts['coarsening.type'] = 'ruge_stuben'
                precon_opts['coarsening.eps_strong'] = 0.5
            else:
                precon_opts['coarsening.type'] = 'smoothed_aggregation'

            precon_opts['relax.type'] = 'damped_jacobi'
            precon_opts['relax.damping'] = 0.53333333
            precon_opts['npre'] = 1
            precon_opts['npost'] = 2

           #precon_opts['relax.type'] = 'gauss_seidel'

            for key in opts:
                precon_opts[k] = opts[key]

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

    def __init__( self, A, name = 'gmres', opts = {} ):

        self.opts    = {}
        self.func    = None
        self.monitor = None
        self.name    = None
        self.package = None
        self.A       = A

        package = 'scipy'
        solver = name

        # Split the package name out of the request if present.
        names = name.split('.')

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

        if params.absolute:
            self.opts['tol'] = params.tolerance / linalg.norm(b)
            print("Reset tolerance to absolution: {}".format( self.opts['tol'] ))

        if self.monitor is not None:
            self.monitor.reset(b=b)

        return self.func( A=self.A, M=M, b=b, callback=self.monitor, **self.opts )

normb = linalg.norm( b )
print("normb = {}".format(normb))

verbose = (params.verbosity > 0)

bs = 9
test_bsr( A, b, bs )
#BlkM = block_diagonal_precond( A, bs )
#sys.exit(0)

if True:
    _D = block_diagonal_precond (A, bs)
    b = _D.Dinv * b
    #A = (_D.Dinv * _D.A).tocsr()
    A = (_D.Dinv * _D.A)

def residual( x ):
    return b - A * x

if True:

    solver = params.solver
    precon = params.precon

    M = create_preconditioner( A, name=precon )
    S = Solver( A, name=solver )

    print(S)

    print("Starting solver test: {}, {}".format(solver, precon) )

    run_time = 0.0
    info = None

    num_tests = 0
    max_tests = params.numtrials

#   S.opts['tol'] = params.tolerance * normb

    while num_tests < max_tests:

        time_start = timestamp()

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
        print("Failed to solve the problem in {} iterations {}".format( info, linalg.norm( residual(x) ) ) )
    elif info < 0:
        print("Error: invalid input parameter")
        sys.exit(10)
    else:
        print("Solved linear system {} times in {} iterations and {:.2f} ms".format( num_tests, niters, 1000.*run_time/num_tests ))

        if S.monitor.residuals is not None:
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

        normr = linalg.norm( residual(x) )
        print("abserr= {}, relerr= {}".format( normr, normr / normb ))
        if xref is not None:
            norm_xref = linalg.norm( residual(xref) )
            print("xref: abserr= {}, relerr= {}".format( norm_xref, norm_xref / normb ) )

# Test out PyAMGCL
if has_pyamgcl and False:
#try:
#    import pyamgcl
#except ImportError:
#    print("PyAMGCL isn't available")
#else:

    time_start = timestamp()

    opts = {}
    opts['coarse_enough'] = 500
    #opts['direct_coarse'] = True
    #opts['coarsening.type'] = 'ruge_stuben'
    #opts['coarsening.eps_strong'] = 0.5
    opts['relax.type'] = 'damped_jacobi'
    opts['relax.damping'] = 0.5333333333
    #opts['relax.type'] = 'gauss_seidel'
    opts['npre'] = 1
    opts['npost'] = 2

    M = pyamgcl.amgcl(A, opts)

    time_end = timestamp()
    print("PyAMGCL build took {} (ms)".format(1000.*(time_end-time_start)))

    print(M)

    time_start = timestamp()

    #monitor = Monitor( takes_residual=True, keep_residuals=True, frequency=int(verbose), normb=normb )
    monitor = Monitor( takes_residual=True, keep_residuals=True, frequency=1, normb=normb )

    def monitor_callback(vec):
        monitor.op(vec)

    args = {}

    #solver_name = ['scipy','gmres']
    solver_name = ['pyamgcl','fgmres']
    solver_name = ['pyamgcl','bicgstab']

    x = None
    info = None
    niters = 0

    if solver_name[0] == 'scipy':

        args['restart'] = params.restart
        args['callback'] = monitor_callback
        args['tol'] = params.tolerance
        args['maxiter'] = params.maxiters

        solver = sp.linalg.gmres
        x, info = solver( A, b, M=M, **args )

        niters = monitor.niters

    elif solver_name[0] == 'pyamgcl':

        args['type']  = solver_name[1]
        if 'gmres' in solver_name[1]:
            args['M']     = params.restart
        args['tol']   = params.tolerance
        args['maxiter'] = params.maxiters

        solver = pyamgcl.solver( P=M, prm=args )

        x = solver( b )

        print("err: {} iters: {}".format(solver.error, solver.iters ))
        info = solver.error < params.tolerance
        niters = solver.iters

    time_end = timestamp()

    print('info: {}, {}'.format( info, niters ))
    print("{}.{} with PyAMGCL solve took {} (ms)".format(solver_name[0], solver_name[1], 1000.*(time_end-time_start)))

    normr = linalg.norm( residual(x) )
    print("abserr: {} relerr: {}".format(normr, normr/normb ))

    if info == 0 and monitor.residuals is not None:
        for i in range(0,len(monitor.residuals),params.restart):
            r = monitor.residuals[i]
            print("iter {}: {:10.5e}, {:10.5e}".format(i,r,r/normb))
