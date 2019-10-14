import sys

delim=','
if len(sys.argv) < 3:
    print("Usage: <program> <ref-file> <test-file>")
    sys.exit(0)

skip = None
if len(sys.argv) > 3:
    skip = int(sys.argv[3])

tol = 1e-12

skip = [0]

def load_data(fl):
    dat = []
    with open(fl) as f:
        for line in f:
           dat.append( line.strip().split(',') )

    return dat

ref = load_data( sys.argv[1] )
dat = load_data( sys.argv[2] )

n = min( len(ref), len(dat) )

if len(ref) != len(dat):
    print("file lengths do not match: {} {}".format(len(ref), len(dat)))
    print("Only comparing the first {} rows".format(n))
    #sys.exit(0)

for i in range(n):
    if len(ref[i]) != len(dat[i]):
        print("record widths do not match on row: {} {} {}".format( i, len(ref[i]), len(dat[i])))
        print(ref[i])
        print(dat[i])
        sys.exit(0)

def isfloat(f):
    try:
        g = float(f)
        return True
    except ValueError:
        return False

def isint(f):
    try:
        g = int(f)
        return True
    except ValueError:
        return False

def isnumeric(f):
    return isfloat(f) or isint(f)

def error(p,q):
    p = float(p)
    q = float(q)
    if abs(p) > 1e-15:
        return abs(p-q) / abs(p)
    else:
        return abs(p-q)

passed = True
maxerr = 0.0
for i in range(n):
    ncols = len(ref[i])
    for j in range(ncols):
        if j in skip:
            continue
        else:
            p = ref[i][j]
            q = dat[i][j]
            if isnumeric(p) and isnumeric(q):
                err = error( p, q )
                if err > tol:
                    print('{},{} do not match {} {} {}'.format(i,j,p,q, err))
                    print('ref: ', ref[i])
                    print('dat: ', dat[i])
                    passed = False
                    break
                else:
                    maxerr = max(maxerr, err)

    if not passed:
        break

if passed:
    print("max error: {}".format(maxerr))
    sys.exit(0)
else:
    sys.exit(1)
