# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 09:44:20 2020

@author: chrst
"""

import sys
import numpy as np
import enum
# import json

from sortedcontainers import SortedKeyList
  
import time
getTimeStamp = time.time


@enum.unique
class GeometricTypes(enum.IntEnum):
    POINT = 0 # 0d
    LINE2 = 1 # 1d
    TRI3  = 4 # 2d
    QUAD4 = 5
    TET4  = 8 # 3d
    HEX8  = 9


def nVerts(gtype: GeometricTypes) -> int:
    return {GeometricTypes.POINT: 1,
            GeometricTypes.LINE2: 2,
            GeometricTypes.TRI3:  3,
            GeometricTypes.QUAD4: 4,
            GeometricTypes.TET4:  4,
            GeometricTypes.HEX8:  8}.get(gtype)


def nDims(gtype: GeometricTypes) -> int:
    return {GeometricTypes.POINT: 0,
            GeometricTypes.LINE2: 1,
            GeometricTypes.TRI3:  2,
            GeometricTypes.QUAD4: 2,
            GeometricTypes.TET4:  3,
            GeometricTypes.HEX8:  3}.get(gtype)

def nFaces(gtype: GeometricTypes) -> int:
    return {GeometricTypes.POINT: 0,
            GeometricTypes.LINE2: 0,
            GeometricTypes.TRI3:  3,
            GeometricTypes.QUAD4: 4,
            GeometricTypes.TET4:  4,
            GeometricTypes.HEX8:  6}.get(gtype)


BndryIds = None
BndryIdMap = None
    
class Element:
    
    def __init__(self, etype, verts = [], eidx = None):
        self.etype  = GeometricTypes(etype)
        self.verts  = verts
        self.eindex = eidx
        self.faces  = []
        self.neighbors = []
        self.isRight = None
        self.bndrys = []
        
    def __str__(self):
        s = '{{idx: {}, type: {}, verts: {}, faces: {}, nghbrs: {}, right: {}, bnd: {}}}'.format(
            self.eindex, self.etype.name, self.verts, self.faces, self.neighbors,
            self.isRight, self.bndrys)
        return s

    # def getBndryFaces(self):
        # return [BndryIds[-k] if k < 0 else BndryIds[0] for k in self.neighbors]


# Load a mesh file in su2 format.
def load_mesh_su2(filename):
    try:
        f = open(filename,"r")
    except:
        print("SU2 mesh file open failed: {}".format(filename))
        sys.exit(1)

    def split_keyword(s, key = None, astype = None):
        if s.find("=") == -1:
            print("SU2 keyword line missing = separator: {}".format(s))
            sys.exit(1)
            
        items = s.split('=')
        if len(items) != 2:
            print("SU2 key/value pair invalid: {}".format(s))
            sys.exit(1)
            
        if key is not None:
            if key != items[0].strip():
                print("SU2 key {} does not match {}".format(key, items[0]))

        if astype is not None:
            return astype(items[1].strip())
        else:
            return items[1].strip()
    
    try:
        # NDIME= #
        ndims = split_keyword(f.readline(), 'NDIME', int)
        print('ndims: {}'.format(ndims))
        assert ndims == 2 or ndims == 3

        # NELEM= #
        nelems = split_keyword(f.readline(), 'NELEM', int)
        print('nelems: {}'.format(nelems))
    except:
        print("SU2 Error: failed to read ndim/nelem header")
        sys.exit(1)

    elems = []
    try:
        for i in range(nelems):
            items = f.readline().split()
            etype = int(items[0])
            n  = len(items[1:-1])
            v  = np.array([int(c) for c in items[1:-1]], dtype='i')
            if etype == 12: # HEX8
                assert len(v) == 8
                assert ndims == 3
                elems.append(Element(GeometricTypes.HEX8, v, eidx=i))
            elif etype == 9: # QUAD's
                # print(items, npts)
                assert ndims == 2
                elems.append(Element(GeometricTypes.QUAD4, v, eidx=i))
            else:
                print("SU2: unknown element type {} {}".format(etype, items))
                sys.exit(1)
    except:
        print("SU2 Error: failed to read the elements")
        sys.exit(1)
        
    try:
        # NPOIN= #
        npts = split_keyword(f.readline(), 'NPOIN', int)
        print('npts: {}'.format(npts))
    except:
        print("SU2 Error: failed to read the # of points")
        sys.exit(1)
        
    coords = np.empty((npts,ndims), dtype='d')
    # print(coords.shape, ndims, nverts)
    try:
        for i in range(npts):
            items = f.readline().split()
            assert len(items) == (ndims+1)
            idx = int(items[-1])
            for j in range(ndims):
                coords[idx][j] = float(items[j])
    except:
        print("SU2 Error: failed to read the vertices")
        sys.exit(1)

    for v in elems[0].verts:
        print(v, coords[v])

    try:
        nmarkers = split_keyword(f.readline(), 'NMARK', int)
        print('nmarkers: {}'.format(nmarkers))
    except:
        print("SU2 Error: failed to read the # of boundary markers")
        sys.exit(1)

    bndry_sections = {}        
    for i in range(nmarkers):
        tag = split_keyword(f.readline(), 'MARKER_TAG', str)
        _nelems = split_keyword(f.readline(), 'MARKER_ELEMS', int)
        print(tag, _nelems)
        assert tag not in bndry_sections
        bndry_sections[tag] = bnd = []
        for j in range(_nelems):
            items = f.readline().split()
            etype = int(items[0])
            n = len(items[1:])
            v = np.array([int(c) for c in items[1:]], dtype='i')
            if etype == 9: # QUAD4
                assert ndims == 3
                assert len(items) == 5
                bnd.append(Element(GeometricTypes.QUAD4, v))
            elif etype == 3: # LINE
                assert ndims == 2
                bnd.append(Element(GeometricTypes.LINE2, v))
            else:
                print("SU2 Error: invalid element in bndry section {} {}".format(tag, items))
                sys.exit(1)

        
    # Verify the input data.
    for i in range(nelems):
        e = elems[i]
        if not all([(e.verts >= 0).all(), (e.verts < npts).all()]):
            print("Element {} has invalid vertices {}".format(i, e.verts[:]))
            sys.exit(1)

    for tag, faces in bndry_sections.items():
        for i, f in enumerate(faces):
            if not all([(f.verts >= 0).all(), (f.verts < npts).all()]):
                print("Bndry {} face {} has invalid verts {}".format(tag, i, f))

    return ndims, elems, coords, bndry_sections


FaceElementMaps = {
    # South, East, North, West, In, Out
    # - Nodes are (0,0,0), (1,0,0), (1,1,0), (0,1,0),
    #             (0,0,1), (1,0,1), (1,1,1), (0,1,1)
    # - All face normals point out from the element center.
    # - First corner is closest to 0,0,0
    GeometricTypes.HEX8:
        np.array([[0,1,5,4],
                  [1,2,6,5],
                  [3,7,6,2],
                  [0,4,7,3],
                  [0,3,2,1],
                  [4,5,6,7]], dtype='i'),
        
    GeometricTypes.QUAD4:
        np.array([[0,1],
                  [1,2],
                  [2,3],
                  [3,0]], dtype='i')
        }
    
def face_normal(coords, fverts, ftype):
    # CGNS definition of normal is V(0->1) x V(0->2)
    x0 = coords[fverts[0]]
    x1 = coords[fverts[1]]
    x2 = coords[fverts[2]]
    v01 = x1-x0
    v02 = x2-x0
    return np.cross(v01, v02)
    

def create_face_list(elems):
    
    """
    Given a list of elements, form the unique list of faces.
    """
    def find_face_in_list(faces, f, nverts):
        
        def _find(faces, idx, start):
            try:
                i = faces.index(idx, start)
                return i
            except ValueError:
                return -1
            
        n = len(faces)
        # print("looking for: ", f, n//nverts)
        assert n % nverts == 0
        start = 0
        while start < n:
            idx = _find(faces, f[0], start)
            if idx == -1:
                break
            i = idx // nverts # Outer index of face.
            j = i * nverts # Inner index into face vertices.
            fj = faces[j:j+nverts]
            # print('idx: ', idx, start, i, j, fj)
            if all([v in fj for v in f[1:]]):
                # print('found: ', i)
                return i, True
            start = j+nverts
            
        # face doesn't exist. append to list.
        faces.extend(f)
        return n // nverts, False
    
  
    class face_list:
        def __init__(self, fidx, verts):
            self.fidx = fidx
            self.verts = verts
            
        def __repr__(self):
            return '({},{})'.format(self.fidx, self.verts)
            
    def comp_verts(it):
        return min(it.verts)
    
    def comp_fidx(it):
        return it.fidx

    all_faces = {}
    for el in elems:

        ftype = None
        if el.etype == GeometricTypes.HEX8:
            ftype = GeometricTypes.QUAD4
        elif el.etype == GeometricTypes.QUAD4:
            ftype = GeometricTypes.LINE2
        else:
            raise NameError('Element type {} supported so far.'.format(el.etype._name_))

        if not ftype in all_faces:
            all_faces[ftype] = SortedKeyList(key=comp_verts)

        faces = all_faces[ftype]
        n = len(faces)

        efaces = el.verts[FaceElementMaps[el.etype]]
        nverts = nVerts(ftype)
        # for i, f in enumerate(efaces.tolist()):
        #     fidx, match = find_face_in_list(faces, f, nverts)
        #     el.faces.append(fidx)
        #     # print('face: ', i, f, fidx, len(elems)*6, match)
        
        for i, f in enumerate(efaces.tolist()):
            #print(i, f, len(faces))
            fit = face_list(-1,f)
            lo = faces.bisect_left(fit)
            hi = faces.bisect_right(fit)
            #print('lo ', lo, 'nil' if lo >= n else faces[lo])
            #print('hi ', hi, 'hil' if hi >= n else faces[hi])
            fs = sorted(f)
            #print(fs)
            fidx = -1
            for j in range(lo, hi):
                if fs == sorted(faces[j].verts):
                    # Match found
                    fidx = faces[j].fidx
                    break

            if fidx == -1:
                # Face doesn't exist yet. Append new face.
                fidx = n
                faces.add(face_list(n, f))
                n += 1
                #print('added')
                
            el.faces.append(fidx)
                
            
    for ftype in all_faces:
        faces = all_faces[ftype]
        n = len(faces)
        nverts = nVerts(ftype)
        # new_faces = np.array(faces, dtype='i').reshape(n//nverts,nverts)
        new_faces0 = SortedKeyList(faces, key=comp_fidx)
        new_faces = np.array([it.verts for it in new_faces0], dtype='i')
        all_faces[ftype] = new_faces
        print(ftype, len(faces), n, nverts, new_faces.shape)
        
    return all_faces


def create_connectivity(elements, all_faces):
    
    # Assuming all faces are the same.
    assert len(all_faces) == 1
    faces = None
    ftype = None
    for k in all_faces.keys():
        assert k == GeometricTypes.QUAD4 or k == GeometricTypes.LINE2

        faces = all_faces[k]
        ftype = k

    # And all elements are the same.
    for e in elements:
        assert e.etype == elements[0].etype
        assert e.etype == GeometricTypes.QUAD4 or e.etype == GeometricTypes.HEX8

    # Build a full list of all face indices across the elements
    # and turn into a 2d np array for quick searching.

    nelems = len(elements)
    nfaces = len(faces)
    etype  = elements[0].etype
    nfaces_per_element = nFaces(etype)
    
    t_start = getTimeStamp()

    # efaces = np.empty((nelems,nfaces_per_element), dtype='i')
    efaces = []
    for i, e in enumerate(elements):
        e.faces = np.array(e.faces, dtype='i')
        # efaces[i,:] = e.faces[:]
        efaces.extend([(fidx, i) for fidx in e.faces])
        e.neighbors = np.full((nfaces_per_element), -1, dtype='i')
        e.isRight = np.full((nfaces_per_element), 0, dtype='i')

    #print("Prep: {}".format(getTimeStamp() - t_start))

    pairs = np.zeros((nfaces,2), dtype='i')

    if False:

        efaces = []
        for i, e in enumerate(elements):
            efaces.extend([(fidx, i) for fidx in e.faces])

        # This will fail if nfaces isn't uniform.
        #print('efaces: ', efaces.shape, nfaces)
        # print('efaces: ', len(efaces), nfaces)
        # print(efaces[:10])
        efaces = sorted(efaces) #sorted(efaces, key=lambda it: it[0])
        # print(efaces[:10])

        import bisect
        
        # Populate an array with left/right element neighbor for each face.
        for i in range(nfaces):
            #els = np.where(efaces == i)[0]
            #assert els.shape[0] == 1 or els.shape[0] == 2
            #pairs[i,0] = els[0]
            #pairs[i,1] = -1 if els.shape[0] == 1 else els[1]
            
            j = bisect.bisect_left(efaces, (i, -1))
            assert j < len(efaces)
            left = efaces[j][1]
            right = -1
            if j != len(efaces)-1:
                if efaces[j+1][0] == i:
                    right = efaces[j+1][1]
            pairs[i,0] = left
            pairs[i,1] = right
            
        print('face pairing (search): {}'.format(getTimeStamp() - t_start))

    else:

        #rowptr = [0]
        rowptr = [nfaces_per_element * i for i in range(nelems+1)]
        colidx = []
        #values = []

        for i, e in enumerate(elements):
            #for j, f in enumerate(e.faces):
            #    colidx.append(f)
            #    #values.append(j)
            colidx.extend([f for f in e.faces])

            #rowptr.append(rowptr[-1] + nfaces_per_element)

        # print(rowptr[:10])
        # print(values[:20])
        # print(colidx[:20])

        import scipy.sparse as sparse

        csr_matrix = sparse.csr_matrix((np.ones(len(colidx), dtype='i'), #np.array(values, dtype='i'),
                                        np.array(colidx, dtype='i'),
                                        np.array(rowptr, dtype='i')),
                                       #shape=(nelems, nfaces),
                                       dtype='i')

        csc_matrix = csr_matrix.tocsc()

        # print('csr:\n', csr_matrix.shape)
        # print('data   :\n', csr_matrix.data.shape, csr_matrix.data[:10])
        # print('indices:\n', csr_matrix.indices.shape, csr_matrix.indices[:10])
        # print('indptr :\n', csr_matrix.indptr.shape, csr_matrix.indptr[:10])

        # print('csc:\n', csc_matrix)
        # print('csc:\n', csc_matrix.has_sorted_indices)
        assert csc_matrix.has_sorted_indices
        # print('data   :\n', csc_matrix.data[:10])
        # print('indices:\n', csc_matrix.indices[:10])
        # print('indptr :\n', csc_matrix.indptr[:10])

        #print("Matrix: {}".format(getTimeStamp() - t_start))
    
        #pairs0 = np.zeros((nfaces,2), dtype='i')

        # Scan over each column (i.e., each face) and look for the # of entries.
        # If 1 entry, this is an unconnected face.
        # If 2, it's connected.
        colptr = csc_matrix.indptr
        rowidx = csc_matrix.indices
        for j in range(nfaces):
            offset = colptr[j]
            n = csc_matrix.indptr[j+1] - offset
            pairs[j,0] = rowidx[offset]
            assert n == 1 or n == 2
            pairs[j,1] = -1 if n == 1 else rowidx[offset+1]

        #print("Face pairing (matrix): {}".format(getTimeStamp() - t_start))

        #for i in range(10):
        #    print(pairs[i,:], pairs0[i,:])

        #print(np.all(pairs == pairs0))


    # print(pairs[:10,:])    
    # Add neighbor data into element list.
    for fidx in range(nfaces):
        eleft, eright = pairs[fidx,0], pairs[fidx,1]
        # Find the local face index for the left element
        e = elements[eleft]
        face_lidx = np.where(e.faces == fidx)[0]
        if face_lidx.shape[0] != 1:
            print('Left local face not found: ', fidx, pairs[fidx,:], face_lidx, e)
        i = face_lidx[0]
        e.neighbors[i] = eright

        # Same from the right side (if it exists).
        if eright != -1:
            e = elements[eright]
            face_lidx = np.where(e.faces == fidx)[0]
            if face_lidx.shape[0] != 1:
                print('Right local face not found: ', fidx, pairs[fidx,:], face_lidx, e)
            i = face_lidx[0]
            e.neighbors[i] = eleft
            e.isRight[i] = True

    #print('element update: {}'.format(getTimeStamp() - t_start))

    return pairs


def set_boundary_faces(elems, all_faces, bndry_sections):
    
    """
    Given an ndarray (nfs,nes) of faces, find the index of the face.
    """
    def find_face_index(faces, f):
        assert faces.shape[1] == f.shape[0]

        rows = np.where(faces == f[0])[0]
        for i in rows:
            row = faces[i]
            if all([j in row for j in f[1:]]):
                return True, i

        return False, -1

    assert len(all_faces) == 1
    faces = None
    faces = []
    for key in all_faces:
        assert key == GeometricTypes.QUAD4 or GeometricTypes.LINE2
        # faces = all_faces[key]
        for i, f in enumerate(all_faces[key]):
            faces.append((sorted(f.tolist()), i))

    # print(faces.shape)
    print(faces[:10])
    
    sfaces = SortedKeyList(faces, key=lambda it: it[0])
    print(sfaces[:10])
    
    for bnd in bndry_sections:
        print(bnd)
        bfs = bndry_sections[bnd]
        bfaces = []
        for bf in bfs:
            #found, fidx = find_face_index(faces, bf.verts)
            sbf = sorted(bf.verts)
            left = sfaces.bisect_left((sbf, -1))
            right = sfaces.bisect_right((sbf,-1))
            # print(left, right)
            found, fidx = (right - left == 1), sfaces[left][1]
            # assert found
            if not found:
                print("Error: bnd face not found {}".format(bf.verts))
            #print(bf.verts, fidx, quads[fidx])
            bfaces.append(fidx)
        bndry_sections[bnd] = np.array(bfaces, dtype='i')
        print(bndry_sections[bnd][:10])


    bndryIdMap = {None: 0}
    bndryIds = [None]
    for bnd in bndry_sections:
        i = len(bndryIds)
        bndryIdMap[bnd] = i
        bndryIds.append(bnd)

    print(bndryIds)
    print(bndryIdMap)

    for ei, e in enumerate(elems):
        nfs = nFaces(e.etype)
        ebnds = np.full((nfs), bndryIdMap[None], dtype='i')
        for i, (f, n) in enumerate(zip(e.faces, e.neighbors)):
            #print(e.eindex, i, f, n)
            if n < 0:
                # Find which boundary this matches with. Or error out.
                bnd = None
                for b in bndry_sections:
                    if f in bndry_sections[b]:
                        # Are any bfaces repeated?
                        assert bnd is None
                        bnd = b
                if bnd is None:
                    print("no bnd face found for face {} in e {}".format(i, e))
                #print(e.eindex, i, f, n, bnd)
                ebnds[i] = bndryIdMap[bnd]
        e.bndrys = ebnds
        if ei < 10:
            print(e.eindex, e.neighbors, [bndryIds[k] for k in e.bndrys])

    return bndryIdMap


def create_conn_matrix_from_vec(vec):
    m = np.zeros((3,3), dtype='i')
    for j, v in enumerate(vec):
        i = abs(v)-1
        m[i,j] = 1 if v >= 0 else -1

    return m
    

def create_neighbor_mappings(elements):
    
    """For each face-neighbor, find the nodal mapping from left to right"""
    def ivec(i,j,k):
        return np.array([i,j,k], dtype='i')

    S, E, N, W, I, O = 0, 1, 2, 3, 4, 5
    Orient = {0: "S", 1: "E", 2: "N", 3: "W", 4: "I", 5: "O"}

    normal_fmaps = {}
    normal_fmaps[E,S] = [ 2, 0, 0]
    normal_fmaps[E,N] = [-2, 0, 0]
    normal_fmaps[E,E] = [-1, 0, 0]
    normal_fmaps[E,W] = [ 1, 0, 0]
    normal_fmaps[W,W] = [-1, 0, 0]
    normal_fmaps[S,N] = [ 0, 2, 0]
    normal_fmaps[S,S] = [ 0,-2, 0]
    normal_fmaps[S,W] = [ 0,-1, 0]
    normal_fmaps[N,W] = [ 0, 1, 0]
    normal_fmaps[N,N] = [ 0,-2, 0]
    normal_fmaps[I,O] = [ 0, 0, 3]
    
    def create_fmap_matrix(vec):
        m = np.zeros((3,3), dtype='i')
        for j, v in enumerate(vec):
            i = abs(v)
            m[i,j] = 1 if v >= 0 else -1

        return m
    
    def create_fmap_vector(m):
        vec = np.zeros((3), dtype='i')
        for j in range(3):
            row = -1
            for i in range(3):
                if m[i,j] != 0:
                    row = i
                    break
            if row != -1:
                i = row
                vec[j] = i if m[i,j] > 0 else -i
        return vec
    
    if True:
        def transpose(v):
            i = [idx for idx, val in enumerate(v) if val != 0][0]
            j = abs(v[i]) - 1
            k = (i+1) if v[i] > 0 else -(i+1)
            vt = [0,0,0]
            vt[j] = k
            return vt
            
        for i in range(6):
            for j in range(6):
                # if j < i and (j,i) in normal_fmaps:
                if not (i,j) in normal_fmaps and (j,i) in normal_fmaps: 
                    normal_fmaps[i,j] = transpose(normal_fmaps[j,i])

    for (i,j) in normal_fmaps:
        print(i, j, normal_fmaps[i,j], Orient[i], Orient[j])
        
    
    ijk2v = {}
    ijk2v[0,0,0] = 0
    ijk2v[1,0,0] = 1
    ijk2v[1,1,0] = 2
    ijk2v[0,1,0] = 3
    ijk2v[0,0,1] = 0+4
    ijk2v[1,0,1] = 1+4
    ijk2v[1,1,1] = 2+4
    ijk2v[0,1,1] = 3+4

    print(ijk2v)
    
    v2ijk = {}
    for (i,j,k), v in ijk2v.items():
        v2ijk[v] = (i,j,k)
        
    print(v2ijk)
    
    def np_find(mask):
        assert np.count_nonzero(mask) == 1
        return np.where(mask)[0][0]
    
    face_map_vec = {}

    for ei, e in enumerate(elements[:10]):
        etype = e.etype
        print('el: {} {}'.format(ei, etype._name_))
        assert etype == GeometricTypes.QUAD4 or etype == GeometricTypes.HEX8
        
        for fi, (n, r) in enumerate(zip(e.neighbors, e.isRight)):
            if n != -1:
                # if a shared face. Find the left/right local face indices.
                fidx = e.faces[fi]
                if r == True:
                    print("Skipping right side faces {}".format(fidx))
                    continue

                assert not fidx in face_map_vec
                #if fidx in face_map_vec:
                #    continue

                fj = np_find(elements[n].faces == fidx)
                if not (fi,fj) in normal_fmaps:
                    print('faces {} {} not in normal_fmaps {} {}'.format(fi, fj, Orient[fi], Orient[fj]))
                    sys.exit(2)
                print(" face: {}, {} neigh: {}, {}, {} dir: {} => {}".format(fi, fidx, n, r, fj, Orient[fi], Orient[fj]))
                # dijk = np.zeros((3), dtype='i')
                dijk = np.array(normal_fmaps[fi,fj], dtype='i')
                if etype == GeometricTypes.QUAD4:
                    dijk[2] = 3
                    
                rorigin = None
                if fi in (E,W):
                    i = 0 if fi == W else 1
                    # dijk[0] = normal_fmaps[fi,fj]
                    
                    lv1 = ijk2v[i,0,0]
                    rv1 = np_find(elements[n].verts == e.verts[lv1]) # local id
                    print(lv1, rv1, e.verts[lv1])
                    rorigin = v2ijk[rv1]
                    
                    # +j
                    lv2 = ijk2v[i,1,0]
                    rv2 = np_find(elements[n].verts == e.verts[lv2])
                    dj = ivec(*v2ijk[rv2]) - ivec(*v2ijk[rv1])
                    idj = np_find(dj != 0)
                    dijk[1] = (idj+1) if dj[idj] > 0 else -(idj+1)
                    print(lv2, rv2, e.verts[lv2])
                    print('j+ ', dj, idj, dj[idj])
                    
                    # +k
                    if etype == GeometricTypes.HEX8:
                        lv3 = ijk2v[i,0,1]
                        rv3 = np_find(elements[n].verts == e.verts[lv3])
                        dk = ivec(*v2ijk[rv3]) - ivec(*v2ijk[rv1])
                        idk = np_find(dk != 0)
                        dijk[2] = (idk+1) if dk[idk] > 0 else -(idk+1)
                        print(lv3, rv3, e.verts[lv3])
                        print('k+ ', dk, idk, dk[idk])

                    # print("{} => {}".format(Orient[fi], Orient[fj]))
                    # print(lv1, lv2, lv3, rv1, rv2, rv3, e.verts[lv1], e.verts[lv2], e.verts[lv3], v2ijk[lv1], v2ijk[rv1])
                    
                    
                    

                elif fi in (S,N):
                    j = 0 if fi == S else 1
                    # dijk[1] = normal_fmaps[fi,fj]
                    
                    lv1 = ijk2v[0,j,0]
                    rv1 = np_find(elements[n].verts == e.verts[lv1]) # local id
                    rorigin = v2ijk[rv1]
                    print(lv1, rv1, e.verts[lv1])
                    
                    # +i
                    lv2 = ijk2v[1,j,0]
                    rv2 = np_find(elements[n].verts == e.verts[lv2])
                    di = ivec(*v2ijk[rv2]) - ivec(*v2ijk[rv1])
                    idi = np_find(di != 0)
                    dijk[0] = (idi+1) if di[idi] > 0 else -(idi+1)
                    print(lv2, rv2, e.verts[lv2])
                    print('i+ ', di, idi, di[idi])
                    
                    if etype == GeometricTypes.HEX8:
                        # +k
                        lv3 = ijk2v[0,j,1]
                        rv3 = np_find(elements[n].verts == e.verts[lv3])
                        dk = ivec(*v2ijk[rv3]) - ivec(*v2ijk[rv1])
                        idk = np_find(dk != 0)
                        dijk[2] = (idk+1) if dk[idk] > 0 else -(idk+1)
                        
                        print(lv3, rv3, e.verts[lv3])
                        print('k+ ', dk, idk, dk[idk])

                    # print("N/S")
                    # print(lv1, rv1, lv2, rv2, dj, idj, dj[idj])
                    # print(lv1, rv1, lv3, rv3, dk, idk, dk[idk])
                    # print(lv1, lv2, lv3, rv1, rv2, rv3, e.verts[lv1], e.verts[lv2], e.verts[lv3])
                    # print(lv1, lv2, lv3, rv1, rv2, rv3, e.verts[lv1], e.verts[lv2], e.verts[lv3], v2ijk[lv1], v2ijk[rv1])
                    

                    
                else:
                    raise 'Not implemented yet'
                    
                print(dijk)
                print(create_conn_matrix_from_vec(dijk))
                if fidx not in face_map_vec:
                    el, er = ((ei, fi), (n, fj)) if r == 0 else ((n, fj), (ei, fi))
                    face_map_vec[fidx] = dijk, el, er, rorigin
                
                # Create a little internal mesh on face and verify that it's
                # the same on both sides.

    # print(face_map_vec)
    return face_map_vec



def test_mesh_hex8(testid):

    ndims = 3
    elems = []
    coords = None
    bndry_sections = {}
    
    print(testid)
    dleft, dright = testid.upper().split(',')
    print(dleft, dright)

    # The input defines the downward face for the left/right elements. 'South'
    # is natural.
    rot2d = {'S': [1,2,3],
             'W': [2,-1,3],
             'N': [-1,-2,3],
             'E': [-2,1,3]}

    left = rot2d[dleft]
    right = rot2d[dright]
    # left = [1,2,3]
    # right = [2,-1,3]
    print(left, right)

    # Create a simple 2-cell mesh.
    nx, ny, nz = 3, 2, 2

    def ijk2n(i,j,k):
        return i + nx*j + nx*ny*k

    coords = []
    for z in range(0,nz):
        for y in range(0,ny):
            for x in range(0,nx):
                coords.append([x,y,z])
    coords = np.array(coords, dtype='f')

    elems = []
    for k in range(0,nz-1):
        for j in range(0,ny-1):
            for i in range(0,nx-1):
                nodes = [ijk2n(i,j,k),
                         ijk2n(i+1,j,k),
                         ijk2n(i+1,j+1,k),
                         ijk2n(i,j+1,k),
                         ijk2n(i,j,k+1),
                         ijk2n(i+1,j,k+1),
                         ijk2n(i+1,j+1,k+1),
                         ijk2n(i,j+1,k+1)]

                verts = np.array(nodes, dtype='i')
                eidx = i + j * (nx-1) + k * (nx-1)*(ny-1)
                elems.append(Element(GeometricTypes.HEX8, verts, eidx))
                
    # create fake bnd's
    if True:
        bndry_sections['empty'] = bndry = []
        for i in [0, nx-1]:
            for k in range(0,nz-1):
                for j in range(0,ny-1):
                    verts = [ijk2n(i,j,k),
                             ijk2n(i,j+1,k),
                             ijk2n(i,j,k+1),
                             ijk2n(i,j+1,k+1)]
                    verts = np.array(verts, dtype='i')
                    bndry.append(Element(GeometricTypes.QUAD4, verts))
                    
        for j in [0, ny-1]:
            for k in range(0,nz-1):
                for i in range(0,nx-1):
                    verts = [ijk2n(i,j,k),
                             ijk2n(i+1,j,k),
                             ijk2n(i,j,k+1),
                             ijk2n(i+1,j,k+1)]
                    verts = np.array(verts, dtype='i')
                    bndry.append(Element(GeometricTypes.QUAD4, verts))
                    
        for k in [0, nz-1]:
            for j in range(0,ny-1):
                for i in range(0,nx-1):
                    verts = [ijk2n(i,j,k),
                             ijk2n(i+1,j,k),
                             ijk2n(i+1,j+1,k),
                             ijk2n(i,j+1,k)]
                    verts = np.array(verts, dtype='i')
                    bndry.append(Element(GeometricTypes.QUAD4, verts))


    def rotate_element(rotvec):

        rotmat = create_conn_matrix_from_vec(rotvec)
 
        ijk2v = {}
        ijk2v[0,0,0] = 0
        ijk2v[1,0,0] = 1
        ijk2v[1,1,0] = 2
        ijk2v[0,1,0] = 3
        ijk2v[0,0,1] = 0+4
        ijk2v[1,0,1] = 1+4
        ijk2v[1,1,1] = 2+4
        ijk2v[0,1,1] = 3+4

        # print(ijk2v)
        
        v2ijk = {}
        for (i,j,k), v in ijk2v.items():
            v2ijk[v] = (i,j,k)
            
        # print(v2ijk)
        
        v = [[0,0,0],
             [1,0,0],
             [1,1,0],
             [0,1,0],
             [0,0,1],
             [1,0,1],
             [1,1,1],
             [0,1,1]]
        v0 = np.array(v, dtype='i')
        v1 = np.zeros_like(v0)
        
        for i, row in enumerate(v):
            v1[i] = rotmat.dot(row)
            
        offset = np.zeros((3), dtype='i')
        for d in range(3):
            offset[d] = np.min(v1[:,d])

        # print(rotvec)
        # print(rotmat)
        # print(v0)
        # print(v1)
        # print(offset)
        v1 -= offset
        # print(v1)
        r = [ijk2v[tuple(v)] for v in v1]

        return r


    for i, rotvec in enumerate([left, right]):
        if rotvec == [1,2,3]:
            pass
        else:
            print("Rotating element {}: {}".format(i, rotvec))
            r = rotate_element(rotvec)
            print(r)
            v = elems[i].verts[:]
            print(v)
            elems[i].verts = v[r] 
            print(elems[i].verts)
    
    # if right == 0:
    #     e = elems[1]
    #     v = [4,1,2,5]
    #     v.extend([i+6 for i in v])
    #     e.verts = np.array(v, dtype='i')
    #     print(e)
    # elif right == 1:
    #     e = elems[1]
    #     v = [5,4,1,2]
    #     v.extend([i+6 for i in v])
    #     e.verts = np.array(v, dtype='i')
    #     print(e)
    # elif right == 2:
    #     e = elems[1]
    #     v = e.verts
    #     p = [1,2,3,0,5,6,7,4]
    #     print(e)
    #     print(p)
    #     print(v[p])
    #     e.verts = v[p]
    #     print(e)
    # elif right == 3:
    #     pass
    # else:
    #     raise "test right id {} not known".format(right)


    return ndims, elems, coords, bndry_sections




def main():

    import getopt

    infile = 'Step_Expansion.su2'
    testid = None
    outfile = None

    if len(sys.argv) > 1:

        try:
            opts, args = getopt.getopt(sys.argv[1:], 'f:t:o:', ['file=', 'test=', 'out='])
        except getopt.GetoptError as err:
            print('getopt error {}'.format(str(err)))
            sys.exit(2)

        print(opts)

        for opt, arg in opts:
            if opt in ('-f', '--file'):
                infile = arg
            elif opt in ('-o', '--out'):
                outfile = arg
            elif opt in ('-t', '--test'):
                testid = arg
            else:
                print('Unknown option {}'.format(opt))
                sys.exit(1)

    if testid is None:
        ndims, elems, coords, bndry_sections = load_mesh_su2(infile)
    else:
        ndims, elems, coords, bndry_sections = test_mesh_hex8(testid)

    original_data = {'ndims': ndims,
                     'elements': [e.verts.copy() for e in elems],
                     'bndfaces': {}
                    }

    for key in bndry_sections:
        bnd = bndry_sections[key]
        original_data['bndfaces'][key] = [e.verts.copy() for e in bnd]
    
    
    report = True
    for i, e in enumerate(elems):
        if not nDims(e.etype) == ndims:
            print("element {} {} is not {}d".format(i, e.etype._name_, ndims))
            sys.exit(1)
            
        # Drop to P1 elements if high-order
        n = len(e.verts)
        nv = nVerts(e.etype)
        if n > nv:
            e.verts = e.verts[:nv]
            if report:
                print("Found high-order element {} ... reducing to P1 {}".format(n, nv))
                report = False
                
    
    report = True
    for bnd_name in bndry_sections:
        print(bnd_name)
        bnd_faces = bndry_sections[bnd_name]
        for i, f in enumerate(bnd_faces):
            if not nDims(f.etype) == (ndims-1):
                print("bnd element {} {} is not {}d-1".format(i, f.etype._name_, ndims))
                sys.exit(1)
                
            # Drop to linear elements if high-oder
            n = len(f.verts)
            nv = nVerts(f.etype)
            if n > nv:
                f.verts = f.verts[:nv]
                if report:
                    print("Found high-order bnd element {} ... reducing to P1 {}".format(n, nv))
                    report = False


    print('midpoints')
    for i, e in enumerate(elems[:10]):
        n = nVerts(e.etype)
        print(i, e.verts, sum(coords[e.verts])/n)


    t_start = getTimeStamp()
    
    all_faces = create_face_list(elems)
    
    print("create_face_list took: {} secs".format(getTimeStamp()-t_start))

    t_start = getTimeStamp()
    
    pairs = create_connectivity(elems, all_faces)
    
    print("create_connectivity took: {} secs".format(getTimeStamp()-t_start))

    
    # for e in elems[:30]:
    #     print(e)

    t_start = getTimeStamp()
    
    bndryIdMap = set_boundary_faces(elems, all_faces, bndry_sections)
    
    print("set_boundary_faces took: {} secs".format(getTimeStamp()-t_start))


    # Export the data into a simple format for reading.
    if outfile:
        print('Writing connectivity to {}'.format(outfile))

        with open(outfile, 'w') as f:
            f.write('ndims:\n{}\n'.format(ndims))
            f.write('nelems:\n{}\n'.format(len(elems)))
            for i, e in enumerate(elems):
                f.write('{} '.format(i))
                # for v in elems.verts:
                for v in original_data['elements'][i]:
                    f.write('{} '.format(v))
                f.write('\n')

            # assert len(all_faces) == 1
            # faces = None
            # for k in all_faces:
            #     faces = all_faces[k]

            # f.write('nfaces:\n{}\n'.format( len(faces) ))
            # #for i, f in enumerate(faces):
            # #    f.write('{} '.format(i)

            f.write('npoints:\n{}\n'.format(len(coords)))
            for i, v in enumerate(coords):
                if ndims == 2:
                    f.write('{} {} {}\n'.format(i, v[0], v[1]))
                elif ndims == 3:
                    f.write('{} {} {} {}\n'.format(i, v[0], v[1], v[2]))

            f.write('connectivity:\n')
            for i, e in enumerate(elems):
                f.write('{} '.format(i))
                for j, (n, b) in enumerate(zip(e.neighbors, e.bndrys)):
                    if n < 0:
                        f.write('{} '.format(-b))
                    else:
                        f.write('{} '.format(n))
                f.write('\n')

            f.write('nbndry:\n{}\n'.format(len(bndryIdMap)-1))
            for k, v in bndryIdMap.items():
                if k:
                    f.write('{} {}\n'.format(v, k))

        return


    print('neighbor cell maps')
    face_map_vecs = create_neighbor_mappings(elems)
    print(len(face_map_vecs))


    if testid is None:
        
        maxlines = 25
        i = 0
        for fidx, (dijk, (ei, fi), (ej, fj), rorigin) in face_map_vecs.items():
            if i > maxlines:
                break
            i += 1
            print(fidx, dijk, ei, ej)
    else:
        
        def trilin(ix, px):
            assert all([x >= 0 and x <= 1.0 for x in ix])
            xi, et, zt = ix
            c00 = px[0] * (1.0 - xi) + px[1] * xi
            c10 = px[3] * (1.0 - xi) + px[2] * xi
            c01 = px[4] * (1.0 - xi) + px[5] * xi
            c11 = px[7] * (1.0 - xi) + px[6] * xi
            c0 = c00 * (1.0 - et) + c10 * et
            c1 = c01 * (1.0 - et) + c11 * et
            return c0 * (1.0 - zt) + c1 * zt
            
        print("reference cell")
        # print(face_map_vecs)
        nxi, net, nzt = 3, 3, 3 if ndims == 2 else 3
        
        xref = np.zeros((3,nxi,net,nzt), dtype='f')
        for i in range(nxi):
            for j in range(net):
                for k in range(nzt):
                    dxi = 1.0 / (nxi-1)
                    det = 1.0 / (net-1)
                    dzt = 1.0 / (nzt-1 if nzt > 1 else 1.0)
                    xi = float(i) * dxi
                    et = float(j) * det
                    zt = float(k) * dzt
                    xref[:,i,j,k] = xi, et, zt
                    # print(xref[:,i,j,k])

        maxiters = 5
        niters = 0
        for _, (dijk, (ei, fi), (ej, fj), rorigin) in face_map_vecs.items():
            if niters > maxiters:
                break
            niters += 1
            
            lxyz = coords[elems[ei].verts]
            rxyz = coords[elems[ej].verts]
            # print('before', lxyz)
            if ndims == 2:
                l = []
                r = []
                for i in range(4):
                    l.append([lxyz[i,0], lxyz[i,1], 0])
                    r.append([rxyz[i,0], rxyz[i,1], 0])
                    
                for i in range(4):
                    l.append([lxyz[i,0], lxyz[i,1], 1])
                    r.append([rxyz[i,0], rxyz[i,1], 1])

                lxyz = np.array(l)
                rxyz = np.array(r)
            print(ei, elems[ei].verts)
            print(ej, elems[ej].verts)
            # print(lxyz)
            # print(rxyz)
            lijk = np.zeros((3,nxi,net,nzt), dtype='f')
            rijk = np.zeros((3,nxi,net,nzt), dtype='f')
            for i in range(nxi):
                for j in range(net):
                    for k in range(nzt):
                        xi, et, zt = xref[:,i,j,k]
                        lijk[:,i,j,k] = [trilin(xref[:,i,j,k], lxyz[:,d]) for d in range(3)]
                        rijk[:,i,j,k] = [trilin(xref[:,i,j,k], rxyz[:,d]) for d in range(3)]
                        # print('l: ', xi, et, zt, lijk[:,i,j,k])
                        # print('r: ', xi, et, zt, rijk[:,i,j,k])
                    
            # print(lijk[:,-1,:])
            # print(rijk[:,-1,:])

            ldat = []
            rdat = []

            assert net == nxi
            ldat = np.full((3,net,nzt), dtype='f', fill_value=-1)
            rdat = np.full((3,net,nzt), dtype='f', fill_value=-1)

            if fi in (1,3): # E/W
                i = nxi-1 if fi == 1 else 0
                for k in range(nzt):
                    for j in range(net):
                        # ldat.append(lijk[:,i,j,k])
                        ldat[:,j,k] = lijk[:,i,j,k]
                        print(i, j, k, lijk[:,i,j,k])
            elif fi in (2,0): # N/Shh
                j = net-1 if fi == 2 else 0
                for k in range(nzt):
                    for i in range(nxi):
                        # ldat.append(lijk[:,i,j,k])
                        ldat[:,i,k] = lijk[:,i,j,k]
                        print(i, j, k, lijk[:,i,j,k])
            else:
                print('Not ready here')
                sys.exit(2)

            cmat = create_conn_matrix_from_vec(dijk)
            print(cmat)
            print(ldat.shape[1:])
            print(cmat.dot(np.array([1,net,nzt])))
            ls = np.array([nxi-1,0    ,0    ], dtype='i')
            le = np.array([nxi-1,net-1,nzt-1], dtype='i')
            rs = np.zeros((3), dtype='i')
            rs[0] = 0 if rorigin[0] == 0 else nxi-1
            rs[1] = 0 if rorigin[1] == 0 else net-1
            rs[2] = 0 if rorigin[2] == 0 else nzt-1
            re = rs + cmat.dot(le-ls)
            print('lorig: ', ls, le)
            print('rorig: ', rs, re)

            if True:
                for k in range(ls[2], le[2]+1):
                    for j in range(ls[1], le[1]+1):
                        for i in range(ls[0], le[0]+1):
                            di = i - ls[0]
                            dj = j - ls[1]
                            dk = k - ls[2]
                            dl = np.array([di,dj,dk])
                            dr = cmat.dot(dl)
                            print(i, j, k, rs + dr)

                
            if fj in (1,3): # E/W
                i = -1 if fj == 1 else 0
                dj = dijk[np.where(abs(dijk) == 2)[0][0]]
                dk = dijk[np.where(abs(dijk) == 3)[0][0]]
                assert dk == 3
                for k in range(nzt):
                    if dj == 2:
                        for j in range(net):
                            # rdat.append(rijk[:,i,j,k])
                            rdat[:,j,k] = rijk[:,i,j,k]
                            print(i, j, k, rijk[:,i,j,k])
                    elif dj == -2:
                        for j in range(net):
                            jj = (net-1) - j
                            # rdat.append(rijk[:,i,jj,k])
                            rdat[:,j,k] = rijk[:,i,jj,k]
                            print(i, jj, k, rijk[:,i,jj,k])
                    else:
                        raise
            elif fj in (0,2): # S/N
                j = -1 if fj == 2 else 0
                di = dijk[np.where(abs(dijk) == 1)[0][0]]
                dk = dijk[np.where(abs(dijk) == 3)[0][0]]
                assert dk == 3
                for k in range(nzt):
                    if di == 1:
                        for i in range(nxi):
                            # rdat.append(rijk[:,i,j,k])
                            rdat[:,i,k] = rijk[:,i,j,k]
                            print(i, j, k, rijk[:,i,j,k])
                    elif di == -1:
                         for i in range(nxi):
                            ii = (nxi-1) - i
                            # rdat.append(rijk[:,ii,j,k])
                            rdat[:,i,k] = rijk[:,ii,j,k]
                            print(ii, j, k, rijk[:,ii,j,k])
                    else:
                        raise
            else:
                raise

            # print(ldat[0,:,:])
            # print(rdat[0,:,:])
            # if not all([np.all(li == ri) for li, ri in zip(ldat, rdat)]):
            if not np.array_equal(ldat, rdat):
                print('Match failed')
                print('--left--')
                print(ldat)
                print('--right--')
                print(rdat)
                raise
            else:
                print("Passed")

    if ndims == 3:
        # Test if the 3d mesh can be cast into 2d easily.
        quads = all_faces[GeometricTypes.QUAD4]
        for j, e in enumerate(elems):
            f = quads[e.faces[4]]
            zs = coords[f,2]
            allzero = np.all(abs(zs) < 1e-12)
            #print(j, coords[f,2], allzero)
            if not allzero:
                print("Face 4 of element {} is not on the z=0 plane".format(j))
                sys.exit(1)



if __name__ == "__main__":
    main()
