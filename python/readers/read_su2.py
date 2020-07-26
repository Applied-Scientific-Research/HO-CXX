# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 09:44:20 2020

@author: chrst
"""

import sys
import numpy as np
import enum
# import json

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
            if etype == 12: # HEX8
                assert len(items) == 10
                v = np.array([int(c) for c in items[1:-1]], dtype='i')
                elems.append(Element(GeometricTypes.HEX8, v, eidx=i))
            else:
                print("SU2: unknown element type {}".format(etype))
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
            if etype == 9: # QUAD4
                assert len(items) == 5
                v = np.array([int(c) for c in items[1:]], dtype='i')
                bnd.append(Element(GeometricTypes.QUAD4, v))
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
                  [4,5,6,7]], dtype='i')
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
    Given a list of elements (3d), form the unique list of faces.
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
    
    all_faces = {}
    for el in elems:
        faces = None
        ftype = None
        if el.etype == GeometricTypes.HEX8:
            ftype = GeometricTypes.QUAD4
            if ftype in all_faces:
                faces = all_faces[ftype]
            else:
                all_faces[ftype] = faces = []
        else:
            raise NameError('Only Hexs supported so far.')
                
        efaces = el.verts[FaceElementMaps[el.etype]]
        nverts = nVerts(ftype)
        for i, f in enumerate(efaces.tolist()):
            fidx, match = find_face_in_list(faces, f, nverts)
            el.faces.append(fidx)
            # print('face: ', i, f, fidx, len(elems)*6, match)
            
    for ftype in all_faces:
        faces = all_faces[ftype]
        n = len(faces)
        nverts = nVerts(ftype)
        new_faces = np.array(faces, dtype='i').reshape(n//nverts,nverts)
        all_faces[ftype] = new_faces
        print(ftype, len(faces), n, nverts, new_faces.shape)
        
    return all_faces


def create_connectivity(elements, all_faces):
    
    # Assuming everything is a HEX elem and QUAD face.
    assert len(all_faces) == 1
    for k in all_faces.keys():
        assert k == GeometricTypes.QUAD4

    faces = all_faces[GeometricTypes.QUAD4]

    # Build a full list of all face indices across the elements
    # and turn into a 2d np array for quick searching.

    nelems = len(elements)
    nfaces = len(faces)
    
    efaces = np.empty((nelems,6), dtype='i')
    for i, e in enumerate(elements):
        assert e.etype == GeometricTypes.HEX8
        e.faces = np.array(e.faces, dtype='i')
        efaces[i,:] = e.faces[:]
        e.neighbors = np.full((6), -1)
        e.isRight = np.full((6), 0, dtype='i')

    
    # This will fail if nfaces isn't uniform.
    print('efaces: ', efaces.shape, nfaces)
    
    # Populate an array with left/right element neighbor for each face.
    pairs = np.zeros((nfaces,2), dtype='i')
    for i in range(nfaces):
        els = np.where(efaces == i)[0]
        assert els.shape[0] == 1 or els.shape[0] == 2
        pairs[i,0] = els[0]
        pairs[i,1] = -1 if els.shape[0] == 1 else els[1]
    
    print(pairs[:10,:])    
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

    return pairs


def set_boundary_faces(elems, faces, bndry_sections):
    
    """
    Given an ndarray (n,4) of quads, find the index of the quad.
    """
    def find_face_index(faces, f):
        assert faces.shape[1] == f.shape[0]

        rows = np.where(faces == f[0])[0]
        for i in rows:
            row = faces[i]
            if all([j in row for j in f[1:]]):
                return True, i

        return False, -1

    quads = faces[GeometricTypes.QUAD4]
    print(quads.shape)
    for bnd in bndry_sections:
        print(bnd)
        bfs = bndry_sections[bnd]
        bfaces = []
        for bf in bfs:
            found, fidx = find_face_index(quads, bf.verts)
            assert found
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
        ebnds = np.full((6), bndryIdMap[None], dtype='i')
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
            print(e.eindex, e.neighbors, i, f, n, [bndryIds[k] for k in e.bndrys])

    return


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
    normal_fmaps[S,N] = [ 0, 2, 0]
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
        print('el: {}'.format(ei))
        for fi, (n, r) in enumerate(zip(e.neighbors, e.isRight)):
            if n != -1:
                # if a shared face. Find the left/right local face indices.
                fidx = e.faces[fi]
                fj = np_find(elements[n].faces == fidx)
                if not (fi,fj) in normal_fmaps:
                    print('faces {} {} not in normal_fmaps {} {}'.format(fi, fj, Orient[fi], Orient[fj]))
                    sys.exit(2)
                print(" face: {}, {} neigh: {}, {}, {} dir: {} => {}".format(fi, fidx, n, r, fj, Orient[fi], Orient[fj]))
                # dijk = np.zeros((3), dtype='i')
                dijk = np.array(normal_fmaps[fi,fj], dtype='i')
                if fi in (E,W):
                    i = 0 if fi == W else 1
                    # dijk[0] = normal_fmaps[fi,fj]
                    
                    lv1 = ijk2v[i,0,0]
                    rv1 = np_find(elements[n].verts == e.verts[lv1]) # local id
                    # +j
                    lv2 = ijk2v[i,1,0]
                    rv2 = np_find(elements[n].verts == e.verts[lv2])
                    dj = ivec(*v2ijk[rv2]) - ivec(*v2ijk[rv1])
                    idj = np_find(dj != 0)
                    dijk[1] = (idj+1) if dj[idj] > 0 else -(idj+1)
                    # +k
                    lv3 = ijk2v[i,0,1]
                    rv3 = np_find(elements[n].verts == e.verts[lv3])
                    dk = ivec(*v2ijk[rv3]) - ivec(*v2ijk[rv1])
                    idk = np_find(dk != 0)
                    dijk[2] = (idk+1) if dk[idk] > 0 else -(idk+1)

                    # print("{} => {}".format(Orient[fi], Orient[fj]))
                    print(lv1, lv2, lv3, rv1, rv2, rv3, e.verts[lv1], e.verts[lv2], e.verts[lv3])
                    print('j+ ', dj, idj, dj[idj])
                    print('k+ ', dk, idk, dk[idk])

                elif fi in (S,N):
                    j = 0 if fi == S else 1
                    # dijk[1] = normal_fmaps[fi,fj]
                    
                    lv1 = ijk2v[0,j,0]
                    rv1 = np_find(elements[n].verts == e.verts[lv1]) # local id
                    # +i
                    lv2 = ijk2v[1,j,0]
                    rv2 = np_find(elements[n].verts == e.verts[lv2])
                    di = ivec(*v2ijk[rv2]) - ivec(*v2ijk[rv1])
                    idi = np_find(di != 0)
                    dijk[0] = (idi+1) if di[idi] > 0 else -(idi+1)
                    # +k
                    lv3 = ijk2v[0,j,1]
                    rv3 = np_find(elements[n].verts == e.verts[lv3])
                    dk = ivec(*v2ijk[rv3]) - ivec(*v2ijk[rv1])
                    idk = np_find(dk != 0)
                    dijk[2] = (idk+1) if dk[idk] > 0 else -(idk+1)

                    # print("N/S")
                    # print(lv1, rv1, lv2, rv2, dj, idj, dj[idj])
                    # print(lv1, rv1, lv3, rv3, dk, idk, dk[idk])
                    print(lv1, lv2, lv3, rv1, rv2, rv3, e.verts[lv1], e.verts[lv2], e.verts[lv3])
                    print('i+ ', di, idi, di[idi])
                    print('k+ ', dk, idk, dk[idk])
                    
                else:
                    raise 'Not implemented yet'
                    
                print(dijk)
                print(create_conn_matrix_from_vec(dijk))
                if fidx not in face_map_vec:
                    el, er = ((ei, fi), (n, fj)) if r == 0 else ((n, fj), (ei, fi))
                    face_map_vec[fidx] = dijk, el, er
                
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


    if left == [1,2,3]:
        pass
    else:
        print("left {} not supported yet".format(left))
        
    if right == [1,2,3]:
        pass
    else:
        rotmat = create_conn_matrix_from_vec(right)
 
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

        print(rotmat)
        print(v0)
        print(v1)
        print(offset)
        v1 -= offset
        print(v1)
        r = [ijk2v[tuple(v)] for v in v1]
        print(r)
        v = elems[1].verts[:]
        print(v)
        elems[1].verts = v[r] 
        print(elems[1].verts)

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

            
if __name__ == "__main__":

    import getopt

    infile = 'Step_Expansion.su2'
    testid = None

    if len(sys.argv) > 1:

        try:
            opts, args = getopt.getopt(sys.argv[1:], 'f:t:', ['file=', 'test='])
        except getopt.GetoptError as err:
            print('getopt error {}'.format(str(err)))
            sys.exit(2)

        print(opts)

        for opt, arg in opts:
            if opt in ('-f', '--file'):
                infile = arg
            elif opt in ('-t', '--test'):
                testid = arg
            else:
                print('Unknown option {}'.format(opt))
                sys.exit(1)

    if testid is None:
        ndims, elems, coords, bndry_sections = load_mesh_su2(infile)
    else:
        ndims, elems, coords, bndry_sections = test_mesh_hex8(testid)


    print('midpoints')
    for i, e in enumerate(elems[:10]):
        print(i, e.verts, sum(coords[e.verts])/8)
    all_faces = create_face_list(elems)

    pairs = create_connectivity(elems, all_faces)

    
    # for e in elems[:30]:
    #     print(e)

    print('boundary faces')
    set_boundary_faces(elems, all_faces, bndry_sections)

    print('neighbor cell maps')
    face_map_vecs = create_neighbor_mappings(elems)

    if testid is not None:
        
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
        print(face_map_vecs)
        nk = 3
        xref = np.zeros((3,nk,nk), dtype='f')
        for i in range(nk):
            for j in range(nk):
                xi = float(i) / (nk-1)
                et = float(j) / (nk-1)
                # zt = float(k) / (nk-1)
                zt = 0.0
                xref[:,i,j] = xi, et, 0.0
                print(xref[:,i,j])
                
        for _, (dijk, (ei, fi), (ej, fj)) in face_map_vecs.items():
            
            lxyz = coords[elems[ei].verts]
            rxyz = coords[elems[ej].verts]
            print(ei, elems[ei].verts)
            print(ej, elems[ej].verts)
            lijk = np.zeros((3,nk,nk), dtype='f')
            rijk = np.zeros((3,nk,nk), dtype='f')
            for i in range(nk):
                for j in range(nk):
                    xi, et, zt = xref[:,i,j]
                    lijk[:,i,j] = [trilin(xref[:,i,j], lxyz[:,d]) for d in range(3)]
                    rijk[:,i,j] = [trilin(xref[:,i,j], rxyz[:,d]) for d in range(3)]
                    print('l: ', xi, et, lijk[:,i,j])
                    print('r: ', xi, et, rijk[:,i,j])
                    
            # print(lijk[:,-1,:])
            # print(rijk[:,-1,:])
            
            if fi == 1: # 'E'
                for j in range(nk):
                    print(j, lijk[:,-1,j])
            else:
                print('Not ready here')
                sys.exit(2)
                
            if fj in (1,3): # E/W
                i = -1 if fj == 1 else 0
                dj = dijk[np.where(abs(dijk) == 2)[0][0]]
                if dj == 2:
                    for j in range(nk):
                        print(j, rijk[:,i,j])
                elif dijk[1] == -2:
                    for j in range(nk):
                        jj = (nk-1) - j
                        print(jj, rijk[:,i,jj])
                else:
                    raise
            elif fj in (0,2): # S/N
                j = -1 if fj == 2 else 0
                di = dijk[np.where(abs(dijk) == 1)[0][0]]
                if di == 1:
                    for i in range(nk):
                        print(i, rijk[:,i,j])
                elif di == -1:
                    for i in range(nk):
                        ii = (nk-1) - i
                        print(ii, rijk[:,ii,j])
            else:
                raise
                sys.exit(2)

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