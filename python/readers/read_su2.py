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

            
            
if __name__ == "__main__":
    f = 'C:/Users/chrst/Downloads/Step_Expansion.su2'
    if len(sys.argv) > 1:
        f = sys.argv[1]

    ndims, elems, coords, bndrys = load_mesh_su2(f)
    for i, e in enumerate(elems[:10]):
        print(i, e.verts, sum(coords[e.verts])/8)
    all_faces = create_face_list(elems)
    
    pairs = create_connectivity(elems, all_faces)

    
    for e in elems[:30]:
        print(e)


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

    quads = all_faces[GeometricTypes.QUAD4]
    print(quads.shape)
    for bnd in bndrys:
        print(bnd)
        bfs = bndrys[bnd]
        faces = []
        for bf in bfs:
            found, fidx = find_face_index(quads, bf.verts)
            assert found
            #print(bf.verts, fidx, quads[fidx])
            faces.append(fidx)
        bndrys[bnd] = np.array(faces, dtype='i')
        print(bndrys[bnd][:10])


    bndryIdMap = {None: 0}
    bndryIds = [None]
    for bnd in bndrys:
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
                for b in bndrys:
                    if f in bndrys[b]:
                        assert bnd is None
                        bnd = b
                if bnd is None:
                    print("no bnd face found for face {} in e {}".format(i, e))
                #print(e.eindex, i, f, n, bnd)
                ebnds[i] = bndryIdMap[bnd]
        e.bndrys = ebnds
        if ei < 10:
            print(e.eindex, e.neighbors, i, f, n, [bndryIds[k] for k in e.bndrys])

    # Test if the 3d mesh can be cast into 2d easily.
    for j, e in enumerate(elems):
        f = quads[e.faces[4]]
        zs = coords[f,2]
        allzero = np.all(abs(zs) < 1e-12)
        #print(j, coords[f,2], allzero)
        if not allzero:
            print("Face 4 of element {} is not on the z=0 plane".format(j))