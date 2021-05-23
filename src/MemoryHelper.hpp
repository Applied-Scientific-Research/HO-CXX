/*
 * MemoryHelper.h - Non-class methods for 2D-5D memory allocation
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

// generalized alloc/dealloc for 2D arrays in contiguous memory

template <class T>
T** allocate_array(const size_t nx, const size_t ny) {

    T** array = new T* [nx];
	array[0] = new T [nx*ny];

	for (size_t i=1; i<nx; ++i) {
		array[i] = array[0] + i * ny;
	}

	return array;
}

template <class T>
void free_array(T** array) {
	if (array == nullptr) return;
	delete[] array[0];
	delete[] array;
	return;
}

// generalized alloc/dealloc for 3D arrays in contiguous memory

template <class T>
T*** allocate_array(const size_t nx, const size_t ny, const size_t nz) {

    T*** array = new T** [nx];
	array[0] = new T* [nx*ny];
	array[0][0] = new T [nx*ny*nz];

	for (size_t i=1; i<nx; ++i) {
		array[i] = array[0] + i * ny;
	}

	for (size_t i=0; i<nx; ++i) {
		if (i != 0) array[i][0] = array[0][0] + i * ny * nz;
		for (size_t j=1; j<ny; ++j) {
			array[i][j] = array[i][0] + j * nz;
		}
	}

	return array;
}

template <class T>
void free_array(T*** array) {
	if (array == nullptr) return;
	delete[] array[0][0];
	delete[] array[0];
	delete[] array;
	return;
}

// generalized alloc/dealloc for 4D arrays in contiguous memory

template <class T>
T**** allocate_array(const size_t nx, const size_t ny, const size_t nz, const size_t nq) {

    T**** array = new T*** [nx];
	array[0] = new T** [nx*ny];
	array[0][0] = new T* [nx*ny*nz];
	array[0][0][0] = new T [nx*ny*nz*nq];

	for (size_t i=1; i<nx; ++i) {
		array[i] = array[0] + i * ny;
	}

	for (size_t i=0; i<nx; ++i) {
		if (i != 0) array[i][0] = array[0][0] + i * ny * nz;
		for (size_t j=1; j<ny; ++j) {
			array[i][j] = array[i][0] + j * nz;
		}
	}

	for (size_t i=0; i<nx; ++i) {
		if (i != 0) array[i][0][0] = array[0][0][0] + i * ny * nz * nq;
		for (size_t j=0; j<ny; ++j) {
			if (j != 0) array[i][j][0] = array[i][0][0] + j * nz * nq;
			for (size_t k=1; k<nz; ++k) {
				array[i][j][k] = array[i][j][0] + k * nq;
			}
		}
	}

	return array;
}

template <class T>
void free_array(T**** array) {
	if (array == nullptr) return;
	delete[] array[0][0][0];
	delete[] array[0][0];
	delete[] array[0];
	delete[] array;
	return;
}

