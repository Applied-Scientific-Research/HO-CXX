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

//
// old routines
//

// allocator helpers
// there is a better way to do this, where all data is held in a contiguous memory chunk
// and we can even do it in a templated function: allocate_2d_array<double>(nx, ny);

double** allocate_2d_array_d(const size_t nx, const size_t ny) {
	double** outptr = nullptr;
	outptr = new double* [nx];
	for (size_t j = 0; j < nx; ++j) {
		outptr[j] = new double [ny];
	}
	return outptr;
}

double*** allocate_3d_array_d(const size_t nx, const size_t ny, const size_t nz) {
	double*** outptr = nullptr;
	outptr = new double** [nx];
	for (size_t j = 0; j < nx; ++j) {
		outptr[j] = new double* [ny];
		for (size_t i = 0; i < ny; ++i) {
			outptr[j][i] = new double[nz];
		}
	}
	return outptr;
}

double**** allocate_4d_array_d(const size_t nx, const size_t ny, const size_t nz, const size_t nq) {
	double**** outptr = nullptr;
	outptr = new double*** [nx];
	for (size_t i = 0; i < nx; ++i) {
		outptr[i] = new double** [ny];
		for (size_t j = 0; j < ny; ++j) {
			outptr[i][j] = new double* [nz];
			for (size_t k = 0; k < nz; ++k) {
				outptr[i][j][k] = new double[nq];
			}
		}
	}
	return outptr;
}

void deallocate_2d_array_d(double** arry, const size_t nx) {
	if (arry == nullptr) return;
	for (size_t k = 0; k < nx; ++k) {
		delete[] arry[k];
	}
	delete[] arry;
}

void deallocate_3d_array_d(double*** arry, const size_t nx, const size_t ny) {
	if (arry == nullptr) return;
	for (size_t k = 0; k < nx; ++k) {
		for (size_t j = 0; j < ny; ++j) {
			delete[] arry[k][j];
		}
		delete[] arry[k];
	}
	delete[] arry;
}

void deallocate_4d_array_d(double**** arry, const size_t nx, const size_t ny, const size_t nz) {
	if (arry == nullptr) return;
	for (size_t k = 0; k < nx; ++k) {
		for (size_t j = 0; j < ny; ++j) {
			for (size_t i = 0; i < nz; ++i) {
				delete[] arry[k][j][i];
			}
			delete[] arry[k][j];
		}
		delete[] arry[k];
	}
	delete[] arry;
}

