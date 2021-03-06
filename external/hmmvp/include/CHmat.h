/* Limited C interface to hmmvp matrix-vector product. */

#ifndef INCLUDE_HMMVP_C_HMAT
#define INCLUDE_HMMVP_C_HMAT

extern "C" {
  typedef void* CHmat;

  /*   All return values are 0 on success, -1 on failure.
       hm is an opaque pointer to the H-matrix. x and y must be allocated by the
     caller.
       nrhs is the maximum number of columns x and y will have.
       The matrix is m by n. */
  int chmat_init(const char* hmat_fn, const int nrhs, const int max_nthreads,
                 CHmat* hm);
  int chmat_get_info(const CHmat hm, int* m, int* n);
  int chmat_mvp(CHmat hm, const double* x, double* y, const int ncol);
  void chmat_cleanup(CHmat hm);
}

#endif
