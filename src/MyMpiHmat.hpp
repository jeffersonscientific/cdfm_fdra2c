/* Routines to parallelize Hmat::Mvp(). We assume that MPI_Init has already been
   called.
     Basic usage:

       MPI_Init(&argc, &argv);
       Hm::MpiHmat<float>* hm;
       try {
         hm = new Hm::MpiHmat<float>(hmat_fn, nrhs);
       } catch (const Exception& e) {
         cout << e.GetMsg() << endl;
         exit(-1);
       }
       hm->Mvp(x, y, nrhs);
       delete hm;
       MPI_Finalize();

     Compile with mpic++. Run with a command of the sort (in Unix):
       mpirun -np 4 ./a.out
*/

/* TODO: In Mvp(x, y, ncol), one could imagine having a class for x and y that
   puts different pieces on each proc automatically. Then all of x and y
   wouldn't end up on every proc.  */

#ifndef INCLUDE_HMMVP_MPI_HMAT
#define INCLUDE_HMMVP_MPI_HMAT

#include "util/include/Defs.hpp"
#include "util/include/Exception.hpp"
#include "util/include/Mpi.hpp"

namespace hmmvp {
  class Hmat;
}

namespace hmmvp {
  using namespace std;
  using namespace util;

  typedef long long Blint;
  
  template<typename real>
  class MpiHmat {
  public:
    MpiHmat(const string& filename, UInt ncol = 1)
      throw (FileException, Exception);
    ~MpiHmat();

    // Number of rows.
    Blint GetM() const;
    // Number of columns.
    Blint GetN() const;
    // Number of nonzeros in the H-matrix approximation.
    Blint GetNnz() const;

    // y = B*x, where x has ncol columns. y is allocated by the caller.
    //   x and y are the full vectors.
    //   If perms are on, then x needs to be non-NULL only on the root cpu.
    //   y must be non-Null and have the size of the full problem. On exit, its
    // values are valid only on the root cpu.
    void Mvp(const real* x, real* y, UInt ncol) const
      throw (Exception);

    void TurnOnPermute()  { do_perms = true;  }
    void TurnOffPermute() { do_perms = false; }

    const Blint* GetQ() const;
    const Blint* GetP() const;

  private:
    hmmvp::Hmat* hm;
    mutable vector<real> workr;
    bool am_root;
    Blint nnz;
    bool do_perms;
  };

}

#include "MyMpiHmat_inl.hpp"

#endif
