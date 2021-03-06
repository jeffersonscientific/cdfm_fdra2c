#include "util/include/Elastostatics.hpp"

// Base class for all Elem-based Okada Green's functions for rectangular meshes
// having no mirror symmetry.
//todo Optionally periodic.
class OkadaRectGf : public ImplGreensFn {
public:
  virtual Hd* ComputeHd();
  virtual bool Call(const vector<UInt>& rs, const vector<UInt>& cs, double* B);
  virtual void DoExtraTasksSerial() throw (Exception);
  
protected:
  string _bc_filename;
  size_t _component;
  vector<util::es::dc3::Elem> _es, _bc_es;
  util::es::LameParms _lp;
  double _disl[3];
};

class OkadaRectTmGf : public OkadaRectGf {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
};

void OkadaRectGf::DoExtraTasksSerial() throw (Exception)
{
  // Compute the BC.
  Matrix<double> B(_es.size());
  double* pB = B.GetPtr();
  for (size_t ie = 0; ie < _es.size(); ie++) {
    pB[ie] = 0.0;
    for (size_t ib = 0; ib < _bc_es.size(); ib++)
      pB[ie] += util::es::GetTractionComp
        (_lp, _bc_es[ib], _disl, _es[ie], _component);
  }

  FILE* fid = fopen(_bc_filename.c_str(), "w");
  if (!fid) throw Exception("Can't write BC file.");
  write(pB, _es.size(), fid);
  fclose(fid);
}

Hd* OkadaRectGf::ComputeHd()
{
  size_t nes = _es.size();
  Matrix<double> ctrs(3, nes);
  double* pc = ctrs.GetPtr();
  for (size_t i = 0; i < nes; i++) {
    memcpy(pc, _es[i].Center(), 3*sizeof(double));
    pc += 3;
  }
  Hd* hd = NewHd(ctrs);
  errpr("nbr blocks = %ld\n", hd->NbrBlocks());
  return hd;
}

bool OkadaRectGf::
Call(const vector<UInt>& rs, const vector<UInt>& cs, double* B)
{
  for (size_t ic = 0, k = 0; ic < cs.size(); ic++)
    for (size_t ir = 0; ir < rs.size(); ir++, k++)
      B[k] = util::es::GetTractionComp
        (_lp, _es[cs[ic] - 1], _disl, _es[rs[ir] - 1], _component);
  return true;
}

void SetOkadaRectBcs
(const double x0, const double xf, const double eta0, const double etaf,
 const double y_min, const double depth_min, const double dipdeg,
 vector<util::es::dc3::Elem>* _bc_es)
{
  using namespace util::es;
  using namespace dc3;
  Matrix<double> xlim(2), etalim(2);
  xlim(1) = x0; xlim(2) = xf;
  etalim(1) = eta0; etalim(2) = etaf;
  vector<Elem> bc_es;
  double L = 100.0*std::max(xlim(2) - xlim(1), etalim(2) - etalim(1));
  // -x
  xlim(1) = x0 - L;
  xlim(2) = x0;
  PlanarTensorMeshToElems(xlim, etalim, y_min, depth_min, dipdeg, bc_es);
  // +x
  xlim(1) = xf;
  xlim(2) = xlim(1) + L;
  PlanarTensorMeshToElems(xlim, etalim, y_min, depth_min, dipdeg, bc_es);
  // -y or y
  xlim(1) = x0 - L;
  xlim(2) = xf + L;
  double dm, ym;
  double deta = etaf - eta0;
  if (dipdeg >= 0.0) {
    etalim(1) = eta0 - L;
    etalim(2) = eta0;
    dm = depth_min + deta*sind(dipdeg);
    ym = y_min - L*cosd(dipdeg);
  } else {
    etalim(1) = etaf;
    etalim(2) = etalim(1) + L;
    dm = depth_min - deta*sind(dipdeg);
    ym = y_min + deta*cosd(dipdeg);
  }
  PlanarTensorMeshToElems(xlim, etalim, ym, dm, dipdeg, bc_es);
  _bc_es->insert(_bc_es->end(), bc_es.begin(), bc_es.end());
}

void OkadaRectTmGf::Init(const KeyValueFile* kvf) throw (Exception)
{
  using namespace util::es;
  using namespace dc3;
  const string* s;
  const Matrix<double>* m;
  double d;

  Matrix<double> x, eta;
  if (!kvf->GetMatd("x", m)) throw Exception("Missing x.");
  x = *m;
  if (!kvf->GetMatd("eta", m)) throw Exception("Missing eta.");
  eta = *m;

  if (!kvf->GetString("bc_filename", s))
    throw Exception("Missing bc_filename.");
  _bc_filename = *s;

#define get(a) if (!kvf->GetDouble(#a, a)) throw Exception("Missing " #a);
  double mu, nu, dipdeg, y_min, depth_min, disl_dip, disl_strike, disl_tensile;
  get(mu); get(nu);
  get(dipdeg); get(depth_min); y_min = 0.0;
  get(disl_dip); get(disl_strike); get(disl_tensile);
#undef get

  _component = 1;
  if (kvf->GetDouble("component", d)) _component = (size_t) d;
  if (_component > 2) throw Exception("component must be 0, 1, or 2");

  _lp.Set(mu, nu);
  _disl[0] = disl_strike; _disl[1] = disl_dip; _disl[2] = disl_tensile;

  // Fault elements.
  PlanarTensorMeshToElems(x, eta, y_min, depth_min, dipdeg, _es);

  // Get these elements for the boundary conditions. The BCs are implemented by
  // huge slabs adjacent to the simulated fault.
  SetOkadaRectBcs(x(1), x(x.Size()), eta(1), eta(eta.Size()), y_min, depth_min,
                  dipdeg, &_bc_es);
}

// -----------------------------------------------------------------------------
// Quad-tree-based mesh
class OkadaRectQmGf : public OkadaRectGf {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
};

void OkadaRectQmGf::Init(const KeyValueFile* kvf) throw (Exception)
{
  using namespace util::es;
  using namespace dc3;
  const string* s;
  double d;

  if (!kvf->GetString("rmesh_filename", s))
    throw Exception("Missing rmesh_filename.");
  string rmesh_filename = *s + ".ser";

  if (!kvf->GetString("bc_filename", s))
    throw Exception("Missing bc_filename.");
  _bc_filename = *s;

#define get(a) if (!kvf->GetDouble(#a, a)) throw Exception("Missing " #a);
  double mu, nu, dipdeg, y_min, depth_min, disl_dip, disl_strike, disl_tensile;
  get(mu); get(nu);
  get(dipdeg); get(depth_min); y_min = 0.0;
  get(disl_dip); get(disl_strike); get(disl_tensile);
#undef get

  _component = 1;
  if (kvf->GetDouble("component", d)) _component = (size_t) d;
  if (_component > 2) throw Exception("component must be 0, 1, or 2");

  _lp.Set(mu, nu);
  _disl[0] = disl_strike; _disl[1] = disl_dip; _disl[2] = disl_tensile;

  rmesh::Rect r;
  { rmesh::RectMeshUnstruct* rm = rmesh::NewRectMeshUnstruct(rmesh_filename);
    r = rm->GetDomain();
    // Fault elements.
    RectMeshUnstructToElems(*rm, y_min, depth_min, dipdeg, _es);
    DeleteRectMeshUnstruct(rm); }

  // Get these elements for the boundary conditions. The BCs are implemented by
  // huge slabs adjacent to the simulated fault.
  Matrix<double> xlim(2), etalim(2);
  xlim(1) = r.x; xlim(2) = r.x + r.dx;
  etalim(1) = r.y; etalim(2) = r.y + r.dy;
  SetOkadaRectBcs(r.x, r.x + r.dx, r.y, r.y + r.dy, y_min, depth_min, dipdeg,
                  &_bc_es);
}

// -----------------------------------------------------------------------------
// Base class for all Elem-based, possibly along-strike periodic, Okada Green's
// functions for rectangular meshes having mirror symmetry across the strike = 0
// line.
class OssrpGf : public ImplGreensFn {
public:
  virtual Hd* ComputeHd();
  virtual void DoExtraTasksSerial() throw (Exception);
  virtual bool Call(const vector<UInt>& rs, const vector<UInt>& cs,
                    double* B);

protected:
  size_t _component;
  string _bc_filename;
  vector<util::es::dc3::Elem> _es, _bc_es;
  double _L; // Along-strike length
  util::es::LameParms _lp;
  double _disl[3];
  size_t _nrepeat; // Number of repeated images in each of -/+x directions
};

class OssrptmGf : public OssrpGf {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  
private:
  size_t _nx, _neta;
};

// Duplicate with -x for the other half of the fault.
static void CreateOtherHalf(vector<util::es::dc3::Elem>& es)
{
  using namespace util::es::dc3;
  for (size_t i = 0, n = es.size(); i < n; i++) {
    // Specify using the center so that we can flip the sign of the x component.
    const Elem& e = es[i];
    double dx = 0.5*(e.al1() + e.al2());
    double de = 0.5*(e.aw1() + e.aw2());
    es.push_back(Elem(-e.Center()[2], e.dipdeg(), dx, dx, de, de,
                      -e.Center()[0], e.Center()[1]));
  }
}

void OssrptmGf::Init(const KeyValueFile* kvf) throw (Exception)
{
  using namespace util::es;
  using namespace dc3;
  const string* s;
  const Matrix<double>* m;
  
  Matrix<double> x, eta;
  if (!kvf->GetMatd("x", m)) throw Exception("Missing x.");
  x = *m;
  if (!kvf->GetMatd("eta", m)) throw Exception("Missing eta.");
  eta = *m;

  _nx = 2*(x.Size() - 1);
  _neta = eta.Size() - 1;
  _L = 2*x(x.Size());

  if (!kvf->GetString("bc_filename", s))
    throw Exception("Missing bc_filename.");
  _bc_filename = *s;

#define get(a) if (!kvf->GetDouble(#a, a)) throw Exception("Missing " #a);
  double mu, nu, dipdeg, y_min, depth_min, nrepeat;
  get(mu); get(nu);
  get(dipdeg); get(depth_min); get(nrepeat); y_min = 0.0;
#undef get

  if (nrepeat < 0.0 || nrepeat > 100.0)
    throw Exception("nrepeat must be between 0 and 100.");
  _nrepeat = (size_t) nrepeat;
  _component = 1;
  _disl[0] = 0.0; _disl[1] = 1.0; _disl[2] = 0.0;
  _lp.Set(mu, nu);

  // Fault elements.
  PlanarTensorMeshToElems(x, eta, y_min, depth_min, dipdeg, _es);
  CreateOtherHalf(_es);

  // Get these elements for the boundary conditions. The BCs are implemented by
  // huge slabs adjacent to the simulated fault.
  Matrix<double> xlim(2), etalim(2);
  double L = std::max(100.0, 10*(1 + nrepeat)) *
    std::max(2*x(x.Size()), eta(eta.Size()) - eta(1));
  // -y or y
  xlim(2) = x(x.Size()) + L;
  xlim(1) = -xlim(2);
  double dm, ym;
  if (dipdeg >= 0.0) {
    etalim(1) = eta(1) - L;
    etalim(2) = eta(1);
    dm = depth_min + (eta(eta.Size()) - eta(1))*sind(dipdeg);
    ym = y_min - L*cosd(dipdeg);
  } else {
    etalim(1) = eta(eta.Size());
    etalim(2) = etalim(1) + L;
    dm = depth_min - (eta(eta.Size()) - eta(1))*sind(dipdeg);
    ym = y_min + (eta(eta.Size()) - eta(1))*cosd(dipdeg);
  }
  PlanarTensorMeshToElems(xlim, etalim, ym, dm, dipdeg, _bc_es);
}

void OssrpGf::DoExtraTasksSerial() throw (Exception)
{
  FILE* fid = fopen(_bc_filename.c_str(), "w");
  if (!fid) throw Exception("OssrpGf: Can't write BC file " + _bc_filename);
  // Compute the BC on the range elements.
  size_t n = _es.size() / 2;
  Matrix<double> B(n);
  double* pB = B.GetPtr();
  for (size_t ie = 0; ie < n; ie++) {
    pB[ie] = 0.0;
    for (size_t ib = 0; ib < _bc_es.size(); ib++)
      pB[ie] += util::es::GetTractionComp
        (_lp, _bc_es[ib], _disl, _es[ie], _component);
  }
  write(pB, n, fid);
  fclose(fid);
}

Hd* OssrpGf::ComputeHd()
{
  size_t nes = _es.size();
  Matrix<double> D(3, nes);   // Domain is the full fault.
  Matrix<double> R(3, nes/2); // Range is the +x side of the fault.

  double* pD = D.GetPtr();
  double* pR = R.GetPtr();
  for (size_t i = 0; i < nes; i++) {
    memcpy(pD, _es[i].Center(), 3*sizeof(double));
    pD += 3;
    if (i < nes/2) {
      memcpy(pR, _es[i].Center(), 3*sizeof(double));
      pR += 3;
    }
  }

  Matrix<double>* ppb = NULL;
  Matrix<double> pb(2, 3);
  if (_nrepeat > 0) {
    // Handle periodicity in x direction.
    ppb = &pb;
    pb.Zero();
    pb(1, 1) = -0.5*_L;
    pb(2, 1) =  0.5*_L;
  }

  Hd* hd = NewHd(D, R, ppb);
  errpr("nbr blocks = %ld\n", hd->NbrBlocks());
  return hd;
}

bool OssrpGf::
Call(const vector<UInt>& rs, const vector<UInt>& cs, double* B)
{
  using namespace util::es;
  using namespace dc3;
  for (size_t ic = 0, k = 0; ic < cs.size(); ic++)
    for (size_t ir = 0; ir < rs.size(); ir++, k++) {
      const Elem& es = _es[cs[ic] - 1];
      B[k] = GetTractionComp(_lp, es, _disl, _es[rs[ir] - 1], _component);
      // Periodic images
      for (size_t ip = 1; ip < _nrepeat + 1; ip++) {
        Elem esr1(es.depth(), es.dipdeg(), es.al1(), es.al2(), es.aw1(),
                  es.aw2(), es.gx() + ip*_L, es.gy());
        B[k] += GetTractionComp(_lp, esr1, _disl, _es[rs[ir] - 1], _component);
        Elem esr2(es.depth(), es.dipdeg(), es.al1(), es.al2(), es.aw1(),
                  es.aw2(), es.gx() - ip*_L, es.gy());
        B[k] += GetTractionComp(_lp, esr2, _disl, _es[rs[ir] - 1], _component);
      }
    }
  return true;
}

// -----------------------------------------------------------------------------
class OssrpqmGf : public OssrpGf {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
};

void OssrpqmGf::Init(const KeyValueFile* kvf) throw (Exception)
{
  using namespace util::es;
  using namespace dc3;
  const string* s;

  if (!kvf->GetString("rmesh_filename", s))
    throw Exception("Missing rmesh_filename.");
  string rmesh_filename = *s + ".ser";

  if (!kvf->GetString("bc_filename", s))
    throw Exception("Missing bc_filename.");
  _bc_filename = *s;

#define get(a) if (!kvf->GetDouble(#a, a)) throw Exception("Missing " #a);
  double mu, nu, dipdeg, y_min, depth_min, nrepeat;
  get(mu); get(nu);
  get(dipdeg); get(depth_min); get(nrepeat); y_min = 0.0;
#undef get

  if (nrepeat < 0.0 || nrepeat > 100.0)
    throw Exception("nrepeat must be between 0 and 100.");
  _nrepeat = (size_t) nrepeat;
  _component = 1;
  _disl[0] = 0.0; _disl[1] = 1.0; _disl[2] = 0.0;
  _lp.Set(mu, nu);

  rmesh::Rect r;
  { rmesh::RectMeshUnstruct* rm = rmesh::NewRectMeshUnstruct(rmesh_filename);
    r = rm->GetDomain();
    // Make sure the user understands the setup.
    if (r.x != 0.0) throw Exception("r.x should be 0 in the rmesh.");
    _L = 2.0*r.dx;

    // Fault elements.
    RectMeshUnstructToElems(*rm, y_min, depth_min, dipdeg, _es);
    DeleteRectMeshUnstruct(rm); }
  CreateOtherHalf(_es);

  // Get these elements for the boundary conditions. The BCs are implemented by
  // huge slabs adjacent to the simulated fault.
  Matrix<double> xlim(2), etalim(2);
  double L = 100.0*(1 + _nrepeat)*std::max(2.0*r.dx, r.dy);
  // -y or y
  xlim(2) = r.dx + L;
  xlim(1) = -xlim(2);
  double dm, ym;
  if (dipdeg >= 0.0) {
    etalim(1) = r.y - L;
    etalim(2) = r.y;
    dm = depth_min + r.dy*sind(dipdeg);
    ym = y_min - L*cosd(dipdeg);
  } else {
    etalim(1) = r.y + r.dy;
    etalim(2) = etalim(1) + L;
    dm = depth_min - r.dy*sind(dipdeg);
    ym = y_min + r.dy*cosd(dipdeg);
  }
  PlanarTensorMeshToElems(xlim, etalim, ym, dm, dipdeg, _bc_es);
}
