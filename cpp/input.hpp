


#include <dolfin.h>

using namespace dolfin;

class Source : public Expression
{
  public: 
  Source(): Expression(2) {}
  void eval(Array<double>& values, const Array<double>& x) const
  {

    values[0] = -((1.0)*( 4.0 * M_PI*M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]) - 2.0 * M_PI*M_PI *cos( M_PI * x[0]) * cos( M_PI * x[1]) ));;
    values[1] = -((1.0)*( 4.0 * M_PI*M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]) - 2.0 * M_PI*M_PI *cos( M_PI * x[0]) * cos( M_PI * x[1]) ));;
  }
};

class ZeroSource : public Expression
{
  public: 
  ZeroSource(): Expression(2) {}
  void eval(Array<double>& values, const Array<double>& x) const
  {

    values[0] = 0.0;
    values[1] = 0.0;
  }
};




// Boundary source for flux boundary condition
class BoundarySource : public Expression
{
public:

  BoundarySource(const Mesh& mesh) : Expression(2), mesh(mesh) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 1.0;
    values[1] = 1.0;
  }

private:

  const Mesh& mesh;

};

class BoundarySource1 : public Expression
{
public:

  BoundarySource1(const Mesh& mesh) : Expression(2), mesh(mesh) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 1.0;
    values[1] = 1.0;
  }

private:

  const Mesh& mesh;

};


// Sub domain for essential boundary condition
class EssentialBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS);
  }
};




