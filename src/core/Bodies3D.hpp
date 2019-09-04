#ifndef BODIES3D_HPP
#define BODIES3D_HPP
#include "../plugin/delaunaymassspringdamper/BodyDelaunayMassSpringDamper.hpp"
#include "../plugin/diffusion/BodyDiffusion.hpp"
#include "Plugins3D.hpp"

using namespace DelaunayMassSpringDamper;
using namespace Diffusion;

template <class cell_t>
class Bodies3D : public BodyDelaunayMassSpringDamper, public BodyDiffusion {

private:

public:
  Bodies3D(cell_t *, MecaCell::Vector3D pos = MecaCell::Vector3D::zero())
      : BodyDelaunayMassSpringDamper(pos), BodyDiffusion(3,nullptr) {}


  double getBaseRadius() const { return MecaCell::Config::DEFAULT_CELL_RADIUS; }

  double getBaseMass() const { return MecaCell::Config::DEFAULT_CELL_MASS; }

  using embedded_plugin_t = Plugins3D<cell_t>;
};
#endif
