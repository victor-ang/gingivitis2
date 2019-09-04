#ifndef PLUGINS3D_HPP
#define PLUGINS3D_HPP
#include "../plugin/delaunaymassspringdamper/PluginDelaunayMassSpringDamper.hpp"
#include "../plugin/diffusion/PluginDiffusion.hpp"
#include <mecacell/genericconnectionplugin.hpp>
#include <mecacell/springconnection.hpp>

using namespace DelaunayMassSpringDamper;
using namespace Diffusion;

template <typename cell_t> struct Plugins3D {

  PluginDelaunayMassSpringDamper<cell_t> physicsPlugin;
  PluginDiffusion diffusionPlugin = PluginDiffusion(20., 0.01);
  Plugins3D() {}

  template <typename world_t> void onRegister(world_t *w) {
    w->registerPlugin(physicsPlugin);
    w->registerPlugin(diffusionPlugin); // commenter cette ligne pour debrancher la diffusion
  }

  template <typename world_t> void preBehaviorUpdate(world_t *w) {}
};

#endif
