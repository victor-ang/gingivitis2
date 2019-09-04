#ifndef ONKO3D_3_0_PLUGINDIFFUSION_HPP
#define ONKO3D_3_0_PLUGINDIFFUSION_HPP


/*!
 * \file PluginDiffusion.hpp
 * \brief grid diffusion based plugin
 *
 * contains the grid diffusion based plugin class
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/intersections.h>
#include <vector>
#include <mecacell/mecacell.h>
#include <math.h>
#include "DiffusionGrid.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;

/*!
 * \namespace Diffusion
 */
namespace Diffusion {

    /*!
     * \class PluginDiffusion
     * \brief grid diffusion based plugin class
     */
    class PluginDiffusion {

    private:
        MecaCell::Vec centroid; /*!< center of gravity*/
        double ro = 0.; /*!< radius of the spheroid*/
        DiffusionGrid grid;

        typedef std::vector<double> Data;

        /*!
        * \brief gaussian probability density function
        *
        * µ = 0 and sigma = 1
        *
        * @param x
        * @return probability density in x
        */
        inline static double Kernel(const double &x) { return 1.f / sqrt(2.f * M_PI) * exp(pow(x, 2) / -2.f); }

        /*!
         * \brief centroid calculator and getter
         * @tparam world_t
         * @param w
         * @return centroid
         */
        template<typename world_t>
        MecaCell::Vec computeCentroid(world_t *w) {
            centroid = MecaCell::Vec(0, 0, 0);
            for (auto &c : w->cells) {
                centroid += c->getBody().getPosition();
            }
            return centroid /= w->cells.size();
        }

    public:

        /*!
         * \brief default constructor
         */
        PluginDiffusion() = default;

        /*!
         * \brief constructor
         * @param dx
         * @param accuracy
         */
        explicit inline PluginDiffusion(double dx, double accuracy): centroid(0.,0.,0.), grid(dx,accuracy) {}

        inline DiffusionGrid *getGrid(){ return &grid; }

        /*!
         * \brief centroid getter
         * @return centroid
         */
        inline MecaCell::Vec getCentroid() const { return centroid; }

        /*!
         * \brief spheroid radius calculator and getter
         *
         * computes and returns ro
         *
         * @tparam world_t
         * @param w
         * @return ro
         */
        template<typename world_t>
        double getSpheroidRadius(world_t *w) {
            computeCentroid(w);
            double r = 0;
            Data v;
            std::vector<Point_3> points;
            double max = -1;
            for (auto &c : w->cells) {
                r = (c->getBody().getPosition() - centroid).length();
                v.push_back(r);
                points.push_back(Point_3(c->getBody().getPosition().x(), c->getBody().getPosition().y(), c->getBody().getPosition().z()));
                if (r > max) {
                    max = r;
                }
            }

            double maxDensity = -1;
            double xMaxDensity = -1;
            double val = 0;
            double lambda = 10.f;
            double sum = 0;
            for (double i = 0.; i <= max; i += 0.2) {
                sum = 0;
                for (size_t j = 0; j < v.size(); ++j) {
                    sum += Kernel((i - v[j]) / lambda);
                }
                val = sum / lambda;
                if (val > maxDensity) {
                    maxDensity = val;
                    xMaxDensity = i;
                }
            }
            ro = (max - (max - xMaxDensity) / 2);
            return ro;
        }

        /*!
         * \brief molecule adder
         * @param m : new molecule to be added
         */
        inline void addMolecule(Molecule m) { grid.addMolecule(m); }


        /*!
         * \brief preBehaviorUpdate MecaCell hook
         *
         * computes and updates the oxygen quantity for each cell
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void preBehaviorUpdate(world_t *w) {
            grid.insertCell(w);
            grid.computeMolecules(w);
            for (auto &c : w->cells) {
                for(unsigned i = 0; i < grid.getMolecules().size(); ++i) {
                    c->getBody().setQuantity(i,grid.getMoleculeRealPos(c->getBody().getPosition(), i));
                }
            }
        }
    };
}

#endif //ONKO3D_3_0_PLUGINDIFFUSION_HPP
