#ifndef PLUGINBODYDELAUNAYMASSSPRINGDAMPER_HPP
#define PLUGINBODYDELAUNAYMASSSPRINGDAMPER_HPP

/*!
 * \file PluginDelaunayMassSpringDamper.hpp
 * \brief Delaunay triangulation & Mass spring based plugin
 *
 * contains the plugin class which works with the DelaunayMassSpringDamper body
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <vector>
#include <mecacell/mecacell.h>
#include "MassSpringDamper.hpp"

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_3<void *, K>;
using Tds =  CGAL::Triangulation_data_structure_3<Vb>;
//Use the Fast_location tag. Default or Compact_location works too.
using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
using Point = Delaunay::Point;

/*!
 * \namespace DelaunayMassSpringDamper
 */
namespace DelaunayMassSpringDamper {
    /*!
     * \class PluginDelaunayMassSpringDamper
     * @tparam cell_t
     * \brief Delaunay triangulation & mass spring based plugin class
     */
    template<typename cell_t>
    class PluginDelaunayMassSpringDamper {

    private:
        vector<MassSpringDamper<cell_t> > springList;  /*!< list of springs representing the connections between cells*/
        Delaunay triangulation;  /*!< Contains and compute the Delaunay triangulation of the cells*/
        bool areForcesNegligible;
        bool areMovementsNegligible;
        double physicsDt = 0.01;  /*!< physical dt is set to 0.01s supposing w->dt is in seconds*/
        double maxSpeed = 10.;  /*!< max velocity of a cell in µm/s*/
        double minForce = 1.;  /*!< minimal force to consider in pN*/
        double neighbourCoeff = 1.1;  /*!< coefficient used to set the maximal distance to consider 2 cells as neighbours*/
        double coherenceCoeff = 0.2;  /*!< coefficient used to to determine the spatial coherence of each cell*/


        /*!
         * \brief triangulation coherence checker
         * @return false if the triangulation has to be recomputed
         */
        bool isDelaunayCoherent(){
            bool coherence = true;

            if(triangulation.number_of_vertices()>0) {

                typename Delaunay::Finite_vertices_iterator vit;

                for (vit = triangulation.finite_vertices_begin(); vit != triangulation.finite_vertices_end(); ++vit) { //for each vertex
                    cell_t *c = static_cast<cell_t *>(vit->info());
                    MecaCell::Vec pos = MecaCell::Vec(vit->point().x(), vit->point().y(), vit->point().z());
                    if ((c->getBody().getPosition() - pos).length()
                    >= c->getBody().getBoundingBoxRadius() * coherenceCoeff) { //verify the coherence between the real position of a cell and its position in the triangulation
                        coherence = false;
                        break;
                    }
                }
            }
            else coherence = false;
            return coherence;
        }

        /*!
         * \brief triangulation points replacer
         *
         * moves the points which positions in the triangulation are not coherent anymore with their real positions
         */
        void moveDelaunay() {
            typename Delaunay::Finite_vertices_iterator vit;
            vit = triangulation.finite_vertices_begin();
            for (int i=0; i < triangulation.number_of_vertices(); ++i) { //for each vertex
                cell_t *c = static_cast<cell_t *>(vit->info());
                c->clearConnectedCells();
                MecaCell::Vec pos = MecaCell::Vec(vit->point().x(), vit->point().y(), vit->point().z());
                if ((c->getBody().getPosition() - pos).length()
                >= c->getBody().getBoundingBoxRadius() * coherenceCoeff){ //verify the coherence between the real position of a cell and its position in the triangulation
                    Point p(c->getBody().getPosition().x(), c->getBody().getPosition().y(), c->getBody().getPosition().z());
                    triangulation.move(vit, p);
                }
                ++vit;
            }
        }

        /*!
         * \brief triangulation calculator
         *
         * computes the triangulation from the list of cells
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void computeDelaunay(world_t *w) {
            for (cell_t *c : w->cells) c->clearConnectedCells();
            if (!isDelaunayCoherent()) {
                std::vector<Point> points;
                std::vector<void *> indices;
                for (cell_t *c : w->cells) {
                    c->clearConnectedCells();

                    points.push_back(Point(c->getBody().getPosition().x(), c->getBody().getPosition().y(),
                                           c->getBody().getPosition().z()));
                    indices.push_back(c);
                }

                triangulation.clear();
                triangulation = Delaunay(boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())),
                             boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));

            }
        }

        /*!
         * \brief springs generator
         *
         * creates the list of springs
         */
        void springCreation(){
            if (!triangulation.is_valid()) {
                std::cout << "UNVALID DELAUNAY !!" << std::endl;
            }

            springList.clear();
            typename Delaunay::Finite_edges_iterator eit;


            for (eit = triangulation.finite_edges_begin(); eit != triangulation.finite_edges_end(); ++eit) { //for each spring
                cell_t *c1 = static_cast<cell_t *>(eit->first->vertex(eit->second)->info());
                cell_t *c2 = static_cast<cell_t *>(eit->first->vertex(eit->third)->info());
                if (c1 != nullptr && c2 != nullptr) {
                    if ((c1->getBody().getPosition() - c2->getBody().getPosition()).length() <=
                    (c1->getBody().getBoundingBoxRadius() + c2->getBody().getBoundingBoxRadius()) * neighbourCoeff ){ //max distance to consider 2 neighbour cells
                        float_t r = 1.; //damping
                        float_t k = 10*(c1->getBody().getMass() + c2->getBody().getMass())/2; //stiffness proportional to the mass of the cells
                        float_t l = c1->getBody().getBoundingBoxRadius() * c1->getBody().getAdhesion() + c2->getBody().getBoundingBoxRadius() * c2->getBody().getAdhesion(); //rest length
                        float_t c = r * 2.0 * sqrt((c1->getBody().getMass() + c2->getBody().getMass()) * k); //damping coefficient

                        springList.push_back(MassSpringDamper<cell_t>(c1, c2, k, c, l));
                        c1->addConnectedCell(c2);
                        c2->addConnectedCell(c1);
                    }
                } else {
                    exit(41);
                }
            }
        }

        /*!
         * \brief new cells inserter
         *
         * add new cells to the triangulation
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void addCells(world_t *w){
            std::vector<Point> points;
            std::vector<void *> indices;
            for (cell_t *c : w->newCells) {
                points.push_back(Point(c->getBody().getPosition().x(), c->getBody().getPosition().y(), c->getBody().getPosition().z()));
                indices.push_back(c);
            }

            triangulation.insert(boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())),
                     boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));

        }

        /*!
         * \brief dead cells remover
         *
         * remove the dead cells from the triangulation
         */
        void removeCells(){
            typename Delaunay::Finite_vertices_iterator vit;

            for (vit = triangulation.finite_vertices_begin(); vit != triangulation.finite_vertices_end(); ++vit) { //for each vertex
                auto *c1 = static_cast<cell_t *>(vit->info());

                if(c1->isDead()){
                    triangulation.remove(vit);
                }
            }

        }

        /*!
         * \brief forces updater
         *
         * updates forces and checks if the forces are too weak to be considered
         */
        void updateForces() {
            areForcesNegligible = true;
            for (auto &msd : springList) {
                msd.computeForces();
                if(areForcesNegligible &&
                (msd.getConnection().first->getBody().getForce().length() > minForce || msd.getConnection().first->getBody().getForce().length() > minForce))
                    areForcesNegligible = false;
            }
        }

        /*!
         * \brief positions updater
         *
         * updates positions and checks if the movements are too weak to be considered
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void updatePositions(world_t *w) {
            areMovementsNegligible = true;
            for (cell_t *c : w->cells) {
                MecaCell::Vector3D vel = c->getBody().getVelocity() + c->getBody().getForce() * physicsDt / c->getBody().getMass();
                if(vel.length() < maxSpeed) c->getBody().setVelocity(vel);
                else c->getBody().setVelocity(maxSpeed * vel/vel.length());
                c->getBody().setPrevposition(c->getBody().getPosition());
                c->getBody().setPosition(c->getBody().getPosition() + c->getBody().getVelocity() * physicsDt);
                if(areMovementsNegligible && vel.length() > coherenceCoeff * c->getBody().getBoundingBoxRadius())
                    areMovementsNegligible = false;
            }
        }

        /*!
         * \brief forces resetter
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        inline void resetForces(world_t *w) {
            for (cell_t *c : w->cells) {
                c->getBody().resetForce();
            }
        }

        /*!
         * \brief physics applyer
         *
         * applies the physics with a coherent physic dt
         * stops when movements or foces are too weak
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void applyPhysics(world_t *w){
            double numberOfSteps = w->dt / physicsDt;
            for (unsigned i = 0; i < numberOfSteps; i++) {
                resetForces(w);
                updateForces();
                if(areForcesNegligible) break;
                updatePositions(w);
                if(areMovementsNegligible) break;
            }
        }

    public:
        /*!
         * \brief constructor
         *
         * initializes the springlist and the triangulation
         */
        inline PluginDelaunayMassSpringDamper()
        :springList(), triangulation()
        {}

        /*!
         * \brief physicsDt setter
         * @param dt
         */
        inline void setPhysicsDt(double dt){ physicsDt = dt; }

        /*!
         * \brief coherenceCoeff setter
         * @param c
         */
        inline void setCoherenceCoeff(double c){ coherenceCoeff = c; }

        /*!
         * \brief maxSpeed setter
         * @param s
         */
        inline void setMaxSpeed(double s){ maxSpeed = s; }

        /*!
         * \brief minForce setter
         * @param f
         */
        inline void setMinForce(double f){ minForce = f; }

        /*!
         * \brief neighbourCoeff setter
         * @param c
         */
        inline void setNeighbourCoeff(double c){ neighbourCoeff = c; }

        /*!
         * \brief endUpdate MecaCell hook
         *
         * recomputes or readjusts the triangulation and regenerates the springs
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void endUpdate(world_t *w){
            //computeDelaunay(w);
            moveDelaunay();
            springCreation();
        }

        /*!
         * \brief onAddCell MecaCell hook
         *
         * adds the new cells to the triangulation
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void onAddCell(world_t *w) {
            addCells(w);
        }

        /*!
         * \brief preBehaviorUpdate MecaCell hook
         *
         * applies the physics to each cell
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void preBehaviorUpdate(world_t *w) {
            applyPhysics(w);
        }

        /*!
         * \brief preDeleteDeadCellsUpdate MecaCell hook
         *
         * removes the dead cells from the triangulation
         *
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void preDeleteDeadCellsUpdate(world_t *w) {
                for (cell_t *c : w->cells) {
                    if (c->isDead()) {
                        removeCells();
                        break;
                    }
                }
        }

    };
}
#endif
