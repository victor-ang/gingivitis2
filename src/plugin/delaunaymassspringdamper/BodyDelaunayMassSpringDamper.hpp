#ifndef BODYMASSSPRINGDAMPER_HPP
#define BODYMASSSPRINGDAMPER_HPP

/*!
 * \file BodyDeLaunayMassSpringDamper.hpp
 * \brief Delaunay triangulation & Mass spring based body
 *
 * contains the body class which works with the DelaunayMassSpringDamper plugin
 */

#include <mecacell/mecacell.h>
#include <mecacell/movable.h>
#include "PluginDelaunayMassSpringDamper.hpp"
#include <math.h>

/*!
 * \namespace DelaunayMassSpringDamper
 */
namespace DelaunayMassSpringDamper {

    /*!
     * \class BodyDelaunayMassSpringDamper
     * \brief Delaunay triangulation & mass spring based body class
     */
    class BodyDelaunayMassSpringDamper : public MecaCell::Movable {

    private:
        double radius; /*!< radius in µm*/
        double density; /*!< density in ng/µm^3*/
        double adhesion; /*!< adhesion coefficient*/

    public:
        /*!
         * \brief Constructor
         *
         * initializes the radius to 40µm and the density to the water density
         *
         * @param pos
         */
        explicit BodyDelaunayMassSpringDamper(const MecaCell::Vector3D &pos = MecaCell::Vector3D::zero())
        {
            // Log-normal distribution for cell size
            std::random_device rd;
            std::mt19937 gen(rd());
            std::lognormal_distribution<> d(0,1);
            float variation = d(gen);
            // radius = 35 + variation;
            if (variation <=5){
                radius = 35+ variation;
            }
            else
            {
                radius = 40;
            }
            //radius = 40; //40 µm radius
            density = 0.001; //in ng/µm^3 (density of water ~ 1000 kg/m^3)
            adhesion = 0.75;
            this->setMass(density *4* M_PI/3 * pow(radius,3)); //mass = density * 4pi/3 * rayon^3 ng // Ajout du *4??
            this->setPosition(pos + MecaCell::Vec::randomUnit() * getBoundingBoxRadius() * 2.0 * adhesion);
        }

        /*!
         * \brief radius getter
         * @return radius
         */
        inline double getBoundingBoxRadius() const { return radius; }

        /*!
         * \brief radius setter
         *
         * changes the mass knowing the current density
         *
         * @param rad
         */
        inline void setRadius(double rad) {
            radius = rad;
            this->setMass(density *4* M_PI/3 * pow(radius,3)); //mass = density * 4pi/3 * rayon^3 ng // Ajout du *4??
        }

        /*!
         * \brief density setter
         *
         * changes the mass knowing the current radius
         *
         * @param d
         */
        inline void setDensity(double d) {
            density = d;
            this->setMass(density *4* M_PI/3 * pow(radius,3)); //mass = density * 4pi/3 * rayon^3 ng //Ajout du *4??
        }

        /*!
         * \brief adhesion getter
         * @return
         */
        inline double getAdhesion() const { return adhesion; }

        /*!
         * \brief adhesion setter
         * @param a
         */
        inline void setAdhesion(double a) { adhesion = a; }

        /*!
         * \brief moves a cell to position v
         * @param v
         */
        void moveTo(const MecaCell::Vec v) {
            this->setPosition(v);
        }
    };
}

#endif
