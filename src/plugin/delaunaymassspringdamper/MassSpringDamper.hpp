#ifndef MASSSPRINGDAMPER_HPP
#define MASSSPRINGDAMPER_HPP

/*!
 * \file MassSpringDamper.hpp
 * \brief mass spring damper physics
 *
 * contains the mass spring class which is used in the DelaunayMassSpringDamper plugin
 */

#include <mecacell/mecacell.h>
#include <mecacell/utilities/logger.hpp>

#define MECACELL_TERMINAL_COLORS

/*!
 * \class MassSpringDamper
 * @tparam cell_t
 * \brief mass spring damper based class
 */
template <typename cell_t>
class MassSpringDamper {

private:
    std::pair<cell_t*,cell_t*> connection; /*!< pair of cells that are connected*/
    float_t stiffness; /*!< spring's stiffness in ng/s^2*/
    float_t dampCoef; /*!< spring's damping coefficient in ng/s*/
    float_t restLength ; /*!< spring's rest length in µm*/
    float_t length; /*!< spring's current length in µm*/
    MecaCell::Vec direction; /*!< spring's current direction from node 0 to 1*/

    /*!
     * \brief updates the current length and direction
     * @param p0
     * @param p1
     */
    void updateLengthDirection(const MecaCell::Vec &p0, const MecaCell::Vec &p1) {
        direction = p1 - p0;
        length = direction.length();
        if (length > 0) direction /= length;
    }

    /*!
     * \brief get the total speed
     * @return the cumulate cells' speeds
     */
    double computeSpeedFromCells(){
        double sp1 = connection.first->getBody().getVelocity().dot(direction);
        double sp2 = -connection.second->getBody().getVelocity().dot(direction);
        return -(sp1+sp2);

    }

public:
    /*!
     * \brief default constructor
     */
    inline MassSpringDamper() = default;

    /*!
     * \brief constructor
     * @param c1 : first cell
     * @param c2 : second cell
     * @param K : stiffness
     * @param C : damping coefficient
     * @param L : rest length
     */
    inline MassSpringDamper(cell_t* c1,cell_t* c2,const float_t &K, const float_t &C, const float_t &L)
        : connection(c1,c2), stiffness(K), dampCoef(C), restLength(L), length(L)
        {}


    /*!
     * \brief rest length setter
     * @param L
     */
    inline void setRestLength(float_t L) { restLength = L; }

    /*!
     * \brief restLength getter
     * @return rest length
     */
    inline float_t getRestLength() const { return restLength; }

    /*!
     * \brief length getter
     * @return current length
     */
    inline float_t getLength() const { return length; }

    /*!
     * \brief connection getter
     * @return connection
     */
    inline pair<cell_t *, cell_t *> getConnection() const { return connection; }

    /*!
     * \brief compute forces for the 2 cells
     */
    void computeForces() {
        // BASIC SPRING
        updateLengthDirection(connection.first->getBody().getPosition(),
                                 connection.second->getBody().getPosition());

        float_t x = length - restLength;  // actual compression or elongation
        bool compression = x < 0;
        double v = computeSpeedFromCells();

        float_t f = (-stiffness * x - dampCoef * v) / 2.0;
        connection.first->getBody().receiveForce(f, -direction, compression);
        connection.second->getBody().receiveForce(f, direction, compression);

    }

    /*!
     * \brief get the connection as a string
     * @return "first cell's id <-> second cell's id"
     */
    inline string toString(){
        return std::to_string(connection.first->id) + " <-> " + std::to_string(connection.second->id);
    }
};

#endif
