#ifndef ONKO3D_3_0_DIFFUSIONGRID_HPP
#define ONKO3D_3_0_DIFFUSIONGRID_HPP
/*!
 * \file DiffusionGrid.hpp
 * \brief molecules grid behaviour
 *
 * contains the molecules grid based class that is used in the diffusion plugin
 */

#include <vector>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <mecacell/utilities/utils.hpp>



/*!
 * \namespace Diffusion
 */
namespace Diffusion {
/*!
 * \struct Molecule
 */
    struct Molecule {
        double diffusion; /*!< diffusion constant in µm^2/s*/
        double defaultQuantity; /*!< quantity of the molecule in the environment without cells in mmHg*/
        double omega; /*!< Henry constant for the molecule at human temperature multiplied by density of a tumor and density of the molecule in mmHg.ng/µm^3*/
        double defaultEvaporation; /*!< consumption constant of the environment µm^3/ng/s*/
        inline Molecule(double d, double dQ, double density, double dE) : diffusion(d), defaultQuantity(dQ),
        omega(30.318 * density / 1.428), defaultEvaporation(dE) {}

    };


/*!
 * \enum CELLTYPE
 */
    enum CELLTYPE {
        CONSUMPTION, EMPTY
    };

/*!
 * \struct GridCell
 */
    struct GridCell {
        vector<double> prevQuantities; /*!< previous quantities of each molecule*/
        vector<double> quantities; /*!< current quantities of each molecule*/
        vector<double> consumptions; /*!< local consumptions of each molecule*/
        CELLTYPE type;

        /*!
         * \brief default constructor
         */
        inline GridCell() = default;

        /*!
         * \brief constructor
         * @param m : existing molecules in the simulation
         */
        explicit GridCell(vector<Molecule> m) {
            for (unsigned i = 0; i < m.size(); ++i) {
                quantities.push_back(m[i].defaultQuantity);
                prevQuantities.push_back(m[i].defaultQuantity);
                type = CELLTYPE::EMPTY;
                consumptions.push_back(0.);
            }
        }
    };

/*!
 * \class DiffusionGrid
 * \brief molecules grid based class
 */
    class DiffusionGrid {

    private:
        double dx; /*!< cell size*/
        double accuracy; /*!< precision of the diffusion*/
        unordered_map<MecaCell::Vec, GridCell> grid; /*!< gridcells hashmap*/
        vector<Molecule> molecules; /*!< existing molecules in the simulation*/

        //Boundaries
        unordered_map<MecaCell::Vec, int> xMin;
        unordered_map<MecaCell::Vec, int> xMax;
        unordered_map<MecaCell::Vec, int> yMin;
        unordered_map<MecaCell::Vec, int> yMax;
        unordered_map<MecaCell::Vec, int> zMin;
        unordered_map<MecaCell::Vec, int> zMax;

        /*!
         * \brief 2nd derivative calculator
         * @param xm1
         * @param x
         * @param xp1
         * @return approximate 2nd derivative in x
         */
        inline double deriv2nd(double xm1, double x, double xp1) { return ((xm1 - (2.0 * x) + xp1) / (dx * dx)); }

        /*!
         * \brief next step cells initializer
         */
        void nextStep(unsigned pos) {
            for (auto &c : grid) {
                c.second.prevQuantities[pos] = c.second.quantities[pos];
            }
        }

        /*!
         * \brief boundaries clearer
         */
        void clearBoundary() {
            xMin.clear();
            yMin.clear();
            zMin.clear();
            xMax.clear();
            yMax.clear();
            zMax.clear();
        }

        /*!
         * \brief cells emptier
         */
        void allCellEmpty() {
            for (auto &c : grid) {
                c.second.type = CELLTYPE::EMPTY;
            }
        }

        /*!
         * \brief out of boundaries cells cleaner
         */
        void cleanCellOutBoundary() {
            unordered_map<MecaCell::Vec, GridCell> newGrid = unordered_map<MecaCell::Vec, GridCell>();
            for (auto &c: grid) {
                if (!isOutBoundary(c.first)) {
                    newGrid.insert(c);
                }
            }
            grid = newGrid;
        }

        /*!
         * \brief cell boundaries verifier
         * @param v
         * @return true if v is out of the boundaries
         */
        bool isOutBoundary(const MecaCell::Vec &v) {
            MecaCell::Vec onX = MecaCell::Vec(v.y(), v.z(), 0);
            MecaCell::Vec onY = MecaCell::Vec(v.x(), v.z(), 0);
            MecaCell::Vec onZ = MecaCell::Vec(v.x(), v.y(), 0);
            return (v.x() < xMin[forward<MecaCell::Vec>(onX)] || v.x() > xMax[forward<MecaCell::Vec>(onX)] ||
                    v.y() < yMin[forward<MecaCell::Vec>(onY)] || v.y() > yMax[forward<MecaCell::Vec>(onY)] ||
                    v.z() < zMin[forward<MecaCell::Vec>(onZ)] || v.z() > zMax[forward<MecaCell::Vec>(onZ)]);
        }

        

        /*!
         * \brief cell inserter
         * @param pos
         */
        void insertPartOfCell(MecaCell::Vec pos, vector<double> consumptions) {
            if (grid.count(forward<MecaCell::Vec>(pos)) <= 0)
                grid[forward<MecaCell::Vec>(pos)] = GridCell(molecules);
            for (unsigned i = 0; i < molecules.size(); ++i) {
                grid[forward<MecaCell::Vec>(pos)].consumptions[i] = consumptions[i];
            }
            grid[forward<MecaCell::Vec>(pos)].type = CELLTYPE::CONSUMPTION;

            boundaryUpdate(pos);
        }

        /*!
         * \brief min boundary updater
         * @param bd
         * @param coordFace
         */
        static void updateMinBoundary(unordered_map<MecaCell::Vec, int> &bd, pair<MecaCell::Vec, int> coordFace) {
            if (bd.count(forward<MecaCell::Vec>(coordFace.first)) <= 0) {
                bd.insert(coordFace);
            } else {
                if (coordFace.second < bd[forward<MecaCell::Vec>(coordFace.first)]) {
                    bd[forward<MecaCell::Vec>(coordFace.first)] = coordFace.second;
                }
            }
        }

        /*!
         * \brief max boundary updater
         * @param bd
         * @param coordFace
         */
        static void updateMaxBoundary(unordered_map<MecaCell::Vec, int> &bd, pair<MecaCell::Vec, int> coordFace) {
            if (bd.count(forward<MecaCell::Vec>(coordFace.first)) <= 0) {
                bd.insert(coordFace);
            } else {
                if (coordFace.second > bd[forward<MecaCell::Vec>(coordFace.first)]) {
                    bd[forward<MecaCell::Vec>(coordFace.first)] = coordFace.second;
                }
            }
        }

        /*!
         * \brief boundaries updater
         * @param v
         */
        void boundaryUpdate(const MecaCell::Vec &v) {
            pair<MecaCell::Vec, int> coordX = make_pair<MecaCell::Vec, int>(MecaCell::Vec(v.y(), v.z(), 0), v.x());
            updateMinBoundary(xMin, coordX);
            updateMaxBoundary(xMax, coordX);
            pair<MecaCell::Vec, int> coordY = make_pair<MecaCell::Vec, int>(MecaCell::Vec(v.x(), v.z(), 0), v.y());
            updateMinBoundary(yMin, coordY);
            updateMaxBoundary(yMax, coordY);
            pair<MecaCell::Vec, int> coordZ = make_pair<MecaCell::Vec, int>(MecaCell::Vec(v.x(), v.y(), 0), v.z());
            updateMinBoundary(zMin, coordZ);
            updateMaxBoundary(zMax, coordZ);
        }




        /*!
         * \brief diffusion dt getter
         * @param D : diffusion constant of the molecule
         * @return diffusion dt
         */
        inline double getDt(double D) const { return (dx * dx) / (6 * D); }

        /*!
         * \brief molecules updater
         *
         * updates the quantity of a molecule for each cell and returns the total amount of this molecule
         *
         * @return total quantity of the molecule
         */
        double computeStep(unsigned pos) {
            Molecule m = molecules[pos]; //considered molecule
            double D = m.diffusion; //Diffusion constant
            double dt = getDt(D); //Diffusion dt
            double total = 0; //total amount of considered molecule
            for (auto &c : grid) { //for each cell in the grid
                GridCell &gc = c.second;
                double laplacien = deriv2nd(getMolecule(c.first - MecaCell::Vec(1, 0, 0), pos), gc.prevQuantities[pos],
                                            getMolecule(c.first + MecaCell::Vec(1, 0, 0), pos)) +
                                   deriv2nd(getMolecule(c.first - MecaCell::Vec(0, 1, 0), pos), gc.prevQuantities[pos],
                                            getMolecule(c.first + MecaCell::Vec(0, 1, 0), pos)) +
                                   deriv2nd(getMolecule(c.first - MecaCell::Vec(0, 0, 1), pos), gc.prevQuantities[pos],
                                            getMolecule(c.first + MecaCell::Vec(0, 0, 1), pos));
                if (gc.type == CELLTYPE::CONSUMPTION) {
                    gc.quantities[pos] =
                            gc.prevQuantities[pos] + (D * laplacien - ((gc.consumptions[pos] + m.defaultEvaporation) * m.omega)) * dt;
                } else {
                    gc.quantities[pos] = gc.prevQuantities[pos] + (D * laplacien - m.defaultEvaporation * m.omega) * dt;
                }
                if (gc.quantities[pos] < 0.) { //can't have a negative amount of molecule
                    gc.quantities[pos] = 0.;
                }
                total += gc.quantities[pos];
            }
            return (total);
        }

    public:

        /*!
         * \brief constructor
         * @param dx
         * @param accuracy
         */
        inline DiffusionGrid(double dx, double accuracy) : dx(dx), accuracy(accuracy), grid(), molecules() {}

        /*!
         * \brief grid's size getter
         * @return grid's size
         */
        inline size_t size() const { return grid.size(); };

        /*!
         * \brief grid getter
         * @return grid
         */
        inline unordered_map<MecaCell::Vec, GridCell> &getGrid() { return grid; }

        /*!
         * \brief dx getter
         * @return dx
         */
        inline double getCellSize() const { return dx; }

        /*!
         * \brief molecules getter
         * @return molecules
         */
        inline vector<Molecule> getMolecules() const { return molecules; }

        /*!
         * \brief molecule adder
         * @param m : molecule to be added
         */
        inline void addMolecule(Molecule m) { molecules.push_back(m); }

        /*!
         * \brief cell inserter
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void insertCell(world_t *w) {
            clearBoundary();
            allCellEmpty();
            for (auto &c : w->cells) { //for each cell
                const MecaCell::Vec center = c->getBody().getPosition();
                const double &radius = c->getBody().getBoundingBoxRadius();
                const double sqRadius = radius * radius;
                MecaCell::Vec minCorner = getIndexFromPosition(center - radius);
                MecaCell::Vec maxCorner = getIndexFromPosition(center + radius);
                for (double i = minCorner.x(); i <= maxCorner.x(); ++i) {
                    for (double j = minCorner.y(); j <= maxCorner.y(); ++j) {
                        for (double k = minCorner.z(); k <= maxCorner.z(); ++k) {
                            double cx = (i + 0.5) * dx;
                            double cy = (j + 0.5) * dx;
                            double cz = (k + 0.5) * dx;
                            MecaCell::Vec cubeCenter(cx, cy, cz);
                            if ((cubeCenter - center).sqlength() < sqRadius) {
                                insertPartOfCell(MecaCell::Vec(i, j, k), c->getBody().getConsumptions());
                            }
                        }
                    }
                }
            }
            cleanCellOutBoundary();
        }

        /*!
         * \brief molecules quantities updater
         * @tparam world_t
         * @param w
         */
        template<typename world_t>
        void computeMolecules(world_t *w) {
            vector<double> lasts; //last step quantities
            vector<double> news; //current step quantities
            for (unsigned i = 0; i < molecules.size(); ++i) { //for each molecule
                lasts.push_back(0);
                news.push_back(INFINITY);
            }
            for (unsigned i = 0; i < molecules.size(); ++i) { //for each molecule
                if (i==2) { // a enlever pour faire diffuser les 3 molecules
                    double t = 0; //elapsed time
                    double dt = getDt(molecules[i].diffusion); //Diffusion dt
                    while (abs(lasts[i] - news[i]) / grid.size() >= accuracy && t < w->dt) {
                        lasts[i] = news[i];
                        nextStep(i);
                        news[i] = computeStep(i);
                        t += dt;
                    }
                } 
                

            }
        }

        /*!
         * \brief molecule quantity getter
         * @param v
         * @return molecule quantity at position v
         */
        inline double getMoleculeRealPos(const MecaCell::Vec &v, unsigned pos) {
            return (getMolecule(getIndexFromPosition(v), pos));
        }


        /*!
         * \brief index getter
         * @param v
         * @return grid index of the v position
         */
        MecaCell::Vec getIndexFromPosition(const MecaCell::Vec &v) const {
            MecaCell::Vec res = v / dx;
            return MecaCell::Vec(floor(res.x()), floor(res.y()), floor(res.z()));
        }
        /*!
         * \brief molecule quantity getter
         * @param v
         * @return molecule_t quantity at position v
         */
        double getMolecule(MecaCell::Vec v, unsigned pos) {
            if (grid.count(v) <= 0)
                return (molecules[pos].defaultQuantity);
            else
                return (grid[forward<MecaCell::Vec>(v)].prevQuantities[pos]);

        }
    };

    
}


#endif //ONKO3D_3_0_DIFFUSIONGRID_HPP
