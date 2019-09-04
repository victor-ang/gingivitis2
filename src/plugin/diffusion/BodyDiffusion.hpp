#ifndef ONKO3D_3_0_BODYDIFFUSION_HPP
#define ONKO3D_3_0_BODYDIFFUSION_HPP

/*!
 * \file BodyDiffusion.hpp
 * \brief molecules diffusion body
 *
 * contains the body class which works with the diffusion plugin
 */

#include "DiffusionGrid.hpp"
#include <vector>
#include <hash_map>

using namespace std;

/*!
 * \namespace Diffusion
 */
namespace Diffusion {

/*!
 * \brief molecules diffusion body class
 */
    class BodyDiffusion {

    private:
        vector<double> quantities; /*!< quantities of each molecules in mmHg*/
        vector<double> consumptions; /*!< consumption of each molecules in µm^3/ng/s*/
        DiffusionGrid *grid;

    public:

        /*!
         * \brief default constructor
         */
        inline BodyDiffusion() : quantities(), consumptions(), grid() {}

        /*!
         * \brief constructor
         * @param size : number of existing molecules
         */
        inline explicit BodyDiffusion(unsigned long size, DiffusionGrid *g)
        : grid(g)
        {
            for(unsigned i = 0; i < size; ++i) {
                quantities.push_back(0);
                consumptions.push_back(0);
            }
        }

        inline DiffusionGrid *getGrid() const { return grid; }

        inline void setGrid(DiffusionGrid *g) { grid = g;}
        /*!
         * \brief quantities getter
         * @return quantities
         */
        inline vector<double> getQuantities() const{ return quantities; }

        /*!
        * \brief consumptions getter
        * @return consumptions
        */
        inline vector<double> getConsumptions() const{ return consumptions; }

        /*!
         * \brief quantity setter
         * @param i : index of the molecule
         * @param value
         */
        inline void setQuantity(unsigned i, double value){ quantities[i] = value;}

        /*!
         * \brief consumption setter
         * @param i : index of the molecule
         * @param value
         */
        inline void setConsumption(unsigned i, double value){ consumptions[i] = value;}

    };

}


#endif //ONKO3D_3_0_BODYDIFFUSION_HPP
