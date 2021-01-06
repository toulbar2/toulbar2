/*
 * \file tb2dgvns.hpp
 * \brief sequential decomposition guided variable neighborhood search method
 *
 *  Created on: 12 December 2016
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#ifndef TB2DGVNS_HPP_
#define TB2DGVNS_HPP_

#include "tb2vns.hpp"

class VNSSolver : public LocalSearch {
public:
    VNSSolver(Cost initUpperBound)
        : LocalSearch(initUpperBound)
    {
    }
    ~VNSSolver() {}
    bool solve(bool first = true);
};

#endif /* TB2DGVNS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
