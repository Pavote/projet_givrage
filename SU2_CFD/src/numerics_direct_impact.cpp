/*!
 * \file numerics_direct_impact.cpp
 * \brief This file contains all the convective term discretization.
 * \author B. Constant, M. Fleurotte, A. Motte, I. Moufid
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>

CSourceDropletDrag::CSourceDropletDrag(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
    
    /*--- Store the pointer to the constant body force vector. ---*/
    
    Body_Force_Vector = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Body_Force_Vector[iDim] = config->GetBody_Force_Vector()[iDim];
    
}

CSourceDropletDrag::~CSourceDropletDrag(void) {
    
    if (Body_Force_Vector != NULL) delete [] Body_Force_Vector;
    
}

void CSourceDropletDrag::ComputeResidual(su2double *val_residual, CConfig *config) {
    
    unsigned short iDim;
    su2double Droplet_LWC = config->GetDroplet_LWC();
    
    /*--- Zero the continuity contribution ---*/
    
    val_residual[0] = 0.0;
    
    /*--- Momentum contribution ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = -Volume * U_i[0] * Body_Force_Vector[iDim] / Droplet_LWC;
    
    /*--- Energy contribution ---*/
    
    val_residual[nDim+1] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
        val_residual[nDim+1] += -Volume * U_i[iDim+1] * Body_Force_Vector[iDim] / Droplet_LWC;
    
}
