/******************************************************************************
 * Copyright (C) Siarhei Uzunbajakau, 2023.
 *
 * This program is free software. You can use, modify, and redistribute it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 or (at your option) any later version.
 * This program is distributed without any warranty.
 *
 * Refer to COPYING.LESSER for more details.
 ******************************************************************************/

#ifndef ProjectThetaToH_H__
#define ProjectThetaToH_H__

#include "project_PSI_to_H.hpp"

namespace StaticScalarSolver {
template<int dim, int stage = 1>
using ProjectTHETAtoH = ProjectPSItoH<dim, stage>;
}

#endif
