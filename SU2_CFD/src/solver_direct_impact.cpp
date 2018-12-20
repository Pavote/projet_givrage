/*!
 * \file solution_direct_impact.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author B. Constant, M. Fleurotte, A. Motte, I. Moufid, F. Morency
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

#include "../include/solver_structure.hpp"

CImpactSolver::CImpactSolver(void) : CSolver() {

  cout << "CImpactSolver*************************************************************************************" << endl;
  /*--- Basic array initialization ---*/

  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt = NULL; CL_Mnt = NULL; CSF_Mnt = NULL;  CEff_Mnt = NULL;
  CMx_Mnt = NULL; CMy_Mnt = NULL; CMz_Mnt = NULL;
  CFx_Mnt = NULL; CFy_Mnt = NULL; CFz_Mnt = NULL;
  CoPx_Mnt = NULL; CoPy_Mnt = NULL; CoPz_Mnt = NULL;

  YPlus = NULL;

  /*--- Surface based array initialization ---*/

  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt = NULL; Surface_CD_Mnt = NULL; Surface_CSF_Mnt = NULL; Surface_CEff_Mnt = NULL;
  Surface_CFx_Mnt = NULL; Surface_CFy_Mnt = NULL; Surface_CFz_Mnt = NULL;
  Surface_CMx_Mnt = NULL; Surface_CMy_Mnt = NULL; Surface_CMz_Mnt = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Rotorcraft simulation array initialization ---*/

  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;

  CMerit_Mnt = NULL;  CT_Mnt = NULL;  CQ_Mnt = NULL;

  /*--- Supersonic simulation array initialization ---*/

  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;

  /*--- Engine simulation array initialization ---*/

  Inflow_MassFlow = NULL;   Inflow_Pressure = NULL;
  Inflow_Mach = NULL;       Inflow_Area = NULL;
  Exhaust_Pressure = NULL;  Exhaust_Temperature = NULL;
  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;

  /*--- Numerical methods array initialization ---*/

  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;

  DonorPrimVar = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT = NULL;

  Smatrix = NULL; Cvector = NULL;

  Secondary = NULL; Secondary_i = NULL; Secondary_j = NULL;

  /*--- Fixed CL mode initialization (cauchy criteria) ---*/

  Cauchy_Value   = 0;
  Cauchy_Func    = 0;
  Old_Func       = 0;
  New_Func       = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;

  AoA_FD_Change = false;

  FluidModel   = NULL;

  SlidingState     = NULL;
  SlidingStateNodes = NULL;

  /*--- Initialize quantities for the average process for internal flow ---*/

  AverageVelocity            = NULL;
  AverageTurboVelocity       = NULL;
  OldAverageTurboVelocity    = NULL;
  ExtAverageTurboVelocity    = NULL;
  AverageFlux                = NULL;
  SpanTotalFlux              = NULL;
  AveragePressure            = NULL;
  OldAveragePressure         = NULL;
  RadialEquilibriumPressure  = NULL;
  ExtAveragePressure         = NULL;
  AverageDensity             = NULL;
  OldAverageDensity          = NULL;
  ExtAverageDensity          = NULL;
  AverageNu                  = NULL;
  AverageKine                = NULL;
  AverageOmega               = NULL;
  ExtAverageNu               = NULL;
  ExtAverageKine             = NULL;
  ExtAverageOmega            = NULL;


  /*--- Initialize primitive quantities for turboperformace ---*/

  DensityIn                     = NULL;
  PressureIn                    = NULL;
  TurboVelocityIn               = NULL;
  DensityOut                    = NULL;
  PressureOut                   = NULL;
  TurboVelocityOut              = NULL;
  KineIn                        = NULL;
  OmegaIn                       = NULL;
  NuIn                          = NULL;
  KineOut                       = NULL;
  OmegaOut                      = NULL;
  NuOut                         = NULL;


  /*--- Initialize quantities for Giles BC---*/

  CkInflow                      = NULL;
  CkOutflow1                    = NULL;
  CkOutflow2                    = NULL;

}

CImpactSolver::CImpactSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

  cout << "CImpactSolver*************************************************************************************" << endl;
  unsigned long iPoint, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  su2double StaticEnergy, Density, Velocity2, Pressure, Temperature;
  unsigned short nZone = geometry->GetnZone();
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));
  unsigned short direct_diff = config->GetDirectDiff();
  int Unst_RestartIter = 0;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  bool low_mach_prec = config->Low_Mach_Preconditioning();

  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  bool fsi     = config->GetFSI_Simulation();
  string filename_ = config->GetSolution_FlowFileName();

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone);

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/

    Read_SU2_Restart_Metadata(geometry, config, false, filename_);

  }

  /*--- Array initialization ---*/

  /*--- Basic array initialization ---*/

  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt= NULL; CL_Mnt= NULL; CSF_Mnt= NULL; CEff_Mnt= NULL;
  CMx_Mnt= NULL;   CMy_Mnt= NULL;   CMz_Mnt= NULL;
  CFx_Mnt= NULL;   CFy_Mnt= NULL;   CFz_Mnt= NULL;
  CoPx_Mnt= NULL;   CoPy_Mnt= NULL;   CoPz_Mnt= NULL;

  YPlus = NULL;

  /*--- Surface based array initialization ---*/

  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt= NULL; Surface_CD_Mnt= NULL; Surface_CSF_Mnt= NULL; Surface_CEff_Mnt= NULL;
  Surface_CFx_Mnt= NULL;   Surface_CFy_Mnt= NULL;   Surface_CFz_Mnt= NULL;
  Surface_CMx_Mnt= NULL;   Surface_CMy_Mnt= NULL;   Surface_CMz_Mnt = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Rotorcraft simulation array initialization ---*/

  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;

  CMerit_Mnt = NULL; CT_Mnt = NULL; CQ_Mnt = NULL;

  /*--- Supersonic simulation array initialization ---*/

  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;

  /*--- Engine simulation array initialization ---*/

  Inflow_MassFlow = NULL;   Inflow_Pressure = NULL;
  Inflow_Mach = NULL;       Inflow_Area = NULL;
  Exhaust_Pressure = NULL;  Exhaust_Temperature = NULL;
  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;

  /*--- Numerical methods array initialization ---*/

  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  DonorPrimVar = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT = NULL;

  Smatrix = NULL; Cvector = NULL;

  Secondary=NULL; Secondary_i=NULL; Secondary_j=NULL;

  /*--- Fixed CL mode initialization (cauchy criteria) ---*/

  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;

  AoA_FD_Change = false;

  FluidModel = NULL;

  /*--- Initialize quantities for the average process for internal flow ---*/

  AverageVelocity                   = NULL;
  AverageTurboVelocity              = NULL;
  OldAverageTurboVelocity           = NULL;
  ExtAverageTurboVelocity           = NULL;
  AverageFlux                       = NULL;
  SpanTotalFlux                     = NULL;
  AveragePressure                   = NULL;
  OldAveragePressure                = NULL;
  RadialEquilibriumPressure         = NULL;
  ExtAveragePressure                = NULL;
  AverageDensity                    = NULL;
  OldAverageDensity                 = NULL;
  ExtAverageDensity                 = NULL;
  AverageNu                         = NULL;
  AverageKine                       = NULL;
  AverageOmega                      = NULL;
  ExtAverageNu                      = NULL;
  ExtAverageKine                    = NULL;
  ExtAverageOmega                   = NULL;


  /*--- Initialize primitive quantities for turboperformace ---*/

  DensityIn                     = NULL;
  PressureIn                    = NULL;
  TurboVelocityIn               = NULL;
  DensityOut                    = NULL;
  PressureOut                   = NULL;
  TurboVelocityOut              = NULL;
  KineIn                        = NULL;
  OmegaIn                       = NULL;
  NuIn                          = NULL;
  KineOut                       = NULL;
  OmegaOut                      = NULL;
  NuOut                         = NULL;


  /*--- Initialize quantities for Giles BC---*/

  CkInflow                      = NULL;
  CkOutflow1                    = NULL;
  CkOutflow2                    = NULL;

  /*--- Set the gamma value ---*/

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure
  /*--- For now the impact solver ( T, vx, vy, (vz), P, rho, h, c, lamMu, EddyMu, ThCond, Cp,
        rhoair, uair, vair, (wair), lamMu_air, T_air ).  We still need to remove some variable like
        EddyMu, Cp, etc. that are meaningless for impact. ---*/

  nDim = geometry->GetnDim();

  nVar = nDim+2;
  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  /*--- Add air flow variable ---*/
  //FM nPrimVar = nPrimVar+nDim+3;
  nSecondaryVar = 2; nSecondaryVarGrad = 2;


  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nPrimVarGrad;

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/

  SetNondimensionalization(geometry, config, iMesh);

  /*--- Allocate the node variables ---*/

  node = new CVariable*[nPoint];

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the geometry ---*/

  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

  /*--- Define some auxiliary vectors related to the primitive solution ---*/

  Primitive   = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the Secondary solution ---*/

  Secondary   = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
  Secondary_i = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
  Secondary_j = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the air solution ---*/

  //Solution_Air   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_Air[iVar]   = 0.0;

  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }

  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/

  if (roe_turkel || low_mach_prec) {
    LowMach_Precontioner = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new su2double[nVar];
  }

  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

  }

  else {
    if (rank == MASTER_NODE) cout << "Explicit scheme. No Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
  }

  /*--- Define some auxiliary vectors for computing flow variable
   gradients by least squares, S matrix := inv(R)*traspose(inv(R)),
   c vector := transpose(WA)*(Wb) ---*/

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];

    Cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      Cvector[iVar] = new su2double [nDim];

  }

  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/

  CharacPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }

  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/

  DonorPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      if (rans) {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar+2];
        for (iVar = 0; iVar < nPrimVar + 2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
        for (iVar = 0; iVar < nPrimVar ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
    }
  }

  /*--- Store the value of the characteristic primitive variables index at the boundaries ---*/

  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
    }
  }

  /*--- Store the value of the Delta P at the Actuator Disk ---*/

  ActDisk_DeltaP = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaP[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaP[iMarker][iVertex] = 0;
    }
  }

  /*--- Store the value of the Delta T at the Actuator Disk ---*/

  ActDisk_DeltaT = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaT[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaT[iMarker][iVertex] = 0;
    }
  }

  /*--- Store the value of the Total Pressure at the inlet BC ---*/

  Inlet_Ttotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Inlet_Ttotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_Ttotal[iMarker][iVertex] = 0;
      }
  }

  /*--- Store the value of the Total Temperature at the inlet BC ---*/

  Inlet_Ptotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Inlet_Ptotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_Ptotal[iMarker][iVertex] = 0;
      }
  }

  /*--- Store the value of the Flow direction at the inlet BC ---*/

  Inlet_FlowDir = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Inlet_FlowDir[iMarker] = new su2double* [geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        Inlet_FlowDir[iMarker][iVertex] = new su2double [nDim];
        for (iDim = 0; iDim < nDim; iDim++) {
          Inlet_FlowDir[iMarker][iVertex][iDim] = 0;
        }
      }
  }

  /*--- Non-dimensional coefficients ---*/

  CD_Inv         = new su2double[nMarker];
  CL_Inv         = new su2double[nMarker];
  CSF_Inv        = new su2double[nMarker];
  CMx_Inv           = new su2double[nMarker];
  CMy_Inv           = new su2double[nMarker];
  CMz_Inv           = new su2double[nMarker];
  CEff_Inv          = new su2double[nMarker];
  CFx_Inv           = new su2double[nMarker];
  CFy_Inv           = new su2double[nMarker];
  CFz_Inv           = new su2double[nMarker];
  CoPx_Inv           = new su2double[nMarker];
  CoPy_Inv           = new su2double[nMarker];
  CoPz_Inv           = new su2double[nMarker];

  CD_Mnt        = new su2double[nMarker];
  CL_Mnt        = new su2double[nMarker];
  CSF_Mnt       = new su2double[nMarker];
  CMx_Mnt          = new su2double[nMarker];
  CMy_Mnt          = new su2double[nMarker];
  CMz_Mnt          = new su2double[nMarker];
  CEff_Mnt         = new su2double[nMarker];
  CFx_Mnt          = new su2double[nMarker];
  CFy_Mnt          = new su2double[nMarker];
  CFz_Mnt          = new su2double[nMarker];
  CoPx_Mnt          = new su2double[nMarker];
  CoPy_Mnt          = new su2double[nMarker];
  CoPz_Mnt          = new su2double[nMarker];

  Surface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Inv     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL_Mnt     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Mnt     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Mnt= new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];

  /*--- Rotorcraft coefficients ---*/

  CT_Inv           = new su2double[nMarker];
  CQ_Inv           = new su2double[nMarker];
  CMerit_Inv       = new su2double[nMarker];

  CT_Mnt           = new su2double[nMarker];
  CQ_Mnt           = new su2double[nMarker];
  CMerit_Mnt       = new su2double[nMarker];

  /*--- Supersonic coefficients ---*/

  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];

  /*--- Engine simulation ---*/

  Inflow_MassFlow     = new su2double[nMarker];
  Inflow_Pressure     = new su2double[nMarker];
  Inflow_Mach         = new su2double[nMarker];
  Inflow_Area         = new su2double[nMarker];

  Exhaust_MassFlow    = new su2double[nMarker];
  Exhaust_Pressure    = new su2double[nMarker];
  Exhaust_Temperature = new su2double[nMarker];
  Exhaust_Area        = new su2double[nMarker];

  /*--- Init total coefficients ---*/

  Total_CD      = 0.0;    Total_CL           = 0.0;    Total_CSF          = 0.0;
  Total_CMx     = 0.0;    Total_CMy          = 0.0;    Total_CMz          = 0.0;
  Total_CoPx    = 0.0;    Total_CoPy         = 0.0;    Total_CoPz         = 0.0;
  Total_CEff    = 0.0;    Total_CEquivArea   = 0.0;    Total_CNearFieldOF = 0.0;
  Total_CFx     = 0.0;    Total_CFy          = 0.0;    Total_CFz          = 0.0;
  Total_CT      = 0.0;    Total_CQ           = 0.0;    Total_CMerit       = 0.0;
  Total_ComboObj     = 0.0;
  Total_CpDiff  = 0.0;    Total_Custom_ObjFunc=0.0;
  Total_NetThrust = 0.0;
  Total_Power = 0.0;      AoA_Prev           = 0.0;
  Total_CL_Prev = 0.0;    Total_CD_Prev      = 0.0;
  Total_CMx_Prev      = 0.0; Total_CMy_Prev      = 0.0; Total_CMz_Prev      = 0.0;
  Total_AeroCD = 0.0;  Total_SolidCD = 0.0;   Total_IDR   = 0.0;    Total_IDC   = 0.0;

  /*--- Read farfield conditions ---*/

  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Mach_Inf        = config->GetMach();

  /*--- Initialize the secondary values for direct derivative approxiations ---*/

  switch(direct_diff) {
    case NO_DERIVATIVE:
      /*--- Default ---*/
      break;
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
      break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
      break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
      break;
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Already done in postprocessing of config ---*/
      break;
    default:
      break;
  }


  /*--- Initialize fan face pressure, fan face mach number, and mass flow rate ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;

    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;
  }

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++){
    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];

        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }

    }
  }


  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CImpactVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/

  counter_local = 0;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    Density = node[iPoint]->GetSolution(0);

    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);

    StaticEnergy= node[iPoint]->GetSolution(nDim+1)/Density - 0.5*Velocity2;

    FluidModel->SetTDState_rhoe(Density, StaticEnergy);
    Pressure= FluidModel->GetPressure();
    Temperature= FluidModel->GetTemperature();

    /*--- Use the values at the infinity ---*/

    if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);
      counter_local++;
    }

  }

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (fsi){
    Residual_BGS      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_Max_BGS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 0.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }
  }

  /*--- Warning message about non-physical points ---*/

  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }

  /*--- Define solver parameters needed for execution of destructor ---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;

  /*--- Perform the MPI communication of the solution ---*/

  Set_MPI_Solution(geometry, config);

}

CImpactSolver::~CImpactSolver(void) {

  unsigned short iVar, iMarker, iSpan;

  unsigned long iVertex;

  /*--- Array deallocation ---*/

  if (CD_Inv != NULL)         delete [] CD_Inv;
  if (CL_Inv != NULL)         delete [] CL_Inv;
  if (CSF_Inv != NULL)    delete [] CSF_Inv;
  if (CMx_Inv != NULL)           delete [] CMx_Inv;
  if (CMy_Inv != NULL)           delete [] CMy_Inv;
  if (CMz_Inv != NULL)           delete [] CMz_Inv;
  if (CFx_Inv != NULL)           delete [] CFx_Inv;
  if (CFy_Inv != NULL)           delete [] CFy_Inv;
  if (CFz_Inv != NULL)           delete [] CFz_Inv;
  if (CoPx_Inv != NULL)           delete [] CoPx_Inv;
  if (CoPy_Inv != NULL)           delete [] CoPy_Inv;
  if (CoPz_Inv != NULL)           delete [] CoPz_Inv;
  if (Surface_CL_Inv != NULL) delete[] Surface_CL_Inv;
  if (Surface_CD_Inv != NULL) delete[] Surface_CD_Inv;
  if (Surface_CSF_Inv != NULL) delete[] Surface_CSF_Inv;
  if (Surface_CEff_Inv != NULL) delete[] Surface_CEff_Inv;
  if (Surface_CFx_Inv != NULL)  delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv != NULL)  delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv != NULL)  delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;

  if (CD_Mnt != NULL)         delete [] CD_Mnt;
  if (CL_Mnt != NULL)         delete [] CL_Mnt;
  if (CSF_Mnt != NULL)    delete [] CSF_Mnt;
  if (CFx_Mnt != NULL)           delete [] CFx_Mnt;
  if (CFy_Mnt != NULL)           delete [] CFy_Mnt;
  if (CFz_Mnt != NULL)           delete [] CFz_Mnt;
  if (CMx_Mnt != NULL)           delete [] CMx_Mnt;
  if (CMy_Mnt != NULL)           delete [] CMy_Mnt;
  if (CMz_Mnt != NULL)           delete [] CMz_Mnt;
  if (CoPx_Mnt != NULL)           delete [] CoPx_Mnt;
  if (CoPy_Mnt != NULL)           delete [] CoPy_Mnt;
  if (CoPz_Mnt != NULL)           delete [] CoPz_Mnt;
  if (Surface_CL_Mnt != NULL) delete[] Surface_CL_Mnt;
  if (Surface_CD_Mnt != NULL) delete[] Surface_CD_Mnt;
  if (Surface_CSF_Mnt != NULL) delete[] Surface_CSF_Mnt;
  if (Surface_CEff_Mnt != NULL) delete[] Surface_CEff_Mnt;
  if (Surface_CFx_Mnt != NULL)  delete [] Surface_CFx_Mnt;
  if (Surface_CFy_Mnt != NULL)  delete [] Surface_CFy_Mnt;
  if (Surface_CFz_Mnt != NULL)  delete [] Surface_CFz_Mnt;
  if (Surface_CMx_Mnt != NULL)  delete [] Surface_CMx_Mnt;
  if (Surface_CMy_Mnt != NULL)  delete [] Surface_CMy_Mnt;
  if (Surface_CMz_Mnt != NULL)  delete [] Surface_CMz_Mnt;

  if (Surface_CL != NULL)    delete [] Surface_CL;
  if (Surface_CD != NULL)    delete [] Surface_CD;
  if (Surface_CSF != NULL) delete [] Surface_CSF;
  if (Surface_CEff != NULL) delete [] Surface_CEff;
  if (Surface_CFx != NULL)      delete [] Surface_CFx;
  if (Surface_CFy != NULL)      delete [] Surface_CFy;
  if (Surface_CFz != NULL)      delete [] Surface_CFz;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)          delete [] CEff_Inv;
  if (CMerit_Inv != NULL)        delete [] CMerit_Inv;
  if (CT_Inv != NULL)            delete [] CT_Inv;
  if (CQ_Inv != NULL)            delete [] CQ_Inv;
  if (CEquivArea_Inv != NULL)    delete [] CEquivArea_Inv;
  if (CNearFieldOF_Inv != NULL)  delete [] CNearFieldOF_Inv;

  if (CEff_Mnt != NULL)          delete [] CEff_Mnt;
  if (CMerit_Mnt != NULL)        delete [] CMerit_Mnt;
  if (CT_Mnt != NULL)            delete [] CT_Mnt;
  if (CQ_Mnt != NULL)            delete [] CQ_Mnt;

  if (Inflow_MassFlow != NULL)  delete [] Inflow_MassFlow;
  if (Exhaust_MassFlow != NULL)  delete [] Exhaust_MassFlow;
  if (Exhaust_Area != NULL)      delete [] Exhaust_Area;
  if (Inflow_Pressure != NULL)  delete [] Inflow_Pressure;
  if (Inflow_Mach != NULL)      delete [] Inflow_Mach;
  if (Inflow_Area != NULL)      delete [] Inflow_Area;

  if (Exhaust_Pressure != NULL)  delete [] Exhaust_Pressure;
  if (Exhaust_Temperature != NULL)      delete [] Exhaust_Temperature;

  if (iPoint_UndLapl != NULL)       delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL)       delete [] jPoint_UndLapl;

  if (Primitive != NULL)        delete [] Primitive;
  if (Primitive_i != NULL)      delete [] Primitive_i;
  if (Primitive_j != NULL)      delete [] Primitive_j;

  if (Secondary != NULL)        delete [] Secondary;
  if (Secondary_i != NULL)      delete [] Secondary_i;
  if (Secondary_j != NULL)      delete [] Secondary_j;

  if (LowMach_Precontioner != NULL) {
    for (iVar = 0; iVar < nVar; iVar ++)
      delete [] LowMach_Precontioner[iVar];
    delete [] LowMach_Precontioner;
  }

  if (CharacPrimVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] CharacPrimVar[iMarker][iVertex];
      delete [] CharacPrimVar[iMarker];
    }
    delete [] CharacPrimVar;
  }

  if (SlidingState != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != NULL ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }

  if ( SlidingStateNodes != NULL ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
        if (SlidingStateNodes[iMarker] != NULL)
            delete [] SlidingStateNodes[iMarker];
    }
    delete [] SlidingStateNodes;
  }

  if (DonorPrimVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] DonorPrimVar[iMarker][iVertex];
      delete [] DonorPrimVar[iMarker];
    }
    delete [] DonorPrimVar;
  }

  if (DonorGlobalIndex != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] DonorGlobalIndex[iMarker];
    delete [] DonorGlobalIndex;
  }

  if (ActDisk_DeltaP != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ActDisk_DeltaP[iMarker];
    delete [] ActDisk_DeltaP;
  }

  if (ActDisk_DeltaT != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ActDisk_DeltaT[iMarker];
    delete [] ActDisk_DeltaT;
  }

  if (Inlet_Ttotal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Inlet_Ttotal[iMarker] != NULL)
        delete [] Inlet_Ttotal[iMarker];
    delete [] Inlet_Ttotal;
  }

  if (Inlet_Ptotal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Inlet_Ptotal[iMarker] != NULL)
        delete [] Inlet_Ptotal[iMarker];
    delete [] Inlet_Ptotal;
  }

  if (Inlet_FlowDir != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Inlet_FlowDir[iMarker] != NULL) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          delete [] Inlet_FlowDir[iMarker][iVertex];
        delete [] Inlet_FlowDir[iMarker];
      }
    }
    delete [] Inlet_FlowDir;
  }

  if (nVertex != NULL)  delete [] nVertex;

  if (YPlus != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] YPlus[iMarker];
    }
    delete [] YPlus;
  }

  if (Cauchy_Serie != NULL)  delete [] Cauchy_Serie;

  if (FluidModel != NULL) delete FluidModel;

  if(AverageVelocity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] AverageVelocity[iMarker][iSpan];
      delete  [] AverageVelocity[iMarker];
    }
    delete [] AverageVelocity;
  }

  if(AverageTurboVelocity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] AverageTurboVelocity[iMarker][iSpan];
      delete  [] AverageTurboVelocity[iMarker];
    }
    delete [] AverageTurboVelocity;
  }

  if(OldAverageTurboVelocity != NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] OldAverageTurboVelocity[iMarker][iSpan];
      delete  [] OldAverageTurboVelocity[iMarker];
    }
    delete [] OldAverageTurboVelocity;
  }

  if(ExtAverageTurboVelocity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] ExtAverageTurboVelocity[iMarker][iSpan];
      delete  [] ExtAverageTurboVelocity[iMarker];
    }
    delete [] ExtAverageTurboVelocity;
  }


  if(AverageFlux !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] AverageFlux[iMarker][iSpan];
      delete  [] AverageFlux[iMarker];
    }
    delete [] AverageFlux;
  }

  if(SpanTotalFlux !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
        delete [] SpanTotalFlux[iMarker][iSpan];
      delete  [] SpanTotalFlux[iMarker];
    }
    delete [] SpanTotalFlux;
  }

  if(AveragePressure !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AveragePressure[iMarker];
    delete [] AveragePressure;
  }

  if(OldAveragePressure !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] OldAveragePressure[iMarker];
    delete [] OldAveragePressure;
  }

  if(RadialEquilibriumPressure !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] RadialEquilibriumPressure[iMarker];
    delete [] RadialEquilibriumPressure;
  }

  if(ExtAveragePressure !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAveragePressure[iMarker];
    delete [] ExtAveragePressure;
  }

  if(AverageDensity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AverageDensity[iMarker];
    delete [] AverageDensity;
  }

  if(OldAverageDensity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] OldAverageDensity[iMarker];
    delete [] OldAverageDensity;
  }

  if(ExtAverageDensity !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAverageDensity[iMarker];
    delete [] ExtAverageDensity;
  }

  if(AverageKine !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AverageKine[iMarker];
    delete [] AverageKine;
  }

  if(AverageOmega !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AverageOmega[iMarker];
    delete [] AverageOmega;
  }

  if(AverageNu !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] AverageNu[iMarker];
    delete [] AverageNu;
  }

  if(ExtAverageKine !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAverageKine[iMarker];
    delete [] ExtAverageKine;
  }

  if(ExtAverageOmega !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAverageOmega[iMarker];
    delete [] ExtAverageOmega;
  }

  if(ExtAverageNu !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ExtAverageNu[iMarker];
    delete [] ExtAverageNu;
  }

  if(TurboVelocityIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
      for (iSpan = 0; iSpan < nSpanMax + 1; iSpan++){
        delete [] TurboVelocityIn[iMarker][iSpan];
      }
      delete [] TurboVelocityIn[iMarker];
    }
    delete [] TurboVelocityIn;
  }

  if(TurboVelocityOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
      for (iSpan = 0; iSpan < nSpanMax + 1; iSpan++){
        delete [] TurboVelocityOut[iMarker][iSpan];
      }
      delete [] TurboVelocityOut[iMarker];
    }
    delete [] TurboVelocityOut;
  }

  if(DensityIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] DensityIn[iMarker];
    delete [] DensityIn;
  }

  if(PressureIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] PressureIn[iMarker];
    delete [] PressureIn;
  }

  if(DensityOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] DensityOut[iMarker];
    delete [] DensityOut;
  }

  if(PressureOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] PressureOut[iMarker];
    delete [] PressureOut;
  }

  if(KineIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] KineIn[iMarker];
    delete [] KineIn;
  }

  if(OmegaIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] OmegaIn[iMarker];
    delete [] OmegaIn;
  }

  if(NuIn !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] NuIn[iMarker];
    delete [] NuIn;
  }

  if(KineOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] KineOut[iMarker];
    delete [] KineOut;
  }

  if(OmegaOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] OmegaOut[iMarker];
    delete [] OmegaOut;
  }

  if(NuOut !=NULL){
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
      delete  [] NuOut[iMarker];
    delete [] NuOut;
  }

  if(CkInflow !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
        delete [] CkInflow[iMarker][iSpan];
      delete  [] CkInflow[iMarker];
    }
    delete [] CkInflow;
  }

  if(CkOutflow1 !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
        delete [] CkOutflow1[iMarker][iSpan];
      delete  [] CkOutflow1[iMarker];
    }
    delete [] CkOutflow1;
  }

  if(CkOutflow2 !=NULL){
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
        delete [] CkOutflow2[iMarker][iSpan];
      delete  [] CkOutflow2[iMarker];
    }
    delete [] CkOutflow2;
  }

}

void CImpactSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

  cout << "CImpactSolver Set_MPI_Solution*************************************************************************************" << endl;
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }

        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Solution[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }

}

void CImpactSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

  cout << "CImpactSolver Set_MPI_Solution_Old*************************************************************************************" << endl;
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }

        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Solution[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }
}

void CImpactSolver::Set_MPI_Solution_Air(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

  cout << "CImpactSolver Set_MPI_Solution_Air*************************************************************************************" << endl;
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Air(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }

        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Air(iVar, Solution[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }
}

void CImpactSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;

  cout << "CImpactSolver Set_MPI_Undivided_Laplacian*************************************************************************************" << endl;
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new su2double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex] = Buffer_Send_Undivided_Laplacian[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Undivided_Laplacian;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
        }

        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetUndivided_Laplacian(iVar, Solution[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Undivided_Laplacian;

    }

  }

}

void CImpactSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR, *Buffer_Receive_Neighbor = NULL, *Buffer_Send_Neighbor = NULL;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;

cout << "CImpactSolver Set_MPI_MaxEigenvalue*************************************************************************************" << endl;
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];
      Buffer_Receive_Neighbor = new unsigned short [nBufferR_Vector];
      Buffer_Send_Neighbor = new unsigned short[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetLambda();
        Buffer_Send_Neighbor[iVertex] = geometry->node[iPoint]->GetnPoint();
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_Neighbor, nBufferS_Vector, MPI_UNSIGNED_SHORT, send_to, 1,
                        Buffer_Receive_Neighbor, nBufferR_Vector, MPI_UNSIGNED_SHORT, receive_from, 1, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
        Buffer_Receive_Neighbor[iVertex] = Buffer_Send_Neighbor[iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      delete [] Buffer_Send_Neighbor;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetLambda(Buffer_Receive_Lambda[iVertex]);
        geometry->node[iPoint]->SetnNeighbor(Buffer_Receive_Neighbor[iVertex]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      delete [] Buffer_Receive_Neighbor;

    }

  }
}

void CImpactSolver::Set_MPI_Sensor(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;

  cout << "CImpactSolver Set_MPI_Sensor*************************************************************************************" << endl;
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetSensor();
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetSensor(Buffer_Receive_Lambda[iVertex]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;

    }

  }
}

void CImpactSolver::Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;

  cout << "CImpactSolver Set_MPI_Primitive_Gradient*************************************************************************************" << endl;
  su2double **Gradient = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Gradient[iVar] = new su2double[nDim];

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad*nDim;        nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Primitive(iVar, iDim);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];

        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
        }

        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient_Primitive(iVar, iDim, Gradient[iVar][iDim]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;

    }

  }

  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;

}

void CImpactSolver::Set_MPI_Primitive_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;

  su2double *Limiter = new su2double [nPrimVarGrad];

  cout << "CImpactSolver Set_MPI_Primitive_Limiter*************************************************************************************" << endl;
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad;        nBufferR_Vector = nVertexR*nPrimVarGrad;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter_Primitive(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }

        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          node[iPoint]->SetLimiter_Primitive(iVar, Limiter[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;

    }

  }

  delete [] Limiter;

}

void CImpactSolver::Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) {

  cout << "CImpactSolver Set_MPI_Nearfield*************************************************************************************" << endl;
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0;
  long iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;

#ifdef HAVE_MPI

  /*--- MPI status and request arrays for non-blocking communications ---*/

  SU2_MPI::Status status, status_;


#endif

  /*--- Define buffer vector interior domain ---*/

  su2double        *Buffer_Send_PrimVar          = NULL;

  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  su2double        *iPrimVar          = new su2double [nPrimVar];

  unsigned long Buffer_Size_PrimVar          = 0;

  unsigned long PointTotal_Counter = 0;

  /*--- Allocate the memory that we only need if we have MPI support ---*/

  su2double        *Buffer_Receive_PrimVar          = NULL;

  /*--- Basic dimensionalization ---*/

  nDomain = size;

  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/

  for (iDomain = 0; iDomain < nDomain; iDomain++) {

    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/

    Buffer_Send_nPointTotal = 0;

    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }

    /*--- Store the counts on a partition by partition basis. ---*/

    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;

    /*--- Total counts for allocating send buffers below ---*/

    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(nPrimVar+3);

  }

  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/

  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];

  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/

  for (iDomain = 0; iDomain < nDomain; iDomain++) {

    /*--- A rank does not communicate with itself through MPI ---*/

    if (rank != iDomain) {

#ifdef HAVE_MPI

      /*--- Communicate the counts to iDomain with non-blocking sends ---*/

      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);

#endif

    } else {

      /*--- If iDomain = rank, we simply copy values into place in memory ---*/

      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];

    }

    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/

    if (rank == iDomain) {

      for (jDomain = 0; jDomain < size; jDomain++) {

        /*--- A rank does not communicate with itself through MPI ---*/

        if (rank != jDomain) {

#ifdef HAVE_MPI

          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/

          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);

#endif

        }
      }

    }
  }

  /*--- Wait for the non-blocking sends to complete. ---*/

#ifdef HAVE_MPI

  SU2_MPI::Barrier(MPI_COMM_WORLD);

#endif

  /*--- Initialize the counters for the larger send buffers (by domain) ---*/

  PointTotal_Counter  = 0;

  for (iDomain = 0; iDomain < nDomain; iDomain++) {

    /*--- Set the value of the interior geometry. Initialize counters. ---*/

    iPointTotal = 0;

    /*--- Load up the actual values into the buffers for sending. ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetPrimitive(iVar);
            }
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+0)]  = su2double(iGlobalIndex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = su2double(jVertex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+2)]  = su2double(jMarker);

            iPointTotal++;

          }

        }

      }

    }

    /*--- Send the buffers with the geometrical information ---*/

    if (iDomain != rank) {

#ifdef HAVE_MPI

      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/

      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)],
                     nPointTotal_s[iDomain]*(nPrimVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);

#endif

    }

    else {

      /*--- Allocate local memory for the local recv of the elements ---*/

      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar+3)];

      for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar+3); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)+iter];

      /*--- Recv the point data from ourselves (same procedure as above) ---*/

      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {

        iGlobal       =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];

        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;

        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;

        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);

        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);

      }

      /*--- Delete memory for recv the point stuff ---*/

      delete [] Buffer_Receive_PrimVar;

    }

    /*--- Increment the counters for the send buffers (iDomain loop) ---*/

    PointTotal_Counter += iPointTotal;

  }

  /*--- Wait for the non-blocking sends to complete. ---*/

#ifdef HAVE_MPI

  SU2_MPI::Barrier(MPI_COMM_WORLD);

#endif

  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/

  for (iDomain = 0; iDomain < size; iDomain++) {

    if (rank != iDomain) {

#ifdef HAVE_MPI

      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/

      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar+3)];

      /*--- Receive the buffers with the coords, global index, and colors ---*/

      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar+3) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);

      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/

      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {

        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];

        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;

        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;

        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar,  iPrimVar[iVar]);

        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);

      }

      /*--- Delete memory for recv the point stuff ---*/

      delete [] Buffer_Receive_PrimVar;

#endif

    }

  }

  /*--- Wait for the non-blocking sends to complete. ---*/

#ifdef HAVE_MPI

  SU2_MPI::Barrier(MPI_COMM_WORLD);

#endif

  /*--- Free all of the memory used for communicating points and elements ---*/

  delete[] Buffer_Send_PrimVar;

  /*--- Release all of the temporary memory ---*/

  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;

}

void CImpactSolver::Set_MPI_Interface(CGeometry *geometry, CConfig *config) {

  cout << "CImpactSolver Set_MPI_Interface*************************************************************************************" << endl;
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;

#ifdef HAVE_MPI

  /*--- MPI status and request arrays for non-blocking communications ---*/

  SU2_MPI::Status status, status_;


#endif

  /*--- Define buffer vector interior domain ---*/

  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [nPrimVar];

  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];

  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;

  /*--- Allocate the memory that we only need if we have MPI support ---*/

  su2double        *Buffer_Receive_PrimVar          = NULL;

  /*--- Basic dimensionalization ---*/

  nDomain = size;

  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/

  for (iDomain = 0; iDomain < nDomain; iDomain++) {

    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/

    Buffer_Send_nPointTotal = 0;

    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }

    /*--- Store the counts on a partition by partition basis. ---*/

    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;

    /*--- Total counts for allocating send buffers below ---*/

    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(nPrimVar+3);

  }

  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/

  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];

  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/

  for (iDomain = 0; iDomain < nDomain; iDomain++) {

    /*--- A rank does not communicate with itself through MPI ---*/

    if (rank != iDomain) {

#ifdef HAVE_MPI

      /*--- Communicate the counts to iDomain with non-blocking sends ---*/

      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);

#endif

    } else {

      /*--- If iDomain = rank, we simply copy values into place in memory ---*/

      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];

    }

    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/

    if (rank == iDomain) {

      for (jDomain = 0; jDomain < size; jDomain++) {

        /*--- A rank does not communicate with itself through MPI ---*/

        if (rank != jDomain) {

#ifdef HAVE_MPI

          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/

          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);

#endif

        }
      }

    }
  }

  /*--- Wait for the non-blocking sends to complete. ---*/

#ifdef HAVE_MPI

  SU2_MPI::Barrier(MPI_COMM_WORLD);

#endif

  /*--- Initialize the counters for the larger send buffers (by domain) ---*/

  PointTotal_Counter  = 0;

  for (iDomain = 0; iDomain < nDomain; iDomain++) {

    /*--- Set the value of the interior geometry. Initialize counters. ---*/

    iPointTotal = 0;

    /*--- Load up the actual values into the buffers for sending. ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();

          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {

            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();

            for (iVar = 0; iVar < nPrimVar; iVar++) {
              Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetPrimitive(iVar);
            }
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+0)]  = su2double(iGlobalIndex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = su2double(jVertex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+2)]  = su2double(jMarker);

            iPointTotal++;

          }

        }

      }

    }

    /*--- Send the buffers with the geometrical information ---*/

    if (iDomain != rank) {

#ifdef HAVE_MPI

      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/

      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)],
                     nPointTotal_s[iDomain]*(nPrimVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);

#endif

    }

    else {

      /*--- Allocate local memory for the local recv of the elements ---*/

      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar+3)];

      for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar+3); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)+iter];

      /*--- Recv the point data from ourselves (same procedure as above) ---*/

      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {

        iGlobal       =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];

        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;

        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);

        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);

      }

      /*--- Delete memory for recv the point stuff ---*/

      delete [] Buffer_Receive_PrimVar;

    }

    /*--- Increment the counters for the send buffers (iDomain loop) ---*/

    PointTotal_Counter += iPointTotal;

  }

  /*--- Wait for the non-blocking sends to complete. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/

  for (iDomain = 0; iDomain < size; iDomain++) {

    if (rank != iDomain) {

#ifdef HAVE_MPI

      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/

      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar+3)];

      /*--- Receive the buffers with the coords, global index, and colors ---*/

      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar+3) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);

      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/

      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {

        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];

        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;

        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;

        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);

        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);

      }

      /*--- Delete memory for recv the point stuff ---*/

      delete [] Buffer_Receive_PrimVar;

#endif

    }

  }

  /*--- Wait for the non-blocking sends to complete. ---*/

#ifdef HAVE_MPI

  SU2_MPI::Barrier(MPI_COMM_WORLD);

#endif

  /*--- Free all of the memory used for communicating points and elements ---*/

  delete[] Buffer_Send_PrimVar;

  /*--- Release all of the temporary memory ---*/

  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;

}

void CImpactSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) {

  cout << "CImpactSolver SetNondimensionalization*************************************************************************************" << endl;
  su2double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0,
  Energy_FreeStream = 0.0, ModVel_FreeStreamND = 0.0, Velocity_Reynolds = 0.0,
  Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
  Temperature_Ref = 0.0, Time_Ref = 0.0, Omega_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Energy_Ref= 0.0,
  Froude = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0, TgammaR = 0.0;

  unsigned short iDim;

  /*--- Local variables ---*/

  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double Reynolds         = config->GetReynolds();
  bool unsteady           = (config->GetUnsteady_Simulation() != NO);
  bool grid_movement      = config->GetGrid_Movement();
  bool gravity            = config->GetGravityForce();
  bool turbulent          = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
  bool tkeNeeded          = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);
  bool aeroelastic        = config->GetAeroelastic_Simulation();

  /*--- Set temperature via the flutter speed index ---*/
  if (aeroelastic) {
    su2double vf             = config->GetAeroelastic_Flutter_Speed_Index();
    su2double w_alpha        = config->GetAeroelastic_Frequency_Pitch();
    su2double b              = config->GetLength_Reynolds()/2.0; // airfoil semichord, Reynolds length is by defaul 1.0
    su2double mu             = config->GetAeroelastic_Airfoil_Mass_Ratio();
    // The temperature times gamma times the gas constant. Depending on the FluidModel temp is calculated below.
    TgammaR = ((vf*vf)*(b*b)*(w_alpha*w_alpha)*mu) / (Mach*Mach);
  }

  /*--- Compressible non dimensionalization ---*/

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);

      FluidModel = new CIdealGas(1.4, config->GetGas_Constant());
      if (free_stream_temp) {
        if (aeroelastic) {
          Temperature_FreeStream = TgammaR / (config->GetGas_Constant()*1.4);
          config->SetTemperature_FreeStream(Temperature_FreeStream);
        }
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case IDEAL_GAS:

      FluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case VW_GAS:

      FluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                       config->GetPressure_Critical(), config->GetTemperature_Critical());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case PR_GAS:

      FluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                     config->GetTemperature_Critical(), config->GetAcentric_Factor());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

  }

  Mach2Vel_FreeStream = FluidModel->GetSoundSpeed();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }

  /*--- Compute the modulus of the free stream velocity ---*/

  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- For inviscid flow, energy is calculated from the specified
  FreeStream quantities using the proper gas law. ---*/

  Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

  /*-- Compute the freestream energy. ---*/

  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute non dimensional quantities. By definition,
     Lref is one because we have converted the grid to meters. ---*/

  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
    Density_Ref       = 1.0;
    Temperature_Ref   = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
    Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDARD_GRAVITY*Length_Ref);         config->SetFroude(Froude);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/

  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);

  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);


  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);

  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);

  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/(Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStream(Omega_FreeStream);

  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/(Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);

  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/

  /*--- Delete the original (dimensional) FluidModel object before replacing. ---*/

  delete FluidModel;

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:
      FluidModel = new CIdealGas(1.4, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case IDEAL_GAS:
      FluidModel = new CIdealGas(Gamma, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case VW_GAS:
      FluidModel = new CVanDerWaalsGas(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                       config->GetTemperature_Critical()/config->GetTemperature_Ref());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case PR_GAS:
      FluidModel = new CPengRobinson(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                     config->GetTemperature_Critical()/config->GetTemperature_Ref(), config->GetAcentric_Factor());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

  }

  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is the master node and first domain ---*/

  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {

    cout.precision(6);

    cout << "Inviscid flow: Computing density based on free-stream" << endl;
    cout << "temperature and pressure using the ideal gas law." << endl;

    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;

    cout <<"-- Input conditions:"<< endl;

    switch (config->GetKind_FluidModel()) {

      case STANDARD_AIR:
        cout << "Fluid Model: STANDARD_AIR "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant();
        if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;

      case IDEAL_GAS:
        cout << "Fluid Model: IDEAL_GAS "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;

      case VW_GAS:
        cout << "Fluid Model: Van der Waals "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;

      case PR_GAS:
        cout << "Fluid Model: Peng-Robinson "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;

    }

    cout << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Free-stream total pressure: " << config->GetPressure_FreeStream() * pow( 1.0+Mach*Mach*0.5*(Gamma-1.0), Gamma/(Gamma-1.0) );
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Free-stream temperature: " << config->GetTemperature_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;

    cout << "Free-stream density: " << config->GetDensity_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;

    if (nDim == 2) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
    }
    if (config->GetSystemMeasurements() == SI) cout << " m/s. ";
    else if (config->GetSystemMeasurements() == US) cout << " ft/s. ";

    cout << "Magnitude: "   << config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m/s (" << config->GetModVel_FreeStream()*1.94384 << " KTS)." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s (" << config->GetModVel_FreeStream()*0.592484 << " KTS)." << endl;

    cout << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;

    if (unsteady) { cout << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << endl; }

    /*--- Print out reference values. ---*/

    cout <<"-- Reference values:"<< endl;

    cout << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;

    cout << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Reference temperature: " << config->GetTemperature_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;

    cout << "Reference density: " << config->GetDensity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;

    cout << "Reference velocity: " << config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;

    cout << "Reference energy per unit mass: " << config->GetEnergy_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;


    if (unsteady) cout << "Reference time: " << config->GetTime_Ref() <<" s." << endl;

    /*--- Print out resulting non-dim values here. ---*/

    cout << "-- Resulting non-dimensional state:" << endl;
    cout << "Mach number (non-dim): " << config->GetMach() << endl;

    if (gravity) {
      cout << "Froude number (non-dim): " << Froude << endl;
      cout << "Lenght of the baseline wave (non-dim): " << 2.0*PI_NUMBER*Froude*Froude << endl;
    }

    cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
    cout << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << endl;

    cout << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << endl;

    cout << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << endl;

    if (nDim == 2) {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    cout << "Magnitude: "    << config->GetModVel_FreeStreamND() << endl;

    cout << "Free-stream total energy per unit mass (non-dim): " << config->GetEnergy_FreeStreamND() << endl;

    if (unsteady) {
      cout << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << endl;
      cout << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << endl;
    }

    cout << endl;

  }

}

void CImpactSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

  cout << "CImpactSolver SetInitialCondition*************************************************************************************" << endl;
  unsigned long iPoint;
  unsigned short iMesh, iDim;
  su2double X0[3] = {0.0,0.0,0.0}, X1[3] = {0.0,0.0,0.0}, X2[3] = {0.0,0.0,0.0},
  X1_X0[3] = {0.0,0.0,0.0}, X2_X0[3] = {0.0,0.0,0.0}, X2_X1[3] = {0.0,0.0,0.0},
  CP[3] = {0.0,0.0,0.0}, Distance, DotCheck, Radius, AirVelocity, AirDensity;

  unsigned short nDim = geometry[MESH_0]->GetnDim();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == RANS) ||
               (config->GetKind_Solver() == ADJ_RANS) ||
               (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool SubsonicEngine = config->GetSubsonicEngine();

//  /*--- If no restart file, read the air solution file and use it to init impact model ---*/
  if (dual_time && (ExtIter == 0 || (restart && (long)ExtIter == config->GetUnst_RestartIter()))) {
    cout << "read air velocity" << endl;
    if (!restart) {
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1), true);
  /*--- Get the air velocity solution as initial velocity for impact solution ---*/
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          AirDensity = solver_container[iMesh][FLOW_SOL]->node[iPoint]->GetSolution(0);
          solver_container[iMesh][IMPACT_SOL]->node[iPoint]->SetSolution(0,AirDensity);
          for (iDim = 0; iDim < nDim; iDim++) {
            AirVelocity = solver_container[iMesh][FLOW_SOL]->node[iPoint]->GetSolution(iDim+1)/AirDensity;
            solver_container[iMesh][IMPACT_SOL]->node[iPoint]->SetSolution(iDim+1,AirVelocity);
          }
        }
      }
    }


  /*--- The primitive air solution variable file information is stored for drag calculation ---*/
      su2double Viscosity, Velocity2, StaticEnergy, Pressure, Temperature;
      FluidModel->SetLaminarViscosityModel(config);
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          AirDensity = solver_container[iMesh][FLOW_SOL]->node[iPoint]->GetSolution(0);
//        AirDensity = 1.0;
          node[iPoint]->SetDensityAir(AirDensity);
          Velocity2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            AirVelocity = solver_container[iMesh][FLOW_SOL]->node[iPoint]->GetSolution(iDim+1)/AirDensity;
//          AirVelocity = 0.0;
            Velocity2 += AirVelocity*AirVelocity;
            node[iPoint]->SetVelocityAir(iDim,AirVelocity);
          }
        StaticEnergy = solver_container[iMesh][FLOW_SOL]->node[iPoint]->GetSolution(nDim+1)/AirDensity - 0.5*Velocity2;
//      StaticEnergy = 200000.0;
        FluidModel->SetTDState_rhoe(AirDensity, StaticEnergy);
        Pressure = FluidModel->GetPressure();
        Temperature = FluidModel->GetTemperature();
        node[iPoint]->SetTemperatureAir(Temperature);
        Viscosity = FluidModel->GetLaminarViscosity();
        node[iPoint]->SetViscosityAir(Viscosity);
      }
    }
  }
  /*--- Make sure that the solution is well initialized for unsteady
   calculations with dual time-stepping (load additional restarts for 2nd-order). ---*/
  /*--- TO DO Unsteady calculation for IMPACT need to read additional airflow solution ---*/

  if (dual_time && (ExtIter == 0 || (restart && (long)ExtIter == config->GetUnst_RestartIter()))) {

    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/

    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        solver_container[iMesh][IMPACT_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][IMPACT_SOL]->node[iPoint]->Set_Solution_time_n1();
      }
    }

    if ((restart && (long)ExtIter == config->GetUnst_RestartIter()) &&
        (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

      /*--- Load an additional restart file for a 2nd-order restart ---*/
      solver_container[MESH_0][IMPACT_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1), true);

      /*--- Push back this new solution to time level N. ---*/

      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          solver_container[iMesh][IMPACT_SOL]->node[iPoint]->Set_Solution_time_n();
          }
        }
      }
    }
  }


void CImpactSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  cout << endl <<"CImpactSolver Preprocessing*************************************************************************************" << endl;
  unsigned long ErrorCounter = 0;

  unsigned long ExtIter = config->GetExtIter();
  bool cont_adjoint     = config->GetContinuous_Adjoint();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
  bool limiter          = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst       = center && (config->GetKind_Centered_Flow() == JST);
  bool engine           = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk    = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool nearfield        = (config->GetnMarker_NearFieldBound() != 0);
  bool interface        = (config->GetnMarker_InterfaceBound() != 0);
  bool fixed_cl         = config->GetFixed_CL_Mode();
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  unsigned short kind_row_dissipation = config->GetKind_RoeLowDiss();
  bool roe_low_dissipation  = (kind_row_dissipation != NO_ROELOWDISS) && (config->GetKind_Upwind_Flow() == ROE);

  /*--- Set the primitive variables ---*/

  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);


  /*--- Compute Interface MPI ---*/

  if (interface) { Set_MPI_Interface(geometry, config); }

  /*--- Compute NearField MPI ---*/

  if (nearfield) { Set_MPI_Nearfield(geometry, config); }


  /*--- Upwind second order reconstruction ---*/

  if ((muscl && !center) && (iMesh == MESH_0) && !Output) {

    /*--- Gradient computation ---*/

    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config);
    }

    /*--- Limiter computation ---*/

    if (limiter && (iMesh == MESH_0)
        && !Output && !van_albada) { SetPrimitive_Limiter(geometry, config); }

  }

  /*--- Artificial dissipation ---*/

  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Roe Low Dissipation Sensor ---*/

  if (roe_low_dissipation){
    SetRoe_Dissipation(geometry, config);
  }

  /*--- Initialize the Jacobian matrices ---*/

  if (implicit && !disc_adjoint) Jacobian.SetValZero();

  /*--- Error message ---*/

  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }

}

void CImpactSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) {cout << "CImpactSolver Postprocessing*************************************************************************************" << endl; }

unsigned long CImpactSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {

  cout << "CImpactSolver SetPrimitive_Variables*************************************************************************************" << endl;
  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the non-physical points vector ---*/

    node[iPoint]->SetNon_Physical(false);

    /*--- impact flow, primitive variables 2*nDim+12, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/

    RightSol = node[iPoint]->SetPrimVar(FluidModel);
    node[iPoint]->SetSecondaryVar(FluidModel);


    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  return ErrorCounter;
}

void CImpactSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {

  cout << "CImpactSolver SetTime_Step*************************************************************************************" << endl;
  su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
    bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetMax_Lambda_Inv(0.0);

  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

    /*--- Mean Values ---*/

    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;

    /*--- Adjustment for grid movement ---*/

    if (grid_movement) {
      su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }

    /*--- Inviscid contribution ---*/

    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

      /*--- Mean Values ---*/

      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

      /*--- Adjustment for grid movement ---*/

      if (grid_movement) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }

      /*--- Inviscid contribution ---*/
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }

    }
  }

  /*--- Each element uses their own speed, steady state simulation ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->node[iPoint]->GetVolume();

    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }

  }


  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;

    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }

  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/

  if (time_steping) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

            /*--- Sets the regular CFL equal to the unsteady CFL ---*/
            config->SetCFL(iMesh,config->GetUnst_CFL());

            /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
             it computes the time step based on the unsteady CFL ---*/
            if (config->GetCFL(iMesh) == 0.0) {
                node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
            } else {
                node[iPoint]->SetDelta_Time(Global_Delta_Time);
            }
        }
  }

  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/

  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }

  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/

  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }

}

void CImpactSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  cout << "CImpactSolver Centered_Residual*************************************************************************************" << endl;
  unsigned long iEdge, iPoint, jPoint;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool jst_scheme = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));
  bool grid_movement = config->GetGrid_Movement();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge, set normal vectors, and number of neighbors ---*/

    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());

    /*--- Set primitive variables w/o reconstruction ---*/

    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());

    /*--- Set the largest convective eigenvalue ---*/

    numerics->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());

    /*--- Set undivided laplacian an pressure based sensor ---*/

    if (jst_scheme) {
      numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
      numerics->SetSensor(node[iPoint]->GetSensor(), node[jPoint]->GetSensor());
    }

    /*--- Grid movement ---*/

    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }

    /*--- Compute residuals, and Jacobians ---*/

    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

    /*--- Update convective and artificial dissipation residuals ---*/

    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);

    /*--- Set implicit computation ---*/
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
  }

}

void CImpactSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                   CConfig *config, unsigned short iMesh) {

  cout << "CImpactSolver Source_Residual*************************************************************************************" << endl;
  unsigned short iVar;
  unsigned long iPoint;
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool rotating_frame   = config->GetRotating_Frame();
  bool axisymmetric     = config->GetAxisymmetric();
  bool gravity          = (config->GetGravityForce() == YES);
  bool harmonic_balance = (config->GetUnsteady_Simulation() == HARMONIC_BALANCE);
  bool windgust         = config->GetWind_Gust();
  bool body_force       = config->GetBody_Force();

  /*--- Initialize the source residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

  /*--- compute the drag force for the droplet equation ---*/
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());

      /*--- Load the air variable ---*/
      numerics->SetPrimitiveAir(node[iPoint]->GetSolution_Air(),
                                node[iPoint]->GetSolution_Air());

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);

    }

  if (body_force) {

    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);

    }
  }

  if (rotating_frame) {

    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

    }
  }

  if (axisymmetric) {

    /*--- Zero out Jacobian structure ---*/
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }

    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Set solution  ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Set y coordinate ---*/
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());

      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);

      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Implicit part ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }

  if (gravity) {

    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Set solution  ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, config);

      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);

    }

  }

  if (harmonic_balance) {

    su2double Volume, Source;

    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Get control volume ---*/
      Volume = geometry->node[iPoint]->GetVolume();

      /*--- Get stored time spectral source term ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Source = node[iPoint]->GetHarmonicBalance_Source(iVar);
        Residual[iVar] = Source*Volume;
      }

      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);

    }
  }

  if (windgust) {

    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the wind gust ---*/
      numerics->SetWindGust(node[iPoint]->GetWindGust(), node[iPoint]->GetWindGust());

      /*--- Load the wind gust derivatives ---*/
      numerics->SetWindGustDer(node[iPoint]->GetWindGustDer(), node[iPoint]->GetWindGustDer());

      /*--- Load the primitive variables ---*/
      numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[iPoint]->GetPrimitive());

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

    }
  }

}

void CImpactSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {

  cout << "CImpactSolver Source_Template*************************************************************************************" << endl;
  /* This method should be used to call any new source terms for a particular problem*/
  /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */

  /* Next we describe how to get access to some important quanties for this method */
  /* Access to all points in the current geometric mesh by saying: nPointDomain */
  /* Get the vector of conservative variables at some point iPoint = node[iPoint]->GetSolution() */
  /* Get the volume (or area in 2D) associated with iPoint = node[iPoint]->GetVolume() */
  /* Get the vector of geometric coordinates of point iPoint = node[iPoint]->GetCoord() */

}

void CImpactSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {

  cout << "CImpactSolver SetMax_Eigenvalue*************************************************************************************" << endl;
  su2double *Normal, Area, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda,
  ProjVel, ProjVel_i, ProjVel_j, *GridVel, *GridVel_i, *GridVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  bool grid_movement = config->GetGrid_Movement();

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetLambda(0.0);
  }

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

    /*--- Mean Values ---*/

    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;

    /*--- Adjustment for grid movement ---*/

    if (grid_movement) {
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }

    /*--- Inviscid contribution ---*/

    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

      /*--- Mean Values ---*/

      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

      /*--- Adjustment for grid movement ---*/

      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }

      /*--- Inviscid contribution ---*/

      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddLambda(Lambda);
      }

    }
  }

  /*--- MPI parallelization ---*/

  Set_MPI_MaxEigenvalue(geometry, config);

}

void CImpactSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {

  cout << "CImpactSolver SetUndivided_Laplacian*************************************************************************************" << endl;
  unsigned long iPoint, jPoint, iEdge;
  su2double Pressure_i = 0, Pressure_j = 0, *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;

  Diff = new su2double[nVar];

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetUnd_LaplZero();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Solution differences ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);

    /*--- Correction for compressible flows which use the enthalpy ---*/

    Pressure_i = node[iPoint]->GetPressure();
    Pressure_j = node[jPoint]->GetPressure();
    Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + Pressure_i) - (node[jPoint]->GetSolution(nVar-1) + Pressure_j);

    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();

    /*--- Both points inside the domain, or both in the boundary ---*/

    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    }

    /*--- iPoint inside the domain, jPoint on the boundary ---*/

    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);

    /*--- jPoint inside the domain, iPoint on the boundary ---*/

    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);

  }

  /*--- MPI parallelization ---*/

  Set_MPI_Undivided_Laplacian(geometry, config);

  delete [] Diff;

}

void CImpactSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) {

  cout << "CImpactSolver SetCentered_Dissipation_Sensor*************************************************************************************" << endl;
  unsigned long iEdge, iPoint, jPoint;
  su2double Pressure_i = 0.0, Pressure_j = 0.0;
  bool boundary_i, boundary_j;

  /*--- Reset variables to store the undivided pressure ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    iPoint_UndLapl[iPoint] = 0.0;
    jPoint_UndLapl[iPoint] = 0.0;
  }

  /*--- Evaluate the pressure sensor ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    Pressure_i = node[iPoint]->GetPressure();
    Pressure_j = node[jPoint]->GetPressure();

    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();

    /*--- Both points inside the domain, or both on the boundary ---*/

    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    }

    /*--- iPoint inside the domain, jPoint on the boundary ---*/

    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }

    /*--- jPoint inside the domain, iPoint on the boundary ---*/

    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }

  }

  /*--- Set pressure switch for each point ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetSensor(fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);

  /*--- MPI parallelization ---*/

  Set_MPI_Sensor(geometry, config);

}

void CImpactSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  cout << "CImpactSolver ExplicitRK_Iteration*************************************************************************************" << endl;
  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;

  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }

  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);


}

void CImpactSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  cout << "CImpactSolver ClassicalRK4_Iteration*************************************************************************************" << endl;
  su2double *Residual, *Res_TruncError, Vol, Delta, Res, tmp_time, tmp_func;
  unsigned short iVar;
  unsigned long iPoint;

  /*--- Hard-coded classical RK4 coefficients. Will be added to config. ---*/
  su2double RK_FuncCoeff[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
  su2double RK_TimeCoeff[4] = {0.5, 0.5, 1.0, 1.0};

  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);

    tmp_time = -1.0*RK_TimeCoeff[iRKStep]*Delta;
    tmp_func = -1.0*RK_FuncCoeff[iRKStep]*Delta;

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        if (iRKStep < 3) {
          /* Base Solution Update */
          node[iPoint]->AddSolution(iVar, tmp_time*Res);

          /* New Solution Update */
          node[iPoint]->AddSolution_New(iVar, tmp_func*Res);
        } else {
          node[iPoint]->SetSolution(iVar, node[iPoint]->GetSolution_New(iVar) + tmp_func*Res);
        }

        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }

  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

void CImpactSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  cout << "CImpactSolver ExplicitEuler_Iteration*************************************************************************************" << endl;
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;

  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;

    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = local_Residual[iVar] + local_Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }

  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

void CImpactSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  cout << "CImpactSolver ImplicitEuler_Iteration*************************************************************************************" << endl;
  unsigned short iVar, jVar;
  unsigned long iPoint, total_index, IterLinSol = 0;
  su2double Delta, *local_Res_TruncError, Vol;

  bool adjoint = config->GetContinuous_Adjoint();
  bool roe_turkel = config->GetKind_Upwind_Flow() == TURKEL;
  bool low_mach_prec = config->Low_Mach_Preconditioning();

  /*--- Set maximum residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the residual ---*/

    local_Res_TruncError = node[iPoint]->GetResTruncError();

    /*--- Read the volume ---*/

    Vol = geometry->node[iPoint]->GetVolume();

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/


    if (node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / node[iPoint]->GetDelta_Time();
      if (roe_turkel || low_mach_prec) {
        SetPreconditioner(config, iPoint);
        for (iVar = 0; iVar < nVar; iVar ++ )
          for (jVar = 0; jVar < nVar; jVar ++ )
            LowMach_Precontioner[iVar][jVar] = Delta*LowMach_Precontioner[iVar][jVar];
        Jacobian.AddBlock(iPoint, iPoint, LowMach_Precontioner);
      }
      else {
        Jacobian.AddVal2Diag(iPoint, Delta);
      }
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- Initialize residual and solution at the ghost points ---*/

  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }

  /*--- Solve or smooth the linear system ---*/

  CSysSolve system;
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  /*--- The the number of iterations of the linear solver ---*/

  SetIterLinSolver(IterLinSol);

  /*--- Update solution (system written in terms of increments) ---*/

  if (!adjoint) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        node[iPoint]->AddSolution(iVar, config->GetRelaxation_Factor_Flow()*LinSysSol[iPoint*nVar+iVar]);
      }
    }
  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

void CImpactSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) {
  cout << "CImpactSolver SetPrimitive_Gradient_GG*************************************************************************************" << endl;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  su2double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
  Partial_Gradient, Partial_Res, *Normal;

  /*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho) ---*/

  PrimVar_Vertex = new su2double [nPrimVarGrad];
  PrimVar_i = new su2double [nPrimVarGrad];
  PrimVar_j = new su2double [nPrimVarGrad];

  /*--- Set Gradient_Primitive to zero ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
      PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);
    }

    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = PrimVar_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
      }
    }
  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()) {

        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          PrimVar_Vertex[iVar] = node[iPoint]->GetPrimitive(iVar);

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++) {
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
            node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
          }
      }
    }
  }

  /*--- Update gradient value ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar, iDim) / (geometry->node[iPoint]->GetVolume());
        node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
      }
    }
  }

  delete [] PrimVar_Vertex;
  delete [] PrimVar_i;
  delete [] PrimVar_j;

  Set_MPI_Primitive_Gradient(geometry, config);

}

void CImpactSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) {
  cout << "CImpactSolver SetPrimitive_Gradient_LS*************************************************************************************" << endl;
  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular;

  /*--- Loop over points of the grid ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Set the value of the singular ---*/
    singular = false;

    /*--- Get coordinates ---*/

    Coord_i = geometry->node[iPoint]->GetCoord();

    /*--- Get primitives from CVariable ---*/

    PrimVar_i = node[iPoint]->GetPrimitive();

    /*--- Inizialization of variables ---*/

    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;

    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;

    AD::StartPreacc();
    AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    AD::SetPreaccIn(Coord_i, nDim);

    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();

      PrimVar_j = node[jPoint]->GetPrimitive();

      AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

      /*--- Sumations for entries of upper triangular matrix R ---*/

      if (weight != 0.0) {

        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;

        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }

        /*--- Entries of c:= transpose(A)*b ---*/

        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;

      }

    }

    /*--- Entries of upper triangular matrix R ---*/

    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;

    if (nDim == 3) {
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }

    /*--- Compute determinant ---*/

    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);

    /*--- Detect singular matrices ---*/

    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }

    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/

    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }

    /*--- Computation of the gradient: S*c ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        }

        node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
      }
    }

    AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
    AD::EndPreacc();
  }

  Set_MPI_Primitive_Gradient(geometry, config);

}

void CImpactSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) {
  cout << "CImpactSolver SetPrimitive_Limiter*************************************************************************************" << endl;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j,
  *Primitive, *Primitive_i, *Primitive_j, *LocalMinPrimitive, *LocalMaxPrimitive,
  *GlobalMinPrimitive, *GlobalMaxPrimitive,
  dave, LimK, eps2, eps1, dm, dp, du, y, limiter;

  dave = config->GetRefElemLength();
  LimK = config->GetVenkat_LimiterCoeff();

  if (config->GetKind_SlopeLimit_Flow() == NO_LIMITER) {

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        node[iPoint]->SetLimiter_Primitive(iVar, 1.0);
      }
    }

  }

  else {

    /*--- Initialize solution max and solution min and the limiter in the entire domain --*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        node[iPoint]->SetSolution_Max(iVar, -EPS);
        node[iPoint]->SetSolution_Min(iVar, EPS);
        node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
      }
    }

    /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);

      /*--- Get the primitive variables ---*/

      Primitive_i = node[iPoint]->GetPrimitive();
      Primitive_j = node[jPoint]->GetPrimitive();

      /*--- Compute the maximum, and minimum values for nodes i & j ---*/

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        du = (Primitive_j[iVar] - Primitive_i[iVar]);
        node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
        node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
        node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
        node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
      }

    }

  }

  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/

  if (config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();

      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {

        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/

        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
          else dp = node[iPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }

        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }

        /*--- Calculate the interface right gradient, delta+ (dp) ---*/

        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
          else dp = node[jPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }

        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }

      }

      AD::EndPreacc();

    }

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        y =  node[iPoint]->GetLimiter_Primitive(iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
      }
    }

  }

  /*--- Venkatakrishnan limiter ---*/

  if ((config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) ||
      (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG)) {

    /*--- Allocate memory for the max and min primitive value --*/

    LocalMinPrimitive = new su2double [nPrimVarGrad]; GlobalMinPrimitive = new su2double [nPrimVarGrad];
    LocalMaxPrimitive = new su2double [nPrimVarGrad]; GlobalMaxPrimitive = new su2double [nPrimVarGrad];

    /*--- Compute the max value and min value of the solution ---*/

    Primitive = node[0]->GetPrimitive();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      LocalMinPrimitive[iVar] = Primitive[iVar];
      LocalMaxPrimitive[iVar] = Primitive[iVar];
    }

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

      /*--- Get the primitive variables ---*/

      Primitive = node[iPoint]->GetPrimitive();

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        LocalMinPrimitive[iVar] = min (LocalMinPrimitive[iVar], Primitive[iVar]);
        LocalMaxPrimitive[iVar] = max (LocalMaxPrimitive[iVar], Primitive[iVar]);
      }

    }

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(LocalMinPrimitive, GlobalMinPrimitive, nPrimVarGrad, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(LocalMaxPrimitive, GlobalMaxPrimitive, nPrimVarGrad, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      GlobalMinPrimitive[iVar] = LocalMinPrimitive[iVar];
      GlobalMaxPrimitive[iVar] = LocalMaxPrimitive[iVar];
    }
#endif

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();

      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(eps2);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {

        if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {
          eps1 = LimK * (GlobalMaxPrimitive[iVar] - GlobalMinPrimitive[iVar]);
          eps2 = eps1*eps1;
        }
        else {
          eps1 = LimK*dave;
          eps2 = eps1*eps1*eps1;
        }

        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/

        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

        /*--- Calculate the interface right gradient, delta+ (dp) ---*/

        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
        else dp = node[iPoint]->GetSolution_Min(iVar);

        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }

        /*-- Repeat for point j on the edge ---*/

        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);

        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }

      }

      AD::EndPreacc();

    }

    delete [] LocalMinPrimitive; delete [] GlobalMinPrimitive;
    delete [] LocalMaxPrimitive; delete [] GlobalMaxPrimitive;

  }

  /*--- Limiter MPI ---*/

  Set_MPI_Primitive_Limiter(geometry, config);

}

void CImpactSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) {
  cout << "CImpactSolver SetPreconditioner*************************************************************************************" << endl;
  unsigned short iDim, jDim, iVar, jVar;
  su2double local_Mach, rho, enthalpy, soundspeed, sq_vel;
  su2double *U_i = NULL;
  //su2double *Uair_i = NULL;
  su2double Beta_max = config->GetmaxTurkelBeta();
  su2double Mach_infty2, Mach_lim2, aux, parameter;

  /*--- Variables to calculate the preconditioner parameter Beta ---*/
  local_Mach = sqrt(node[iPoint]->GetVelocity2())/node[iPoint]->GetSoundSpeed();

  /*--- Weiss and Smith Preconditioning---*/
  Mach_infty2 = pow(config->GetMach(),2.0);
  Mach_lim2 = pow(0.00001,2.0);
  aux = max(pow(local_Mach,2.0),Mach_lim2);
  parameter = min(1.0, max(aux,Beta_max*Mach_infty2));

  U_i = node[iPoint]->GetSolution();
  //Uair_i = node[iPoint]->GetSolutionAir();

  rho = U_i[0];
  enthalpy = node[iPoint]->GetEnthalpy();
  soundspeed = node[iPoint]->GetSoundSpeed();
  sq_vel = node[iPoint]->GetVelocity2();

  /*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
  LowMach_Precontioner[0][0] = 0.5*sq_vel;
  LowMach_Precontioner[0][nVar-1] = 1.0;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[0][1+iDim] = -1.0*U_i[iDim+1]/rho;

  for (iDim = 0; iDim < nDim; iDim ++) {
    LowMach_Precontioner[iDim+1][0] = 0.5*sq_vel*U_i[iDim+1]/rho;
    LowMach_Precontioner[iDim+1][nVar-1] = U_i[iDim+1]/rho;
    for (jDim = 0; jDim < nDim; jDim ++) {
      LowMach_Precontioner[iDim+1][1+jDim] = -1.0*U_i[jDim+1]/rho*U_i[iDim+1]/rho;
    }
  }

  LowMach_Precontioner[nVar-1][0] = 0.5*sq_vel*enthalpy;
  LowMach_Precontioner[nVar-1][nVar-1] = enthalpy;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[nVar-1][1+iDim] = -1.0*U_i[iDim+1]/rho*enthalpy;


  for (iVar = 0; iVar < nVar; iVar ++ ) {
    for (jVar = 0; jVar < nVar; jVar ++ ) {
      LowMach_Precontioner[iVar][jVar] = (parameter - 1.0) * ((Gamma-1.0)/(soundspeed*soundspeed))*LowMach_Precontioner[iVar][jVar];
      if (iVar == jVar)
        LowMach_Precontioner[iVar][iVar] += 1.0;
    }
  }

}


void CImpactSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  cout << "CImpactSolver BC_Euler_Wall*************************************************************************************" << endl;
  unsigned short iDim, iVar, jVar, kVar, jDim;
  unsigned long iPoint, iVertex;
  su2double *Normal = NULL, *GridVel = NULL, Area, UnitNormal[3], *NormalArea,
  ProjGridVel = 0.0, turb_ke;
  su2double Density_b, StaticEnergy_b, Enthalpy_b, *Velocity_b, Kappa_b, Chi_b, Energy_b, VelMagnitude2_b, Pressure_b;
  su2double Density_i, *Velocity_i, ProjVelocity_i = 0.0, Energy_i, VelMagnitude2_i;
  su2double **Jacobian_b, **DubDu;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));

  Normal = new su2double[nDim];
  NormalArea = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_i = new su2double[nDim];
  Jacobian_b = new su2double*[nVar];
  DubDu = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_b[iVar] = new su2double[nVar];
    DubDu[iVar] = new su2double[nVar];
  }

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negative for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++) {
        NormalArea[iDim] = -Normal[iDim];
        UnitNormal[iDim] = -Normal[iDim]/Area;
      }

      /*--- Get the state i ---*/

      VelMagnitude2_i = 0.0; ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
        VelMagnitude2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      Density_i = node[iPoint]->GetDensity();
      Energy_i = node[iPoint]->GetEnergy();

      /*--- Compute the boundary state b ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        Velocity_b[iDim] = Velocity_i[iDim] - ProjVelocity_i * UnitNormal[iDim]; //Force the velocity to be tangential to the surface.

      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        for (iDim = 0; iDim < nDim; iDim++) Velocity_b[iDim] += GridVel[iDim] - ProjGridVel * UnitNormal[iDim];
      }

      VelMagnitude2_b = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        VelMagnitude2_b += Velocity_b[iDim] * Velocity_b[iDim];

      /*--- Compute the residual ---*/

      turb_ke = 0.0;
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);

      Density_b = Density_i;
      StaticEnergy_b = Energy_i - 0.5 * VelMagnitude2_i - turb_ke;
      Energy_b = StaticEnergy_b + 0.5 * VelMagnitude2_b + turb_ke;

      FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);
      Kappa_b = FluidModel->GetdPde_rho() / Density_b;
      Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
      Pressure_b = FluidModel->GetPressure();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;

      numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, NormalArea, Residual);

      /*--- Grid velocity correction to the energy term ---*/
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        Residual[nVar-1] += Pressure_b*ProjGridVel*Area;
      }

      /*--- Add the Reynolds stress tensor contribution ---*/

      if (tkeNeeded) {
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] += (2.0/3.0)*Density_b*turb_ke*NormalArea[iDim];
      }

      /*--- Add value to the residual ---*/

      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Form Jacobians for implicit computations ---*/

      if (implicit) {

        /*--- Initialize Jacobian ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }

        /*--- Compute DubDu ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0.0;
          DubDu[iVar][iVar]= 1.0;
        }

        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim<nDim; jDim++)
            DubDu[iDim+1][jDim+1] -= UnitNormal[iDim]*UnitNormal[jDim];
        DubDu[nVar-1][0] += 0.5*ProjVelocity_i*ProjVelocity_i;
        for (iDim = 0; iDim < nDim; iDim++) {
          DubDu[nVar-1][iDim+1] -= ProjVelocity_i*UnitNormal[iDim];
        }

        /*--- Compute flux Jacobian in state b ---*/

        numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, NormalArea, 1, Jacobian_b);

        // Check for grid movement, should be already considered since Jacobian b is computed from u_b
        // if (grid_movement) {
        // Jacobian_b[nVar-1][0] += 0.5*ProjGridVel*ProjGridVel;
        // for (iDim = 0; iDim < nDim; iDim++)
        // Jacobian_b[nVar-1][iDim+1] -= ProjGridVel * UnitNormal[iDim];
        // }

        /*--- Compute numerical flux Jacobian at node i ---*/

        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            for (kVar = 0; kVar < nVar; kVar++)
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];

        /*--- Add the Jacobian to the sparse matrix ---*/

        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      }
    }
  }

  delete [] Normal;
  delete [] NormalArea;
  delete [] Velocity_b;
  delete [] Velocity_i;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_b[iVar];
    delete [] DubDu[iVar];
  }
  delete [] Jacobian_b;
  delete [] DubDu;

}

void CImpactSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  cout << "CImpactSolver BC_Far_Field*************************************************************************************" << endl;
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double *GridVel;
  su2double Area, UnitNormal[3] = {0.0,0.0,0.0};
  su2double Density, Pressure, Energy,  Velocity[3] = {0.0,0.0,0.0};
  su2double Density_Bound, Pressure_Bound, Vel_Bound[3] = {0.0,0.0,0.0};
  su2double Density_Infty, Pressure_Infty, Vel_Infty[3] = {0.0,0.0,0.0};
  su2double SoundSpeed, Entropy, Velocity2, Vn;
  su2double SoundSpeed_Bound, Entropy_Bound, Vel2_Bound, Vn_Bound;
  su2double SoundSpeed_Infty, Entropy_Infty, Vel2_Infty, Vn_Infty, Qn_Infty;
  su2double RiemannPlus, RiemannMinus;
  su2double *V_infty, *V_domain;

  su2double Gas_Constant     = config->GetGas_ConstantND();

  bool implicit       = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool grid_movement  = config->GetGrid_Movement();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS ) ||
                     (config->GetKind_Solver() == DISC_ADJ_RANS))
                    && (config->GetKind_Turb_Model() == SST));

  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Allocate the value at the infinity ---*/
    V_infty = GetCharacPrimVar(val_marker, iVertex);

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Construct solution state at infinity for compressible flow by
         using Riemann invariants, and then impose a weak boundary condition
         by computing the flux using this new state for U. See CFD texts by
         Hirsch or Blazek for more detail. Adapted from an original
         implementation in the Stanford University multi-block (SUmb) solver
         in the routine bcFarfield.f90 written by Edwin van der Weide,
         last modified 06-12-2005. First, compute the unit normal at the
         boundary nodes. ---*/

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Store primitive variables (density, velocities, velocity squared,
         energy, pressure, and sound speed) at the boundary node, and set some
         other quantities for clarity. Project the current flow velocity vector
         at this boundary node into the local normal direction, i.e. compute
         v_bound.n.  ---*/

      Density_Bound = V_domain[nDim+2];
      Vel2_Bound = 0.0; Vn_Bound = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Bound[iDim] = V_domain[iDim+1];
        Vel2_Bound     += Vel_Bound[iDim]*Vel_Bound[iDim];
        Vn_Bound       += Vel_Bound[iDim]*UnitNormal[iDim];
      }
      Pressure_Bound   = node[iPoint]->GetPressure();
      SoundSpeed_Bound = sqrt(Gamma*Pressure_Bound/Density_Bound);
      Entropy_Bound    = pow(Density_Bound, Gamma)/Pressure_Bound;

      /*--- Store the primitive variable state for the freestream. Project
         the freestream velocity vector into the local normal direction,
         i.e. compute v_infty.n. ---*/

      Density_Infty = GetDensity_Inf();
      Vel2_Infty = 0.0; Vn_Infty = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Infty[iDim] = GetVelocity_Inf(iDim);
        Vel2_Infty     += Vel_Infty[iDim]*Vel_Infty[iDim];
        Vn_Infty       += Vel_Infty[iDim]*UnitNormal[iDim];
      }
      Pressure_Infty   = GetPressure_Inf();
      SoundSpeed_Infty = sqrt(Gamma*Pressure_Infty/Density_Infty);
      Entropy_Infty    = pow(Density_Infty, Gamma)/Pressure_Infty;

      /*--- Adjust the normal freestream velocity for grid movement ---*/

      Qn_Infty = Vn_Infty;
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          Qn_Infty -= GridVel[iDim]*UnitNormal[iDim];
      }

      /*--- Compute acoustic Riemann invariants: R = u.n +/- 2c/(gamma-1).
         These correspond with the eigenvalues (u+c) and (u-c), respectively,
         which represent the acoustic waves. Positive characteristics are
         incoming, and a physical boundary condition is imposed (freestream
         state). This occurs when either (u.n+c) > 0 or (u.n-c) > 0. Negative
         characteristics are leaving the domain, and numerical boundary
         conditions are required by extrapolating from the interior state
         using the Riemann invariants. This occurs when (u.n+c) < 0 or
         (u.n-c) < 0. Note that grid movement is taken into account when
         checking the sign of the eigenvalue. ---*/

      /*--- Check whether (u.n+c) is greater or less than zero ---*/

      if (Qn_Infty > -SoundSpeed_Infty) {
        /*--- Subsonic inflow or outflow ---*/
        RiemannPlus = Vn_Bound + 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Supersonic inflow ---*/
        RiemannPlus = Vn_Infty + 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }

      /*--- Check whether (u.n-c) is greater or less than zero ---*/

      if (Qn_Infty > SoundSpeed_Infty) {
        /*--- Supersonic outflow ---*/
        RiemannMinus = Vn_Bound - 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Subsonic outflow ---*/
        RiemannMinus = Vn_Infty - 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }

      /*--- Compute a new value for the local normal velocity and speed of
         sound from the Riemann invariants. ---*/

      Vn = 0.5 * (RiemannPlus + RiemannMinus);
      SoundSpeed = 0.25 * (RiemannPlus - RiemannMinus)*Gamma_Minus_One;

      /*--- Construct the primitive variable state at the boundary for
         computing the flux for the weak boundary condition. The values
         that we choose to construct the solution (boundary or freestream)
         depend on whether we are at an inflow or outflow. At an outflow, we
         choose boundary information (at most one characteristic is incoming),
         while at an inflow, we choose infinity values (at most one
         characteristic is outgoing). ---*/

      if (Qn_Infty > 0.0)   {
        /*--- Outflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Bound[iDim] + (Vn-Vn_Bound)*UnitNormal[iDim];
        Entropy = Entropy_Bound;
      } else  {
        /*--- Inflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Infty[iDim] + (Vn-Vn_Infty)*UnitNormal[iDim];
        Entropy = Entropy_Infty;
      }

      /*--- Recompute the primitive variables. ---*/

      Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
      Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      if (tkeNeeded) Energy += GetTke_Inf();

      /*--- Store new primitive state for computing the flux. ---*/

      V_infty[0] = Pressure/(Gas_Constant*Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_infty[iDim+1] = Velocity[iDim];
      V_infty[nDim+1] = Pressure;
      V_infty[nDim+2] = Density;
      V_infty[nDim+3] = Energy + Pressure/Density;



      /*--- Set various quantities in the numerics class ---*/

      conv_numerics->SetPrimitive(V_domain, V_infty);

      if (grid_movement) {
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      }

      /*--- Compute the convective residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Convective Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Roe Turkel preconditioning, set the value of beta ---*/

      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());


    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}


void CImpactSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  cout << "CImpactSolver BC_Inlet*************************************************************************************" << endl;
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
  Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
  alpha, aa, bb, cc, dd, Area, UnitNormal[3];
  su2double *V_inlet, *V_domain;

  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  su2double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
  su2double Gas_Constant       = config->GetGas_ConstantND();
  unsigned short Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the inlet ---*/

    V_inlet = GetCharacPrimVar(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Retrieve solution at this boundary node ---*/

      V_domain = node[iPoint]->GetPrimitive();

      /*--- Build the fictitious intlet state based on characteristics ---*/


      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/

      switch (Kind_Inlet) {

        /*--- Total properties have been specified at the inlet. ---*/

        case TOTAL_CONDITIONS:

          /*--- Retrieve the specified total conditions for this inlet. ---*/

          P_Total  = Inlet_Ptotal[val_marker][iVertex];
          T_Total  = Inlet_Ttotal[val_marker][iVertex];
          Flow_Dir = Inlet_FlowDir[val_marker][iVertex];

          /*--- Non-dim. the inputs if necessary. ---*/

          P_Total /= config->GetPressure_Ref();
          T_Total /= config->GetTemperature_Ref();

          /*--- Store primitives and set some variables for clarity. ---*/

          Density = V_domain[nDim+2];
          Velocity2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = V_domain[iDim+1];
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
          Pressure    = V_domain[nDim+1];
          H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
          SoundSpeed2 = Gamma*Pressure/Density;

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Total speed of sound ---*/

          SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

          /*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/

          alpha = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            alpha += UnitNormal[iDim]*Flow_Dir[iDim];

          /*--- Coefficients in the quadratic equation for the velocity ---*/

          aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
          bb = -1.0*Gamma_Minus_One*alpha*Riemann;
          cc =  0.5*Gamma_Minus_One*Riemann*Riemann
              -2.0*SoundSpeed_Total2/Gamma_Minus_One;

          /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/

          dd = bb*bb - 4.0*aa*cc;
          dd = sqrt(max(0.0, dd));
          Vel_Mag   = (-bb + dd)/(2.0*aa);
          Vel_Mag   = max(0.0, Vel_Mag);
          Velocity2 = Vel_Mag*Vel_Mag;

          /*--- Compute speed of sound from total speed of sound eqn. ---*/

          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/

          Mach2 = Velocity2/SoundSpeed2;
          Mach2 = min(1.0, Mach2);
          Velocity2   = Mach2*SoundSpeed2;
          Vel_Mag     = sqrt(Velocity2);
          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Compute new velocity vector at the inlet ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

          /*--- Static temperature from the speed of sound relation ---*/

          Temperature = SoundSpeed2/(Gamma*Gas_Constant);

          /*--- Static pressure using isentropic relation at a point ---*/

          Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);

          /*--- Density at the inlet from the gas law ---*/

          Density = Pressure/(Gas_Constant*Temperature);

          /*--- Using pressure, density, & velocity, compute the energy ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
          if (tkeNeeded) Energy += GetTke_Inf();

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Temperature;
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Velocity[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;

          /*--- Mass flow has been specified at the inlet. ---*/

        case MASS_FLOW:

          /*--- Retrieve the specified mass flow for the inlet. ---*/

          Density  = Inlet_Ttotal[val_marker][iVertex];
          Vel_Mag  = Inlet_Ptotal[val_marker][iVertex];
          Flow_Dir = Inlet_FlowDir[val_marker][iVertex];

          /*--- Non-dim. the inputs if necessary. ---*/

          Density /= config->GetDensity_Ref();
          Vel_Mag /= config->GetVelocity_Ref();

          /*--- Get primitives from current inlet state. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
          Pressure    = node[iPoint]->GetPressure();
          SoundSpeed2 = Gamma*Pressure/V_domain[nDim+2];

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Speed of sound squared for fictitious inlet state ---*/

          SoundSpeed2 = Riemann;
          for (iDim = 0; iDim < nDim; iDim++)
            SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitNormal[iDim];

          SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
          SoundSpeed2 = SoundSpeed2*SoundSpeed2;

          /*--- Pressure for the fictitious inlet state ---*/

          Pressure = SoundSpeed2*Density/Gamma;

          /*--- Energy for the fictitious inlet state ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Vel_Mag*Vel_Mag;
          if (tkeNeeded) Energy += GetTke_Inf();

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Pressure / ( Gas_Constant * Density);
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;
      }

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Roe Turkel preconditioning, set the value of beta ---*/

      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());

    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;

}

void CImpactSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  cout << "CImpactSolver BC_Outlet*************************************************************************************" << endl;
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint;
  su2double Pressure, P_Exit, Velocity[3],
  Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
  Area, UnitNormal[3];
  su2double *V_outlet, *V_domain;

  bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  su2double Gas_Constant     = config->GetGas_ConstantND();
  bool grid_movement      = config->GetGrid_Movement();
  string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the outlet ---*/
    V_outlet = GetCharacPrimVar(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Current solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Build the fictitious intlet state based on characteristics ---*/

      /*--- Retrieve the specified back pressure for this outlet. ---*/
      if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDARD_GRAVITY;
      else P_Exit = config->GetOutlet_Pressure(Marker_Tag);

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit = P_Exit/config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
         of boundary update depends on this. ---*/
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Mach_Exit  = sqrt(Velocity2)/SoundSpeed;

      if (Mach_Exit >= 1.0) {

        /*--- Supersonic exit flow: there are no incoming characteristics,
           so no boundary condition is necessary. Set outlet state to current
           state so that upwinding handles the direction of propagation. ---*/
        for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = V_domain[iVar];

      } else {

        /*--- Subsonic exit flow: there is one incoming characteristic,
           therefore one variable can be specified (back pressure) and is used
           to update the conservative variables. Compute the entropy and the
           acoustic Riemann variable. These invariants, as well as the
           tangential velocity components, are extrapolated. Adapted from an
           original implementation in the Stanford University multi-block
           (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
           der Weide, last modified 09-10-2007. ---*/

        Entropy = Pressure*pow(1.0/Density, Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Compute the new fictious state at the outlet ---*/
        Density    = pow(P_Exit/Entropy,1.0/Gamma);
        Pressure   = P_Exit;
        SoundSpeed = sqrt(Gamma*P_Exit/Density);
        Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();

        /*--- Conservative variables, using the derived quantities ---*/
        V_outlet[0] = Pressure / ( Gas_Constant * Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_outlet[iDim+1] = Velocity[iDim];
        V_outlet[nDim+1] = Pressure;
        V_outlet[nDim+2] = Density;
        V_outlet[nDim+3] = Energy + Pressure/Density;

      }

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }

      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());

    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}




void CImpactSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {
  cout << "CImpactSolver BC_Sym_Plane*************************************************************************************" << endl;
  /*--- Call the Euler residual ---*/

  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);

}

void CImpactSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config) {
  cout << "CImpactSolver BC_Fluid_Interface*************************************************************************************" << endl;
  unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
  unsigned short iDim, iVar, iMarker, nDonorVertex;

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();

  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  su2double *tmp_residual = new su2double[nVar];

  su2double weight;
  su2double P_static, rho_static;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          nDonorVertex = GetnSlidingStates(iMarker, iVertex);

          /*--- Initialize Residual, this will serve to accumulate the average ---*/

          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = 0.0;

          /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

          for (jVertex = 0; jVertex < nDonorVertex; jVertex++){

            Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = GetSlidingState(iMarker, iVertex, iVar, jVertex);
            }

            /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

            weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

            /*--- Set primitive variables ---*/

            conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );

            if( !( config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS ) ) {
              Secondary_i = node[iPoint]->GetSecondary();

              P_static   = PrimVar_j[nDim+1];
              rho_static = PrimVar_j[nDim+2];
              FluidModel->SetTDState_Prho(P_static, rho_static);

              Secondary_j[0] = FluidModel->GetdPdrho_e();
              Secondary_j[1] = FluidModel->GetdPde_rho();

              conv_numerics->SetSecondary(Secondary_i, Secondary_j);
            }

            /*--- Set the normal vector ---*/

            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++)
              Normal[iDim] = -Normal[iDim];

            conv_numerics->SetNormal(Normal);

            if (grid_movement)
              conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

            /*--- Compute the convective residual using an upwind scheme ---*/

            conv_numerics->ComputeResidual(tmp_residual, Jacobian_i, Jacobian_j, config);

            /*--- Accumulate the residuals to compute the average ---*/

            for (iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] += weight*tmp_residual[iVar];

          }

          /*--- Add Residuals and Jacobians ---*/

          LinSysRes.AddBlock(iPoint, Residual);
          if (implicit)
            Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] tmp_residual;
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
}

void CImpactSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container,
                                CConfig *config, unsigned short val_marker) {cout << "CImpactSolver BC_Dirichlet*************************************************************************************" << endl; }

void CImpactSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {cout << "CImpactSolver BC_Custom*************************************************************************************" << endl; }

void CImpactSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {
  cout << "CImpactSolver LoadRestart*************************************************************************************" << endl;
  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  bool grid_movement  = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool static_fsi = ((config->GetUnsteady_Simulation() == STEADY) &&
                     (config->GetFSI_Simulation()));
  bool steady_restart = config->GetSteadyRestart();
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  string UnstExt, text_line;
  ifstream restart_file;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string restart_filename = config->GetSolution_ImpactFileName();

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time || time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      node[iPoint_Local]->SetSolution(Solution);
      iPoint_Global_Local++;

      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/

      if (grid_movement && val_update_geo) {

        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/

        if (turb_model == SA || turb_model == SA_NEG) {
          index++;
        } else if (turb_model == SST) {
          index+=2;
        }

        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/

        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar;
          for (iDim = 0; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[MESH_0]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
      }

      if (static_fsi && val_update_geo) {
       /*--- Rewind the index to retrieve the Coords. ---*/
        index = counter*Restart_Vars[1];
        for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim];}

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
        }
      }

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We alo call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  solver[MESH_0][IMPACT_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][IMPACT_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][IMPACT_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][IMPACT_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][IMPACT_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][IMPACT_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  }

  /*--- Update the geometry for flows on dynamic meshes ---*/

  if (grid_movement && val_update_geo) {

    /*--- Communicate the new coordinates and grid velocities at the halos ---*/

    geometry[MESH_0]->Set_MPI_Coord(config);
    geometry[MESH_0]->Set_MPI_GridVel(config);

    /*--- Recompute the edges and dual mesh control volumes in the
     domain and on the boundaries. ---*/

    geometry[MESH_0]->SetCoord_CG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);

    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine], config);
      }
    }

  /*--- Update the geometry for flows on static FSI problems with moving meshes ---*/

  if (static_fsi && val_update_geo) {

    /*--- Communicate the new coordinates and grid velocities at the halos ---*/

    geometry[MESH_0]->Set_MPI_Coord(config);

    /*--- Recompute the edges and  dual mesh control volumes in the
     domain and on the boundaries. ---*/

    geometry[MESH_0]->SetCoord_CG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);

    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
    }
  }

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}


void CImpactSolver::SetFreeStream_Solution(CConfig *config) {
  cout << "CImpactSolver SetFreeStream_Solution*************************************************************************************" << endl;
  /*--- not used now ---*/
  unsigned long iPoint;
  unsigned short iDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetSolution(0, Density_Inf);
    for (iDim = 0; iDim < nDim; iDim++) {
      node[iPoint]->SetSolution(iDim+1, Density_Inf*Velocity_Inf[iDim]);
    }
    node[iPoint]->SetSolution(nVar-1, Density_Inf*Energy_Inf);
  }
}
