/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

License
    This file is part of ITHACA-FV

    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the reducedUnsteadyNS class


#include "ReducedUnsteadyNS.H"
#include "unsteadyNS.H"
//#include <list>
//#include <iostream>

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor initialization
// original
/*reducedUnsteadyNS::reducedUnsteadyNS()
{
}

reducedUnsteadyNS::reducedUnsteadyNS(unsteadyNS& FOMproblem)
    :
    problem(&FOMproblem)
{
    N_BC = problem->inletIndex.rows();
    //N_BC_P = problem->outletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();

    // Create locally the velocity modes
    for (int k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append((problem->liftfield[k]).clone());
    }

    for (int k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append((problem->Umodes[k]).clone());
    }

    for (int k = 0; k < problem->NSUPmodes; k++)
    {
        Umodes.append((problem->supmodes[k]).clone());
    }

    // Create locally the pressure modes
    for (int k = 0; k < problem->NPmodes; k++)
    {
        Pmodes.append((problem->Pmodes[k]).clone());
    }

    newton_object_sup = newton_unsteadyNS_sup(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                        FOMproblem);
    newton_object_PPE = newton_unsteadyNS_PPE(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                        FOMproblem);
}

// * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //

// Operator to evaluate the residual for the Supremizer approach
int newton_unsteadyNS_sup::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);

    // Choose the order of the numerical difference scheme for approximating the time derivative
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(
                     Nphi_u)) / dt;
    }

    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = tauU(l,
                                   0) * (BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                                         a_tmp);
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_tmp;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i);

        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC; l++)
            {
                fvec(i) += penaltyU(i, l);
            }
        }
    }

    for (int j = 0; j < Nphi_p; j++)
    {
        int k = j + Nphi_u;
        fvec(k) = M3(j);
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);
        }
    }

    return 0;
}*/
reducedUnsteadyNS::reducedUnsteadyNS()
{
}

reducedUnsteadyNS::reducedUnsteadyNS(unsteadyNS& FOMproblem)
    :
    problem(&FOMproblem)
{
    N_BC = problem->inletIndex.rows();
    N_BC_P = problem->outletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols(); 
    cout << "assegnazione Nphi_u:" << Nphi_u << "\n";
    cout << "assegnazione Nphi_p:" << Nphi_p << "\n";
    cout << "nu: " << nu << "\n";
    // Create locally the velocity modes
    for (int k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append((problem->liftfield[k]).clone());
    }

    for (int k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append((problem->Umodes[k]).clone());
    }

    for (int k = 0; k < problem->NSUPmodes; k++)
    {
        Umodes.append((problem->supmodes[k]).clone());
    }

    // Create locally the pressure modes
    //INSERITO DA ME
    Add_liftfieldP = problem->para->ITHACAdict->lookupOrDefault<word>("Add_liftfieldP", "no");
    kernel_stabilization = problem->para->ITHACAdict->lookupOrDefault<word>("kernel_stabilization", "no");
    //M_Assert(Add_liftfieldP == "yes" || Add_liftfieldP == "no",
      //       "The BC method can be set to yes or no");
    std::cout << "Add_liftfieldP" << Add_liftfieldP <<"\n";
    std::cout << "kernel_stabilization" << kernel_stabilization <<"\n";
    if (Add_liftfieldP == "yes"){
    	for (int k = 0; k < problem->liftfieldP.size(); k++)
    	{
		Pmodes.append((problem->liftfieldP[k]).clone());
		//std::cout << "liftfieldP aggiunta in Pmodes " <<"\n";
    	}
	for (int k = 0; k < problem->NPmodes; k++)
    	{
        	Pmodes.append((problem->Pmodes[k]).clone());
    	}
    } 
    else{
	    for (int k = 0; k < problem->NPmodes; k++)
    	    {
        	Pmodes.append((problem->Pmodes[k]).clone());
    	    }
    }
    cout << "Nphi_p:" << Nphi_p << "\n";
    cout << "Nphi_u:" << Nphi_u << "\n";
    cout << "nu: " << nu << "\n";
    newton_object_sup = newton_unsteadyNS_sup(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                        FOMproblem);
    newton_object_PPE = newton_unsteadyNS_PPE(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                        FOMproblem);
    cout << "dopo la def di neton object sup e PPE" << "\n";
}

// * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //

// Operator to evaluate the residual for the Supremizer approach
int newton_unsteadyNS_sup::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    
    //cout << "newton_unsteadyNS_sup Nphi_p:" << Nphi_p << "\n";
    //cout << "newton_unsteadyNS_sup Nphi_u:" << Nphi_u << "\n";
    //cout << "nu: " << nu << "\n";
    // Choose the order of the numerical difference scheme for approximating the time derivative
    //a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(
                     Nphi_u)) / dt;
    }

    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;

    // qui mettiamo un flag FDA
    if (problem->kernel_stabilization == "yes"){
	 cout << "Nphi_u per stabilizzazione: " << Nphi_u << "\n";
         Eigen::VectorXd M1_app = problem->B_matrix * a_tmp * nu;
         //M1 = problem->B_matrix * a_tmp * nu;
         for (int i = 0; i < Nphi_u; i++){
                   M1(i) = (1 + nu_a * i / Nphi_u) * M1_app(i);
                   //cout << "nu_a: " << nu_a << "\n";
         }
         cout << "kernel fatto" << nu_a << "\n";
    }

    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = tauU(l,
                                   0) * (BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                                         a_tmp);
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_tmp;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i);

        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC; l++)
            {
                fvec(i) += penaltyU(i, l);
            }
        }
    }

    for (int j = 0; j < Nphi_p; j++)
    {
        int k = j + Nphi_u;
        fvec(k) = M3(j);
    }

    //cout << "fvec:" << fvec << "\n";

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);
        }
    }
    if (problem->Add_liftfieldP == "yes"){//ho aggiunto problem-> se non funziona e' questo
        for (int j = 0; j < N_BC_P; j++)
        {
            fvec(j + Nphi_u) = x(j + Nphi_u) - BC_P(j);
       	}
    }
//	cout << "fvec:" << fvec << "\n";


    

    //ITHACAstream::exportMatrix(fvec, "fvec", "python",
      //                         "./ITHACAoutput/fvec");

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyNS_sup::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNS_sup> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

// Operator to evaluate the residual for the Pressure Poisson Equation (PPE) approach
int newton_unsteadyNS_PPE::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    //cout << "QUI FATTO3" << "\n";
    // Choose the order of the numerical difference scheme for approximating the time derivative
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(
                     Nphi_u)) / dt;
    }
    
    //kernel_stabilization = problem->para->ITHACAdict->lookupOrDefault<word>("kernel_stabilization", "no");

    // Convective terms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;

    // qui mettiamo un flag FDA
    if (problem->kernel_stabilization == "yes"){ 
         Eigen::VectorXd M1_app = problem->B_matrix * a_tmp * nu;
	 //M1 = problem->B_matrix * a_tmp * nu;
	 for (int i = 0; i < Nphi_u; i++){
 	 	   M1(i) = (1 + nu_a * (pow(i,1))/(Nphi_u)) * M1_app(i);//(1 + (i/Nphi_u)*nu_a) * M1_app(i);
		   //cout << "nu_a: " << nu_a << "\n";
         }
    	 cout << "kernel fatto" << nu_a << "\n";
    }

    /*if (kernel_stabilization == "no"){
	 M1 = problem->B_matrix * a_tmp * nu;
    }*/

    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->D_matrix * b_tmp;
    // BC PPE
    Eigen::VectorXd M7 = problem->BC3_matrix * a_tmp * nu;// * (1 + nu_a);

    //cout << "PRIMA DI KERNEL_STABILIZATION:" << problem->kernel_stabilization << "\n";
    // qui mettiamo un flag FDA
    /*if (problem->kernel_stabilization == "yes"){  
         Eigen::VectorXd M7_app = problem->BC3_matrix * a_tmp * nu;
         //M1 = problem->B_matrix * a_tmp * nu;
         for (int i = 0; i < Nphi_p; i++){
                   M7(i) = (1 + nu_a * ((pow(i,1))/(Nphi_p))) * M7_app(i);//(1 + (i/Nphi_u)*nu_a) * M1_app(i);
                   //cout << "nu_a: " << nu_a << "\n";
         }
	 cout << "kernel fatto" << problem->kernel_stabilization << "\n";
    }*/
    
    // BC PPE time-dependents BCs
    Eigen::VectorXd M8 = problem->BC4_matrix * a_dot;
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = tauU(l,
                                   0) * (BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                                         a_tmp);
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_tmp;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i);

        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC; l++)
            {
                fvec(i) += penaltyU(i, l);
            }
        }
    }

    for (int j = 0; j < Nphi_p; j++)
    {
        int k = j + Nphi_u;
        gg = a_tmp.transpose() * Eigen::SliceFromTensor(problem->gTensor, 0,
                j) * a_tmp;
        fvec(k) = M3(j, 0) + gg(0, 0) - M7(j, 0);

        if (problem->timedepbcMethod == "yes")
        {
            fvec(k) += M8(j, 0);
        }
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);  //di questo non sono sicuro
        }
	/*for (int j = 0; j< N_BC_P; j++)
	{
	    fvec(Nphi_u + j) = x(Nphi_u + j) - BC_P(j);
	}*/
    }
    if (Add_liftfieldP == "yes"){//ho aggiunto problem se non funziona e' colpa sua
        for (int j = 0; j < N_BC_P; j++)
        {
            fvec(j + Nphi_u) = x(j + Nphi_u) - BC_P(j);
        }
    }


    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyNS_PPE::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNS_PPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * Solve Functions supremizer * * * * * * * * * * * //

//original
/*void reducedUnsteadyNS::solveOnline_sup(Eigen::MatrixXd vel,
                                        int startSnap)
{
    M_Assert(exportEvery >= dt,
             "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt,
             "The time step dt must be smaller than storeEvery.");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt) == true,
             "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt) == true,
             "The variable exportEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery) == true,
             "The variable exportEvery must be an integer multiple of the variable storeEvery.");
    int numberOfStores = round(storeEvery / dt);

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[startSnap],
                     Umodes);
    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[startSnap],
                     Pmodes);
    int nextStore = 0;
    int counter2 = 0;

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
    }

    // Set some properties of the newton object
    newton_object_sup.nu = nu;
    newton_object_sup.y_old = y;
    newton_object_sup.yOldOld = newton_object_sup.y_old;
    newton_object_sup.dt = dt;
    newton_object_sup.BC.resize(N_BC);
    newton_object_sup.tauU = tauU;

    for (int j = 0; j < N_BC; j++)
    {
        newton_object_sup.BC(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_solution[counter] = tmp_sol;
    counter ++;
    counter2++;
    nextStore += numberOfStores;
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_sup> hnls(newton_object_sup);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    while (time < finalTime)
    {
        time = time + dt;

        // Set time-dependent BCs
        if (problem->timedepbcMethod == "yes" )
        {
            for (int j = 0; j < N_BC; j++)
            {
                newton_object_sup.BC(j) = vel_now(j, counter);
            }
        }

        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);

        if (problem->bcMethod == "lift")
        {
            for (int j = 0; j < N_BC; j++)
            {
                if (problem->timedepbcMethod == "no" )
                {
                    y(j) = vel_now(j, 0);
                }
                else if (problem->timedepbcMethod == "yes" )
                {
                    y(j) = vel_now(j, counter);
                }
            }
        }

        newton_object_sup.operator()(y, res);
        newton_object_sup.yOldOld = newton_object_sup.y_old;
        newton_object_sup.y_old = y;
        std::cout << "################## Online solve N° " << counter <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;

        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;

        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size())
            {
                online_solution.append(tmp_sol);
            }
            else
            {
                online_solution[counter2] = tmp_sol;
            }

            nextStore += numberOfStores;
            counter2 ++;
        }

        counter ++;
    }

    // Export the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
}*/


void reducedUnsteadyNS::solveOnline_sup(Eigen::MatrixXd vel, Eigen::MatrixXd pressure,
                                        int startSnap)
{
    M_Assert(exportEvery >= dt,
             "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt,
             "The time step dt must be smaller than storeEvery.");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt) == true,
             "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt) == true,
             "The variable exportEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery) == true,
             "The variable exportEvery must be an integer multiple of the variable storeEvery.");

    int numberOfStores = round(storeEvery / dt);

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
	if (Add_liftfieldP == "yes"){
		P_now = setOnlinePressure(pressure);
		cout << "P_now set" << "\n";
	}
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    cout << "Nphi_u:" << Nphi_u << "\n";
    cout << "reducedUnsteadyNS Nphi_p:" << Nphi_p << "\n";
    cout << "y size:" << y.size() << "\n";
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[startSnap],
                     Umodes);
    cout << "y.head(Nphi_u)" << y << "\n";
    cout << "Pmodes size" << Pmodes.size() << "\n";
    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[startSnap],
                     Pmodes);
    cout << "y.tail(Nphi_p)" << y << "\n";
    int nextStore = 0;
    int counter2 = 0;

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
	if (Add_liftfieldP == "yes"){
        	for (int j = 0; j < N_BC_P; j++)
        	{
            		y(j + Nphi_u) = P_now(j, 0);
        	}
	}
    }
    //cout << "FINO A QUI SI" << "\n";
    // Set some properties of the newton object
    newton_object_sup.nu = nu;
    newton_object_sup.nu_a = nu_a;
    newton_object_sup.y_old = y;
    //cout << "Nphi_p prima"<< Nphi_p << "\n";
    newton_object_sup.yOldOld = newton_object_sup.y_old;
    //cout << "Nphi_p dopo" << Nphi_p << "\n";
    newton_object_sup.dt = dt;
    newton_object_sup.BC.resize(N_BC);
    if (Add_liftfieldP == "yes"){newton_object_sup.BC_P.resize(N_BC_P);};
    //newton_object_sup.BC_P.resize(N_BC_P);
    newton_object_sup.tauU = tauU;

    for (int j = 0; j < N_BC; j++)
    {
        newton_object_sup.BC(j) = vel_now(j, 0);
    }
    cout << "Add_liftfieldP" << Add_liftfieldP << "\n";
    if (Add_liftfieldP == "yes"){
	//cout << "PROVA" << "\n";
    	for (int j = 0; j < N_BC_P; j++)
    	{
        	newton_object_sup.BC_P(j) = P_now(j, 0);
		//cout << "PROVA" << "\n";
    	}
    }

    
    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    //cout << "PROVA" << "\n";
    online_solution.resize(onlineSize);
    //cout << "PROVA" << "\n";
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    //cout << "PROVA" << "\n";
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_solution[counter] = tmp_sol;
    counter ++;
    counter2++;
    nextStore += numberOfStores;
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_sup> hnls(newton_object_sup);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    
    while (time < finalTime && time != finalTime)
    {
        cout << "time:" << time << "\n";
        cout << "finalTime:" << finalTime << "\n";
        //cout << "SCEMO" << "\n";
        time = time + dt;

        // Set time-dependent BCs
        if (problem->timedepbcMethod == "yes" )
        {
            for (int j = 0; j < N_BC; j++)
            {
                newton_object_sup.BC(j) = vel_now(j, counter);
            }
	    if (Add_liftfieldP == "yes"){
            	for (int j = 0; j < N_BC_P; j++)
            	{
                	newton_object_sup.BC_P(j) = P_now(j, counter);
            	}
	    }

        }

        Eigen::VectorXd res(y);
        res.setZero();
	//cout << "CIAONE" << "\n";
        hnls.solve(y);
	//cout << "CIAONE" << "\n";
	//cout << "y:" << y << "\n";

        if (problem->bcMethod == "lift")
        {
            for (int j = 0; j < N_BC; j++)
            {
                if (problem->timedepbcMethod == "no" )
                {
                    y(j) = vel_now(j, 0);
                }
                else if (problem->timedepbcMethod == "yes" )
                {
                    y(j) = vel_now(j, counter);
                }
            }
	    if (Add_liftfieldP == "yes"){
            	for (int j = 0; j < N_BC_P; j++)
            	{
                	if (problem->timedepbcMethod == "no" )
                	{
                    		y(j + Nphi_u) = P_now(j, 0);
                	}
                	else if (problem->timedepbcMethod == "yes" )
                	{
                    		//y(j) = vel_now(j, counter);
		   		y(j + Nphi_u) = P_now(j, counter);
                	}
            	}
	    }

        }

        newton_object_sup.operator()(y, res);
        newton_object_sup.yOldOld = newton_object_sup.y_old;
        newton_object_sup.y_old = y;
        std::cout << "################## Online solve N° " << counter <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;

        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }

        tmp_sol(0) = time;
	//cout << "CIAONE" << "\n";
        tmp_sol.col(0).tail(y.rows()) = y;
	//cout << "CIAONE" << "\n";
	//cout << "counter" << counter << "\n";
	//cout << "nextStore" << nextStore << "\n";
        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size())
            {
                online_solution.append(tmp_sol);
		//cout << "CIAONE" << "\n";
            }
            else
            {
                //cout << "counter2" << counter2 << "\n";
		//cout << "tmp_sol size" << tmp_sol.size() << "\n";
		//cout << "online_solution size" << online_solution.size() << "\n";
		online_solution[counter2] = tmp_sol;
		//cout << "CIAONE" << "\n";
            }
	    //cout << "nextstore" << nextStore << "\n";
	    //cout << "numberOfStores" << numberOfStores << "\n";
            nextStore += numberOfStores;
            counter2 ++;
        }
	//cout << "CIAONE" << "\n";
        counter ++;
	cout << "counter" << counter << "\n";
    }
    //cout << "CIAONE" << "\n";
    // Export the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
}

// * * * * * * * * * * * * * * * Solve Functions PPE * * * * * * * * * * * * * //

void reducedUnsteadyNS::solveOnline_PPE(Eigen::MatrixXd vel, Eigen::MatrixXd pressure,
                                        int startSnap)
{
    M_Assert(exportEvery >= dt,
             "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt,
             "The time step dt must be smaller than storeEvery.");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt) == true,
             "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt) == true,
             "The variable exportEvery must be an integer multiple of the time step dt.");

    int numberOfStores = round(storeEvery / dt);

    if (problem->bcMethod == "lift")
    {
	//cout << "vel size:" << vel.size() << "\n";
	vel_now = setOnlineVelocity(vel);
	//cout << "vel_now size:" << vel_now.size() << "\n";
	//P_now = setOnlinePressure(pressure);
        if (Add_liftfieldP == "yes"){
                P_now = setOnlinePressure(pressure);
                //cout << "P_now set" << "\n";
        }

    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }
    
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    // Set Initial Conditions
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[startSnap],
                     Umodes);
    cout << "Nphi_u:" << Nphi_u << "\n";
    cout << "Nphi_p:" << Nphi_p << "\n";
    //cout << "nu_a: " << nu_a << "\n";
    cout << "nu: " << nu << "\n";
    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[startSnap],
                     Pmodes); 
    //cout << "Pmodes:" << Pmodes << "\n"; 
    int nextStore = 0;
    int counter2 = 0;
    //cout << "prima di Change condition for the lift fatto" << "\n";
    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
	/*for (int j = 0; j < N_BC_P; j++)
	{   
	    //y(Nphi_u + j) = P_now(j, 0);
	}*/

        if (Add_liftfieldP == "yes"){
                for (int j = 0; j < N_BC_P; j++)
                {
                        y(j + Nphi_u) = P_now(j, 0);
                }
        }

    }
    // Set some properties of the newton object
    newton_object_PPE.nu = nu;// * (1 + nu_a);
    //if (problem->kernel_stabilization == "yes"){newton_object_PPE.nu_a = nu_a;}else{newton_object_PPE.nu_a = 0;}
    newton_object_PPE.nu_a = nu_a;
    newton_object_PPE.y_old = y;
    newton_object_PPE.yOldOld = newton_object_PPE.y_old;
    newton_object_PPE.dt = dt;
    newton_object_PPE.BC.resize(N_BC);
    newton_object_PPE.BC_P.resize(N_BC_P);
    newton_object_PPE.tauU = tauU;

    for (int j = 0; j < N_BC; j++)
    {
        newton_object_PPE.BC(j) = vel_now(j, 0);
    }
 
    /*for (int j = 0; j< N_BC_P; j++)
    {
	//newton_object_PPE.BC_P(j) = P_now(j, 0);
    }*/

    if (Add_liftfieldP == "yes"){
        for (int j = 0; j < N_BC_P; j++)
        {
                newton_object_PPE.BC_P(j) = P_now(j, 0);
        }
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_solution[counter] = tmp_sol;
    counter ++;
    counter2++;
    nextStore += numberOfStores;
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_PPE> hnls(newton_object_PPE);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // Start the time loop
    while (time < finalTime)
    {
	//cout << "time:" << time << "\n";
	//cout << "finalTime:" << finalTime << "\n";
	//cout << "SCEMO" << "\n";
        time = time + dt;

        // Set time-dependent BCs
        if (problem->timedepbcMethod == "yes" )
        {
            for (int j = 0; j < N_BC; j++)
            {
                newton_object_PPE.BC(j) = vel_now(j, counter);
            }
	    /*for (int j = 0; j < N_BC_P; j++)
	    {
		//newton_object_PPE.BC_P(j) = P_now(j, counter);
	    }*/
            if (Add_liftfieldP == "yes"){
                for (int j = 0; j < N_BC_P; j++)
                {
                        newton_object_PPE.BC_P(j) = P_now(j, counter);
                }
            }
        }
	//cout << "prima di res(y) fatto" << "\n";
        Eigen::VectorXd res(y);
        res.setZero();
	//cout << "prima di solve(y) fatto" << "\n";
	//cout << "y size:" << y.size() <<"\n";
        hnls.solve(y);
	//cout << "dopo di solve(y) fatto" << "\n";
        if (problem->bcMethod == "lift")
        {
            for (int j = 0; j < N_BC; j++)
            {
                if (problem->timedepbcMethod == "no" )
                {
                    y(j) = vel_now(j, 0);
		    //cout << "y(j): " << y(j) << "\n";
                }
                else if (problem->timedepbcMethod == "yes" )
                {
                    y(j) = vel_now(j, counter);
                }
            }
	    /*for (int j = 0; j < N_BC_P; j++)
	    {
		if (problem->timedepbcMethod == "no")
		{
		    //y(j) = P_now(j, 0);
		}
		else if (problem->timedepbcMethod == "yes" )
		{ 
	            //y(Nphi_u + j) = P_now(j, counter);
		}
	    }*/
            if (Add_liftfieldP == "yes"){
                for (int j = 0; j < N_BC_P; j++)
                {
                        if (problem->timedepbcMethod == "no" )
                        {
                                y(j + Nphi_u) = P_now(j, 0);
                        }
                        else if (problem->timedepbcMethod == "yes" )
                        {
                                //y(j) = vel_now(j, counter);
                                y(j + Nphi_u) = P_now(j, counter);
                        }
                }
            }
        }

        newton_object_PPE.operator()(y, res);
        newton_object_PPE.yOldOld = newton_object_PPE.y_old;
        newton_object_PPE.y_old = y;
        std::cout << "################## Online solve N° " << counter <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;

        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;
	//tmp_sol(1) = 1;

        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size())
            {
                online_solution.append(tmp_sol);
            }
            else
            {
                online_solution[counter2] = tmp_sol;
            }

            nextStore += numberOfStores;
            counter2 ++;
        }

        counter ++;
    }

    // Export the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(vel, "vel", "python",
                               "./ITHACAoutput/BC_vel"); 
    ITHACAstream::exportMatrix(vel_now, "vel_now", "python",
                               "./ITHACAoutput/BC_vel_now");
    //ITHACAstream::exportMatrix(vel_now(0,0), "vel_now_0", "python",
    //                           "./ITHACAoutput/BC_vel_now_0");
    //cout << "vel now (0,0) : " << vel_now(0,0) << "\n";
    //cnpy::save(online_solution, "./ITHACAoutput/POD/red_coeff.npy");
}

Eigen::MatrixXd reducedUnsteadyNS::penalty_sup(Eigen::MatrixXd& vel_now,
        Eigen::MatrixXd& tauIter,
        int startSnap)
{
    // Initialize new value on boundaries
    Eigen::MatrixXd valBC = Eigen::MatrixXd::Zero(N_BC, timeStepPenalty);
    // Initialize old values on boundaries
    Eigen::MatrixXd valBC0 = Eigen::MatrixXd::Zero(N_BC, timeStepPenalty);
    int Iter = 0;
    Eigen::VectorXd diffvel =  (vel_now.col(timeStepPenalty - 1) - valBC.col(
                                    timeStepPenalty - 1));
    diffvel = diffvel.cwiseAbs();

    while (diffvel.maxCoeff() > tolerancePenalty && Iter < maxIterPenalty)
    {
        if ((valBC.col(timeStepPenalty - 1) - valBC0.col(timeStepPenalty - 1)).sum() !=
                0)
        {
            for (int j = 0; j < N_BC; j++)
            {
                tauIter(j, 0) = tauIter(j, 0) * diffvel(j) / tolerancePenalty;
            }
        }

        std::cout << "Solving for penalty factor(s): " << tauIter << std::endl;
        std::cout << "number of iterations: " << Iter << std::endl;
        //  Set the old boundary value to the current value
        valBC0  = valBC;
        y.resize(Nphi_u + Nphi_p, 1);
        y.setZero();
        y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[startSnap],
                         Umodes);
        y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[startSnap],
                         Pmodes);
        // Set some properties of the newton object
        newton_object_sup.nu = nu;
        newton_object_sup.y_old = y;
        newton_object_sup.dt = dt;
        newton_object_sup.BC.resize(N_BC);
        newton_object_sup.tauU = tauIter;

        // Set boundary conditions
        for (int j = 0; j < N_BC; j++)
        {
            newton_object_sup.BC(j) = vel_now(j, 0);
        }

        // Create nonlinear solver object
        Eigen::HybridNonLinearSolver<newton_unsteadyNS_sup> hnls(newton_object_sup);
        // Set output colors for fancy output
        Color::Modifier red(Color::FG_RED);
        Color::Modifier green(Color::FG_GREEN);
        Color::Modifier def(Color::FG_DEFAULT);
        // Set initially for convergence check
        Eigen::VectorXd res(y);
        res.setZero();

        // Start the time loop
        for (int i = 1; i < timeStepPenalty; i++)
        {
            // Set boundary conditions
            for (int j = 0; j < N_BC; j++)
            {
                if (problem->timedepbcMethod == "yes" )
                {
                    newton_object_sup.BC(j) = vel_now(j, i);
                }
                else
                {
                    newton_object_sup.BC(j) = vel_now(j, 0);
                }
            }

            Eigen::VectorXd res(y);
            res.setZero();
            hnls.solve(y);
            newton_object_sup.operator()(y, res);
            newton_object_sup.y_old = y;

            if (res.norm() < 1e-5)
            {
                std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                          hnls.iter << " iterations " << def << std::endl << std::endl;
            }
            else
            {
                std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                          hnls.iter << " iterations " << def << std::endl << std::endl;
            }

            volVectorField U_rec("U_rec", Umodes[0] * 0);

            for (int j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * y(j);
            }

            for (int k = 0; k < problem->inletIndex.rows(); k++)
            {
                int BCind = problem->inletIndex(k, 0);
                int BCcomp = problem->inletIndex(k, 1);
                valBC(k, i) = U_rec.boundaryFieldRef()[BCind][0].component(BCcomp);
            }
        }

        for (int j = 0; j < N_BC; j++)
        {
            diffvel(j) = abs(abs(vel_now(j, timeStepPenalty - 1)) - abs(valBC(j,
                             timeStepPenalty - 1)));
        }

        std::cout << "max error: " << diffvel.maxCoeff() << std::endl;
        // Count the number of iterations
        Iter ++;
    }

    std::cout << "Final penalty factor(s): " << tauIter << std::endl;
    std::cout << "Iterations: " << Iter - 1 << std::endl;
    return tauIter;
}

Eigen::MatrixXd reducedUnsteadyNS::penalty_PPE(Eigen::MatrixXd& vel_now,
        Eigen::MatrixXd& tauIter,
        int startSnap)
{
    // Initialize new value on boundaries
    Eigen::MatrixXd valBC = Eigen::MatrixXd::Zero(N_BC, timeStepPenalty);
    // Initialize old values on boundaries
    Eigen::MatrixXd valBC0 = Eigen::MatrixXd::Zero(N_BC, timeStepPenalty);
    int Iter = 0;
    Eigen::VectorXd diffvel =  (vel_now.col(timeStepPenalty - 1) - valBC.col(
                                    timeStepPenalty - 1));
    diffvel = diffvel.cwiseAbs();

    while (diffvel.maxCoeff() > tolerancePenalty && Iter < maxIterPenalty)
    {
        if ((valBC.col(timeStepPenalty - 1) - valBC0.col(timeStepPenalty - 1)).sum() !=
                0)
        {
            for (int j = 0; j < N_BC; j++)
            {
                tauIter(j, 0) = tauIter(j, 0) * diffvel(j) / tolerancePenalty;
            }
        }

        std::cout << "Solving for penalty factor(s): " << tauIter << std::endl;
        std::cout << "number of iterations: " << Iter << std::endl;
        //  Set the old boundary value to the current value
        valBC0  = valBC;
        y.resize(Nphi_u + Nphi_p, 1);
        y.setZero();
        y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[startSnap],
                         Umodes);
        y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[startSnap],
                         Pmodes);
        // Set some properties of the newton object
        newton_object_PPE.nu = nu;
        newton_object_PPE.y_old = y;
        newton_object_PPE.yOldOld = newton_object_PPE.y_old;
        newton_object_PPE.dt = dt;
        newton_object_PPE.BC.resize(N_BC);
        newton_object_PPE.tauU = tauIter;

        // Set boundary conditions
        for (int j = 0; j < N_BC; j++)
        {
            newton_object_PPE.BC(j) = vel_now(j, 0);
        }

        // Create nonlinear solver object
        Eigen::HybridNonLinearSolver<newton_unsteadyNS_PPE> hnls(newton_object_PPE);
        // Set output colors for fancy output
        Color::Modifier red(Color::FG_RED);
        Color::Modifier green(Color::FG_GREEN);
        Color::Modifier def(Color::FG_DEFAULT);
        // Set initially for convergence check
        Eigen::VectorXd res(y);
        res.setZero();

        // Start the time loop
        for (int i = 1; i < timeStepPenalty; i++)
        {
            // Set boundary conditions
            for (int j = 0; j < N_BC; j++)
            {
                if (problem->timedepbcMethod == "yes" )
                {
                    newton_object_PPE.BC(j) = vel_now(j, i);
                }
                else
                {
                    newton_object_PPE.BC(j) = vel_now(j, 0);
                }
            }

            Eigen::VectorXd res(y);
            res.setZero();
            hnls.solve(y);
            newton_object_PPE.operator()(y, res);
            newton_object_PPE.yOldOld = newton_object_PPE.y_old;
            newton_object_PPE.y_old = y;

            if (res.norm() < 1e-5)
            {
                std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                          hnls.iter << " iterations " << def << std::endl << std::endl;
            }
            else
            {
                std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                          hnls.iter << " iterations " << def << std::endl << std::endl;
            }

            volVectorField U_rec("U_rec", Umodes[0] * 0);

            for (int j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * y(j);
            }

            for (int k = 0; k < problem->inletIndex.rows(); k++)
            {
                int BCind = problem->inletIndex(k, 0);
                int BCcomp = problem->inletIndex(k, 1);
                valBC(k, i) = U_rec.boundaryFieldRef()[BCind][0].component(BCcomp);
            }
        }

        for (int j = 0; j < N_BC; j++)
        {
            diffvel(j) = abs(abs(vel_now(j, timeStepPenalty - 1)) - abs(valBC(j,
                             timeStepPenalty - 1)));
        }

        std::cout << "max error: " << diffvel.maxCoeff() << std::endl;
        // Count the number of iterations
        Iter ++;
    }

    std::cout << "Final penalty factor(s): " << tauIter << std::endl;
    std::cout << "Iterations: " << Iter - 1 << std::endl;
    return tauIter;
}

void reducedUnsteadyNS::reconstruct(bool exportFields, fileName folder)
{
    if (exportFields)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    int counter = 0;
    int nextwrite = 0;
    List < Eigen::MatrixXd> CoeffU;
    List < Eigen::MatrixXd> CoeffP;
    List <double> tValues;
    CoeffU.resize(0);
    CoeffP.resize(0);
    tValues.resize(0);
    int exportEveryIndex = round(exportEvery / storeEvery);

    std::cout << "Nphi_u in reconstruct: " << Nphi_u << std::endl;
    std::cout << "Nphi_p in recostruct: " << Nphi_p << std::endl;

    for (int i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            Eigen::MatrixXd currentUCoeff;
            Eigen::MatrixXd currentPCoeff;
            currentUCoeff = online_solution[i].block(1, 0, Nphi_u, 1);
	    //std::cout << "currentUcoeff: " << currentUCoeff.size() << std::endl;
	    //currentUCoeff(0) = 1;
            currentPCoeff = online_solution[i].bottomRows(Nphi_p);
	    //std::cout << "Nphi_p in recostruct: " << Nphi_p << std::endl;
            CoeffU.append(currentUCoeff);
            CoeffP.append(currentPCoeff);
            nextwrite += exportEveryIndex;
            double timeNow = online_solution[i](0, 0);
            tValues.append(timeNow);
        }

        counter++;
    }

    volVectorField uRec("uRec", Umodes[0] * 0);
    volScalarField pRec("pRec", Pmodes[0] * 0);

////////////////////////////////////////////////////////////////
    int counter2 = 1;
   
    cout << "Pmodes size dentro reconstruc: " << Pmodes.size() << "\n";
    cout << "Umodes size dentro reconstruc: " << Umodes.size() << "\n";
    cout << "LUSUPmodes size dentro reconstruc: " << problem->L_U_SUPmodes.size() << "\n";
    
    volVectorField U_rec("U_rec", problem->L_U_SUPmodes[0] * 0);
    volScalarField P_rec("P_rec", Pmodes[0] * 0);

    for (label i = 0; i < online_solution.size(); i++)
    {
	//cout << "counter: " << counter << "\n";
	//cout << "nextwrite: " << nextwrite << "\n";
        if (counter == nextwrite)
        {
    		volScalarField P_rec("P_rec", Pmodes[0] * 0);
		volScalarField P_modes0("P_modes0", Pmodes[0] * 0);
		volScalarField P_modes1("P_modes1", Pmodes[0] * 0);
		volScalarField P_modes2("P_modes2", Pmodes[0] * 0);
            	//Info << "PROVAAAA    " << P_rec.size();
		//cout << "P_rec.size: " << P_rec.size();
            	for (label j = 0; j < Nphi_p; j++)
            	{
                	label k = j + Nphi_u;
                	//Info<< "PROVAAAAA    "<<Pmodes.size();
                	P_rec += Pmodes[j]*online_solution[i](k + 1, 0);
            	}
         	ITHACAstream::exportSolution(P_rec, name(counter2), folder);
		//cout << "EXPORT P_rec " << "\n";
		//P_modes0 = Pmodes[0];
                //P_modes1 = Pmodes[1];
                //P_modes2 = Pmodes[2];                
                //ITHACAstream::exportSolution(P_modes0,  name(counter2), folder);
                //ITHACAstream::exportSolution(P_modes1,  name(counter2), folder);
                //ITHACAstream::exportSolution(P_modes2,  name(counter2), folder);



		//ITHACAstream::exportFields(P_rec, folder,"P_rec");
		//////////////////////////////////////////////////
		volVectorField U_rec("U_rec", problem->L_U_SUPmodes[0] * 0);
		volVectorField U_modes0("U_modes0", problem->L_U_SUPmodes[0] * 0);
    		volVectorField U_modes1("U_modes1", problem->L_U_SUPmodes[0] * 0);
		volVectorField U_modes2("U_modes2", problem->L_U_SUPmodes[0] * 0);
	        volVectorField U_sup1("U_sup1", problem->L_U_SUPmodes[0] * 0);
	        volVectorField U_sup2("U_sup2", problem->L_U_SUPmodes[0] * 0);	
		for (label j = 0; j < Nphi_u; j++)
            	{
                	//U_rec += Umodes[j] * online_solution[i](j + 1, 0);
                	//label k = j;
                 	//label k = j + Nphi_u;// + Nphi_p;
                 	//Uevolve_rec += (1-5e-3)* Uevolvemodes[j] * online_solution[i](j + 1, 0) + 5e-3*Umodes[j] * online_solution[i](k + 1, 0);
                 	U_rec += problem->L_U_SUPmodes[j] * online_solution[i](j + 1, 0);
			//cout << "time:" << i << "\n";
			//cout << "red_coeff recostruct:" << online_solution[i](j + 1, 0) << "\n";
                 	//U_modes += Umodes[j];
            	}
		//U_modes0 = problem->L_U_SUPmodes[0];
		//U_modes1 = problem->L_U_SUPmodes[1];
		//U_modes2 = problem->L_U_SUPmodes[2];
		//U_sup1 = problem->L_U_SUPmodes[3];
		//U_sup2 = problem->L_U_SUPmodes[4];
    		ITHACAstream::exportSolution(U_rec,  name(counter2), folder);
		//ITHACAstream::exportSolution(U_modes0,  name(counter2), folder);
		//ITHACAstream::exportSolution(U_modes1,  name(counter2), folder);
		//ITHACAstream::exportSolution(U_modes2,  name(counter2), folder);
		//ITHACAstream::exportSolution(U_sup1,  name(counter2), folder);
		//ITHACAstream::exportSolution(U_sup2,  name(counter2), folder);

	
		pRecFields.append((P_rec).clone());
		uRecFields.append((U_rec).clone());
	}
	counter2++;
    }
 
 //   for (label i = 0; i < problem->Ufield.size(); i++)
//    {
//    	Eigen::VectorXd Uproj1 = Foam2Eigen::projectField(problem->Ufield[i], problem->L_U_SUPmodes, Nphi_u);
//	cout << "Uproj1" << Uproj1 << "\n";	
    //}
    ///////////////////////////////////////////////////////////////

    //uRecFields = problem->L_U_SUPmodes.reconstruct(uRec, CoeffU, "uRec");
    //pRecFields = problem->Pmodes.reconstruct(pRec, CoeffP, "pRec");
   
    //std::cout << "Pmodes[0][5]: " << Pmodes[0][5] << std::endl;

    if (exportFields)
    {
        ITHACAstream::exportFields(uRecFields, folder,
                                   "uRec");
        ITHACAstream::exportFields(pRecFields, folder,
                                   "pRec");
	ITHACAstream::exportFields(problem->L_U_SUPmodes, folder,
                                   "L_U_SUP_modes");
        ITHACAstream::exportFields(problem->Pmodes, folder,
                                   "problem_P_modes");
        ITHACAstream::exportFields(Pmodes, folder,
                                   "Pmodes");
    }
}

//originale
/*void reducedUnsteadyNS::reconstruct(bool exportFields, fileName folder)
{
    if (exportFields)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    int counter = 0;
    int nextwrite = 0;
    List < Eigen::MatrixXd> CoeffU;
    List < Eigen::MatrixXd> CoeffP;
    List <double> tValues;
    CoeffU.resize(0);
    CoeffP.resize(0);
    tValues.resize(0);
    int exportEveryIndex = round(exportEvery / storeEvery);

    for (int i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            Eigen::MatrixXd currentUCoeff;
            Eigen::MatrixXd currentPCoeff;
            currentUCoeff = online_solution[i].block(1, 0, Nphi_u, 1);
            currentPCoeff = online_solution[i].bottomRows(Nphi_p);
            CoeffU.append(currentUCoeff);
            CoeffP.append(currentPCoeff);
            nextwrite += exportEveryIndex;
            double timeNow = online_solution[i](0, 0);
            tValues.append(timeNow);
        }

        counter++;
    }

    volVectorField uRec("uRec", Umodes[0] * 0);
    volScalarField pRec("pRec", problem->Pmodes[0] * 0);
    uRecFields = problem->L_U_SUPmodes.reconstruct(uRec, CoeffU, "uRec");
    pRecFields = problem->Pmodes.reconstruct(pRec, CoeffP, "pRec");

    if (exportFields)
    {
        ITHACAstream::exportFields(uRecFields, folder,
                                   "uRec");
        ITHACAstream::exportFields(pRecFields, folder,
                                   "pRec");
    }
}*/


Eigen::MatrixXd reducedUnsteadyNS::setOnlineVelocity(Eigen::MatrixXd vel)
{
    assert(problem->inletIndex.rows() == vel.rows()
           && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
    Eigen::MatrixXd vel_scal;
    vel_scal.resize(vel.rows(), vel.cols());

    for (int k = 0; k < problem->inletIndex.rows(); k++)
    {
        int p = problem->inletIndex(k, 0);
        int l = problem->inletIndex(k, 1);
        scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
        scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
                           problem->liftfield[k].boundaryField()[p]).component(l) / area;
	//cout << "u_lf velocity:" << u_lf << "\n";
	for (int i = 0; i < vel.cols(); i++)
        {
            vel_scal(k, i) = vel(k, i)/ u_lf;
        }
    }

    return vel_scal;
}

Eigen::MatrixXd reducedUnsteadyNS::setOnlinePressure(Eigen::MatrixXd pressure)
{
    assert(problem->outletIndex.rows() == pressure.rows()
            && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
    Eigen::MatrixXd pressure_scal;
    pressure_scal.resize(pressure.rows(), pressure.cols());
    for (int k = 0; k < problem->outletIndex.rows(); k++)
    {
            int p = problem->outletIndex(k, 0);
            //int l = problem->outletIndex(k, 1);
            scalar area = gSum(problem->liftfieldP[0].mesh().magSf().boundaryField()[p]);
            scalar u_lf = gSum(problem->liftfieldP[k].mesh().magSf().boundaryField()[p] *
                                problem->liftfieldP[k].boundaryField()[p]) / area;
	    //cout << "u_lf pressure:" << u_lf << "\n";
            for (int i = 0; i < pressure.cols(); i++)
            {
                pressure_scal(k, i) = pressure(k, i) / u_lf;
            }
    }

   return pressure_scal;
}
//************************************************************************* //
