/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2021 by the ITHACA-FV authors
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
Description
    Example of an unsteady NS Reduction Problem with time-dependent boundary
    conditions.
SourceFiles
    17YJunction.C
\*---------------------------------------------------------------------------*/
#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>
#include "Modes.H"


class tutorial17: public unsteadyNS
{
    public:
        explicit tutorial17(int argc, char* argv[])
            :
            unsteadyNS(argc, argv),
            U(_U()),
            p(_p())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;

        void offlineSolve()
        {
            List<scalar> mu_now(1);
            volVectorField U0 = U;
            volScalarField P0 = p;

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }
            else
            {
                U = U0;
                p = P0;
                mu_now[0] = mu(0, 0);
                truthSolve(mu_now);
            }
        }

        // Method to compute the lifting functions for this tutorial
        void liftSolve(Eigen::MatrixXd BCs, Eigen::MatrixXd BCs_P)
        {
	    volVectorField U = _U();
	    //cout << "U in liftfield:" << U << "\n";
	    
	    volVectorField Ulift("Ulift", U);
	    //surfaceScalarField& phi = _phi();
            for (label k = 0; k < inletPatch.rows(); k++)
            {
                Time& runTime = _runTime();
                surfaceScalarField& phi = _phi();
                fvMesh& mesh = _mesh();
                volScalarField p = _p();
                //volVectorField U = _U();
		volScalarField Pippo_p1("Pippo_p1_nuovo", p);
                volVectorField Pippo_u1("Pippo_u1_nuovo", U);
		IOMRFZoneList& MRF = _MRF();
                label BCind = inletPatch(k, 0);
                //volVectorField Ulift("Ulift_nuovo", U);
                instantList Times = runTime.times();
                runTime.setTime(Times[1], 1);
                pisoControl potentialFlow(mesh, "potentialFlow");
                Info << "Solving a lifting Problem" << endl;
                Vector<double> v1(0, 0, 0);
                v1[2] = BCs(0, k);
                //v1[1] = BCs(1, k);
                Vector<double> v0(0, 0, 0);

                for (label j = 0; j < U.boundaryField().size(); j++)
                {
                    if (j == BCind)
                    {
                        assignBC(Ulift, j, v1);
                    }
                    else if (U.boundaryField()[BCind].type() == "fixedValue")
                    {
                        assignBC(Ulift, j, v0);
                    }
                    else
                    {
                    }

                    assignIF(Ulift, v0);
                    phi = linearInterpolate(Ulift) & mesh.Sf();
                }

                Info << "Constructing velocity potential field Phi\n" << endl;
                volScalarField Phi
                (
                    IOobject
                    (
                        "Phi",
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("Phi", dimLength * dimVelocity, 0),
                    p.boundaryField().types()
		    //cout << "vediamo che e':" << p.boundaryField().types() << "\n";
                );
		//cout << "vediamo che e':" << p.boundaryField().types() << "\n";
                label PhiRefCell = 0;
                scalar PhiRefValue = 0;
                setRefCell
                (
                    Phi,
                    potentialFlow.dict(),
                    PhiRefCell,
                    PhiRefValue
                );
                mesh.setFluxRequired(Phi.name());
                runTime.functionObjects().start();
                MRF.makeRelative(phi);
                adjustPhi(phi, Ulift, p); //SEMBRA NON FACCIA NULLA


                while (potentialFlow.correctNonOrthogonal())
                {
                    fvScalarMatrix PhiEqn
                    (
                        fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
                        ==
                        fvc::div(phi)
                    );
                    PhiEqn.setReference(PhiRefCell, PhiRefValue);
                    PhiEqn.solve();

                    if (potentialFlow.finalNonOrthogonalIter())
                    {
                        phi -= PhiEqn.flux();
                    }
                }

                MRF.makeAbsolute(phi);
                Info << "Continuity error = "
                     << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
                     << endl;
                Ulift = fvc::reconstruct(phi);
                Ulift.correctBoundaryConditions();
                Info << "Interpolated velocity error = "
                     << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
                         / sum(mesh.magSf())).value()
                     << endl;
                Ulift.write();
		Pippo_u1.write();
		Pippo_p1.write();
		Phi.write(); //proviamo cosi
                liftfield.append((Ulift).clone());
		//std::cout << "liftfieldU:" << Ulift[0][10] << "\n";
            }
	    
	    
	    //Per calcolare la plift 
	    for (label k = 0; k < outletPatch.rows(); k++)
	    {
		Info<< nl << "Calculating Plift field" << endl;

                Time& runTime = _runTime();
                surfaceScalarField& phi = _phi();
                fvMesh& mesh = _mesh();
                volScalarField p = _p();
                //volVectorField U = _U();
		//volVectorField Ulift("Ulift" + name(k), U);
                IOMRFZoneList& MRF = _MRF();
                label BC_oind = outletPatch(k, 0);
                volScalarField Plift("Plift", p);
                instantList Times = runTime.times();
                runTime.setTime(Times[1], 1);
                pisoControl potentialFlow(mesh, "potentialFlow");
                //Info << "Solving a lifting Problem" << endl;
                double p1 = 0; //v1(0, 0, 0);
                p1 = BCs_P(0, k); //v1[2] = BCs(0, k);
                //std::cout << "p1:\n" << p1 << "\n";
                double p0 = 0;//v0(0, 0, 0);

                for (label j = 0; j < p.boundaryField().size(); j++)
                {
                    //std::cout << "j:\n" << j << "\n";
                    if (j == BC_oind)
                    {
                        assignBC(Plift, j, p1);
                    }
                    else if (p.boundaryField()[BC_oind].type() == "fixedValue")
                    {
                        assignBC(Plift, j, p0);
                    }
                    else
                    {
                    }

                    assignIF(Plift, p0);
                    phi = linearInterpolate(Ulift) & mesh.Sf();
                }


		label pliftRefCell = 0;
                scalar pliftRefValue = 0.0;
        	setRefCell
       	 	(
            		Plift,
            		potentialFlow.dict(),
            		pliftRefCell,
            		pliftRefValue
        	);
        	
		//std::cout << "Ulift:" << Ulift[0][10] << "\n";	
		// Calculate the flow-direction filter tensor
        	volScalarField magSqrU(magSqr(Ulift));
        	volSymmTensorField F(sqr(Ulift)/(magSqrU + SMALL*average(magSqrU)));
		
        	// Calculate the divergence of the flow-direction filtered div(U*U)
        	// Filtering with the flow-direction generates a more reasonable
        	// pressure distribution in regions of high velocity gradient in the
        	// direction of the flow
        	volScalarField divDivUUlift
        	(
            		fvc::div
           		(
                		F & fvc::div(phi, Ulift),
                		"div(div(phi,Ulift))"
            		)
        	);
		
		// Solve a Poisson equation for the approximate pressure lift
        	while (potentialFlow.correctNonOrthogonal())
        	{
            		fvScalarMatrix pliftEqn
            		(
                		fvm::laplacian(Plift) + divDivUUlift
            		);

            		pliftEqn.setReference(pliftRefCell, pliftRefValue);
            		pliftEqn.solve();
        	}

        	Plift.write();
		liftfieldP.append((Plift).clone());
		
		Info<< nl << "Plift field computed" << endl;

	    }
	}
        	
};


/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial17 object
    tutorial17 example(argc, argv);

    //////////////////////////// 
    double t_i = 0.0;
    double t_f = 0.5;
    double dt = 0.001;
    int nb_time_steps = (t_f - t_i)/dt + 1.0; //(reduced.finalTime - reduced.tstart) / reduced.dt + 1;

    int counter = 0;
    double time = t_i;
    Eigen::MatrixXd provaBC(1, nb_time_steps);
    Eigen::MatrixXd provaBC_P(1, nb_time_steps);
    while (time < t_f + t_f * dt) // - 0.5 * dt)
    {
	provaBC(0, counter) = 0.007957747154594767*std::sin(6.283185*time)*std::sin(6.283185*time);
	provaBC_P(0,counter) = 2;//6*std::sin(6.283185*time)*std::sin(6.283185*time)+2.5*(1-std::exp(-12.566371*time)-std::sin(12.566371*time));
	//std::cout << "time:\n" << time << "\n";
	time = time + dt;
        counter ++;	
    }
    
  
    ITHACAstream::exportMatrix(provaBC_P, "BC_P", "matlab",
		               "./ITHACAoutput/BC_P/");
    cnpy::save(provaBC_P, "./ITHACAoutput/POD/prova_BC_P.npy");
    ITHACAstream::exportMatrix(provaBC, "BC", "matlab",
		                "./ITHACAoutput/BC/");
    cnpy::save(provaBC, "./ITHACAoutput/POD/prova_BC.npy"); 
    ///////////////////////////

    // the offline samples for the boundary conditions
    //word par_offline_BC("./timeBCoff"); 
    // Convert boundary values x-direction and y-direction to x,y,z
    //Eigen::MatrixXd timeBCoff2D = ITHACAstream::readMatrix(par_offline_BC);
    //cout << "BCoff : " << timeBCoff2D.col(0);
    Eigen::MatrixXd timeBCoff3D = Eigen::MatrixXd::Zero(3, nb_time_steps);//timeBCoff2D.cols());
    timeBCoff3D.row(2) = provaBC.row(0); 
    //timeBCoff3D.row(0) = timeBCoff2D.row(0); //cambia
    //timeBCoff3D.row(1) = timeBCoff2D.row(1);
    //timeBCoff3D.row(3) = timeBCoff2D.row(2);
    //timeBCoff3D.row(4) = timeBCoff2D.row(3);
    example.timeBCoff = timeBCoff3D;
    example.timeBCoff_P = provaBC_P;
    Eigen::MatrixXd par_on_BC;
    std::cout << "par_on_BC:\n" << par_on_BC << "\n";
    //word par_online_BC("./timeBCon");
    par_on_BC = provaBC.row(0); //ITHACAstream::readMatrix(par_online_BC);
    
    // Read parameters from ITHACAdict file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 10);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 10);
    int NmodesSUPout = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 10);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 2);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 2);
    int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 2);
    int NmodesOut = para->ITHACAdict->lookupOrDefault<int>("NmodesOut", 10);
    
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.004;//0.004;//0.01;
    example.mu_range(0, 1) = 0.004;//0.004;//0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 2;
    //example.inletIndex(1, 0) = 2;
    //example.inletIndex(1, 1) = 2;
    
    example.outletIndex.resize(1, 2);
    example.outletIndex(0, 0) = 2;
    example.outletIndex(0, 1) = 2;
    
    //example.inletIndex.resize(4, 2); // rows: total number of patches
    //example.inletIndex(0, 0) = 1;  // Patch inlet 1
    //example.inletIndex(0, 1) = 0;  // Patch inlet 1: x-direction
    //example.inletIndex(1, 0) = 1;  // Patch inlet 1: y-direction
    //example.inletIndex(1, 1) = 1;  // Patch inlet 2
    //example.inletIndex(2, 0) = 2;  // Patch inlet 2: x-direction
    //example.inletIndex(2, 1) = 0;  // Patch inlet 2: y-direction
    //example.inletIndex(3, 0) = 2;  // Patch inlet 2: x-direction
    //example.inletIndex(3, 1) = 1;  // Patch inlet 2: y-direction
    example.inletPatch.resize(1, 1);
    example.inletPatch(0, 0) = example.inletIndex(0, 0);  // Patch inlet 1

    example.outletPatch.resize(1, 1);
    example.outletPatch(0, 0) = example.outletIndex(0, 0);//outletIndex(0, 0);
    
    //example.inletPatch(1, 0) = example.inletIndex(1, 0);  // Patch inlet 2
    // Time parameters
    example.startTime = 0.0;
    example.finalTime = 0.5; //12;
    example.timeStep = dt;//0.0005;
    example.writeEvery = 0.01;//0.01;//0.03;
    // Perform The Offline Solve;
    example.offlineSolve();
    example.solvesupremizer("snapshots"); //modes or snapshots

    // Perform POD
    if (example.bcMethod == "lift")
    {
        // Search the lift function (da cambiare???)
        Eigen::MatrixXd BCs;
	Eigen::MatrixXd BCs_P;
        BCs.resize(1, 1);
	BCs_P.resize(1, 1);
        BCs(0, 0) = 1;  // Patch inlet 1
	BCs_P(0, 0) = 1;//2; //Patch outlet
        
	//BCs(1, 0) = -1;  // Patch inlet 2
        //BCs(0, 1) = 0;  // Patch inlet 1
        //BCs(1, 1) = 1;  // Patch inlet 2
        example.liftSolve(BCs, BCs_P);
	//example.liftSolveP(BCs_P);

	//example.solvesupremizer("snapshots");

        // Normalize the lifting function
        //ITHACAutilities::normalizeFields(example.liftfield);
	
        //ITHACAutilities::normalizeFields(example.liftfieldP);	
	
	// Create homogeneous basis functions for velocity
        example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
	example.computeLiftP(example.Pfield, example.liftfieldP, example.Pomfield);
        
	//example.solvesupremizer("snapshots");

        ITHACAstream::exportFields(example.Uomfield, "./ITHACAoutput/Uomfield",
                               "Uomfield");
        ITHACAstream::exportFields(example.Pomfield, "./ITHACAoutput/Pomfield",
		                               "Pomfield");	
	
        // Perform a POD decomposition for velocity and pressure
        ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                            example.podex, 0, 0,
                            NmodesUout);
        ITHACAPOD::getModes(example.Pomfield, example.Pmodes, example._p().name(),
                            example.podex, 0, 0,
                            NmodesPout);
	
	/*for (int i = 0; i < example.Pmodes.size(); i++)
	{
		ITHACAutilities::changeBCtype(example.Pmodes[i], "fixedValue", 2);
		ITHACAutilities::assignBC(example.Pmodes[i], 2, 0);
	}
	Vector<double> v(0, 0, 0);
        for (int i = 0; i < example.Umodes.size(); i++)
        {
                ITHACAutilities::changeBCtype(example.Umodes[i], "fixedValue", 0);
		ITHACAutilities::assignBC(example.Umodes[i], 0, v);
        }*/
	
	//example.solvesupremizer("modes");
	ITHACAPOD::getModes(example.supfield, example.supmodes, example._U().name(),
                            example.podex,
                            example.supex, 1, 
			    NmodesSUPout);
    }

    // Reduced Matrices
    //example.projectPPE("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj); //CREDO CHE QUESTO LO FACCIA MALE
    
    cout << "NmodesUproj" << NmodesUproj << "\n";
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    cout << "projectSUP fatto" << "\n";

    reducedUnsteadyNS reduced(example); 
    // Set values of the online phase
    reduced.nu = 0.004;//0.004; //0.01;
    reduced.tstart = 0.0;
    reduced.finalTime = 0.5;//18;
    reduced.dt = dt;//0.0005;
    reduced.storeEvery = 0.01;//0.01;//0.0005;
    reduced.exportEvery = 0.01;//0.01; //0.03;
    // Set values velocity boundary conditions of the online phase
    Eigen::MatrixXd vel_now = par_on_BC;
    Eigen::MatrixXd P_now = provaBC_P;

    if (example.bcMethod == "penalty")
    {
        // Set values of the iterative penalty method
        reduced.maxIterPenalty = 100;
        reduced.tolerancePenalty = 1e-5;
        reduced.timeStepPenalty = 5;
        // Set initial quess for penalty factors
        reduced.tauIter = Eigen::MatrixXd::Zero(4, 1);
        reduced.tauIter <<  1e-6, 1e-6, 1e-6, 1e-6;
        // Solve for the penalty factors with the iterative solver
        reduced.tauU = reduced.penalty_PPE(vel_now, reduced.tauIter);
    }
    // Set the online temperature BC and solve reduced model
    //reduced.solveOnline_PPE(vel_now, P_now, 1);
    reduced.solveOnline_sup(vel_now, P_now, 1);
    cout << "solveOnline_sup fatto" << "\n";
    reduced.reconstruct(true, "./ITHACAoutput/Reconstruction/");

    // Reduced projection coefficient
    cout << "Pmodes size projection" << example.Pmodes.size() << "\n";
    cout << "Umodes size projection" << example.Umodes.size() << "\n";
    Eigen::MatrixXd Proj_coeff = ITHACAutilities::getCoeffs(example.Pfield, example.Pmodes, NmodesPproj, false);
    cnpy::save(Proj_coeff, "./ITHACAoutput/POD/Proj_coeff_P.npy");
    Eigen::MatrixXd Proj_coeff_U = ITHACAutilities::getCoeffs(example.Ufield, example.Umodes, NmodesUproj, false);
    cnpy::save(Proj_coeff_U, "./ITHACAoutput/POD/Proj_coeff_U.npy");
 
    // Error
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield,
                             reduced.uRecFields);
    Eigen::MatrixXd errL2P =  ITHACAutilities::errorL2Rel(example.Pfield,
                              reduced.pRecFields);
    ITHACAstream::exportMatrix(errL2U, "errL2U", "python",
                               "./ITHACAoutput/ErrorsL2U/");
    ITHACAstream::exportMatrix(errL2P, "errL2P", "python",
                               "./ITHACAoutput/ErrorsL2P/");
    
    //Projection error per la pressione
    /*Eigen::MatrixXd Pressure = Foam2Eigen::PtrList2Eigen(example.Pfield);*/
    //Eigen::VectorXd Pproj1 = Foam2Eigen::projectField(example.Ufield[0], example.Umodes, NmodesUproj);

    // Project the full order solution onto the POD space
    PtrList<volVectorField> Umodes_proj_err;
    PtrList<volScalarField> Pmodes_proj_err;

    ITHACAPOD::getModes(example.Uomfield, Umodes_proj_err, example._U().name(),
                        example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pomfield, Pmodes_proj_err, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);

    example.Umodes[0] = example.liftfield[0];
    cout << "lift appesa" << "\n";
    for (int i = 0; i < Umodes_proj_err.size()-1; i++)
    {
        example.Umodes[i+1] = Umodes_proj_err[i];
    }
    cout << "modi appesi" << "\n";
    
    example.Pmodes[0] = example.liftfieldP[0];
    for (int i = 0; i < Pmodes_proj_err.size()-1; i++)
    {
        example.Pmodes[i+1] = Pmodes_proj_err[i];
    }
    
    PtrList<volScalarField> Pproj;

    example.Pmodes.projectSnapshots(example.Pfield, Pproj, NmodesPproj);
    //example.Pmodes.projectSnapshot(example.Pfield[0], NmodesPproj);
    Eigen::MatrixXd errL2P_proj = ITHACAutilities::errorL2Rel(example.Pfield,
                              Pproj);
    ITHACAstream::exportMatrix(errL2P_proj, "errL2P_proj", "python",
                               "./ITHACAoutput/ErrorsL2P_proj/");
    
    // e per la velocita
    PtrList<volVectorField> Uproj;
    example.L_U_SUPmodes.projectSnapshots(example.Ufield, Uproj, NmodesUproj);
    cout << "projectsnap fatto" << "\n";
    Eigen::MatrixXd errL2U_proj = ITHACAutilities::errorL2Rel(example.Ufield,
                              Uproj);
    cout << "err calcolato" << "\n";
    ITHACAstream::exportMatrix(errL2U_proj, "errL2U_proj", "python",
                               "./ITHACAoutput/ErrorsL2U_proj/");

    exit(0);
}

/// \dir 17YJunction Folder of the turorial 17
/// \file
/// \brief Implementation of tutorial 17 for an unsteady Navier-Stokes problem with time-dependent inlet boundary conditions
///
/// \example 17YJunction.C
/// \section intro_UnsteadyNS Introduction to tutorial 17
/// In this tutorial, we contruct a reduced order model for a Y-junction flow problem.
/// The Y-junction consists of two inlets and one outlet channel whose time-dependent
/// inlet boundary conditions are controlled. The angle between each inlet and the horizontal axis is 60 degrees.
/// The length of the channels is 2 m. The two inlets, \f$\Gamma_{i1}\f$ and \f$\Gamma_{i2}\f$, have a width of 0.5 m, while the outlet,
/// \f$\Gamma_{o}\f$, has a width of 1 m. The kinematic viscosity is equal to \f$\nu\f$ = 0.01 m\f$^2\f$/s and the flow is considered laminar.
/// As initial conditions the steady state solution, obtained with the simpleFOAM-solver, for a velocity magnitude of 1 m/s at both inlets is chosen.
/// In this tutorial, new values of the velocity magnitude of the flow at the inlets are imposed in the reduced order model with either an
/// iterative penalty method or a lifting function method.
/// The reduced order model is constructed using Galerkin projection approach together with an exploitation of a pressure Poisson equation
/// during the projection stage.
///
///
/// The following image depicts a sketch of the geometry of the two-dimensional Y-junction.
/// \image html YJunction.png
///
/// \section code17 A detailed look into the code
///
/// In this section we explain the main steps necessary to construct the tutorial N°17.
///
/// \subsection header ITHACA-FV header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
///
/// The header files of ITHACA-FV necessary for this tutorial are: <unsteadyNS.H> for the full order unsteady NS problem,
/// <ITHACAPOD.H> for the POD decomposition, <reducedUnsteadyNS.H> for the construction of the reduced order problem,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
///
/// \dontinclude 17YJunction.C
/// \skip unsteadyNS
/// \until ITHACAstream
///
/// \subsection classtutorial17 Definition of the tutorial17 class
///
/// We define the tutorial17 class as a child of the unsteadyNS class.
/// \skipline tutorial17
/// \until volScalarField& p;
///
/// Inside the tutorial17 class we define the offlineSolve method according to the
/// specific problem that needs to be solved. If the offline solve has
/// been previously performed then the method just reads the existing snapshots from the Offline directory.
/// Otherwise it performs the offline solve.
///
/// \skipline offlineSolve
/// \until }
/// \skipline else
/// \until }
/// \skipline }
///
/// We also define a liftSolve method, which will be needed to construct the lifting functions
/// if the lifting function method is used:
///
/// \skipline liftSolve
/// \until liftfield.append((Ulifty).clone());
/// \skipline }
/// \skipline }
///
/// \subsection main Definition of the main function
///
/// In this section we show the definition of the main function.
/// First we construct the object "example" of type tutorial17:
///
/// \skipline example
///
/// Next, we specify the name of the file that contains the values at the boundary for the full order solver for each time step.
/// The inlet velocity magnitude of, alternately, inlet 1 or 2 is increased or decreased linearly over time between 1.0 m/s to 0.5 m/s.
///
/// \skipline par_offline_BC
///
/// The input file only contains the boundary values at inlet 1 and 2 in the x-direction and y-direction.
/// We are storing these values in an matrix together with zero velocity in the z-direction since it is required to specify the
/// velocity in the x,y and z direction even though the problem is two-dimensional.
///
/// \skipline timeBCoff2D
/// \until example.timeBCoff
///
/// We also load the new inlet velocities for the reduced order solver from a text file and store the values in a matrix.
///
/// \skipline  Eigen::MatrixXd par_on_BC
/// \until par_on_BC = ITHACAstream::readMatrix(par_online_BC)
///
/// Then we parse the ITHACAdict file to determine the number of modes
/// to be written out and also the ones to be used for projection of the velocity and pressure:
/// \skipline ITHACAparameters
/// \until NmodesSUPproj
///
/// we note that a default value can be assigned in case the parser did
/// not find the corresponding string in the ITHACAdict file.
///
/// In our implementation, the viscocity needs to be defined by specifying that
/// Nparameters=1, Nsamples=1, and the parameter ranges from 0.01 to 0.01 equispaced, i.e.
///
/// \skipline example.Pnumber
/// \until example.genEquiPar()
///
/// After that we set the inlet boundaries where we have the non homogeneous BC
///
/// \skipline example.inlet
/// \until example.inletIndex(3, 1) = 1;
///
/// as well as the patch number of inlet 1 and the patch number of inlet 2.
/// \skipline example.inletPatch.resize(2, 1)
/// \until example.inletPatch(1, 0) = example.inletIndex(2, 0)
///
/// We set the parameters for the time integration, so as to simulate 12.0 seconds of simulation time,
/// with a step size = 0.0005 seconds, and the data are dumped every 0.03 seconds, i.e.
///
/// \skipline example.startTime
/// \until example.writeEvery
///
/// Now we are ready to perform the offline stage:
///
/// \skipline Solve()
///
/// After that, the modes for the velocity and pressure are obtained.
/// The lifting function method requires the velocity modes to be homogeneous.
/// Therefore we differentiate here between the lifting function method and the penalty method.
/// Which bcMethod to be used can be specified in the ITHACAdict file.
///
/// In the case of the lifting function method, we specify first a unit value onto the inlet patches in the x- and y-direction
/// \skipline if
/// \until BCs(1, 1) = 1
///
/// Then we solve for the lifting functions; one lifting function for each inlet boundary condition and direction
/// \skipline liftSolve
///
/// The lifting functions are normalized
/// \skipline normalizeFields
///
/// and the velocity snapshots are homogenized by subtracting the normalized lifting functions from the velocity snapshots.
/// \skipline computeLift
///
/// Finally, we obtain the homogeneous velocity modes
/// \skipline getModes
/// \until NmodesUout
///
/// and the pressure modes.
/// \skipline getModes
/// \until }
///
/// If the penalty method is used, the velocity and pressure modes are computed as follows
/// \skipline bcMethod
/// \until }
///
/// Then the projection onto the POD modes is performed using the Pressure Poisson Equation (PPE) approach.
/// \skipline projectPPE
///
/// Now that we obtained all the necessary information from the POD decomposition and the reduced matrices,
/// we are ready to construct the dynamical system for the reduced order model (ROM). We proceed
/// by constructing the object "reduced" of type reducedUnsteadyNS:
///
/// \skipline reducedUnsteadyNS
///
/// And then we can use the constructed ROM to perform the online procedure, from which we can simulate the
/// problem for new values of the inlet velocities that vary linearly over time.
/// We are keeping the time stepping the same as for the offline stage,
/// but the ROM simulation (online stage) is performed for a longer time period.
///
/// \skipline reduced.nu
/// \until exportEvery
///
/// We have to specify a new value for the inlet velocities, which were already loaded from a text file
///
/// \skipline Eigen::
///
/// If the penalty method is used, we also have to set the parameters for the iterative penalty method to determine a suitable penalty factor
///
/// \skipline if (example.bcMethod == "penalty")
/// \until }
///
/// and then the online solve is performed.
///
/// \skipline solveOnline
///
/// Finally the ROM solution is reconstructed.
/// In the case the solution should be exported and exported, put true instead of false in the function:
///
/// \skipline reconstruct
///
/// \section plaincode The plain program
/// Here there's the plain code
///




