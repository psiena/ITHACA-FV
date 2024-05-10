    //Label for outlet patches
    const label outlet = mesh.boundaryMesh().findPatchID("outlet");
    //const label da = mesh.boundaryMesh().findPatchID("outlet_aorta");
    //const label sa = mesh.boundaryMesh().findPatchID("outlet_carotide_sx");
    //const label ia = mesh.boundaryMesh().findPatchID("inlet_aorta");
    //const label wa = mesh.boundaryMesh().findPatchID("aorta");
    //const label ffa = mesh.boundaryMesh().findPatchID("outlet_succlavia_dx");
    //const label ta = mesh.boundaryMesh().findPatchID("outlet_succlavia_sx");
    //const label fsa = mesh.boundaryMesh().findPatchID("outlet_carotide_dx");

    //Pressures initialization
    scalar p_outlet = 0.0;
    //scalar p_da = 0.0;
    //scalar p_sa = 0.0;
    //scalar p_ffa = 0.0;
    //scalar p_ta = 0.0;
    //scalar p_fsa = 0.0;

    scalar p_outlet_2 = 0.0;
    //scalar p_da_2 = 0.0;
    //scalar p_sa_2 = 0.0;
    //scalar p_ffa_2 = 0.0;
    //scalar p_ta_2 = 0.0;
    //scalar p_fsa_2 = 0.0;

    //Calculation of outlets areas
    const surfaceScalarField& magSf = mesh.magSf(); 

    scalar A_outlet = gSum(magSf.boundaryField()[outlet]);
    //scalar A_da = gSum(magSf.boundaryField()[da]);
    //scalar A_sa = gSum(magSf.boundaryField()[sa]);
    //scalar A_ia = gSum(magSf.boundaryField()[ia]);
    //scalar A_ffa = gSum(magSf.boundaryField()[ffa]);
    //scalar A_ta = gSum(magSf.boundaryField()[ta]);
    //scalar A_fsa = gSum(magSf.boundaryField()[fsa]);
    
    scalar A_tot = A_outlet;//A_da + A_sa + A_ffa + A_ta + A_fsa; 

    //Calculation of Rp, Rd e C for each RCR outlet
    scalar R_outlet = (A_tot/A_outlet)*Rtot;
    //scalar R_da = (A_tot/A_da)*Rtot;
    //scalar R_sa = (A_tot/A_sa)*Rtot;
    //scalar R_ffa = (A_tot/A_ffa)*Rtot;
    //scalar R_ta = (A_tot/A_ta)*Rtot;
    //scalar R_fsa = (A_tot/A_fsa)*Rtot;

    scalar Rp_outlet = 0.1e5;
    //scalar Rp_da = R_da/17.85;
    //scalar Rp_sa = R_sa/17.85;    
    //scalar Rp_ffa = R_ffa/17.85;
    //scalar Rp_ta = R_ta/17.85;   
    //scalar Rp_fsa = R_fsa/17.85;
    
    scalar Rd_outlet = R_outlet - Rp_outlet;
    //scalar Rd_da = R_da - Rp_da;
    //scalar Rd_sa = R_sa - Rp_sa;
    //scalar Rd_ffa = R_ffa - Rp_ffa;    
    //scalar Rd_ta = R_ta - Rp_ta;
    //scalar Rd_fsa = R_fsa - Rp_fsa;
    
    scalar C_outlet = (A_outlet/A_tot)*Ctot;
    //scalar C_da = (A_da/A_tot)*Ctot;
    //scalar C_sa = (A_sa/A_tot)*Ctot;    
    //scalar C_ffa = (A_ffa/A_tot)*Ctot;   
    //scalar C_ta = (A_ta/A_tot)*Ctot;    
    //scalar C_fsa = (A_fsa/A_tot)*Ctot;
    
    //Cardiac Cycle period
    scalar Tau_ = Tau;

    //Save data
    std::ofstream file;
    file.open ("results.txt", std::ofstream::out | std::ofstream::app);

    //int indice = 0;
    //int divid = 15;

 

    /*if (Pstream::master())
    {
        file << "t" << "\t" << "p"<< "\t" << "p_outlet" << "\t" << "U" << std::endl << "\n";
    }*/
 
    
