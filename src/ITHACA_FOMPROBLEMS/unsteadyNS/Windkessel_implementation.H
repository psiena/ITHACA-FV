        // time step
        scalar delta_t = runTime.deltaT().value();
	cout << "Delta t Windkessel: " << delta_t << "\n";	

        // Calculation of the outlet flowrates
        scalar Q_outlet = gSum(phi.boundaryField()[outlet]);
	//scalar Q_da = gSum(phi.boundaryField()[da]);
        //scalar Q_sa = gSum(phi.boundaryField()[sa]);
        //scalar Q_ffa = gSum(phi.boundaryField()[ffa]);      
        //scalar Q_ta = gSum(phi.boundaryField()[ta]);
        //scalar Q_fsa = gSum(phi.boundaryField()[fsa]);

        // Calculation of the outlet pressures
        p_outlet_2 = (1 - delta_t/Rd_outlet/C_outlet)*p_outlet_2 + delta_t*Q_outlet/C_outlet;
	//p_da_2 = (1 - delta_t/Rd_da/C_da)*p_da_2 + delta_t*Q_da/C_da;
        //p_sa_2 = (1 - delta_t/Rd_sa/C_sa)*p_sa_2 + delta_t*Q_sa/C_sa;
        //p_ffa_2 = (1 - delta_t/Rd_ffa/C_ffa)*p_ffa_2 + delta_t*Q_ffa/C_ffa;
        //p_ta_2 = (1 - delta_t/Rd_ta/C_ta)*p_ta_2 + delta_t*Q_ta/C_ta;
        //p_fsa_2 = (1 - delta_t/Rd_fsa/C_fsa)*p_fsa_2 + delta_t*Q_fsa/C_fsa;

        p_outlet = p_outlet_2 + Rp_outlet*Q_outlet;
	//p_da = p_da_2 + Rp_da*Q_da;
        //p_sa = p_sa_2 + Rp_sa*Q_sa;
        //p_ffa = p_ffa_2 + Rp_ffa*Q_ffa;
        //p_ta = p_ta_2 + Rp_ta*Q_ta;
        //p_fsa = p_fsa_2 + Rp_fsa*Q_fsa;

	cout << "p_outlet: " << p_outlet/density << "\n";

        //pRCR imposition
        
	forAll (p.boundaryField()[outlet], face)
        {
           p.boundaryFieldRef()[outlet][face] = p_outlet/density; 
        }
	
	//salva in file txt riga per riga
	Pfield_out.push_back(p_outlet/density);
	//for (label i = 0; i < outletPatch.rows(); i++)
	//{
         //  timeBCoff_P(i, counter2) = p_outlet/density;
        //}

	/*forAll (p.boundaryField()[da], face)
        {
           p.boundaryFieldRef()[da][face] = p_da/density; 
        }

        forAll (p.boundaryField()[sa], face)
        {
           p.boundaryFieldRef()[sa][face] = p_sa/density;
        }

        forAll (p.boundaryField()[ffa], face)
        {
           p.boundaryFieldRef()[ffa][face] = p_ffa/density;
        }
        forAll (p.boundaryField()[ta], face)
        {
           p.boundaryFieldRef()[ta][face] = p_ta/density;
        }
        forAll (p.boundaryField()[fsa], face)
        {
           p.boundaryFieldRef()[fsa][face] = p_fsa/density; 
        }*/
