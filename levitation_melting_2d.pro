Include "coil.geo";
Include "stale.geo";

Group {
  Charge = Region[CHARGE_IN];
  Air = Region[AIR_IN];
  AirInf = Region[{}];

  // three-phase current
  // Ind1 = Region[{(100:100+number_of_turns-1:3)}];
  // Ind2 = Region[{(101:100+number_of_turns-1:3)}];
  // Ind3 = Region[{(102:100+number_of_turns-1:3)}];

  //one-phase current
  Ind1 = Region[{(IND_IN:IND_IN+number_of_turns-1)}];
  
  If (number_of_turns2>0)
    Ind2 = Region[{(IND_IN+number_of_turns:IND_IN+number_of_turns+number_of_turns2-1)}];
  Else
    Ind2=Region[{}];
  EndIf

  If (number_of_turns3>0)
    Ind3 = Region[{(IND_IN+number_of_turns+number_of_turns2:IND_IN+number_of_turns+number_of_turns2+number_of_turns3-1)}];
  Else
    Ind3=Region[{}];
  EndIf

  If (number_of_turns4>0)
    Ind4 = Region[{(IND_IN+number_of_turns+number_of_turns2+number_of_turns3:IND_IN+number_of_turns+number_of_turns2+number_of_turns3+number_of_turns4-1)}];
  Else
    Ind4=Region[{}];
  EndIf


  Ind = Region[{Ind1, Ind2, Ind3,Ind4}];

  Surface_ht0 = Region[{}];
  Surface_bn0 = Region[{AXIS}];
  Surface_Inf = Region[{INF}];
  
  //rint = coil_x + 6 * r; // internal radius of infinite region
  //rext = 1.2 * rint; // external radius of infinite region

  Vol_Mag = Region[{Air, Charge, Ind, AirInf}]; // full magnetic domain
  Vol_C_Mag = Region[{Charge, Ind}]; // massive conductors
  Vol_Inf_Mag = Region[AirInf]; // ring-shaped shell for infinite transformation
  //Val_Rint = rint; Val_Rext = rext; // interior and exterior radii of ring
}

Function {

  Flag_FrequencyDomain = 1;
  Freq = 10000;

  mu0 = 4.e-7 * Pi;
  nu[]  = 1. / mu0;
  
  If (distribution==0)
    sigma[ Charge ] = AluminiumConductivity; 
    density[Charge] = AluminiumDensity;
  ElseIf (distribution==1)
    sigma[ Charge ] = (X[]-charge_R)*(X[]-charge_R)+Y[]*Y[]<=(charge_r+0.0005/2)*(charge_r+0.0005/2)?AluminiumConductivity:AluminiumConductivity*0.00001; 
    density[Charge] = (X[]-charge_R)*(X[]-charge_R)+Y[]*Y[]<=(charge_r+0.0005/2)*(charge_r+0.0005/2)?AluminiumDensity:0;    
  Else   
    sigma[Charge]=ScalarField[XYZ[]];
    density[Charge]=ScalarField[XYZ[]]/AluminiumConductivity*AluminiumDensity;
  EndIf
  
    sigma[ Ind ] =CopperConductivity;
  

  



  // reading node-value from a file when parsing the .pro file:
  

  
  
  
  
  //fluentsigma[] = 3.5e7;

  // For a correct definition of the voltage:
  // Axi-symmetric case:
  Flag_Axi = 1;
  CoefGeos[] = 2*Pi;

  // 2D case:
  //Flag_Axi = 0;
  //CoefGeos[] = 1;
}

Constraint {
  { Name MagneticVectorPotential_2D;
    Case {
      { Region Surface_bn0; Value 0; }
      { Region Surface_Inf; Value 0; }
    }
  }
  { Name Current_2D;
    Case {
      { Region Ind1; Value Current1 * Complex[Sin[0], Cos[0]]; }
      { Region Ind2; Value Current2 * Complex[Sin[0], Cos[0]]; }
      { Region Ind3; Value Current3 * Complex[Sin[0], Cos[0]]; }
      { Region Ind4; Value Current4 * Complex[Sin[0], Cos[0]]; }
    }
  }
  { Name Voltage_2D;
    Case {
      { Region Charge; Value 0; }
    }
  }

}

Include "Lib_Magnetodynamics2D_av_Cir.pro";

Resolution {
  { Name Melting;
    System {
      { Name A; NameOfFormulation Magnetodynamics2D_av;
        Type ComplexValue; Frequency Freq;
      }
    }
    Operation {
//      CreateDirectory[resPath];
      /*
      Generate[A]; Solve[A]; SaveSolution[A];
      PostOperation[dyn];
      */

      SetGlobalSolverOptions["-ksp_type gmres -pc_type lu -ksp_rtol 1e-12"];
      
      If (distribution==2)
      GmshRead["sigma.pos"];
      EndIf

      Generate[A]; Solve[A]; SaveSolution[A];
      
      

     
    }
  }
}

PostOperation {
  { Name small; NameOfPostProcessing Magnetodynamics2D_av;
    Operation {
      Print[ f_DC, OnElementsOf Charge ,  Format SimpleTable, Smoothing, File "fAV.st"] ;
      }
  }


  { Name full; NameOfPostProcessing Magnetodynamics2D_av;
    Operation {
      //Print[ az, OnElementsOf Vol_Mag, File "az.pos" ];
      Print[ b, OnElementsOf Vol_Mag, File "b.pos" ];
      Print[ j, OnElementsOf Vol_Mag, File "j.pos" ];
      
      Print[ f_DC, OnElementsOf Charge ,  File "fAV.pos"] ;
      Print[ f_DC, OnElementsOf Charge ,  File "fAV.st", Smoothing, Format SimpleTable] ;

      //Print[ cond, OnElementsOf Charge ,  File "sigma.pos"] ;
      //Print[ cond, OnElementsOf Charge ,  File "sigma.org"] ;


     
      Print[int_f_DC[Charge], OnRegion Charge, Format Table,Smoothing, File "GlobalForce.txt"];
      Print[int_f_DC_r[Charge], OnRegion Charge, Format Table,Smoothing, File "GlobalForce_r.txt"];
      Print[int_mass[Charge], OnRegion Charge, Format Table, Smoothing,File "GlobalGravity.txt"];
      Print[JouleLosses[Vol_C_Mag], OnRegion Vol_C_Mag, Smoothing,  Format Table, File "GlobalJoule.txt"];
     
    }
  }
}

// Choose resolution, computation and post-operation in interactive mode
DefineConstant [
  R_ = {"Melting", Name "GetDP/1ResolutionChoices", Visible 0},
  C_ = {"-solve -pos -v2", Name "GetDP/9ComputeCommand", Visible 0},
  P_ = {"full", Name "GetDP/2PostOperationChoices", Visible 0}
];
