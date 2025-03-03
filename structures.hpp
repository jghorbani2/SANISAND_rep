//
//  structures.hpp
//  umat_for_falcon
//
//  Created by Javad Ghorbani on 14/11/2024.
//

#ifndef structures_hpp
#define structures_hpp

#include <stdio.h>
enum StateVariable {
    StressXX,
    StressYY,
    StressZZ,
    StressZY,
    StressZX,
    StressXY,
    DStrainXX,
    DStrainYY,
    DStrainZZ,
    DStrainZY,
    DStrainZX,
    DStrainXY,
    PoreWaterPressure,
    PoreAirPressure,
    InitialPoreWaterPressure,
    InitialPoreAirPressure,
    VoidRatio,
    DegreeOfSaturation,
    PermW_XX,
    PermW_YY,
    PermW_ZZ,
    PermW_ZY,
    PermW_ZX,
    PermW_XY,
    PermA_XX,
    PermA_YY,
    PermA_ZZ,
    PermA_ZY,
    PermA_ZX,
    PermA_XY,
    Xi,
    Damping,
    NumVariables // Always keep this last to represent the count of variables
};
#endif /* structures_hpp */
