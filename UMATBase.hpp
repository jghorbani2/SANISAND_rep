//
//  UMATBase.hpp
//
//  Created by Javad Ghorbani on 14/11/2024.
//

// UMATBase.hpp
#ifndef UMAT_BASE_HPP
#define UMAT_BASE_HPP

#include <vector>
#include <unordered_map>
#include <string>
#include "structures.hpp" // for StateVariable enum
class UMATBase {
public:
    virtual ~UMATBase() = default;

    struct InputData {
        double strainIncrement[6];
        double velocityIncrement[6];
        double poreWaterPressureIncrement;
        double poreAirPressureIncrement;
        double saturationIncrement;

        // Fixed-size vector for standard state variables
        std::vector<double> stateVariables;
        // Flexible container for custom user-defined state variables
        std::unordered_map<std::string, double> customStateVariables;
        // New fields for element number and Gauss coordinates
        int elementNumber;                 // Element number
        std::vector<double> gaussCoords;  // Gauss point coordinates (e.g., [x, y, z])

        InputData() : stateVariables(StateVariable::NumVariables, 0.0) {}
    };

    struct OutputData {
        double stressIncrement[6];
        double stressStrainMatrix[6][6];

        // Fixed-size vector for standard state variables
        std::vector<double> updatedStateVariables;
        // Flexible container for custom user-defined state variables
        std::unordered_map<std::string, double> updatedCustomStateVariables;

        OutputData() : updatedStateVariables(StateVariable::NumVariables, 0.0) {}
    };

    virtual void calculateStressIncrement(const InputData& inputData, OutputData& outputData) = 0;
    virtual void computeStressStrainMatrix(const InputData& inputData, OutputData& outputData) = 0;
};

#endif // UMAT_BASE_HPP
