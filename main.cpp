//
//  Created by Javad Ghorbani on 14/11/2024.
//
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <cerrno>    // for errno
#include <cstring>   // for strerror
#include <unistd.h>  // for access()
#include "SANISANDModelUMAT.hpp"
#include "UMATBase.hpp"
// Validate analysis parameters
void validateAnalysisParameters(int analysisType, int numSteps, double epsilonYY) {
    if (numSteps <= 0) {
        throw std::runtime_error("Invalid number of steps (" + std::to_string(numSteps) + "). NumSteps must be > 0.");
    }
    // For example, only allow AnalysisType of 2, 3, or 4
    if (analysisType != 2 && analysisType != 3 && analysisType != 4) {
        throw std::runtime_error("Invalid AnalysisType (" + std::to_string(analysisType) + "). Allowed values are 2, 3, or 4.");
    }
    if (epsilonYY == 0.0) {
        throw std::runtime_error("Invalid EpsilonYY (" + std::to_string(epsilonYY) + "). EpsilonYY must be nonzero.");
    }
}
// Function to load environment variables from .env file
std::unordered_map<std::string, std::string> loadEnvFile() {
    std::unordered_map<std::string, std::string> env;
    std::ifstream envFile(".env");
    if (!envFile) {
        std::cerr << "Warning: .env file not found" << std::endl;
        return env;
    }

    std::string line;
    while (std::getline(envFile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        // Find the equals sign
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos + 1);
            env[key] = value;
        }
    }
    return env;
}

int main() {
    try {
        // Load environment variables
        auto env = loadEnvFile();
        
        // Get file paths from environment or use defaults
        std::string inputFilePath = env.count("INPUT_FILE_PATH") ? 
            env["INPUT_FILE_PATH"] : "input_data.txt";
        std::string outputFilePath = env.count("OUTPUT_FILE_PATH") ? 
            env["OUTPUT_FILE_PATH"] : "q_vs_epsilon2.csv";
        
        std::cout << "Using input file: " << inputFilePath << std::endl;
        std::cout << "Using output file: " << outputFilePath << std::endl;

        // Try to open the file
        std::ifstream inputFile(inputFilePath);
        if (!inputFile.is_open()) {
            std::cerr << "Error opening file: " << std::strerror(errno) << std::endl;
            throw std::runtime_error("Could not open input file: " + inputFilePath);
        }

        std::cout << "Successfully opened input file" << std::endl;

        // Store properties, state variables, and other inputs
        std::unordered_map<std::string, double> fileProperties;
        UMATBase::InputData inputData;  // For state variables and custom variables
        int Analysistype = 1;          // Default analysis type
        int numSteps = 1;              // Default number of steps
        double EpsilonYY = -1e-4;      // Default EpsilonYY
        double QUP = 50.0, QD = -50.0; // Default QUP and QD

        // Read from the input file
        std::string line;
        while (std::getline(inputFile, line)) {
            if (line.empty()) {
                continue;  // Skip empty lines
            }
            std::istringstream iss(line);
            std::string label;
            iss >> label;

            if (label == "Property:") {
                // Format: Property: KEY VALUE
                std::string key;
                double val;
                iss >> key >> val;
                fileProperties[key] = val;
            } else if (label == "StateVariable:") {
                // Format: StateVariable: KEY VALUE
                std::string key;
                double val;
                iss >> key >> val;
                if      (key == "VoidRatio")  { inputData.stateVariables[StateVariable::VoidRatio] = val; }
                else if (key == "StressXX")   { inputData.stateVariables[StateVariable::StressXX]  = val; }
                else if (key == "StressYY")   { inputData.stateVariables[StateVariable::StressYY]  = val; }
                else if (key == "StressZZ")   { inputData.stateVariables[StateVariable::StressZZ]  = val; }
                else if (key == "StressZY")   { inputData.stateVariables[StateVariable::StressZY]  = val; }
                else if (key == "StressZX")   { inputData.stateVariables[StateVariable::StressZX]  = val; }
                else if (key == "StressXY")   { inputData.stateVariables[StateVariable::StressXY]  = val; }
            } else if (label == "CustomVariable:") {
                // Format: CustomVariable: KEY VALUE
                std::string key;
                double val;
                iss >> key >> val;
                inputData.customStateVariables[key] = val;
            } else if (label == "AnalysisType:") {
                // Format: AnalysisType: <int>
                iss >> Analysistype;
            } else if (label == "NumSteps:") {
                // Format: NumSteps: <int>
                iss >> numSteps;
            } else if (label == "EpsilonYY:") {
                // Format: EpsilonYY: <double>
                iss >> EpsilonYY;
            } else if (label == "QUP:") {
                iss >> QUP;
            } else if (label == "QD:") {
                iss >> QD;
            }
            // Validate QUP and QD values
            if (QUP == 0.0 || QD == 0.0 || QUP <= QD) {
                std::cerr << "Error: Invalid QUP (" << QUP << ") and/or QD (" << QD 
                        << "). QUP must be nonzero and greater than QD." << std::endl;
                // Throw an exception to stop execution
                throw std::runtime_error("Invalid QUP and QD values.");
            }
        }
        inputFile.close();
        // Validate analysis parameters
        validateAnalysisParameters(Analysistype, numSteps, EpsilonYY);

        // --------------------------------------------------------------------
        // 2. Populate `properties` map
        // --------------------------------------------------------------------
        std::unordered_map<std::string, double> properties = fileProperties;

        // --------------------------------------------------------------------
        // 3. Convert `properties` map to array for UMAT initialization
        // --------------------------------------------------------------------
        static const char* variableNames[] = {
            "G0", "K0", "Mc", "Lambda", "N_c", "alpha_c", "n_b", "ch", "n_d",
            "h0", "A0", "Me", "cz", "zmax", "m_iso", "Patm", "P_min", "STOL", "FTOL", "LTOL"
        };
        const int numProps = sizeof(variableNames) / sizeof(variableNames[0]);
        std::vector<double> variableValues(numProps, 0.0);

        for (int i = 0; i < numProps; ++i) {
            if (properties.find(variableNames[i]) != properties.end()) {
                variableValues[i] = properties[variableNames[i]];
            }
        }

        // --------------------------------------------------------------------
        // 4. Initialize UMAT properties
        // --------------------------------------------------------------------
        initializeUMATProperties(variableNames, variableValues.data(), numProps);

        // --------------------------------------------------------------------
        // 5. Prepare default input data
        // --------------------------------------------------------------------
        inputData.strainIncrement[0] = -0.5 * EpsilonYY;      // Example strain increment
        inputData.strainIncrement[1] = EpsilonYY; // Read from file
        inputData.strainIncrement[2] = -0.5 * EpsilonYY;
        inputData.strainIncrement[3] = 0.0;
        inputData.strainIncrement[4] = 0.0;
        inputData.strainIncrement[5] = 0.0;

        // --------------------------------------------------------------------
        // 6. Create SANISANDModel and run loop
        // --------------------------------------------------------------------
        SANISANDModel model(properties);
        UMATBase::OutputData outputData;

        std::vector<double> q_values;
        std::vector<double> epsilon2_values;
        std::vector<double> sigm_values;
        std::vector<double> specificVoid_values;
        double specificVoid = 0.0, dLambda = 0.0, DPara = 0.0;
        Eigen::VectorXd NMAT;
        double DLAM = 0.0, DDD = 0.0;
        NMAT.setZero(6);
        double stressXX = 0.0, stressYY = 0.0, stressZZ = 0.0;
        double SIGM = 0.0;
        double G = 0.0, K = 0.0;
        double N11 = 0.0, N33 = 0.0;
        // Loop over the simulation steps
        for (int step = 0; step < numSteps; ++step) {
            // Example: Check analysis type & do some adjustments
            if (Analysistype == 2) {
                // Get stress components
                stressXX = inputData.stateVariables[StateVariable::StressXX];
                stressYY = inputData.stateVariables[StateVariable::StressYY];
                stressZZ = inputData.stateVariables[StateVariable::StressZZ];

                // Ensure consistent sign management
                SIGM = - (stressXX + stressYY + stressZZ) / 3.0; // Mean stress
                specificVoid = inputData.stateVariables[StateVariable::VoidRatio] + 1.0;
                // Compute shear modulus G
                G = std::pow(std::max(properties["P_min"], SIGM) / properties["Patm"], 0.5) *  properties["G0"] *  properties["Patm"] * std::pow((2.97 - (-1.0 + specificVoid)), 2.0)/ (1.0 + (specificVoid - 1.0));

                // Update K (bulk modulus) formula
                K = std::pow(std::max(properties["P_min"], SIGM) / properties["Patm"], 2.0 / 3.0) * properties["K0"] * properties["Patm"] * specificVoid/ (specificVoid - 1.0);

                
                N11 = NMAT[1];  // Example, should be updated based on the model
                N33 = NMAT[2];  // Example, should be updated based on the model
                double DEPSYY = -inputData.strainIncrement[1];
                
                double DEPSXX = ((3.0 * K - 2.0 * G) * DEPSYY - 3.0 * std::max(DLAM, 0.0) * K * DDD +
                                std::max(DLAM, 0.0) * 2.0 * G * (N11 - N33)) /
                               (-(6.0 * K + 2.0 * G));

                // Assign DEPSZZ to DEPSXX as per Fortran logic
                double DEPSZZ = DEPSXX;
                inputData.strainIncrement[0] =  - DEPSXX;
                inputData.strainIncrement[2] =  - DEPSZZ;
            }else if (Analysistype == 3) {

                specificVoid = inputData.stateVariables[StateVariable::VoidRatio] + 1.0;
                // Compute strain increments
                double DEPSYY = -inputData.strainIncrement[1];

                // Compute DEPSXX using the given equation
                double DEPSXX = 0.5 * (std::max(DLAM, 0.0) * DDD - DEPSYY);
                // Assign DEPSZZ to DEPSXX as per the provided Fortran logic
                double DEPSZZ = DEPSXX;

                // Update strain increments in the input data
                inputData.strainIncrement[0] = -DEPSXX;
                inputData.strainIncrement[2] = -DEPSZZ;
            } else if (Analysistype == 4) {
                // Get stress components
                double stressXX = - inputData.stateVariables[StateVariable::StressXX];
                double stressYY = - inputData.stateVariables[StateVariable::StressYY];

                // Calculate deviatoric stress difference (QSIG)
                double QSIG = stressYY - stressXX;

                // Define MID and DIS (parameters based on QUP and QD)
                double MID = 0.5 * (QUP + QD);
                double DIS = 0.5 * (QUP - QD);

                // Check condition for QSIG close to MID
                if (not (std::abs(QSIG - MID) < DIS)) {
                    // Update strain increments
                    double DEPSYY = -inputData.strainIncrement[1];
                    double DEPSXX = -0.5 * DEPSYY;
                    double DEPSZZ = DEPSXX;

                    // Assign the computed strain increments back to inputData
                    inputData.strainIncrement[0] = DEPSXX;
                    inputData.strainIncrement[1] = DEPSYY;
                    inputData.strainIncrement[2] = DEPSZZ;
                }
            }

            // Perform stress integration
            calculateStressIncrement(inputData, outputData);

            // Update stress state
            for (int i = 0; i < 6; ++i) {
                inputData.stateVariables[StateVariable::StressXX + i] += outputData.stressIncrement[i];
            }

            // Update void ratio
            specificVoid = inputData.stateVariables[StateVariable::VoidRatio] + 1.0;
            specificVoid *= (1.0 + inputData.strainIncrement[0]
                                  + inputData.strainIncrement[1]
                                  + inputData.strainIncrement[2]);
            inputData.stateVariables[StateVariable::VoidRatio] = specificVoid - 1.0;
            // Store specificVoid for the current step
            specificVoid_values.push_back(specificVoid);

            // Evaluate "q" :
            double stressXX = inputData.stateVariables[StateVariable::StressXX];
            double stressYY = inputData.stateVariables[StateVariable::StressYY];
            double stressZZ = inputData.stateVariables[StateVariable::StressZZ];
            double sigma_m  = -(stressXX + stressYY + stressZZ) / 3.0;
            double q        = -(stressYY - stressXX); // or any other definition

            // Accumulate results
            q_values.push_back(q);
            epsilon2_values.push_back(-inputData.strainIncrement[1] * (step + 1));
            sigm_values.push_back(sigma_m);

            // Update custom state variables from output
            for (const auto& var : {"AlphaXX", "AlphaYY", "AlphaZZ",
                                    "AlphaZY", "AlphaZX", "AlphaXY"}) {
                inputData.customStateVariables[var] = outputData.updatedCustomStateVariables[var];
            }
            for (const auto& var : {"AlphaInitialXX", "AlphaInitialYY",
                                    "AlphaInitialZZ", "AlphaInitialZY",
                                    "AlphaInitialZX", "AlphaInitialXY"}) {
                inputData.customStateVariables[var] = outputData.updatedCustomStateVariables[var];
            }
            for (const auto& var : {"ZXX", "ZYY", "ZZZ", "ZZY", "ZZX", "ZXY"}) {
                inputData.customStateVariables[var] = outputData.updatedCustomStateVariables[var];
            }
                Eigen::VectorXd Sigma0(6);
                Sigma0 << inputData.stateVariables[StateVariable::StressXX],
                          inputData.stateVariables[StateVariable::StressYY],
                          inputData.stateVariables[StateVariable::StressZZ],
                          inputData.stateVariables[StateVariable::StressZY],
                          inputData.stateVariables[StateVariable::StressZX],
                          inputData.stateVariables[StateVariable::StressXY];
                // Update kinematic hardening vector from custom state variables
                Eigen::VectorXd kinematicHardeningVector(6);
                kinematicHardeningVector << inputData.customStateVariables["AlphaXX"],
                                            inputData.customStateVariables["AlphaYY"],
                                            inputData.customStateVariables["AlphaZZ"],
                                            inputData.customStateVariables["AlphaZY"],
                                            inputData.customStateVariables["AlphaZX"],
                                            inputData.customStateVariables["AlphaXY"];

                // Update fabric vector Zvec from custom state variables
                Eigen::VectorXd Zvec(6);
                Zvec << inputData.customStateVariables["ZXX"],
                        inputData.customStateVariables["ZYY"],
                        inputData.customStateVariables["ZZZ"],
                        inputData.customStateVariables["ZZY"],
                        inputData.customStateVariables["ZZX"],
                        inputData.customStateVariables["ZXY"];

                // Update initial backstress vector alphaInitial from custom state variables
                Eigen::VectorXd alphaInitial(6);
                alphaInitial << inputData.customStateVariables["AlphaInitialXX"],
                                inputData.customStateVariables["AlphaInitialYY"],
                                inputData.customStateVariables["AlphaInitialZZ"],
                                inputData.customStateVariables["AlphaInitialZY"],
                                inputData.customStateVariables["AlphaInitialZX"],
                                inputData.customStateVariables["AlphaInitialXY"];
                specificVoid = inputData.stateVariables[StateVariable::VoidRatio] + 1.0;

                Sigma0 *=-1.0;
                // --- Compute the elastic stiffness matrix ---
                Eigen::MatrixXd D = model.computeModuli(Sigma0, kinematicHardeningVector, specificVoid);
                Sigma0 *=-1.0;
                Eigen::Map<Eigen::VectorXd> strainIncrement(inputData.strainIncrement, 6);
                Eigen::VectorXd DStress_ela = D * strainIncrement;

                // --- Prepare deltaStrain (often similar to strainIncrement with sign change) ---
                Eigen::VectorXd deltaStrain(6);
                for (int i = 0; i < 6; ++i) {
                    deltaStrain(i) = -inputData.strainIncrement[i];
                }
                // --- Define the time increment DT -- this is dummy variable  and must be set to 1.0 here ---
                double DT = 1.0;

                // --- Call ComputeDlambda to obtain dLambda and DPara (here used as DDD) ---
                auto [computed_dLambda, computed_DPara, computed_NMAT] = model.ComputeDlambda(
                    Sigma0, DStress_ela, DT, kinematicHardeningVector, Zvec, alphaInitial, specificVoid, deltaStrain);

                DLAM = computed_dLambda;  // Updated consistency parameter
                DDD  = computed_DPara;     // Updated hardening modulus (or DPara)
                NMAT = computed_NMAT;
        }

        // Write results to CSV
        std::ofstream file(outputFilePath);
        file << "Step,Epsilon_2,p,q,specificVoid\n";
        for (size_t i = 0; i < q_values.size(); ++i) {
            file << (i + 1) << ","
                << epsilon2_values[i] << ","
                << sigm_values[i]     << ","
                << q_values[i]        << ","
                << specificVoid_values[i] << "\n";
        }
        file.close();

        std::cout << "Simulation complete. Results written to " << outputFilePath << "." << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
