//  SANISANDModelUMAT.hpp
//  SANISAND
//
//  Created by Javad Ghorbani on 22/12/2024.
//

#ifndef SANISANDMODELUMAT_HPP
#define SANISANDMODELUMAT_HPP
#include "UMATBase.hpp"
#include <unordered_map>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

// Struct to hold stress invariants and related quantities
struct StressInvariants {
    double sigma_m;          // Mean stress
    double J2;               // Second invariant of deviatoric stress
    double J3;               // Third invariant of deviatoric stress
    double SBAR;             // Deviatoric stress magnitude
    double S3TA;             // Lode angle parameter
    double SIGQ;             // Deviatoric stress norm
    double COSTTETA;         // Cosine of the Lode angle
    Eigen::VectorXd devStress; // Deviatoric stress vector
    Eigen::VectorXd RMAT;    // Normalized deviatoric stress tensor
    Eigen::VectorXd NMAT;    // Normal to yield surface
    double J3Q;              // Scaled third invariant
    double THETA;
};

// Main class for SANISAND model
class SANISANDModel : public UMATBase {
public:
    // Constructor
    explicit SANISANDModel(const std::unordered_map<std::string, double>& properties);

    // Public member functions
    StressInvariants computeStressInvariants(const Eigen::VectorXd& stress, const Eigen::VectorXd& backstress) const;
    // Compute stress increment
    void calculateStressIncrement(const InputData& inputData, OutputData& outputData) override;

    // Compute stress-strain matrix
    void computeStressStrainMatrix(const InputData& inputData, OutputData& outputData) override;
    void setCustomVariable(UMATBase::InputData& inputData);
    Eigen::MatrixXd computeModuli(const Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardeningVector, const double specificVoid) const;
    std::tuple<double, double, Eigen::VectorXd> ComputeDlambda(
        const Eigen::VectorXd& stress,
        const Eigen::VectorXd& DStress_ela,
        const double DT,
        const Eigen::VectorXd& kinematicHardeningVector,
        const Eigen::VectorXd& Zvec,
        const Eigen::VectorXd& alpha_init,
        const double specificVoid,
        const Eigen::VectorXd& deltaStrain) const;
private:
    // Material properties
    double G0;       // Shear modulus constant
    double K0;       // Bulk modulus constant
    double Mc;       // Critical state slope in compression
    double lambda;   // Compression index
    double N_c;      // Reference void ratio
    double alpha_c;  // Dilatancy parameter
    double n_b;      // Bounding surface parameter
    double ch;       // Hardening parameter
    double n_d;      // Dilatancy
    double h0;       // Hardening modulus
    double A0;       // Plastic potential parameter
    double Me;       // Critical state slope in extension
    double cz;       // Critical state parameter
    double zmax;     // Maximum fabric tensor
    double m_iso;    // Isotropic hardening parameter
    double patm;     //Atmospheric pressure
    double P_min;
    double STOL;     // Convergence tolerance
    double FTOL;     // Yield tolerance
    double LTOL;     // Load/unload detection tolerance

    Eigen::VectorXd GradientOfYieldToSigma(const Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardeningVector) const;

    Eigen::VectorXd GradientOfPlasticToSigma(
        const Eigen::VectorXd& stress,
        const Eigen::VectorXd& kinematicHardeningVector,
        const Eigen::VectorXd& fabricVector,
        const double specificVoid) const;

    std::pair<Eigen::VectorXd, double> PlasticModulus(
        const Eigen::VectorXd& stress,
        const Eigen::VectorXd& kinematicHardeningVector,
        const Eigen::VectorXd& Zvec,
        const Eigen::VectorXd& alpha_init,
                                                      const double specificVoid) const;
    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, double> ComputeDlambdaAndUpdateStressInc(const Eigen::VectorXd& stress, const Eigen::VectorXd& DStress_ela, const double DT, const Eigen::VectorXd& kinematicHardeningVector,const Eigen::VectorXd& Zvec, const Eigen::VectorXd& alpha_init, const double specificVoid) const;
//    Eigen::VectorXd computeElasticStressIncrement(const Eigen::VectorXd& deltaStrain, Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardening, double specificVoid) const;
    void updateAInit(Eigen::VectorXd& AINIT, const Eigen::VectorXd& SIGMAT, const Eigen::VectorXd& kinHardening);
    double computeCSLVoidRatio(double SIGM) const;
    double calculateRelativeErrors(
        const Eigen::VectorXd& DSIG1, const Eigen::VectorXd& DSIG2,
        const Eigen::VectorXd& Sigma1, const Eigen::VectorXd& Sigma2,
        const Eigen::VectorXd& dKinematicHardening1, const Eigen::VectorXd& dKinematicHardening2,
        const Eigen::VectorXd& kinHardening2,
        const Eigen::VectorXd& dz1, const Eigen::VectorXd& dz2,
                                                  const Eigen::VectorXd& Zvec2);
    double ComputeDPara(const Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardeningVector, const Eigen::VectorXd& Zvec, const double specificVoid) const;
};

// C-style interface functions
extern "C" int getNumRequiredVariables();
extern "C" const char* getRequiredVariableName(int index);
extern "C" void initializeUMATProperties(const char** variableNames, const double* variableValues, int variableCount);
extern "C" void calculateStressIncrement(const SANISANDModel::InputData& inputData, SANISANDModel::OutputData& outputData);
extern "C" void computeStressStrainMatrix(const SANISANDModel::InputData& inputData, SANISANDModel::OutputData& outputData);
extern "C" const double* getStressIncrement(const SANISANDModel::OutputData& outputData);
extern "C" const double* getStressStrainMatrixRow(const SANISANDModel::OutputData& outputData, int row);
extern "C" void initializeCustomVariable(UMATBase::InputData& inputData);
#endif // SANISANDMODELUMAT_HPP
