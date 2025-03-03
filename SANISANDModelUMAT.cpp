//
//  SANISANDModelUMAT.cpp
//  SANISAND
//
//  Created by Javad Ghorbani on 22/12/2024.
//

// SANISANDModelUMAT.cpp
#include "SANISANDModelUMAT.hpp"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iomanip> // for std::setprecision
#include <unordered_map>

// Global instance of SANISANDModel
static SANISANDModel* umatInstance = nullptr;

// Required variable names for SANISANDModel
static const char* requiredVariables[] = {
    "G0", "K0", "Mc", "Lambda", "N_c", "alpha_c", "n_b", "ch", "n_d",
    "h0", "A0", "Me", "cz", "zmax", "m_iso", "Patm", "P_min", "STOL", "FTOL", "LTOL"
};

extern "C" int getNumRequiredVariables() {
    return sizeof(requiredVariables) / sizeof(requiredVariables[0]);
}

extern "C" const char* getRequiredVariableName(int index) {
    if (index < 0 || index >= getNumRequiredVariables()) {
        return nullptr; // Invalid index
    }
    return requiredVariables[index];
}

// C-style function for initializing the SANISANDModel with variable names and values
extern "C" void initializeUMATProperties(const char** variableNames, const double* variableValues, int variableCount) {
    std::cout << "Initializing properties for SANISANDModel..." << std::endl;

    // Ensure the variable count matches the number of required variables
    if (variableCount != getNumRequiredVariables()) {
        std::cerr << "Mismatch in number of required variables for SANISANDModel." << std::endl;
        throw std::invalid_argument("Variable count mismatch.");
    }

    // Create properties map
    std::unordered_map<std::string, double> properties;
    for (int i = 0; i < variableCount; ++i) {
        properties[variableNames[i]] = variableValues[i];
    }

    // Clean up an existing instance if necessary
    if (umatInstance != nullptr) {
        delete umatInstance;
        umatInstance = nullptr;
    }

    // Initialize modelInstance with SANISANDModel
    umatInstance = new SANISANDModel(properties);
}

// Constructor
SANISANDModel::SANISANDModel(const std::unordered_map<std::string, double>& properties) {
    auto findOrThrow = [&properties](const std::string& key) -> double {
        auto it = properties.find(key);
        if (it == properties.end()) {
            throw std::invalid_argument("Missing required material property: " + key);
        }
        return it->second;
    };

    // Retrieve required material properties
    G0 = findOrThrow("G0");               // Shear modulus constant
    K0 = findOrThrow("K0");               // Bulk modulus constant
    Mc = findOrThrow("Mc");               // Critical state slope in compression
    lambda = findOrThrow("Lambda");       // Compression index
    N_c = findOrThrow("N_c");             // Reference void ratio
    alpha_c = findOrThrow("alpha_c");     // Dilatancy parameter
    n_b = findOrThrow("n_b");             // Bounding surface parameter
    ch = findOrThrow("ch");               // Hardening parameter
    n_d = findOrThrow("n_d");             // Dilatancy
    h0 = findOrThrow("h0");               // Hardening modulus
    A0 = findOrThrow("A0");               // Plastic potential parameter
    Me = findOrThrow("Me");               // Critical state slope in extension
    cz = findOrThrow("cz");               // Critical state parameter
    m_iso = findOrThrow("m_iso");         // Isotropic hardening parameter
    zmax = findOrThrow("zmax");           // Maximum fabric tensor
    patm = findOrThrow("Patm");           // Atmospheric reference pressure
    P_min = findOrThrow("P_min");          // minimum pressure for
    STOL = findOrThrow("STOL");           // Convergence tolerance
    FTOL = findOrThrow("FTOL");           // Yield tolerance
    LTOL = findOrThrow("LTOL");           // Load/unload detection tolerance

    // Validate parameters
    if (lambda <= 0.0) {
        throw std::invalid_argument("Lambda must be greater than 0.");
    }
    if (n_b <= 0.0 || n_d <= 0.0) {
        throw std::invalid_argument("n_b and n_d must be greater than 0.");
    }
    if (ch <= 0.0 || h0 <= 0.0) {
        throw std::invalid_argument("ch and h0 must be greater than 0.");
    }
    if (A0 <= 0.0) {
        throw std::invalid_argument("A0 must be greater than 0.");
    }
    if (cz <= 0.0) {
        throw std::invalid_argument("cz must be greater than 0.");
    }
    if (zmax <= 0.0) {
        throw std::invalid_argument("zmax must be greater than 0.");
    }
}


// Compute stress invariants
StressInvariants SANISANDModel::computeStressInvariants(const Eigen::VectorXd& stress, const Eigen::VectorXd& backstress) const {
    StressInvariants invariants;
    const double TINY = 1.0e-12;
    // Mean stress
    invariants.sigma_m = (stress[0] + stress[1] + stress[2]) / 3.0;

    // Deviatoric stress components
    invariants.devStress = Eigen::VectorXd(6);
    invariants.devStress[0] = stress[0] - invariants.sigma_m;
    invariants.devStress[1] = stress[1] - invariants.sigma_m;
    invariants.devStress[2] = stress[2] - invariants.sigma_m;
    invariants.devStress[3] = stress[3]; // Shear components remain the same
    invariants.devStress[4] = stress[4];
    invariants.devStress[5] = stress[5];

    // Second deviatoric invariant J2
    invariants.J2 = 0.5 * (invariants.devStress[0] * invariants.devStress[0] +
                           invariants.devStress[1] * invariants.devStress[1] +
                           invariants.devStress[2] * invariants.devStress[2] +
                           2.0 * (invariants.devStress[3] * invariants.devStress[3] +
                                  invariants.devStress[4] * invariants.devStress[4] +
                                  invariants.devStress[5] * invariants.devStress[5]));

    Eigen::Matrix3d S;
    S << invariants.devStress[0], invariants.devStress[5], invariants.devStress[4],
         invariants.devStress[5], invariants.devStress[1], invariants.devStress[3],
         invariants.devStress[4], invariants.devStress[3], invariants.devStress[2];

    // J3 = det(S)
    invariants.J3 = S.determinant();

    // Guard against numerical round-off
    if (std::fabs(invariants.J3) < TINY) {
        invariants.J3 = 0.0;
    }

    // Deviatoric stress magnitude SBAR and Lode angle parameters
    if (invariants.J2 > TINY) {
        invariants.SBAR = std::sqrt(invariants.J2);
        invariants.S3TA = -4.5 * invariants.J3 / (std::sqrt(3.0) * invariants.SBAR * invariants.J2);
        invariants.S3TA = std::clamp(invariants.S3TA, -1.0, 1.0);
        // Check and adjust S3TA if it is out of bounds
        if (invariants.S3TA < -1.0) {
            invariants.S3TA = -1.0;
        } else if (invariants.S3TA > 1.0) {
            invariants.S3TA = 1.0;
        }
        invariants.THETA = (1.0 / 3.0) * std::asin(invariants.S3TA);
        invariants.SIGQ = std::sqrt(3.0 * invariants.J2);

    } else {
        // Special case of zero deviatoric stress
        invariants.J2 = TINY * TINY;
        invariants.SBAR = TINY;
        invariants.THETA = 0.0;
        invariants.S3TA = 0.0;
        invariants.SIGQ = TINY;
    }
    // Compute J3Q

    invariants.J3Q = -(27.0 * invariants.J3) / (2.0 * std::pow(invariants.SIGQ, 3));
    // Normalized deviatoric stress tensor RMAT(I,1) = dev_stress(I,1)/SIGM
    invariants.RMAT = Eigen::VectorXd(6);
    if (std::fabs(invariants.sigma_m) > TINY) {
        for (int i = 0; i < 6; ++i) {
            invariants.RMAT[i] = invariants.devStress[i] / invariants.sigma_m;
        }
    } else {
        invariants.RMAT.setZero(); // Avoid division by zero
    }
    // Compute normal to yield surface (NMAT = RMAT - Malpha)
    Eigen::VectorXd NMAT(6);
    NMAT = invariants.RMAT - backstress;

    // Compute norm of NMAT
    double ABSNMAT = NMAT.norm();

    // Handle cases where J2 and ABSNMAT are very small
    if (invariants.J2 <= TINY && ABSNMAT < TINY) {
        NMAT << 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0; // Set the first three components to 1 and the rest to 0
    } else {
        NMAT /= ABSNMAT; // Normalize NMAT
    }

    // Assign NMAT to invariants
    invariants.NMAT = NMAT;
    // Compute COSTTETA using NMAT
    double COSTTETA = std::sqrt(6.0) * (std::pow(NMAT[0], 3.0) +
                                        std::pow(NMAT[1], 3.0) +
                                        std::pow(NMAT[2], 3.0));
//                                        std::pow(NMAT[3], 3.0) +
//                                        std::pow(NMAT[4], 3.0) +
//                                        std::pow(NMAT[5], 3.0));

    // Check for isotropic case and set COSTTETA to 1
    if (NMAT[0] == NMAT[1] && NMAT[1] == NMAT[2] &&
        NMAT[3] == 0.0 && NMAT[4] == 0.0 && NMAT[5] == 0.0) {
        COSTTETA = 1.0; // Isotropic case: no shear, normal components equal
    } else {
        // Ensure COSTTETA is within bounds [-1.0, 1.0]
        if (std::abs(COSTTETA) > 1.0) {
            COSTTETA = COSTTETA / std::abs(COSTTETA);
        }
    }

    // Assign COSTTETA to invariants (add COSTTETA to the StressInvariants struct)
    invariants.COSTTETA = COSTTETA;
    
    return invariants;
}
Eigen::VectorXd SANISANDModel::GradientOfYieldToSigma(const Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardeningVector) const {
    // Compute stress invariants
    StressInvariants invariants = computeStressInvariants(stress, kinematicHardeningVector);
    Eigen::VectorXd devStress = invariants.devStress; // Deviatoric stress components
    // Include RMAT and NMAT from invariants
    Eigen::VectorXd RMAT = invariants.RMAT;
    Eigen::VectorXd NMAT = invariants.NMAT;
    // Compute TNMAT
    Eigen::RowVectorXd TNMAT = NMAT.transpose(); // NMAT is now a row vector
    double N2 = TNMAT.dot(kinematicHardeningVector);
    // Initialize GradientYield
    Eigen::VectorXd GradientYield(6);
    
    // Compute GY terms
    for (int i = 0; i < 6; ++i) {
        if (i < 3) { // Normal stress terms (1st, 2nd, 3rd components)
            GradientYield[i] = NMAT[i] - (1.0 / 3.0) * N2;
        } else { // Shear stress terms (4th, 5th, 6th components)
            GradientYield[i] = NMAT[i];
        }
    }
    return GradientYield;
}

Eigen::VectorXd SANISANDModel::GradientOfPlasticToSigma(const Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardeningVector, const Eigen::VectorXd& Zvec, const double specificVoid) const {
    
    // Compute stress invariants
    StressInvariants invariants = computeStressInvariants(stress, kinematicHardeningVector);
    Eigen::VectorXd devStress = invariants.devStress; // Deviatoric stress components
    double SIGM = invariants.sigma_m;                 // Mean stress
    double COSTTETA = invariants.COSTTETA;
    
    // Include RMAT and NMAT from invariants
    Eigen::VectorXd RMAT = invariants.RMAT;
    Eigen::VectorXd NMAT = invariants.NMAT;
    
    
    if (std::fabs(Mc) < 1.0e-12) {
        throw std::runtime_error(" (M_c) is too small, causing division by zero in CPARA.");
    }
    double CPARA = Me / Mc;
    
    // Compute GTETC (as one-liner)
    double GTETC = std::pow((2.0 * std::pow(CPARA, 4.0)) / ((1.0 + std::pow(CPARA, 4.0)) - (1.0 - std::pow(CPARA, 4.0)) * COSTTETA), 1.0 / 4.0);
    
    // Compute critical state void ratio and related terms
    double ECPARA = computeCSLVoidRatio(SIGM); // Critical state void ratio
    double EPARA = specificVoid - 1.0;        // Specific void
    double SI = -ECPARA + EPARA;              // Difference between current and critical void ratios
    
    // Compute ATETD
    double ATETD = GTETC * Mc * std::exp(n_d * SI);
    
    Eigen::VectorXd NMAT2 = NMAT;                   // Normalize NMAT2
    
    // Compute ATETDMAT
    Eigen::VectorXd ATETDMAT = (std::sqrt(2.0 / 3.0)) * ATETD * NMAT2;
    // Compute Z * NMAT (dot product with transpose)
    double ZN = Zvec.transpose().dot(NMAT); // Dot product
    
    double MDAA = A0 * (1.0 + std::max(std::sqrt(3.0 / 2.0) * ZN, 0.0));
    
    // Compute TNMAT2
    Eigen::RowVectorXd TNMAT2 = NMAT2.transpose(); // NMAT is now a row vector
    
    // Compute DPARA (Dot Product)
    double DPARA = MDAA * TNMAT2.dot(ATETDMAT - kinematicHardeningVector);
    DPARA *= std::sqrt(3.0 / 2.0); // Modify by multiplying with sqrt(3.0/2.0)
    
    Eigen::VectorXd GP(6);
    // Compute GP terms
    for (int i = 0; i < 6; ++i) {
        if (i < 3) { // Normal stress terms (1st, 2nd, 3rd components)
            GP[i] = NMAT[i] + (1.0 / 3.0) * DPARA;
        } else { // Shear stress terms (4th, 5th, 6th components)
            GP[i] = NMAT[i];
        }
    }
    return GP;
}

std::pair<Eigen::VectorXd, double> SANISANDModel::PlasticModulus(
    const Eigen::VectorXd& stress,
    const Eigen::VectorXd& kinematicHardeningVector,
    const Eigen::VectorXd& Zvec,
    const Eigen::VectorXd& alpha_init,
    const double specificVoid) const {
    
    // Compute stress invariants
    StressInvariants invariants = computeStressInvariants(stress, kinematicHardeningVector);
    Eigen::VectorXd devStress = invariants.devStress; // Deviatoric stress components
    double SIGM = invariants.sigma_m;                 // Mean stress
    double COSTTETA = invariants.COSTTETA;

    // Include RMAT and NMAT from invariants
    Eigen::VectorXd RMAT = invariants.RMAT;
    Eigen::VectorXd NMAT = invariants.NMAT;

    double CPARA = Me / Mc;
    // Compute GTETC
    double GTETC = std::pow((2.0 * std::pow(CPARA, 4.0)) / ((1.0 + std::pow(CPARA, 4.0)) - (1.0 - std::pow(CPARA, 4.0)) * COSTTETA), 1.0 / 4.0);

    // Compute critical state void ratio and related terms
    double ECPARA = computeCSLVoidRatio(SIGM); // Critical state void ratio
    double EPARA = specificVoid - 1.0;        // Specific void
    double SI = -ECPARA + EPARA;              // Difference between current and critical void ratios

    // Compute ATETD
    double ATETD = GTETC * Mc * std::exp(n_d * SI);

    Eigen::VectorXd NMAT2 = NMAT;                   // Normalize NMAT2

    // Compute ATETDMAT
    Eigen::VectorXd ATETDMAT = (std::sqrt(2.0 / 3.0)) * ATETD * NMAT2;
    double ZN = Zvec.transpose().dot(NMAT);         // Compute Z * NMAT (dot product)
    double MDAA = A0 * (1.0 + std::max(std::sqrt(3.0 / 2.0) * ZN, 0.0));

    // Compute DPARA
    double DPARA = MDAA * (NMAT2.transpose().dot(ATETDMAT - kinematicHardeningVector));
    DPARA *= std::sqrt(3.0 / 2.0);

    // Other computations
    double DGDP = DPARA;
    double DGDQ = 1.0;

    double ATETB = GTETC * Mc * std::max(std::exp(-n_b * SI), 0.0);
    Eigen::VectorXd ATETBMAT = std::sqrt(2.0 / 3.0) * ATETB * NMAT;
    Eigen::VectorXd BPARAMAT = ATETBMAT - kinematicHardeningVector;
    double BREF = G0 * h0 * std::max(1.0 - ch * EPARA, 1.0e-9) * std::pow(SIGM / patm, -0.5);
    Eigen::VectorXd AMINIT = kinematicHardeningVector - alpha_init;
    double NAMANIT = AMINIT.transpose().dot(NMAT);
    double HPARA = BREF / std::max(1.0e-15, NAMANIT);
    Eigen::VectorXd ASTAR = HPARA * BPARAMAT;

    // Compute Kp and B_kin
    double NA = ASTAR.transpose().dot(NMAT);
    double Kp = (2.0 / 3.0) * SIGM * NA * DGDQ;
    Eigen::VectorXd B_kin = (2.0 / 3.0) * ASTAR;

    // Return both as a pair
    return std::make_pair(B_kin, Kp);
}

double SANISANDModel::computeCSLVoidRatio(double SIGM) const{
    // lnSigma: log of sigma
    // Nc: parameter N_c
    // lambda: parameter lambda
    // alpha_c: parameter alpha_c
    
    double lne = std::log(N_c) - lambda * std::log(SIGM + alpha_c); // lne computation
    return std::exp(lne); // e = exp(lne)
}

// Compute stress-strain matrix (elastic stiffness matrix)
Eigen::MatrixXd SANISANDModel::computeModuli(const Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardeningVector, const double specificVoid) const {
    Eigen::MatrixXd D(6, 6);
    D.setZero();
    StressInvariants invariants = computeStressInvariants(stress, kinematicHardeningVector);
    double SIGM = invariants.sigma_m;                 // Mean stress
    double MDA = 0.5;
    // Update G (shear modulus) formula
    double G = std::pow(std::max(P_min, SIGM) / patm, MDA) * G0 * patm * std::pow((2.97 - (-1.0 + specificVoid)), 2.0)/ (1.0 + (specificVoid - 1.0));

    // Update K (bulk modulus) formula
    double K = std::pow(std::max(P_min, SIGM) / patm, 2.0 / 3.0) * K0 * patm * specificVoid/ (specificVoid - 1.0);

    double lambda = K - (2.0 / 3.0) * G;
    double mu = G;
    // Print parameters
//    std::cout << "Input Parameters:" << std::endl;
//    std::cout << "Stress vector: " << stress.transpose() << std::endl;
//    std::cout << "Specific void: " << specificVoid << std::endl;
//    std::cout << "kappa: " << kappa << std::endl;
//
//    std::cout << "\nComputed Values:" << std::endl;
//    std::cout << "Mean stress (SIGM): " << SIGM << std::endl;
//    std::cout << "Bulk modulus (K): " << K << std::endl;
//    std::cout << "Shear modulus (G): " << G << std::endl;
//    std::cout << "Lambda: " << lambda << std::endl;
//    std::cout << "Mu: " << mu << std::endl;

    // Populate stiffness matrix
    D(0, 0) = lambda + 2.0 * mu;
    D(0, 1) = lambda;
    D(0, 2) = lambda;
    D(1, 0) = lambda;
    D(1, 1) = lambda + 2.0 * mu;
    D(1, 2) = lambda;
    D(2, 0) = lambda;
    D(2, 1) = lambda;
    D(2, 2) = lambda + 2.0 * mu;
    D(3, 3) = mu; // Shear modulus for zy
    D(4, 4) = mu; // Shear modulus for zx
    D(5, 5) = mu; // Shear modulus for xy

    return D;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, double> SANISANDModel::ComputeDlambdaAndUpdateStressInc(const Eigen::VectorXd& stress, const Eigen::VectorXd& DStress_ela, const double DT, const Eigen::VectorXd& kinematicHardeningVector,const Eigen::VectorXd& Zvec, const Eigen::VectorXd& alpha_init, const double specificVoid) const {
    // Obtain the gradients of the yield function and plastic potential
    Eigen::VectorXd df_dsigma = GradientOfYieldToSigma(stress, kinematicHardeningVector);
    Eigen::VectorXd dg_dsigma = GradientOfPlasticToSigma(stress, kinematicHardeningVector, Zvec, specificVoid);

    // Get the elastic matrix (De)
    Eigen::MatrixXd De = computeModuli(stress, kinematicHardeningVector, specificVoid);

    // Compute the numerator: df_dsigma.transpose() * De * deltaStrain
    double numerator = DT * df_dsigma.transpose() * DStress_ela;

    // Retrieve both B_iso and Kp from PlasticModulus
    auto plasticModulusResult = PlasticModulus(stress, kinematicHardeningVector, Zvec, alpha_init, specificVoid);
    Eigen::VectorXd B_kin = plasticModulusResult.first;
    double Kp = plasticModulusResult.second;
    
    // Compute the denominator: Kp + df_dsigma.transpose() * De * dg_dsigma
    double denominator = Kp + df_dsigma.transpose() * De * dg_dsigma;

    // Compute the consistency parameter (dLambda)
    double dLambda = numerator / denominator;

    // Compute norms for numerical checks
    double ANORM = df_dsigma.squaredNorm();
    double SNORM = DStress_ela.squaredNorm();
    double ATSIG = df_dsigma.transpose() * DStress_ela;
    double COSALP = ATSIG / std::sqrt(ANORM * SNORM);

    // Check if COSALP is below tolerance
    if (std::abs(COSALP) < -LTOL) {
        std::cerr << "COS ALPHA < ZERO DSTRS1: " << COSALP << ", dLambda: " << dLambda << std::endl;
        dLambda = 0.0;
    }

    // Check if ATDB is negative
    double ATDB = df_dsigma.transpose() * De * dg_dsigma;
    if (ATDB < 0) {
        std::cerr << "*** NEGATIVE ATDB ***: " << ATDB << std::endl;
        if (dLambda < 0) {
            throw std::runtime_error("Negative dLambda encountered during stress integration.");
        }
    }

    // Ensure dLambda is non-negative
    dLambda = std::max(dLambda, 0.0);

    // Compute DB = De * dg_dsigma
    Eigen::VectorXd DB = De * dg_dsigma;

    // Compute the stress increment (DeltaSigma)
    Eigen::VectorXd DeltaSigma = DStress_ela * DT - dLambda * DB;

    // Compute the isotropic hardening increment
    Eigen::VectorXd dkinematicHardening = dLambda * B_kin;
    
    StressInvariants invariants = computeStressInvariants(stress, kinematicHardeningVector);
    Eigen::VectorXd  NMAT = invariants.NMAT;
    // Initialize dz as a vector of size 6
    Eigen::VectorXd dz(6);
    // Define the identity vector for deviatoric strain
    Eigen::VectorXd I(6);
    I << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0; // 1 for normal components, 0 for shear
    double DEPSPV = dLambda * (dg_dsigma[0] + dg_dsigma[1] + dg_dsigma[2]);
    // Compute the deviatoric plastic strain tensor
    Eigen::VectorXd deviatoricPlasticStrain(6);
    deviatoricPlasticStrain = dLambda * dg_dsigma - (DEPSPV * I);

    // Compute equivalent deviatoric plastic strain
    double equivalentDeviatoricPlasticStrain = std::sqrt((2.0 / 3.0) * deviatoricPlasticStrain.dot(deviatoricPlasticStrain));

    // Coefficients
    double factor = -cz * std::max(-DEPSPV, 0.0);
    double scaledZmax = std::sqrt(2.0 / 3.0) * zmax;

    // Compute each component of dz based on the input vectors and parameters
    dz[0] = factor * (scaledZmax * NMAT[0] + Zvec[0]); // DZ11
    dz[1] = factor * (scaledZmax * NMAT[1] + Zvec[1]); // DZ22
    dz[2] = factor * (scaledZmax * NMAT[2] + Zvec[2]); // DZ33
    dz[3] = factor * (scaledZmax * NMAT[3] + Zvec[3]); // DZ12
    dz[4] = factor * (scaledZmax * NMAT[4] + Zvec[4]); // DZ23
    dz[5] = factor * (scaledZmax * NMAT[5] + Zvec[5]); // DZ31

    // Return both DeltaSigma and dIsotropicHardening
    return {DeltaSigma, dkinematicHardening, dz, equivalentDeviatoricPlasticStrain};
}


// Compute stress-strain matrix (wrapper function)
void SANISANDModel::computeStressStrainMatrix(const InputData& inputData, OutputData& outputData) {
    // Define Sigma0 as a 6-element Eigen vector containing the state variables for stress components
    Eigen::VectorXd Sigma0(6);
    Sigma0 << inputData.stateVariables[StateVariable::StressXX],
    inputData.stateVariables[StateVariable::StressYY],
    inputData.stateVariables[StateVariable::StressZZ],
    inputData.stateVariables[StateVariable::StressZY],
    inputData.stateVariables[StateVariable::StressZX],
    inputData.stateVariables[StateVariable::StressXY];
    Sigma0 *= -1.0; // reverse its direction
    // Form the kinematic hardening vector from custom state variables
    Eigen::VectorXd kinematicHardeningVector(6);
    kinematicHardeningVector << inputData.customStateVariables.at("AlphaXX"),
                                inputData.customStateVariables.at("AlphaYY"),
                                inputData.customStateVariables.at("AlphaZZ"),
                                inputData.customStateVariables.at("AlphaZY"),
                                inputData.customStateVariables.at("AlphaZX"),
                                inputData.customStateVariables.at("AlphaXY");
    // Form the fabric vector (Zvec) from custom state variables
    Eigen::VectorXd Zvec(6);
    Zvec << inputData.customStateVariables.at("ZXX"),
            inputData.customStateVariables.at("ZYY"),
            inputData.customStateVariables.at("ZZZ"),
            inputData.customStateVariables.at("ZZY"),
            inputData.customStateVariables.at("ZZX"),
            inputData.customStateVariables.at("ZXY");
    // Form the initial backstress vector (alpha_initial) from custom state variables
    Eigen::VectorXd alphaInitial(6);
    alphaInitial << inputData.customStateVariables.at("AlphaInitialXX"),
                    inputData.customStateVariables.at("AlphaInitialYY"),
                    inputData.customStateVariables.at("AlphaInitialZZ"),
                    inputData.customStateVariables.at("AlphaInitialZY"),
                    inputData.customStateVariables.at("AlphaInitialZX"),
                    inputData.customStateVariables.at("AlphaInitialXY");


    double specificVoid = inputData.stateVariables[StateVariable::VoidRatio] + 1.0;
    Eigen::MatrixXd D = computeModuli(Sigma0, kinematicHardeningVector, specificVoid);
    // Compute stress invariants
    StressInvariants invariants = computeStressInvariants(Sigma0, kinematicHardeningVector);

    // Check if the mean stress (SIGM) is below the threshold
    if (invariants.sigma_m < P_min) {
        // If SIGM is below the threshold, assume elastic behavior
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                outputData.stressStrainMatrix[i][j] = D(i, j); // Populate the elastic stiffness matrix
            }
        }
        return; // Exit the function as no further calculations are needed
    }
    Eigen::VectorXd df_dsigma = GradientOfYieldToSigma(Sigma0, kinematicHardeningVector);
    Eigen::VectorXd dg_dsigma = GradientOfPlasticToSigma(Sigma0, kinematicHardeningVector, Zvec, specificVoid);
    auto plasticModulusResult = PlasticModulus(Sigma0, kinematicHardeningVector, Zvec, alphaInitial, specificVoid);
//    double B_iso = plasticModulusResult.first;
    double Kp = plasticModulusResult.second;
    // Compute the numerator: D * df_dsigma * dg_dsigma.transpose() * D
    Eigen::MatrixXd numerator = D * dg_dsigma * df_dsigma.transpose() * D;

    // Compute the denominator: Kp + df_dsigma.transpose() * De * dg_dsigma
    double denominator = Kp + df_dsigma.transpose() * D * dg_dsigma;
    Eigen::MatrixXd Dep = D - (1.0/denominator) * numerator;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            outputData.stressStrainMatrix[i][j] = Dep(i, j);
        }
    }
    // Optional: Explicit return for clarity
    return;
}

//Eigen::VectorXd SANISANDModel::computeElasticStressIncrement(const Eigen::VectorXd& deltaStrain, Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardening, double specificVoid) const {
//    
//    double T = 0.0;
//    double DT = 1.0;
//    bool FAIL = false;
//    const double TINY = 1.0e-12;
//    const double DTMIN = 0.001;
//
//    while (std::abs(1.0 - T) > TINY) {
//        Eigen::MatrixXd D = computeModuli(stress, kinematicHardening, specificVoid); // Elastic stiffness matrix
//
//        // Compute first estimate of stress increments
//        Eigen::VectorXd DS1 = D * deltaStrain * DT;
//        Eigen::VectorXd stress1 = stress + DS1;
//
////        double specificvoidRatio1 = specificVoid * (1.0 - DT * (deltaStrain[0] + deltaStrain[1] + deltaStrain[2]));
//
//        Eigen::MatrixXd D_mid = computeModuli(stress1, kinematicHardening, specificVoid);
//
//        // Compute second estimate of stress increments
//        Eigen::VectorXd DS2 = D_mid * deltaStrain * DT;
//        // Update stress at the midpoint
//        Eigen::VectorXd stress2 = stress + 0.5 * (DS1 + DS2);
//
//        // Compute error estimates
//        Eigen::VectorXd ERR = stress2 - stress1;
//        double ENORM = ERR.squaredNorm();
//        double SNORM = stress2.squaredNorm();
//        double RELERS = std::sqrt(ENORM) / (std::sqrt(SNORM) + 1.0);
//
//        double EPS = TINY * SNORM / (SNORM + 1.0);
//        double RELERR = std::max(EPS, RELERS);
//
//        if ((RELERR > STOL) && (DT - DTMIN > TINY)) {
//            // Error too large, reduce time step
//            FAIL = true;
//            double Q = 0.9 * std::sqrt(STOL / RELERR);
//            Q = std::max(Q, std::max(0.1, DTMIN / DT));
//            DT *= Q;
//        } else {
//            // Accept the step
//            T += DT;
//            stress = stress2;
//
//            if (std::abs(1.0 - T) > TINY) {
//                // Update VOID ratio
//                specificVoid *= (1.0 - DT * (deltaStrain[0] + deltaStrain[1] + deltaStrain[2]));
//                double Q = 0.9 * std::sqrt(STOL / RELERR);
//
//                // Adjust Q based on FAIL flag
//                if (FAIL) {
//                    Q = std::min(Q, 1.0);
//                    FAIL = false; // Reset FAIL flag after successful step
//                } else {
//                    Q = std::min(Q, 1.1);
//                }
//
//                DT *= Q;
//                DT = std::max(DT, DTMIN);
//                DT = std::min(DT, 1.0 - T);
//            }
//        }
//    }
//
//    return stress; // Return the final computed stress vector
//}

void SANISANDModel::updateAInit(Eigen::VectorXd& AINIT, const Eigen::VectorXd& SIGMAT, const Eigen::VectorXd& kinHardening) {
    if (AINIT.size() != 6 || SIGMAT.size() != 6 || kinHardening.size() != 6) {
        throw std::invalid_argument("updateAInit: All input vectors must have a size of 6.");
    }

    // Compute stress invariants
    StressInvariants invariants = computeStressInvariants(SIGMAT, kinHardening);
    Eigen::VectorXd devStress = invariants.devStress; // Deviatoric stress components

    // Include RMAT and NMAT from invariants
    Eigen::VectorXd RMAT = invariants.RMAT;
    Eigen::VectorXd NMAT = invariants.NMAT;

    // Compute AAINIT
    Eigen::VectorXd AAINIT = -AINIT + kinHardening;

    // Calculate NAAINIT = NMAT^T * AAINIT
    double NAAINIT2 = NMAT.transpose().dot(AAINIT);

    // Check and update A-INIT if necessary
    if (NAAINIT2 < 0.0) {
        AINIT = kinHardening;
    }
}

double SANISANDModel::calculateRelativeErrors(
    const Eigen::VectorXd& DSIG1, const Eigen::VectorXd& DSIG2,
    const Eigen::VectorXd& Sigma1, const Eigen::VectorXd& Sigma2,
    const Eigen::VectorXd& dKinematicHardening1, const Eigen::VectorXd& dKinematicHardening2,
    const Eigen::VectorXd& kinHardening2,
    const Eigen::VectorXd& dz1, const Eigen::VectorXd& dz2,
    const Eigen::VectorXd& Zvec2) {
    // Define a small value for TINY
    const double TINY = 1e-12;

    // Compute the absolute errors in stress components
    Eigen::VectorXd stressError = (Sigma2 - Sigma1).cwiseAbs();
    Eigen::VectorXd dSigmaError = (DSIG2 - DSIG1).cwiseAbs();

    // Compute the absolute errors in kinematic hardening
    Eigen::VectorXd kinematicHardeningError = (dKinematicHardening2 - dKinematicHardening1).cwiseAbs();

    // Compute the absolute errors in fabric evolution (dz)
    Eigen::VectorXd dzError = (dz2 - dz1).cwiseAbs();

    // Calculate the relative error using the max norm for stress
    double ENORM = stressError.maxCoeff(); // Maximum error component
    double SNORM = Sigma2.cwiseAbs().maxCoeff(); // Maximum stress component
    double RELERR_Sigma = ENORM / (SNORM + TINY);

    // Calculate the relative error for kinematic hardening using the max norm
    double ENORM_kin = kinematicHardeningError.maxCoeff();
    double SNORM_kin = kinHardening2.cwiseAbs().maxCoeff();
    double RELERR_kin = ENORM_kin / (SNORM_kin + TINY);

    // Calculate the relative error for dz using the max norm
    double ENORM_dz = dzError.maxCoeff();
    double SNORM_dz = Zvec2.cwiseAbs().maxCoeff();
    double RELERR_dz = ENORM_dz / (SNORM_dz + TINY);

    // Calculate the relative error using the L2 norm for stress
    double ENORM_L2_Sigma = stressError.squaredNorm(); // Sum of squared errors for stress
    double SNORM_L2_Sigma = Sigma2.squaredNorm();      // Sum of squared stress components
    double RELERS_Sigma = std::sqrt(ENORM_L2_Sigma) / (std::sqrt(SNORM_L2_Sigma) + 1.0);

    // Calculate the relative error using the L2 norm for dz
    double ENORM_L2_dz = dzError.squaredNorm(); // Sum of squared errors for dz
    double SNORM_L2_dz = Zvec2.squaredNorm();   // Sum of squared fabric components
    double RELERS_dz = std::sqrt(ENORM_L2_dz) / (std::sqrt(SNORM_L2_dz) + 1.0);

    // Calculate EPS to avoid division by zero
    double EPS_Sigma = TINY * SNORM_L2_Sigma / (SNORM_L2_Sigma + 1.0);

    // Combine all relative errors into a final value
    double RELERR = std::max({EPS_Sigma, RELERR_Sigma, RELERR_kin, RELERR_dz, RELERS_Sigma, RELERS_dz});

    return RELERR;
}

// Stress increment calculation
void SANISANDModel::calculateStressIncrement(const InputData& inputData, OutputData& outputData) {
    int NS=0;
    int NF=0;
    //    int NC=0;
    //    bool foundElaPortion = false; //80
    const double DTMIN = 0.001;
    //    const int MAXITC = 100;
    const double TINY = 1e-12;
    double  specificVoid0, specificVoid1, specificVoid2;
    Eigen::VectorXd  kinHardening0(6), kinHardening1(6), kinHardening2(6), Zvec0(6), Zvec1(6), Zvec2(6);
    Eigen::VectorXd deltaStrain_elastic_correctted(6);
    Eigen::VectorXd deltaStrain_elastic(6);
    Eigen::VectorXd Sigma2(6);
    // Define Sigma0 as a 6-element Eigen vector containing the state variables for stress components
    // Reverse the direction of deltaStrain
    Eigen::VectorXd deltaStrain(6);
    for (int i = 0; i < 6; ++i) {
        deltaStrain(i) = -inputData.strainIncrement[i]; // Reverse the direction
    }

    double eps_p_q0 = inputData.customStateVariables.at("eps_p_q");
    double eps_p_q1, eps_p_q2;

    // Define Sigma0 as a 6-element Eigen vector containing the state variables for stress components
    Eigen::VectorXd Sigma0(6);
    Sigma0 << inputData.stateVariables[StateVariable::StressXX],
              inputData.stateVariables[StateVariable::StressYY],
              inputData.stateVariables[StateVariable::StressZZ],
              inputData.stateVariables[StateVariable::StressZY],
              inputData.stateVariables[StateVariable::StressZX],
              inputData.stateVariables[StateVariable::StressXY];

    // Reverse the direction of Sigma0
    Sigma0 *= -1; // Reverse the direction
    
    // Form the kinematic hardening vector from custom state variables
    kinHardening0 << inputData.customStateVariables.at("AlphaXX"),
    inputData.customStateVariables.at("AlphaYY"),
    inputData.customStateVariables.at("AlphaZZ"),
    inputData.customStateVariables.at("AlphaZY"),
    inputData.customStateVariables.at("AlphaZX"),
    inputData.customStateVariables.at("AlphaXY");
    // Form the fabric vector (Zvec) from custom state variables
    Zvec0 << inputData.customStateVariables.at("ZXX"),
    inputData.customStateVariables.at("ZYY"),
    inputData.customStateVariables.at("ZZZ"),
    inputData.customStateVariables.at("ZZY"),
    inputData.customStateVariables.at("ZZX"),
    inputData.customStateVariables.at("ZXY");
    // Form the initial backstress vector (alpha_initial) from custom state variables
    Eigen::VectorXd alphaInitial(6);
    alphaInitial << inputData.customStateVariables.at("AlphaInitialXX"),
    inputData.customStateVariables.at("AlphaInitialYY"),
    inputData.customStateVariables.at("AlphaInitialZZ"),
    inputData.customStateVariables.at("AlphaInitialZY"),
    inputData.customStateVariables.at("AlphaInitialZX"),
    inputData.customStateVariables.at("AlphaInitialXY");
    specificVoid0 = inputData.stateVariables[StateVariable::VoidRatio] + 1.0;
    
    updateAInit(alphaInitial, Sigma0, kinHardening0);
    
    // Compute plastic strain increments
    Eigen::VectorXd deltaStrainPlastic = deltaStrain;
        
    // Initialize variables
    double T = 0.0;
    double DT = 1.0;
    NS = 0;
    NF = 0;
    bool FAIL = false;
    while(true){
        Eigen::MatrixXd D = computeModuli(Sigma0, kinHardening0, specificVoid0);
        // Compute the elastic stress increment
        Eigen::VectorXd dSigma = D * deltaStrainPlastic;
        auto result = ComputeDlambdaAndUpdateStressInc(Sigma0, dSigma, DT, kinHardening0, Zvec0, alphaInitial, specificVoid0);
        // Call ComputeDlambdaAndUpdateStressInc and unpack the results
        auto [dSigma_1, dKinematicHardening1, dz1, dDeviatoricPlasticStrain1] =
        ComputeDlambdaAndUpdateStressInc(Sigma0, dSigma, DT, kinHardening0, Zvec0, alphaInitial, specificVoid0);
        Eigen::VectorXd Sigma1 = Sigma0 + dSigma_1;
        kinHardening1 = kinHardening0 + dKinematicHardening1;
        specificVoid1 = specificVoid0;
        Zvec1 = Zvec0 + dz1;
        eps_p_q1 = eps_p_q0 + dDeviatoricPlasticStrain1;
        //hardening and void not needed
        D = computeModuli(Sigma1, kinHardening1, specificVoid1);
        // Compute the elastic stress increment
        dSigma = D * deltaStrainPlastic; // these are not essentially needed but broght here to avoid confusion in future devs
        auto [dSigma_2, dKinematicHardening2, dz2, dDeviatoricPlasticStrain2]  = ComputeDlambdaAndUpdateStressInc(Sigma1, dSigma, DT, kinHardening1, Zvec1, alphaInitial, specificVoid1);
        
        Sigma2 = Sigma0 + 0.5 * (dSigma_1 + dSigma_2);//
        kinHardening2 = kinHardening0 + 0.5 * (dKinematicHardening1 + dKinematicHardening2);
        specificVoid2 = specificVoid0;
        Zvec2 = Zvec0 + 0.5 * (dz1 + dz2);
        eps_p_q2 = eps_p_q0 + 0.5 * (dDeviatoricPlasticStrain1 + dDeviatoricPlasticStrain2);
        
        double RELERR = calculateRelativeErrors(dSigma_1, dSigma_2, Sigma1, Sigma2, dKinematicHardening1, dKinematicHardening2, kinHardening2, dz1, dz2, Zvec2);
        
        if ((RELERR > STOL) && (DT - DTMIN > TINY)) {
            // Failure case
            FAIL = true;
            NF += 1;
            
            // Adjust step size
            double Q = 0.9 * std::sqrt(STOL / RELERR);
            Q = std::max({Q, 0.1, DTMIN / DT});
            DT *= Q;
            
            // Restart substepping
            continue; // Equivalent to GOTO 90
        } else {
            // Success case
            NS += 1;
            T += DT;
            
            // Update stresses and hardening parameters
            Sigma0 = Sigma2;
            kinHardening0 = kinHardening2;
            Zvec0 = Zvec2;
            eps_p_q0 = eps_p_q0;
            if (std::abs(1.0 - T) > TINY) {
                // Adjust VOID0 and step size
                specificVoid0 *= (1.0 + DT * (deltaStrainPlastic(0) + deltaStrainPlastic(1) + deltaStrainPlastic(2)));
                double Q = 0.9 * std::sqrt(STOL / RELERR);
                
                if (FAIL) {
                    Q = std::min(Q, 1.0);
                    FAIL = false;
                } else {
                    Q = std::min(Q, 1.1);
                }
                
                DT *= Q;
                DT = std::max(DT, DTMIN);
                
                DT = std::min(DT, 1.0 - T);
                continue; // Restart substepping
            }else{
                break; //exit the while loop upon completion
            }
        }
    }
    // Set final stress increment
    outputData.stressIncrement[0] = - Sigma0[0] - inputData.stateVariables[StateVariable::StressXX];
    outputData.stressIncrement[1] = - Sigma0[1] - inputData.stateVariables[StateVariable::StressYY];
    outputData.stressIncrement[2] = - Sigma0[2] - inputData.stateVariables[StateVariable::StressZZ];
    outputData.stressIncrement[3] = - Sigma0[3] - inputData.stateVariables[StateVariable::StressZY];
    outputData.stressIncrement[4] = - Sigma0[4] - inputData.stateVariables[StateVariable::StressZX];
    outputData.stressIncrement[5] = - Sigma0[5] - inputData.stateVariables[StateVariable::StressXY];

    // Update final custom variables in output
    // Kinematic Hardening
    outputData.updatedCustomStateVariables["AlphaXX"] = kinHardening0[0];
    outputData.updatedCustomStateVariables["AlphaYY"] = kinHardening0[1];
    outputData.updatedCustomStateVariables["AlphaZZ"] = kinHardening0[2];
    outputData.updatedCustomStateVariables["AlphaZY"] = kinHardening0[3];
    outputData.updatedCustomStateVariables["AlphaZX"] = kinHardening0[4];
    outputData.updatedCustomStateVariables["AlphaXY"] = kinHardening0[5];

    // Fabric Tensor (Zvec)
    outputData.updatedCustomStateVariables["ZXX"] = Zvec0[0];
    outputData.updatedCustomStateVariables["ZYY"] = Zvec0[1];
    outputData.updatedCustomStateVariables["ZZZ"] = Zvec0[2];
    outputData.updatedCustomStateVariables["ZZY"] = Zvec0[3];
    outputData.updatedCustomStateVariables["ZZX"] = Zvec0[4];
    outputData.updatedCustomStateVariables["ZXY"] = Zvec0[5];

    // Plastic Deviatoric Strain Measure
    outputData.updatedCustomStateVariables["eps_p_q"] = eps_p_q0;

    // Set custom state variables (AlphaInitial) to updated alphaInitial values
    outputData.updatedCustomStateVariables["AlphaInitialXX"] = alphaInitial[0];
    outputData.updatedCustomStateVariables["AlphaInitialYY"] = alphaInitial[1];
    outputData.updatedCustomStateVariables["AlphaInitialZZ"] = alphaInitial[2];
    outputData.updatedCustomStateVariables["AlphaInitialZY"] = alphaInitial[3];
    outputData.updatedCustomStateVariables["AlphaInitialZX"] = alphaInitial[4];
    outputData.updatedCustomStateVariables["AlphaInitialXY"] = alphaInitial[5];
    
    return;
}

void SANISANDModel::setCustomVariable(UMATBase::InputData& inputData) {
    // Retrieve stress components from stateVariables
    Eigen::VectorXd stress(6); // Assuming 6 components for stress (3 normal, 3 shear)
    stress[0] = inputData.stateVariables[StateVariable::StressXX];
    stress[1] = inputData.stateVariables[StateVariable::StressYY];
    stress[2] = inputData.stateVariables[StateVariable::StressZZ];
    stress[3] = inputData.stateVariables[StateVariable::StressZY]; // Shear stress
    stress[4] = inputData.stateVariables[StateVariable::StressZX];
    stress[5] = inputData.stateVariables[StateVariable::StressXY];
    stress *= -1.0; // Reverse direction

    // Initialize a dummy kinematic hardening vector (all zeros)
    Eigen::VectorXd kinHardening0 = Eigen::VectorXd::Zero(6);

    // Compute RMAT from stress and the dummy kinHardening0
    StressInvariants invariants = computeStressInvariants(stress, kinHardening0);
    Eigen::VectorXd RMAT = invariants.RMAT;

    // Set custom state variables (Alpha) to RMAT values
    inputData.customStateVariables["AlphaXX"] = RMAT[0];
    inputData.customStateVariables["AlphaYY"] = RMAT[1];
    inputData.customStateVariables["AlphaZZ"] = RMAT[2];
    inputData.customStateVariables["AlphaZY"] = RMAT[3];
    inputData.customStateVariables["AlphaZX"] = RMAT[4];
    inputData.customStateVariables["AlphaXY"] = RMAT[5];

    // Set custom state variables (AlphaInitial) to RMAT values
    inputData.customStateVariables["AlphaInitialXX"] = RMAT[0];
    inputData.customStateVariables["AlphaInitialYY"] = RMAT[1];
    inputData.customStateVariables["AlphaInitialZZ"] = RMAT[2];
    inputData.customStateVariables["AlphaInitialZY"] = RMAT[3];
    inputData.customStateVariables["AlphaInitialZX"] = RMAT[4];
    inputData.customStateVariables["AlphaInitialXY"] = RMAT[5];
}

std::tuple<double, double, Eigen::VectorXd> SANISANDModel::ComputeDlambda( const Eigen::VectorXd& stress, const Eigen::VectorXd& DStress_ela, const double DT, const Eigen::VectorXd& kinematicHardeningVector, const Eigen::VectorXd& Zvec, const Eigen::VectorXd& alpha_init, const double specificVoid, const Eigen::VectorXd& deltaStrain) const {
    //this function is only used in element test and does not have any role in FEM
    // Change the direction of stress and DStress_ela
    Eigen::VectorXd reversedStress = -stress;         // Reverse the direction of stress
    Eigen::VectorXd reversedDStress_ela = -DStress_ela; // Reverse the direction of DStress_ela
    
    // Obtain the gradients of the yield function and plastic potential
    Eigen::VectorXd df_dsigma = GradientOfYieldToSigma(reversedStress, kinematicHardeningVector);
    Eigen::VectorXd dg_dsigma = GradientOfPlasticToSigma(reversedStress, kinematicHardeningVector, Zvec, specificVoid);
    
    // Get the elastic matrix (De)
    Eigen::MatrixXd De = computeModuli(reversedStress, kinematicHardeningVector, specificVoid);
    
    // Compute norms for numerical checks
    double ATDB = df_dsigma.transpose() * De * dg_dsigma;
    
    // Retrieve both B_iso and Kp from PlasticModulus
    auto plasticModulusResult = PlasticModulus(reversedStress, kinematicHardeningVector, Zvec, alpha_init, specificVoid);
    Eigen::VectorXd B_kin = plasticModulusResult.first;
    double Kp = plasticModulusResult.second;
    
    
    // Compute volumetric strain increment (DV)
    double DEPSXX = deltaStrain[0];
    double DEPSYY = deltaStrain[1];
    double DEPSZZ = deltaStrain[2];
    double DEPSXY = deltaStrain[5]; // Assuming standard Voigt notation
    double DEPSZY = deltaStrain[3]; // Shear strain
    double DEPSZX = deltaStrain[4]; // Shear strain
    double DV = DEPSXX + DEPSYY + DEPSZZ;
    
    // Deviatoric strain components
    double SDXX = DEPSXX - DV / 3.0;
    double SDYY = DEPSYY - DV / 3.0;
    double SDZZ = DEPSZZ - DV / 3.0;
    double SDXY = DEPSXY; // Shear strain remains as is
    double SDZY = DEPSZY; // Shear strain remains as is
    double SDZX = DEPSZX; // Shear strain remains as is
    
    // Compute NSD using NMAT components
    StressInvariants invariants = computeStressInvariants(reversedStress, kinematicHardeningVector);
    Eigen::VectorXd NMAT = invariants.NMAT; // Normal to the yield surface
    double SIGM = invariants.sigma_m; // Normal to the yield surface
    double NSD = NMAT[0] * SDXX + NMAT[1] * SDYY + NMAT[2] * SDZZ +
    NMAT[3] * SDXY + NMAT[4] * SDZY + NMAT[5] * SDZX;
    
    // Compute shear modulus G
    double G = std::pow(std::max(P_min, SIGM) / patm, 0.5) * G0 * patm * std::pow((2.97 - (-1.0 + specificVoid)), 2.0)/ (1.0 + (specificVoid - 1.0));
    
    // Update K (bulk modulus) formula
    double K = std::pow(std::max(P_min, SIGM) / patm, 2.0 / 3.0) * K0 * patm * specificVoid/ (specificVoid - 1.0);
    // Compute N22 as dot product of NMAT transpose and kinematicHardeningVector
    double N22 = NMAT.transpose().dot(kinematicHardeningVector);
    
    double dLambda = (2.0 * G * NSD - K * N22 * DV) / (ATDB + Kp);
    // Ensure dLambda is non-negative
    dLambda = std::max(dLambda, 0.0);
    
    // Call ComputeDPara to calculate DPara
    double DPara = ComputeDPara(reversedStress, kinematicHardeningVector, Zvec, specificVoid);
    
    // Return dLambda, DPara, and NMAT
    return {dLambda, DPara, NMAT};
}

double SANISANDModel::ComputeDPara(const Eigen::VectorXd& stress, const Eigen::VectorXd& kinematicHardeningVector, const Eigen::VectorXd& Zvec, const double specificVoid) const {
    // Compute stress invariants
    StressInvariants invariants = computeStressInvariants(stress, kinematicHardeningVector);
    Eigen::VectorXd devStress = invariants.devStress; // Deviatoric stress components
    double SIGM = invariants.sigma_m;                 // Mean stress
    double COSTTETA = invariants.COSTTETA;

    // Include RMAT and NMAT from invariants
    Eigen::VectorXd RMAT = invariants.RMAT;
    Eigen::VectorXd NMAT = invariants.NMAT;

    if (std::fabs(Mc) < 1.0e-12) {
        throw std::runtime_error(" (M_c) is too small, causing division by zero in CPARA.");
    }
    double CPARA = Me / Mc;

    // Compute GTETC (as one-liner)
    double GTETC = std::pow((2.0 * std::pow(CPARA, 4.0)) / ((1.0 + std::pow(CPARA, 4.0)) - (1.0 - std::pow(CPARA, 4.0)) * COSTTETA), 1.0 / 4.0);

    // Compute critical state void ratio and related terms
    double ECPARA = computeCSLVoidRatio(SIGM); // Critical state void ratio
    double EPARA = specificVoid - 1.0;        // Specific void
    double SI = -ECPARA + EPARA;              // Difference between current and critical void ratios

    // Compute ATETD
    double ATETD = GTETC * Mc * std::exp(n_d * SI);

    Eigen::VectorXd NMAT2 = NMAT;                   // Normalize NMAT2

    // Compute ATETDMAT
    Eigen::VectorXd ATETDMAT = (std::sqrt(2.0 / 3.0)) * ATETD * NMAT2;
    // Compute Z * NMAT (dot product with transpose)
    double ZN = Zvec.transpose().dot(NMAT); // Dot product

    double MDAA = A0 * (1.0 + std::max(std::sqrt(3.0 / 2.0) * ZN, 0.0));

    // Compute TNMAT2
    Eigen::RowVectorXd TNMAT2 = NMAT2.transpose(); // NMAT is now a row vector

    // Compute DPARA (Dot Product)
    double DPARA = MDAA * TNMAT2.dot(ATETDMAT - kinematicHardeningVector);
    DPARA *= std::sqrt(3.0 / 2.0); // Modify by multiplying with sqrt(3.0/2.0)

    // Return DPARA and NMAT as a pair
    return DPARA;
}


// The rest of the C-interface functions would remain similar, calling the methods of `umatInstance`.
extern "C" void calculateStressIncrement(const UMATBase::InputData& inputData, UMATBase::OutputData& outputData) {
    if (umatInstance) {
        umatInstance->calculateStressIncrement(inputData, outputData);
    } else {
        std::cerr << "Error: CapModel instance not initialized." << std::endl;
    }
}

extern "C" void computeStressStrainMatrix(const UMATBase::InputData& inputData, UMATBase::OutputData& outputData) {
    if (umatInstance) {
        umatInstance->computeStressStrainMatrix(inputData, outputData);
    } else {
        std::cerr << "Error: CapModel instance not initialized." << std::endl;
    }
}

extern "C" const double* getStressIncrement(const UMATBase::OutputData& outputData) {
    return outputData.stressIncrement;
}

extern "C" const double* getStressStrainMatrixRow(const UMATBase::OutputData& outputData, int row) {
    return outputData.stressStrainMatrix[row];
}
// Extern C function that can be called by the mother code to set a custom variable
extern "C" void initializeCustomVariable(UMATBase::InputData& inputData) {
    if (!umatInstance) {
        std::cerr << "Error: CapModelUMAT instance not initialized." << std::endl;
        return;
    }

    // Call the setCustomVariable method to set the custom variable within the inputData object
    umatInstance->setCustomVariable(inputData);

//    std::cout << "Custom variable 'IsotropicHardening' initialized to: "
//              << inputData.customStateVariables.at("IsotropicHardening") << std::endl;
}
