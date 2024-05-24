#include "mex.h"

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

#include <iostream>
#include <fstream>
#include <chrono>

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own directory here.
#ifndef PINOCCHIO_MODEL_DIR
    #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

void copyPointerToEigenVector(Eigen::VectorXd& a, double* b, size_t n, size_t offset = 0) {
    for (size_t i = 0; i < n; i++) {
        a(i) = b[i + offset];
    }
}

void copyEigenVectorToPointer(double* a, Eigen::VectorXd& b, size_t n, size_t offset = 0) {
    for (size_t i = 0; i < n; i++) {
        a[i + offset] = b(i);
    }
}

void copyEigenMatrixToPointer(double* a, Eigen::MatrixXd& b, size_t m, size_t n, size_t offset = 0) {
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            a[j * m + i + offset] = b(i,j);
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    using namespace pinocchio;

    if (nrhs != 4) {
        mexErrMsgTxt("Expect 4 inputs: urdf_filename, q, qd, u!\n");
    }
    
    const std::string urdf_filename = mxArrayToString(prhs[0]);
    
    // Load the URDF model
    Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    model.gravity.linear(Eigen::Vector3d(0,0,0));    
    
    // // Build a data related to model
    Data data(model);

    size_t num_samples = 0;
    size_t q_m = mxGetM(prhs[1]);
    size_t q_n = mxGetN(prhs[1]);
    size_t v_m = mxGetM(prhs[2]);
    size_t v_n = mxGetN(prhs[2]);
    size_t tau_m = mxGetM(prhs[3]);
    size_t tau_n = mxGetN(prhs[3]);

    if (q_m != model.nv || v_m != model.nv || tau_m != model.nv) {
        mexPrintf("Expect matrix with %d rows!\n", model.nv);
        mexErrMsgTxt("\n");
    }

    if (q_n != v_n || q_n != tau_n || v_n != tau_n) {
        mexErrMsgTxt("Expect 3 matrix with same number of columns!\n");
    }    
    else {
        num_samples = q_n;
    }

    plhs[0] = mxCreateNumericMatrix(model.nv, num_samples, mxDOUBLE_CLASS, mxREAL);

    if (nlhs > 1) {
        mwSize dims[3] = {(mwSize)model.nv, (mwSize)model.nv, num_samples};
        plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    }

    size_t i = 0;

    #pragma omp parallel for shared(prhs, plhs, model, data) private(i) schedule(dynamic)
    for (i = 0; i < num_samples; i++) {
        // Allocate joint configuration as well as joint velocity and torque
        Eigen::VectorXd q = randomConfiguration(model);
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
        Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);
    
        // Allocate result container
        Eigen::MatrixXd djoint_acc_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
        Eigen::MatrixXd djoint_acc_dv = Eigen::MatrixXd::Zero(model.nv, model.nv);
        Eigen::MatrixXd djoint_acc_dtau = Eigen::MatrixXd::Zero(model.nv, model.nv);

        // Setup joint configuration as well as joint velocity and torque
        copyPointerToEigenVector(q, (double*)mxGetData(prhs[1]), model.nv, i * model.nv);
        copyPointerToEigenVector(v, (double*)mxGetData(prhs[2]), model.nv, i * model.nv);
        copyPointerToEigenVector(tau, (double*)mxGetData(prhs[3]), model.nv, i * model.nv);

        // Computes the forward dynamics (ABA) derivatives for all the joints of the robot
        computeABADerivatives(model, data, q, v, tau, djoint_acc_dq, djoint_acc_dv, djoint_acc_dtau);

        // Copy result to output array
        copyEigenVectorToPointer((double*)mxGetData(plhs[0]), data.ddq, model.nv, i * model.nv);
        if (nlhs > 1) {
            copyEigenMatrixToPointer((double*)mxGetData(plhs[1]), djoint_acc_dq, model.nv, model.nv, i * model.nv * model.nv);
            copyEigenMatrixToPointer((double*)mxGetData(plhs[2]), djoint_acc_dv, model.nv, model.nv, i * model.nv * model.nv);
            copyEigenMatrixToPointer((double*)mxGetData(plhs[3]), djoint_acc_dtau, model.nv, model.nv, i * model.nv * model.nv);
        }
    }
}