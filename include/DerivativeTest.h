#ifndef DERIVATIVE_TEST_H
#define DERIVATIVE_TEST_H

// forward declaration
class Objective;

// Check derivatives via central finite difference
class DerivativeTest
{    
public:
    DerivativeTest() = delete;
    static void checkGradient(Objective &obj);
    static void checkHessian(Objective &obj);

};

#endif