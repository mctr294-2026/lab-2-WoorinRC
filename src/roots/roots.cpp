#include "roots.hpp"
#include <functional>
#include <cmath>
#include <limits>
#include <stdexcept>

/*This method is a Bracketing Method which means it 
* Narrows down the interval that has the root, starting with IVT 
* Kind of like the squeeze theorem. 
* Plan: start with two points a and b where the y is opposite signs 
* this means there is a root in between them.
* find the midpoint, see the midpoint's y value, then choose new bounds
* depending on sign or accuracy.
*/
bool bisection(std::function<double(double)> f, // f is the function.
               double a, double b, // Initial interval values.
               double *root) { // root is like storage.

                // Plan: Evaluate end points with if else statements.
                double fa = f(a);
                double fb = f(b);
                // Case 1: If fb or fa happen to = 0 then it is stored as the root.
                if (fa == 0.0)
                {
                    *root = a;
                    return true;
                }
                if (fb == 0.0)
                {
                    *root = b;
                    return true;
                }
                // Case 2: If fa and fb have the same sign, then it is invalid.
                if (fa * fb > 0.0) {
                    return false; // invalid interval.
                }
                // Case 3: If all passes, calculate midpoint and solve to a tolerance.
                const double tolerance = 1e-6;
                const int max_cap = 500;// Incase of infinite loop ruining the code.

                for (int i = 0; i < max_cap; i++) {
                    double midpoint = (a + b) / 2.0;
                    double fmid = f(midpoint);
                    // Check if midpoint is root or within tolerance.
                    if (fabs(fmid) < tolerance) {
                        *root = midpoint;
                        return true;
                    }
                    // Narrow down the interval.
                    if (fa * fmid < 0.0) {
                        b = midpoint;
                        fb = fmid;
                    } else {
                        a = midpoint;
                        fa = fmid;
                    }
                }
                *root = 0.5 * (a + b);
return true;


    return false;
}

/* False position method is similar to the bisection
* It finds two points a and b that have opposite signs 
* You can calculate midpoint with a different formula
* c = a - [f(a).(b - a)]/[f(b) - f(a)]
* From my understanding, it should be the same except the mid
* calculation is different.
* Because there is an issue with one side sticking, this is also a change 
* we need to include.
*/
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {// Plan: Evaluate end points with if else statements.
                
    double fa = f(a);
    double fb = f(b);
    // Case 1: If fb or fa happen to = 0 then it is stored as the root.
    if (fa == 0.0)
    {
        *root = a;
        return true;
    }
    if (fb == 0.0)
    {
        *root = b;
        return true;
    }
    // Case 2: If fa and fb have the same sign, then it is invalid.
    if (fa * fb > 0.0) {
        return false; // invalid interval.
    }
    // Case 3: If all passes, calculate midpoint and solve to a tolerance.
    const double tolerance = 1e-6;
    const int max_cap = 500;// Incase of infinite loop ruining the code.
    int count_a = 0;
    int count_b = 0;

    for (int i = 0; i < max_cap; i++) {
        double midpoint = a - (fa * (b - a)) / (fb - fa);
        double fmid = f(midpoint);
        // Check if midpoint is root or within tolerance.
        if (std::fabs(fmid) < tolerance) {
            *root = midpoint;
            return true;
        }
        // Narrow down the interval.
        // A change for regula falsi is that prevented the stuck endpoint
        // by using the Illinois fix.
        // This essentially makes the flatter slopes behave more normally
        if (fa * fmid < 0.0) {
            b = midpoint;
            fb = fmid;
            count_b = 0;
            count_a++;
            if (count_a > 1) fa /= 2.0;
        } else {
            a = midpoint;
            fa = fmid;
            count_a = 0;
            count_b++;
            if (count_b > 1) fb /= 2.0;
        }
    }
    *root = a - fa * (b - a) / (fb - fa);
    return true;
                 

    return false;
}

/* This is an Open Method which means they don't use a interval.
* Instead, this one in particular uses a functions derivative with an 
* initial guess. Then if the value of x at a point is similar to a future point
* minus the ratio of their value and derivative.
*/
bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {
                        // Plan: start with an intial guess, find derivative at this point.
                        // Use the formula to find the next approximation.
                        // At a point, the approximation becomes more and more accurate until it is the guess.
    const double tol = 1e-6;
    const int max_iter = 100;
 
    // Chooses the x value it gives us as our starting point.
    double x = c;

    for (int i = 0; i < max_iter; ++i)
    {
        double fx = f(x);
        double gx = g(x);

        if (gx == 0.0)
        {
            return false;
        }

        double x_new = x - fx / gx;

        if (x_new < a || x_new > b)
        {
            return false;
        }

        if (std::fabs(x_new - x) < tol || std::fabs(f(x_new)) < tol)
        {
            *root = x_new;
            return true;
        }

        x = x_new;
    }
    return false;
}

/* The Secant Method is similar but instead it uses a secant line to approximate.
* you start with initial guesses of x0 and x1 
* It substitutes the derivative with a secant
*/
bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
               const double tol = 1e-6;
    const int max_iter = 100;

    double x_prev = c;
    double x_curr = c + 1e-4;

    for (int i = 0; i < max_iter; ++i)
    {
        double f_prev = f(x_prev);
        double f_curr = f(x_curr);

        double denom = f_curr - f_prev;
        if (denom == 0.0)
        {
            return false;
        }

        double x_next = x_curr
                      - f_curr * (x_curr - x_prev) / denom;

        if (x_next < a || x_next > b)
        {
            return false;
        }

        if (std::fabs(x_next - x_curr) < tol || std::fabs(f(x_next)) < tol)
        {
            *root = x_next;
            return true;
        }

        x_prev = x_curr;
        x_curr = x_next;
    }
  
    return false;

}

