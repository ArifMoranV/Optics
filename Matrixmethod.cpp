#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <numeric>  // for std::partial_sum

int main() {
    int N = 2;                      // number of thin films (intermediate layers)
    int pol = 0;                    // polarization: 0 = TM, 1 = TE
    double lam_ini = 500.0;         // initial wavelength [nm]
    double lam_fin = 830.0;         // final wavelength [nm]
    double ang_ini = 0.0;           // initial incidence angle [deg]
    double ang_fin = 80.0;          // final incidence angle [deg]
    double nsub = 1.5;              // substrate refractive index
    double nsup = 1.0;              // superstrate refractive index
    std::vector<double> d = {50.0}; // thickness of each thin film [nm]

    // -------------------------------
    // Wavelength and angle grids
    // -------------------------------
    std::vector<std::complex<double>> lam;
    for (double l = lam_ini; l <= lam_fin; l += 1.0)
        lam.emplace_back(l, 0.0); // ComplexF64 in Julia → std::complex<double>

    std::vector<std::complex<double>> k0;
    const double pi = M_PI;
    for (auto& l : lam)
        k0.push_back(2.0 * pi / l);  // elementwise 2π/λ

    std::vector<double> ang_inc;
    for (double a = ang_ini; a <= ang_fin; a += 0.1)
        ang_inc.push_back(a);

    std::vector<double> ang;
    for (auto& a : ang_inc)
        ang.push_back(a * pi / 180.0); // deg2rad

    // -------------------------------
    // Film position vector z
    // -------------------------------
    double z0 = 0.0;
    std::vector<double> z(d.size());
    std::partial_sum(d.begin(), d.end(), z.begin());

    // prepend z0
    z.insert(z.begin(), z0);

    // Display a few results
    std::cout << "First few λ values: ";
    for (size_t i = 0; i < 5 && i < lam.size(); ++i)
        std::cout << lam[i] << " ";
    std::cout << "\n";

    std::cout << "z vector: ";
    for (auto val : z) std::cout << val << " ";
    std::cout << "\n";

    return 0;
}
