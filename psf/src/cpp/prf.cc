#include "include/prf.hh"
#include "include/utils.hh"

PRF::PRF(double wavelength, double focal_ratio, 
            std::pair<double,double> variance,
            std::pair<double,double> standard_deviation,
            std::pair<std::uint32_t, std::uint32_t> subpixel_res, 
            std::pair<std::uint32_t, std::uint32_t> map_res
        ) : intensity(wavelength, focal_ratio, standard_deviation, subpixel_res, map_res),
            x_vals(map_res.first), y_vals(map_res.second) { }

void PRF::initialize(std::pair<double, double> max_wavenumber) {
    intensity.execute(max_wavenumber);

    intensity.fill_x_wavenumbers(x_vals.data());
    intensity.fill_y_wavenumbers(y_vals.data());

    auto wavenumber_grid = meshgrid(x_vals, y_vals);
    std::vector<double> const& x_wavenumbers = wavenumber_grid.first;
    std::vector<double> const& y_wavenumbers = wavenumber_grid.second;

    alglib::real_1d_array x_vals_ag;
    alglib::real_1d_array y_vals_ag;
    alglib::real_1d_array i_vals_ag;

    std::vector<double> const& intensities = intensity.get_all_intensities();

    x_vals_ag.setcontent(x_wavenumbers.size(), x_wavenumbers.data());
    y_vals_ag.setcontent(y_wavenumbers.size(), y_wavenumbers.data());
    i_vals_ag.setcontent(intensities.size(), intensities.data());

    alglib::spline2dbuildbicubicv(
        x_vals_ag,
        x_wavenumbers.size(),
        y_vals_ag,
        y_wavenumbers.size(),
        i_vals_ag,
        intensities.size(),
        spline
    );
}

