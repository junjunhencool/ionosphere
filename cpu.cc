/* This implementation is based on the following articles:
 *
 * 1. D. A. Randall, T. D. Ringler, R. P. Heikes, P. Jones, and J. Baumgardner,
 * "Climate Modeling with Spherical Geodesic Grids," Computing in Science & 
 * Engineering, vol. 4, no. 5, pp. 32-41, Sep-2002.
 *
 * 2. J. J. Simpson, R. P. Heikes, and A. Taflove, "FDTD modeling of a novel 
 * ELF Radar for major oil deposits using a three-dimensional geodesic grid of
 * the Earth-ionosphere waveguide,” IEEE Trans. Antennas Propag., vol. 54, 
 * no. 6, pp. 1734–1741, Jun. 2006.
 *
 * Kyungwon Chun <ruddyscent@gmail.com>
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "hdf5.h"

#define er(i,j) er[(i)*er_jm+(j)]
#define h1(i,j) h1[(i)*h1_jm+(j)]
#define h2(i,j) h2[(i)*h2_jm+(j)]
#define h3(i,j) h3[(i)*h3_jm+(j)]

using namespace std;

class Pole {
private:
  double er;
  array<double, 5> h;
  
public:
  void
  write(const string& name) {
    hid_t file, space, dset;
    hsize_t dims[1];
    
    file = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    dims[0] = 1;
    space = H5Screate_simple(1, dims, NULL);
    dset = H5Dcreate(file, "er", H5T_NATIVE_DOUBLE, space, 
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      &er);
    H5Dclose(dset);
    H5Sclose(space);

    dims[0] = 5;
    space = H5Screate_simple(1, dims, NULL);
    dset = H5Dcreate(file, "h", H5T_NATIVE_DOUBLE, space, 
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, h.data());
    H5Dclose(dset);
    H5Sclose(space);
    H5Fclose(file);
  }

  void 
  init(double val) {
    er = val;
    h.fill(val);
  }

  void
  set_h(const array<double, 5>& h_in) {
    for (auto i = 0; i < 5; i++) {
      h[i] = h_in[i];
    }
  }
  // void 
  // syn_north(const Panel& p1, const Panel& p2, 
  //           const Panel& p3, const Panel& p4, const Panel& p5) {
  //   h[0] = p1.get_pole_h1();
  //   h[1] = p2.get_pole_h1();
  //   h[2] = p3.get_pole_h1();
  //   h[3] = p4.get_pole_h1();
  //   h[4] = p5.get_pole_h1();
  //  }

  // void 
  // sync_south(const Panel& p1, const Panel& p2, 
  //            const Panel& p3, const Panel& p4, const Panel& p5) {
  //   h[0] = p1.get_pole_h3();
  //   h[1] = p2.get_pole_h3();
  //   h[2] = p3.get_pole_h3();
  //   h[3] = p4.get_pole_h3();
  //   h[4] = p5.get_pole_h3();
  // }

  void 
  update_er(double ds, double dt) {
    double coeff = 2. * sqrt(15. - 6 * sqrt(5)) * dt / (5 * ds);
    er += coeff * (h[0] + h[1] + h[2] + h[3] + h[4]);
  }

  double
  get_er() const {
    return er;
  }
}; // class Pole


class Panel 
{    
private:
  size_t er_im, h1_im, h2_im, h3_im;
  size_t er_jm, h1_jm, h2_jm, h3_jm;
  size_t er_size, h1_size, h2_size, h3_size;

public:
  double *er, *h1, *h2, *h3;

  void
  write(const string& name) {
    hid_t file, space, dset;
    hsize_t dims[2];
    
    file = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = er_im;
    dims[1] = er_jm;
    space = H5Screate_simple(2, dims, NULL);
    dset = H5Dcreate(file, "er", H5T_NATIVE_DOUBLE, space, 
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, er);
    H5Dclose(dset);
    H5Sclose(space);

    dims[0] = h1_im;
    dims[1] = h1_jm;
    space = H5Screate_simple(2, dims, NULL);
    dset = H5Dcreate(file, "h1", H5T_NATIVE_DOUBLE, space, 
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, h1);
    H5Dclose(dset);
    H5Sclose(space);

    dims[0] = h2_im;
    dims[1] = h2_jm;
    space = H5Screate_simple(2, dims, NULL);
    dset = H5Dcreate(file, "h2", H5T_NATIVE_DOUBLE, space, 
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, h2);
    H5Dclose(dset);
    H5Sclose(space);

    dims[0] = h3_im;
    dims[1] = h3_jm;
    space = H5Screate_simple(2, dims, NULL);
    dset = H5Dcreate(file, "h3", H5T_NATIVE_DOUBLE, space, 
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, h3);
    H5Dclose(dset);
    H5Sclose(space);

    H5Fclose(file);
  }

  Panel(int im, int jm) 
  {
    er_im = im;
    er_jm = jm;
    er_size = er_im * er_jm;

    h1_im = im - 1;
    h1_jm = jm - 2;
    h1_size = h1_im * h1_jm;

    h2_im = im - 1;
    h2_jm = jm - 1;
    h2_size = h2_im * h2_jm;

    h3_im = im - 2;
    h3_jm = jm - 1;
    h3_size = h3_im * h3_jm;
    er = new double[er_size];
    h1 = new double[h1_size];
    h2 = new double[h2_size];
    h3 = new double[h3_size];
  }
  
  ~Panel() {
    delete[] er;
    delete[] h1;
    delete[] h2;
    delete[] h3;
  }
  
  void
  init(double val) {
    for (size_t i = 0; i < er_size; i++) {
      er[i] = val;
    }
    for (size_t i = 0; i < h1_size; i++) {
      h1[i] = val;
    }
    for (size_t i = 0; i < h2_size; i++) {
      h2[i] = val;
    }
    for (size_t i = 0; i < h3_size; i++) {
      h3[i] = val;
    }
  }
  
  void
  sync(const Panel& west, const Panel& east, const Pole& north, const Pole& south)
  {
    er(0,er_jm-2) = north.get_er();

    for (size_t j = 1; j < er_jm - 1; j++) {
      er(j,er_jm-1) = east.er(1,er_jm-1-j);
      er(j+er_jm-2,er_jm-1) = east.er(j,1);
      er(er_im-1,er_jm-1-j) = east.er(er_jm-2+j,1);

      er(0,er_jm-2-j) = west.er(j,er_jm-2);
      er(j,0) = west.er(er_jm-2+j,er_jm-1);
      er(er_jm-3+j,0) = west.er(er_im-1,er_jm-1-j);
    }
    
    er(er_im-2,0) = south.get_er();
  }

  void
  update_er(double ds, double dt)
  {
    const double coeff1 =  2. * dt / (3. * ds);
    for (size_t i = 1; i < er_im - 1; i++) {
      for (size_t j = 1; j < er_jm - 2; j++) {
        er(i,j) += coeff1 * (h1(i,j-1) - h1(i-1,j-1) + h2(i,j) 
                             - h2(i-1,j-1) - h3(i-1,j-1) + h3(i-1,j));
        if (er(i,j) != 0) {
          cout << "i: " << i << " j: " << j << endl;
          cout << "h1(i,j-1)" << h1(i,j-1) << endl;
          cout << "h1(i-1,j-1)" << h1(i-1,j-1) << endl;
          cout << "h2(i,j)" << h2(i,j) << endl;
          cout << "h2(i-1,j-1)" <<  h2(i-1,j-1) << endl;
          cout << "h3(i-1,j-1)" <<  h3(i-1,j-1) << endl;
          cout << "h3(i-1,j)" <<  h3(i-1,j) << endl;
        }
      }
    }
    
    // for (size_t i = 1, j = er_jm-1; i < er_jm - 2; i++) {
    //   er(i,j) += coeff1 * (h1(i,j-1) - h1(i-1,j-1) + h2(i,j) 
    //                        - h2(i-1,j-1) - h3(i-1,j-1) + h3(i-1,j));
    // }

    // for (size_t i = 9, j = er_jm-1; i < er_im - 1; i++) {
    //   er(i,j) += coeff1 * (h1(i,j-1) - h1(i-1,j-1) + h2(i,j) 
    //                        - h2(i-1,j-1) - h3(i-1,j-1) + h3(i-1,j));
    // }

    // const double coeff2 = 2. * sqrt(15. - 6 * sqrt(5)) * dt / (5 * ds);
    // for (size_t i = er_jm-2, j = er_jm-2; i < er_im - 1; i += er_jm-2) {
    //   er(i,j) += coeff2 * (h1(i,j-1) - h1(i-1,j-1) 
    //                        + h2(i,j) - h2(i-1,j-1) - h3(i-1,j-1));
    // }
  }

  void
  update_h1(double ds, double dt) {
    for (size_t i = 0; i < h1_im; i++) {
      for (size_t j = 0; j < h1_jm; j++) {
        h1(i,j) += dt / ds * (er(i+1,j) - er(i,j));
      }
    }
  }           

  void
  update_h2(double ds, double dt) {
    for (size_t i = 0; i < h2_im; i++) {
      for (size_t j = 0; j < h2_jm; j++) {
        h2(i,j) += dt / ds * (er(i+1,j+1) - er(i,j));
      }
    }
  }           

  void
  update_h3(double ds, double dt) {
    for (size_t i = 0; i < h3_im; i++) {
      for (size_t j = 0; j < h3_jm; j++) {
        h3(i,j) += dt / ds * (er(i,j+1) - er(i,j));
      }
    }
  }
  
  double
  get_pole_h1() const {
    return h1(0,h1_jm-1);
  }

  double
  get_pole_h3() const {
    return h3(h3_im-1,0);
  }
}; // class Panel

int 
main(int argc, char* argv[]) {
  double ds = 4e-2;
  double dt = 2e-2;

  // for 642 cells
  array<Panel*, 5> panel;
  for (size_t i = 0; i < panel.size(); i++) {
    panel[i] = new Panel(18, 10);
    panel[i]->init(0);
  }

  Pole north = Pole();
  Pole south = Pole();
  north.init(0);
  south.init(0);

  // syncronization of er
  panel[0]->sync(*(panel[panel.size()-1]), *(panel[1]), north, south);
  for (size_t i = 1; i < panel.size() - 1; i++) {
    panel[i]->sync(*(panel[i-1]), *(panel[i+1]), north, south);
  }
  panel[panel.size()-1]->sync(*(panel[panel.size()-2]), *(panel[0]), north, south);

  // update h
  for (size_t i = 0; i < panel.size(); i++) {
    panel[i]->update_h1(ds, dt);
    panel[i]->update_h2(ds, dt);
    panel[i]->update_h3(ds, dt);
  }
  
  // syncronization of h
  array<double, 5> north_buf, south_buf;
  for (size_t i = 0; i < panel.size(); i++) {
       north_buf[i] = panel[i]->get_pole_h1();
       south_buf[i] = panel[i]->get_pole_h3();
  }
  north.set_h(north_buf);
  south.set_h(south_buf);

  // update er
  for (size_t i = 0; i < panel.size(); i++) {
    panel[i]->update_er(ds, dt);
  }
  north.update_er(ds, dt);
  south.update_er(ds, dt);

  // output ot file
  array<string, 5> fname;
  for (size_t i = 0; i < panel.size(); i++) {
    fname[i] = string("panel") + to_string(i) + string(".h5");
    panel[i]->write(fname[i]);
  }

  string fname_np("north.h5");
  north.write(fname_np);
  string fname_sp("south.h5");
  north.write(fname_sp);

  return EXIT_SUCCESS;
}
