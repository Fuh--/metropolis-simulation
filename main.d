import std.stdio;
import std.random;
import std.math;
import std.conv;
import std.range;
import std.parallelism;


void main() {
  //Input parameters
  //Exchange interaction energy
  immutable int coef = 1;
  //Length of each size
  immutable uint size = 10;
  //Number of Cycles
  immutable uint cycle = 100_000;
  //Minimum temperature (kbT)
  immutable min_temp = 0.1;
  //Maximum temperature (kbT)
  immutable max_temp = 10.0;
  //Temperature step size
  immutable step_size = 0.1;
  //Number of trials
  immutable uint trial = 1000;

  //Check input parameters
  assert(coef != 0);
  assert(size > 3);
  assert(cycle > 0);
  assert(min_temp > 0.0);
  assert(step_size > 0.0);
  assert(trial > 0);

  foreach (t; iota(min_temp, max_temp + step_size, step_size).parallel) {
    long squared_energy_sum, energy_sum;
    foreach (_i; 0..trial) {
      //Set initial states
      auto ising = new int[][](size, size);
      foreach (e1; ising) {
        foreach (ref e2; e1) {
          immutable r = uniform(0, 2);
          if (r == 0) {
            e2 = -1;
          } else {
            e2 = 1;
          }
        }
      }

      long prev_energy;
      auto prev_ising = new int[][](size, size);
      uint ty, tx;
      foreach (c; 0..cycle) {
        //Except first cycle
        if (c != 0) {
          //Invert an element
          ty = uniform(0, ising.length);
          tx = uniform(0, ising[0].length);
          ising[ty][tx] *= -1;
        }

        //Energy calculation
        long energy;
        foreach (i, e1; ising) {
          foreach (j, e2; e1) {
            foreach (dy; iota(-1, 2, 2)) {
              foreach (dx; iota(-1, 2, 2)) {
                auto x = j + dx;
                auto y = i + dy;
                //Periodic boundary condition
                if (y == -1) y = ising.length.to!int - 1;
                else if (y == ising.length) y = 0;
                if (x == -1) x = e1.length.to!int - 1;
                else if (x == e1.length) x = 0;

                energy += ising[i][j] * ising[y][x];
              }
            }
          }
        }
        energy *= -coef;

        //Except first cycle
        if (c != 0) {
          //Metropolis method
          bool flag;
          auto dE = prev_energy - energy;
          if (dE <= 0) {
            flag = true;
          } else {
            immutable r = uniform!"[]"(0.0, 1.0);
            immutable f = exp(-dE / t);
            if (r <= f) {
              flag = true;
            }
          }

          if (flag) {
            //Update process
            prev_energy = energy;
            foreach (i, e; ising) {
              prev_ising[i] = e.dup();
            }
          } else {
            //Rehection process
            ising[ty][tx] *= -1;
          }
        }

        //process of first cycle
        if (c == 0) {
          prev_energy = energy;
          foreach (i, e; ising) {
            prev_ising[i] = e.dup();
          }
          continue;
        }
      }

      squared_energy_sum += prev_energy ^^ 2;
      energy_sum += prev_energy;
    }
    //Specific heat calculation
    immutable mean_squared_energy = squared_energy_sum.to!double / trial;
    immutable squared_mean_energy = (energy_sum.to!double / trial) ^^ 2;
    immutable specific_heat = (mean_squared_energy - squared_mean_energy) / (t ^^ 2);

    writefln("%.1f\t%.2f", t, specific_heat);
  }
}
