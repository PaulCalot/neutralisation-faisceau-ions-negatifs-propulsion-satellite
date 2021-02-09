
int add_flux_particles(fmatrix & p, int & n_active_particles, const double temperature, const double v_drift, const double mass, const double n_add){
    // add flux of particles from a boundary of lenght 'y_thruster' at the left of a rectangular domain
    // 'p' is an array containing all position and velocities of the particles and it's accessed using row major order

    double f_n_add = floor(n_add);
    int n_new = (r_unif() <= (n_add - f_n_add) ?  f_n_add + 1 : f_n_add);

    double v_temperature = sqrt(q * temperature / mass); //temperature in eV
    
    for (int i = n_active_particles; i < n_active_particles + n_new; i++)
    {   
        p.val[i * 6 + 3] = v_temperature * sqrt(- 2 *  log(r_unif())) + v_drift;
        p.val[i * 6 + 4] = r_norm(0.0, v_temperature); //r_norm(double mean, double sigma): random samples from a normal (Gaussian) distribution
        p.val[i * 6 + 5] = r_norm(0.0, v_temperature);

        p.val[i * 6 + 0] = dt * p.val[i * 6 + 3] * r_unif(); //r_unif(): random samples from a uniform distribution; dt is the time step of the simulation
        p.val[i * 6 + 1] = y_thruster * r_unif(); 
        p.val[i * 6 + 2] = 0.0;
    }
    n_active_particles += n_new;
    
    return n_new;
}

