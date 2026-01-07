#ifndef OUTPUT_H
#define OUTPUT_H


void open_output(char filepath[], FILE **file);
void close_output(FILE *file);
void print_output(FILE *file, struct stellar_properties *star_props, double R_e_BM, double R_p_BM, double Mass_BM, double Mass_0_BM, struct stellar_properties *DM_props, double R_e_DM, double R_p_DM, double Mass_DM, double Mass_0_DM, double fdm, double Omega_BM, double Omega_DM, double J_BM, double J_DM, double r_ratio_BM, double r_ratio_DM, double Omega_K_BM, double Omega_K_DM, double T_BM, double T_DM, double W_BM, double W_DM);

char *getDiffRotTypeString(enum DiffRotType *type);
int getDiffRotTypeFromString(const char *str);

void hdf5_save_var(double(*s_gp), int sdiv, double(*mu), int mdiv, char filename[], struct evolution_variables *ev, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, double *r_ratio_BM, double *r_ratio_DM, double *Omega_BM, double *Omega_DM, double *r_e_BM, double *r_e_DM, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, char *input_parfile);

void hdf5_read_var(double *s_gp, int *sdiv, double *mu, int *mdiv, char filename[], struct evolution_variables *ev, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, double *r_ratio_BM, double *r_ratio_DM, double *Omega_BM, double *Omega_DM, double *r_e_BM, double *r_e_DM, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM);

void write_profile(double **variable, const char *filename);

#endif