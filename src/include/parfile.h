#ifndef PARFILE_H
#define PARFILE_H


#include "consts.h"
#include <stdbool.h>

void update_int(const char *value, void *var);
void update_double(const char *value, void *var);
void update_string(const char *value, void *var);
void update_bool(const char *value, void *var);

int compare_config_entries(const void *a, const void *b);
int binary_search(const char *key, ConfigEntry *config_entries, int size);

void readConfig(const char *filename, struct EOS *eosBM, struct stellar_properties *BM_prop, double *accuracy, double *r_ratio_BM, bool *OUTPUT1D, bool *OUTPUT2D, struct EOS *eosDM, struct stellar_properties *DM_prop, double *r_ratio_DM, double *fdm_target, double *fdm_tol, double *MDM, int *verbose, bool *counter, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, double *r_step, char *id_file, char *output_name, bool *load_id);

int readToken(char *pair, char **key, char **value);

int writeConfig(const char *filename, struct EOS *eosBM, struct stellar_properties *BM_prop, double accuracy, double r_ratio_BM, bool OUTPUT1D, bool OUTPUT2D, struct EOS *eosDM, struct stellar_properties *DM_prop, double r_ratio_DM, double fdm_target, double fdm_tol, double MDM, int verbose, bool counter, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, double r_step);

#endif