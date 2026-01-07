#include "parfile.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "consts.h"

void update_int(const char* value, void* var) {
    *(int*)var = atoi(value);
}

void update_double(const char* value, void* var) {
    *(double*)var = atof(value);
}

void update_string(const char* value, void* var) {
    strcpy((char*)var, value);
}

void update_bool(const char* value, void* var) {
    if (strcasecmp(value, "true") == 0 || strcmp(value, "1") == 0) {
        *(bool*)var = true;
    } else if (strcasecmp(value, "false") == 0 || strcmp(value, "0") == 0) {
        *(bool*)var = false;
    } else {
        printf("Error: Invalid boolean value '%s'.\n", value);
        exit(1);
    }
}

int readToken(char* pair, char** key, char** value) {
    char* trimmed_line = pair;
    while (*trimmed_line == ' ' || *trimmed_line == '\t') {
        trimmed_line++;
    }

    // Skip comments and empty lines
    if (trimmed_line[0] == '#' || trimmed_line[0] == '\n' || trimmed_line[0] == '\0') {
        return 100;
    }

    char delimiters[] = " =#\t\n";
    *key = strtok(trimmed_line, delimiters);
    *value = strtok(NULL, delimiters);
    return 1;
}

// Comparison function for sorting configuration entries
int compare_config_entries(const void* a, const void* b) {
    ConfigEntry* entry_a = (ConfigEntry*)a;
    ConfigEntry* entry_b = (ConfigEntry*)b;
    return strcasecmp(entry_a->key, entry_b->key);
}

void readConfig(const char* filename, struct EOS* eosBM, struct stellar_properties* BM_prop, double* accuracy, double* r_ratio_BM, bool* OUTPUT1D, bool* OUTPUT2D, struct EOS* eosDM, struct stellar_properties* DM_prop, double* r_ratio_DM, double* fdm_target, double* fdm_tol, double* MDM, int* verbose, bool* counter, struct DiffRotParams* DiffRotBM, struct DiffRotParams* DiffRotDM, double* r_step, char* id_file, char* output_name, bool* load_id) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Could not open config file %s\n", filename);
        exit(1);
    }

    // Array of configuration entries
    int num_entries = 35;
    ConfigEntry config_entries[] = {
        {"id_file", update_string, id_file},
        {"EoS_BM_type", update_int, &eosBM->type},
        {"EoS_BM_file", update_string, &eosBM->filepath},
        {"EoS_BM_N", update_double, &eosBM->n_P},
        {"EoS_DM_type", update_int, &eosDM->type},
        {"EoS_DM_file", update_string, &eosDM->filepath},
        {"EoS_DM_N", update_double, &eosDM->n_P},
        {"accuracy", update_double, accuracy},
        {"BM_central_energy", update_double, &BM_prop->e_center},
        {"DM_central_energy", update_double, &DM_prop->e_center},
        {"DM_particle_mass", update_double, MDM},
        {"DM_fraction", update_double, fdm_target},
        {"DM_fraction_tol", update_double, fdm_tol},
        {"BM_rotation_type", update_int, &BM_prop->RotType},
        {"BM_A", update_double, &DiffRotBM->A},
        {"BM_B", update_double, &DiffRotBM->B},
        {"BM_lambda1", update_double, &DiffRotBM->lambda1},
        {"BM_lambda2", update_double, &DiffRotBM->lambda2},
        {"BM_p", update_double, &DiffRotBM->p},
        {"BM_csi", update_double, &DiffRotBM->csi},
        {"DM_rotation_type", update_int, &DM_prop->RotType},
        {"DM_A", update_double, &DiffRotDM->A},
        {"DM_B", update_double, &DiffRotDM->B},
        {"DM_lambda1", update_double, &DiffRotDM->lambda1},
        {"DM_lambda2", update_double, &DiffRotDM->lambda2},
        {"DM_p", update_double, &DiffRotDM->p},
        {"DM_csi", update_double, &DiffRotDM->csi},
        {"counter_rotation", update_bool, counter},
        {"r_ratio_BM", update_double, r_ratio_BM},
        {"r_ratio_DM", update_double, r_ratio_DM},
        {"r_step", update_double, r_step},
        {"1Doutput", update_bool, OUTPUT1D},
        {"2Doutput", update_bool, OUTPUT2D},
        {"output_name", update_string, output_name},
        {"verbose", update_int, verbose}};

    qsort(config_entries, num_entries, sizeof(ConfigEntry), compare_config_entries);

    // Read all parameters from the config file
    char line[256];
    int index;
    while (fgets(line, sizeof(line), file)) {
        // Read the key-value pair
        char *key, *value;
        if (readToken(line, &key, &value) != 1) {
            continue;  // Skip comments and empty lines
        }

        // Perform a binary search to find the key
        index = binary_search(key, config_entries, num_entries);
        if (index != -1) {
            config_entries[index].update_func(value, config_entries[index].variable);
        } else {
            printf("Parameter '%s' not recognized.\n", key);
            exit(1);
        }
    }

    fclose(file);

    // Derived parameters
    *MDM *= MeVc2_to_grams;
    eosBM->Gamma_P = 1.0 + 1.0 / eosBM->n_P;
    eosDM->Gamma_P = 1.0 + 1.0 / eosDM->n_P;

    *load_id = true;
    if (strcasecmp(id_file, "0") == 0) {
        *load_id = false;
    }
}

// Function to perform binary search on the configuration entries
int binary_search(const char* key, ConfigEntry* config_entries, int size) {
    int low = 0;
    int high = size - 1;

    while (low <= high) {
        int mid = (low + high) / 2;
        int cmp = strcasecmp(key, config_entries[mid].key);

        if (cmp == 0) {
            return mid;  // Key found
        } else if (cmp < 0) {
            high = mid - 1;  // Search left
        } else {
            low = mid + 1;  // Search right
        }
    }

    return -1;  // Key not found
}

int writeConfig(const char* filename, struct EOS* eosBM, struct stellar_properties* BM_prop, double accuracy, double r_ratio_BM, bool OUTPUT1D, bool OUTPUT2D, struct EOS* eosDM, struct stellar_properties* DM_prop, double r_ratio_DM, double fdm_target, double fdm_tol, double MDM, int verbose, bool counter, struct DiffRotParams* DiffRotBM, struct DiffRotParams* DiffRotDM, double r_step) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not create config file %s\n", filename);
        return -1;
    }

    fprintf(file, "# Configuration file for RNS\n");
    fprintf(file, "# WARNING: The parameters are not checked for consistency, be careful!\n");
    fprintf(file, "\n");
    fprintf(file, "# SETTINGS\n");
    fprintf(file, "id_file            \t%-10s\t# Path to initial data file, set to 0 or remove to ignore", "0");
    fprintf(file, "EoS_BM_type        \t%-10d\t# {0 -> tab, 1 -> poly}\n", eosBM->type);
    fprintf(file, "EoS_BM_file        \t%-10s\t# Required if EoS_type = tab\n", eosBM->filepath);
    fprintf(file, "EoS_BM_N           \t%-10.2f\t# P \\propto \\rho^(1 + 1 / N) Required if EoS_type = poly\n", eosBM->n_P);
    fprintf(file, "\n");
    fprintf(file, "EoS_DM_type        \t%-10d\t# {0 -> tab, 1 -> poly}\n", eosDM->type);
    fprintf(file, "EoS_DM_file        \t%-10s\t# Required if EoS_type = tab\n", eosDM->filepath);
    fprintf(file, "EoS_DM_N           \t%-10.2f\t# P \\propto \\rho^(1 + 1 / N) Required if EoS_type = poly\n", eosDM->n_P);
    fprintf(file, "accuracy           \t%-10.6g# Accuracy goal\n", accuracy);
    fprintf(file, "\n");
    fprintf(file, "# PARAMETERS\n");
    fprintf(file, "BM_central_energy  \t%-10.6g# Central energy density of the baryonic matter in 10^15 g / cm^3 \n", BM_prop->e_center);
    fprintf(file, "DM_central_energy  \t%-10.6g# Central energy density of the dark matter in 10^15 g / cm^3\n", DM_prop->e_center);
    fprintf(file, "\n");
    fprintf(file, "DM_particle_mass   \t%-10.6g# Mass of the DM particle in MeV \n", MDM);
    fprintf(file, "DM_fraction        \t%-10.6g# Target dark matter fraction \n", fdm_target);
    fprintf(file, "DM_fraction_tol    \t%-10.6g# Tolerance in determination of target dark matter fraction \n", fdm_tol);
    fprintf(file, "\n");
    fprintf(file, "BM_rotation_type   \t%-10d# Type of rotation for BM {0 -> Uniform, 1 -> J constant, 2 -> Uryu 8, 3 -> Uryu Extended}\n", BM_prop->RotType);
    fprintf(file, "BM_A               \t%-10.6g# Uryu parameter A\n", DiffRotBM->A);
    fprintf(file, "BM_B               \t%-10.6g# Uryu parameter B\n", DiffRotBM->B);
    fprintf(file, "BM_lambda1         \t%-10.6g# Omega_equator / Omega_central\n", DiffRotBM->lambda1);
    fprintf(file, "BM_lambda2         \t%-10.6g# Omega_max / Omega_central\n", DiffRotBM->lambda2);
    fprintf(file, "BM_p               \t%-10.6g# Uryu parameter p, Uryu Extended only\n", DiffRotBM->p);
    fprintf(file, "BM_csi             \t%-10.6g# Uryu parameter csi, Uryu Extended only\n", DiffRotBM->csi);
    fprintf(file, "\n");
    fprintf(file, "DM_rotation_type   \t%-10d# Type of rotation for DM {0 -> Uniform, 1 -> J constant, 2 -> Uryu 8, 3 -> Uryu Extended}\n", DM_prop->RotType);
    fprintf(file, "DM_A               \t%-10.6g# Uryu parameter A\n", DiffRotDM->A);
    fprintf(file, "DM_B               \t%-10.6g# Uryu parameter B\n", DiffRotDM->B);
    fprintf(file, "DM_lambda1         \t%-10.6g# Omega_max / Omega_central\n", DiffRotDM->lambda1);
    fprintf(file, "DM_lambda2         \t%-10.6g# Omega_equator / Omega_central\n", DiffRotDM->lambda2);
    fprintf(file, "DM_p               \t%-10.6g# Uryu parameter p, Uryu Extended only\n", DiffRotDM->p);
    fprintf(file, "DM_csi             \t%-10.6g# Uryu parameter csi, Uryu Extended only\n", DiffRotDM->csi);
    fprintf(file, "\n");
    fprintf(file, "counter_rotation   \t%-10d# Set wheter DM counter-rotates with respect to the BM  \n", counter ? 1 : 0);
    fprintf(file, "\n");
    fprintf(file, "r_ratio_BM         \t%-10.6g# Ratio of polar and equatorial radius for baryonic matter \n", r_ratio_BM);
    fprintf(file, "r_ratio_DM         \t%-10.6g# Ratio of polar and equatorial radius for dark matter \n", r_ratio_DM);
    fprintf(file, "r_step             \t%-10.6g# Step size for equatiorial ratios\n", r_step);
    fprintf(file, "\n");
    fprintf(file, "1Doutput           \t%-10d# Enable 1D output {0 -> no, 1 -> yes}\n", 0);
    fprintf(file, "2Doutput           \t%-10d# Enable 2D output {0 -> no, 1 -> yes}\n", 0);
    fprintf(file, "output_name        \t%-10s# Name of the h5 2D output. Only usefull if 2Doutput is 1\n", "auto");
    fprintf(file, "verbose            \t%-10d# Enable different levels of verbosity {0-4}\n", verbose);

    return 1;
}