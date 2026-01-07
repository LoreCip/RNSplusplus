#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef H5_USE_18_API
#define H5_USE_18_API
#endif
#define H5Acreate_vers 2
#define H5Dcreate_vers 2
#define H5Dopen_vers 2

#include "consts.h"
#include "equil.h"
#include "output.h"

void open_output(char filepath[], FILE **file) {
    *file = fopen(filepath, "w");
    if (*file == NULL) {
        printf("File not opened.\n");
        exit(1);
    }

    fprintf(*file,
            "# %16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\n",
            "BM_e_center", "R_e_BM (km)", "R_p_BM (km)", "Mass_BM (M_sun)", "Mass_0_BM (M_sun)", "DM_e_center", "R_e_DM (km)", "R_p_DM (km)", "Mass_DM (M_sun)", "Mass_0_DM (M_sun)", "Mtot (M_sun)", "fdm", "Omega_BM (Hz)", "Omega_DM (Hz)", "J_BM/M_BM^2", "J_DM/M_DM^2", "r_ratio_BM", "r_ratio_DM", "Omega_K_BM (Hz)", "Omega_K_DM (Hz)", "T_BM", "T_DM", "W_BM", "W_DM", "beta_BM", "beta_DM");
}

void close_output(FILE *file) {
    fclose(file);
}

void print_output(FILE *file, struct stellar_properties *star_props, double R_e_BM, double R_p_BM, double Mass_BM, double Mass_0_BM, struct stellar_properties *DM_props, double R_e_DM, double R_p_DM, double Mass_DM, double Mass_0_DM, double fdm, double Omega_BM, double Omega_DM, double J_BM, double J_DM, double r_ratio_BM, double r_ratio_DM, double Omega_K_BM, double Omega_K_DM, double T_BM, double T_DM, double W_BM, double W_DM) {
    fprintf(file, "%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\t%16.8g\n",
            star_props->e_center, R_e_BM / 1.0e5, R_p_BM / 1e5, Mass_BM / MSUN, Mass_0_BM / MSUN, DM_props->e_center, R_e_DM / 1.0e5, R_p_DM / 1e5, Mass_DM / MSUN, Mass_0_DM / MSUN, Mass_BM / MSUN + Mass_DM / MSUN, fdm, Omega_BM, Omega_DM, C * J_BM / (G * SQ(Mass_BM)), C * J_DM / (G * SQ(Mass_DM) + DBL_EPSILON), r_ratio_BM, r_ratio_DM, Omega_K_BM, Omega_K_DM, T_BM, T_DM, W_BM, W_DM, T_BM / W_BM, T_DM / W_DM);
    fflush(file);
}

char *getDiffRotTypeString(enum DiffRotType *type) {
    switch (*type) {
        case UNIFORM:
            return "Uniform";
        case DIFF:
            return "J Constant";
        case U8:
            return "Uryu 8";
        case U9EXT:
            return "Uryu 9 Extended";
        default:
            return "UNKNOWN_ROT";
    }
}
int getDiffRotTypeFromString(const char *str) {
    if (strcmp(str, "Uniform") == 0)
        return UNIFORM;
    else if (strcmp(str, "J Constant") == 0)
        return DIFF;
    else if (strcmp(str, "Uryu 8") == 0)
        return U8;
    else if (strcmp(str, "Uryu 9 Extended") == 0)
        return U9EXT;
    else
        return UNKNOWN_ROT;
}

void hdf5_save_var(double(*s_gp), int sdiv, double(*mu), int mdiv, char filename[], struct evolution_variables *ev, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, double *r_ratio_BM, double *r_ratio_DM, double *Omega_BM, double *Omega_DM, double *r_e_BM, double *r_e_DM, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM, char *input_parfile) {
    hid_t file_id; /* file identifier */
    herr_t status;

    hid_t attribute_id; /* identifiers for attributes*/
    hid_t attributeH5type, attributeH5c128;

    hid_t dataset_id, dataspace_id; /* identifiers for dsets*/
    hsize_t dims[3];

    hid_t dcpl; /* Compression parameter */

    int i, j;
    double *dset_data;
    double **var;
    int varIndex;
    /* =========================================== */
    /* Create a new file using default properties. */
    /* =========================================== */
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* ============================================= */
    /*   Create and save the data attributes         */
    /* ============================================= */
    int attrblistLEN = 32;
    struct
    {
        char name[128];
        hid_t H5type;
        void *data;
    } attrblist[] = {
        {"BM poly Gamma", H5T_NATIVE_DOUBLE, &eosBM->Gamma_P},
        {"BM n_P", H5T_NATIVE_DOUBLE, &eosBM->n_P},
        {"axis_ratio_BM", H5T_NATIVE_DOUBLE, r_ratio_BM},
        {"axis_ratio_DM", H5T_NATIVE_DOUBLE, r_ratio_DM},
        {"re_BM", H5T_NATIVE_DOUBLE, r_e_BM},
        {"re_DM", H5T_NATIVE_DOUBLE, r_e_DM},
        {"BM ec", H5T_NATIVE_DOUBLE, &star_props->e_center},
        {"BM Omega", H5T_NATIVE_DOUBLE, Omega_BM},
        {"SDIV", H5T_NATIVE_INT, &sdiv},
        {"MDIV", H5T_NATIVE_INT, &mdiv},
        {"BM eos_type", H5T_NATIVE_INT, &eosBM->type},
        {"BM eos_file", H5T_C_S1, &eosBM->filepath},
        {"DM poly Gamma", H5T_NATIVE_DOUBLE, &eosDM->Gamma_P},
        {"DM n_P", H5T_NATIVE_DOUBLE, &eosDM->n_P},
        {"DM ec", H5T_NATIVE_DOUBLE, &DM_props->e_center},
        {"DM Omega", H5T_NATIVE_DOUBLE, Omega_DM},
        {"DM eos_type", H5T_NATIVE_INT, &eosDM->type},
        {"DM eos_file", H5T_C_S1, &eosDM->filepath},
        {"BM RotType", H5T_C_S1, getDiffRotTypeString(&star_props->RotType)},
        {"DM RotType", H5T_C_S1, getDiffRotTypeString(&DM_props->RotType)},
        {"BM A", H5T_NATIVE_DOUBLE, &DiffRotBM->A},
        {"BM B", H5T_NATIVE_DOUBLE, &DiffRotBM->B},
        {"BM lambda_1", H5T_NATIVE_DOUBLE, &DiffRotBM->lambda1},
        {"BM lambda_2", H5T_NATIVE_DOUBLE, &DiffRotBM->lambda2},
        {"BM p", H5T_NATIVE_DOUBLE, &DiffRotBM->p},
        {"BM csi", H5T_NATIVE_DOUBLE, &DiffRotBM->csi},
        {"DM A", H5T_NATIVE_DOUBLE, &DiffRotDM->A},
        {"DM B", H5T_NATIVE_DOUBLE, &DiffRotDM->B},
        {"DM lambda_1", H5T_NATIVE_DOUBLE, &DiffRotDM->lambda1},
        {"DM lambda_2", H5T_NATIVE_DOUBLE, &DiffRotDM->lambda2},
        {"DM p", H5T_NATIVE_DOUBLE, &DiffRotDM->p},
        {"DM csi", H5T_NATIVE_DOUBLE, &DiffRotDM->csi},
    };

    /* ============================================= */
    /*  Read the data set for variables and save it  */
    /* ============================================= */
    int varlistLEN = 14;
    struct
    {
        char name[128];
        double **data;
        char description[128];
    } varlist[] = {
        {"/rho_potential", ev->rho, "Values for RNSID variable rho_potential"},
        {"/gama", ev->gama, "Values for RNSID variable gama"},
        {"/alpha", ev->alpha, "Values for RNSID variable alpha"},
        {"/omega", ev->omega, "Values for RNSID variable omega"},
        {"/energy_BM", ev->energyBM, "Values for RNSID variable energyBM"},
        {"/pressure_BM", ev->pressureBM, "Values for RNSID variable pressureBM"},
        {"/enthalpy_BM", ev->enthalpyBM, "Values for RNSID variable enthalpyBM"},
        {"/velocity_sq_BM", ev->velocity_sqBM, "Values for RNSID variable velocity_sqBM"},
        {"/energy_DM", ev->energyDM, "Values for RNSID variable energyDM"},
        {"/pressure_DM", ev->pressureDM, "Values for RNSID variable pressureDM"},
        {"/enthalpy_DM", ev->enthalpyDM, "Values for RNSID variable enthalpyDM"},
        {"/velocity_sq_DM", ev->velocity_sqDM, "Values for RNSID variable velocity_sqDM"},
        {"/Omega_diffBM", ev->Omega_diffBM, "Values for RNSID variable Omega_diffBM"},
        {"/Omega_diffDM", ev->Omega_diffDM, "Values for RNSID variable Omega_diffDM"}};

    /* ============================================= */
    /*   Create and save the data attributes         */
    /* ============================================= */
    for (varIndex = 0; varIndex < attrblistLEN; varIndex++) {
        if (attrblist[varIndex].H5type == H5T_C_S1) {
            attributeH5type = H5Tcopy(H5T_C_S1);
            H5Tset_size(attributeH5type, strlen(attrblist[varIndex].data) + 1);
            dims[0] = 1;
        } else {
            dims[0] = 1;
            attributeH5type = attrblist[varIndex].H5type;
        }
        dataspace_id = H5Screate_simple(1, dims, NULL);
        attribute_id = H5Acreate(file_id,
                                 attrblist[varIndex].name, attributeH5type,
                                 dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite(attribute_id, attributeH5type, attrblist[varIndex].data);
        status = H5Aclose(attribute_id);
        status = H5Sclose(dataspace_id);

        if (attrblist[varIndex].H5type == H5T_C_S1) {
            status = H5Tclose(attributeH5type);
        }
    }
    /* ============================================= */
    /*   Create and save the data set GRID points    */
    /* ============================================= */
    double s_tosave[SDIV], mu_tosave[MDIV];
    for (int s = 0; s < SDIV; s++) {
        s_tosave[s] = s_gp[s + 1];
    }
    for (int m = 0; m < MDIV; m++) {
        mu_tosave[m] = mu[m + 1];
    }

    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_shuffle(dcpl);
    H5Pset_deflate(dcpl, GZIP_COMPRESSION_LEVEL);

    dims[0] = (sdiv);
    dataspace_id = H5Screate_simple(1, dims, NULL);
    H5Pset_chunk(dcpl, 1, dims);
    dataset_id = H5Dcreate(file_id, "/s_qp",
                           H5T_NATIVE_DOUBLE, dataspace_id,
                           H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, s_tosave);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);

    dims[0] = (mdiv);
    dataspace_id = H5Screate_simple(1, dims, NULL);
    H5Pset_chunk(dcpl, 1, dims);
    dataset_id = H5Dcreate(file_id, "/mu",
                           H5T_NATIVE_DOUBLE, dataspace_id,
                           H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mu_tosave);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);

    /* ============================================= */
    /*   Create and save the data set for variables  */
    /* ============================================= */
    attributeH5c128 = H5Tcopy(H5T_C_S1);
    H5Tset_size(attributeH5c128, 128);
    for (varIndex = 0; varIndex < varlistLEN; varIndex++) {
        dset_data = malloc(sdiv * mdiv * sizeof(double));
        dims[0] = sdiv;
        dims[1] = mdiv;
        var = varlist[varIndex].data;

        for (i = 0; i < sdiv; i++) {
            for (j = 0; j < mdiv; j++) {
                dset_data[i * mdiv + j] = var[i + 1][j + 1];
            }
        }

        dataspace_id = H5Screate_simple(2, dims, NULL);
        hsize_t chunk_dims[2] = {200, 200};
        H5Pset_chunk(dcpl, 2, chunk_dims);
        dataset_id = H5Dcreate(file_id, varlist[varIndex].name,
                               H5T_NATIVE_DOUBLE, dataspace_id,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
        status = H5Sclose(dataspace_id);
        status = H5Dclose(dataset_id);
        free(dset_data);
        /* ========================================================= */
        /* ---- Now add the attributes to the Just created data sets */
        /* ========================================================= */
        dataset_id = H5Dopen(file_id, varlist[varIndex].name, H5P_DEFAULT);
        dims[0] = 1;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        attribute_id = H5Acreate(dataset_id, "Description", attributeH5c128,
                                 dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite(attribute_id, attributeH5c128, varlist[varIndex].description);
        status = H5Aclose(attribute_id);
        status = H5Dclose(dataset_id);
    }

    status = H5Pclose(dcpl);
    status = H5Tclose(attributeH5c128);

    // === Save the parfile buffer as a dataset '/input_parfile_contents' ===
    FILE *fp = fopen(input_parfile, "rb");
    if (!fp) {
        fprintf(stderr, "Error: Could not open parameter file %s\n", input_parfile);
        return;
    }

    fseek(fp, 0, SEEK_END);
    size_t parfile_size = ftell(fp);
    rewind(fp);

    // Allocate buffer
    char *parfile_buf = malloc(parfile_size);
    if (!parfile_buf) {
        fprintf(stderr, "Error: Failed to allocate memory for parameter file buffer\n");
        fclose(fp);
        return;
    }

    // Read entire file into buffer
    size_t read_bytes = fread(parfile_buf, 1, parfile_size, fp);
    if (read_bytes != parfile_size) {
        fprintf(stderr, "Error: Failed to read entire parameter file\n");
        free(parfile_buf);
        fclose(fp);
        return;
    }
    fclose(fp);

    hsize_t dims_par[1] = {parfile_size};
    dataspace_id = H5Screate_simple(1, dims_par, NULL);

    // Fixed length string datatype for the entire buffer
    hid_t strdatatype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(strdatatype, parfile_size);
    if (status < 0) {
        fprintf(stderr, "Error setting string datatype size\n");
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        return;
    }

    dataset_id = H5Dcreate(file_id, "/parfile", strdatatype, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error creating dataset /parfile\n");
        H5Tclose(strdatatype);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        return;
    }

    status = H5Dwrite(dataset_id, strdatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, parfile_buf);
    if (status < 0) {
        fprintf(stderr, "Error writing dataset /parfile\n");
    }

    // Cleanup
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(strdatatype);
    free(parfile_buf);

    /* Terminate access to the file. */
    status = H5Fclose(file_id);
}

void hdf5_read_var(double *s_gp, int *sdiv, double *mu, int *mdiv, char filename[], struct evolution_variables *ev, struct stellar_properties *star_props, struct stellar_properties *DM_props, struct EOS *eosBM, struct EOS *eosDM, double *r_ratio_BM, double *r_ratio_DM, double *Omega_BM, double *Omega_DM, double *r_e_BM, double *r_e_DM, struct DiffRotParams *DiffRotBM, struct DiffRotParams *DiffRotDM) {
    hid_t file_id; /* file identifier */
    herr_t status;

    hid_t attribute_id; /* identifiers for attributes*/
    hid_t attributeH5type, attributeH5c128;
    hid_t datasetH5type;
    hid_t dataset_id, dataspace_id; /* identifiers for dsets*/
    hsize_t dims[3];

    int i, j;
    double *dset_data;
    double **var;
    int varIndex;
    /* =========================================== */
    /* Create a new file using default properties. */
    /* =========================================== */
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    /* ============================================= */
    /*   Create and save the data attributes         */
    /* ============================================= */
    char BM_RotType_str[128];
    char DM_RotType_str[128];
    int attrblistLEN = 32;
    struct
    {
        char name[128];
        hid_t H5type;
        void *data;
    } attrblist[] = {
        {"BM poly Gamma", H5T_NATIVE_DOUBLE, &eosBM->Gamma_P},
        {"BM n_P", H5T_NATIVE_DOUBLE, &eosBM->n_P},
        {"axis_ratio_BM", H5T_NATIVE_DOUBLE, r_ratio_BM},
        {"axis_ratio_DM", H5T_NATIVE_DOUBLE, r_ratio_DM},
        {"re_BM", H5T_NATIVE_DOUBLE, r_e_BM},
        {"re_DM", H5T_NATIVE_DOUBLE, r_e_DM},
        {"BM ec", H5T_NATIVE_DOUBLE, &star_props->e_center},
        {"BM Omega", H5T_NATIVE_DOUBLE, Omega_BM},
        {"SDIV", H5T_NATIVE_INT, sdiv},
        {"MDIV", H5T_NATIVE_INT, mdiv},
        {"BM eos_type", H5T_NATIVE_INT, &eosBM->type},
        {"BM eos_file", H5T_C_S1, &eosBM->filepath},
        {"DM poly Gamma", H5T_NATIVE_DOUBLE, &eosDM->Gamma_P},
        {"DM n_P", H5T_NATIVE_DOUBLE, &eosDM->n_P},
        {"DM ec", H5T_NATIVE_DOUBLE, &DM_props->e_center},
        {"DM Omega", H5T_NATIVE_DOUBLE, Omega_DM},
        {"DM eos_type", H5T_NATIVE_INT, &eosDM->type},
        {"DM eos_file", H5T_C_S1, &eosDM->filepath},
        {"BM RotType", H5T_C_S1, BM_RotType_str},
        {"DM RotType", H5T_C_S1, DM_RotType_str},
        {"BM A", H5T_NATIVE_DOUBLE, &DiffRotBM->A},
        {"BM B", H5T_NATIVE_DOUBLE, &DiffRotBM->B},
        {"BM lambda_1", H5T_NATIVE_DOUBLE, &DiffRotBM->lambda1},
        {"BM lambda_2", H5T_NATIVE_DOUBLE, &DiffRotBM->lambda2},
        {"BM p", H5T_NATIVE_DOUBLE, &DiffRotBM->p},
        {"BM csi", H5T_NATIVE_DOUBLE, &DiffRotBM->csi},
        {"DM A", H5T_NATIVE_DOUBLE, &DiffRotDM->A},
        {"DM B", H5T_NATIVE_DOUBLE, &DiffRotDM->B},
        {"DM lambda_1", H5T_NATIVE_DOUBLE, &DiffRotDM->lambda1},
        {"DM lambda_2", H5T_NATIVE_DOUBLE, &DiffRotDM->lambda2},
        {"DM p", H5T_NATIVE_DOUBLE, &DiffRotDM->p},
        {"DM csi", H5T_NATIVE_DOUBLE, &DiffRotDM->csi},
    };

    /* ============================================= */
    /*  Read the data set for variables and save it  */
    /* ============================================= */
    int varlistLEN = 14;
    struct
    {
        char name[128];
        double **data;
        char description[128];
    } varlist[] = {
        {"/rho_potential", ev->rho, "Values for RNSID variable rho_potential"},
        {"/gama", ev->gama, "Values for RNSID variable gama         "},
        {"/alpha", ev->alpha, "Values for RNSID variable alpha        "},
        {"/omega", ev->omega, "Values for RNSID variable omega        "},
        {"/energy_BM", ev->energyBM, "Values for RNSID variable energy       "},
        {"/pressure_BM", ev->pressureBM, "Values for RNSID variable pressure     "},
        {"/enthalpy_BM", ev->enthalpyBM, "Values for RNSID variable enthalpy     "},
        {"/velocity_sq_BM", ev->velocity_sqBM, "Values for RNSID variable velocity_sq  "},
        {"/energy_DM", ev->energyDM, "Values for RNSID variable energy       "},
        {"/pressure_DM", ev->pressureDM, "Values for RNSID variable pressure     "},
        {"/enthalpy_DM", ev->enthalpyDM, "Values for RNSID variable enthalpy     "},
        {"/velocity_sq_DM", ev->velocity_sqDM, "Values for RNSID variable velocity_sq  "},
        {"/Omega_diffBM", ev->Omega_diffBM, "Values for RNSID variable Omega_diffBM   "},
        {"/Omega_diffDM", ev->Omega_diffDM, "Values for RNSID variable Omega_diffDM   "}};

    /* ============================================= */
    /*   Read the data attributes         */
    /* ============================================= */
    for (varIndex = 0; varIndex < attrblistLEN; varIndex++) {
        if (attrblist[varIndex].H5type == H5T_C_S1) {
            attributeH5type = H5Tcopy(H5T_C_S1);
            H5Tset_size(attributeH5type, strlen(attrblist[varIndex].data) + 1);
            dims[0] = 1;
        } else {
            dims[0] = 1;
            attributeH5type = attrblist[varIndex].H5type;
        }

        dataspace_id = H5Screate_simple(1, dims, NULL);
        attribute_id = H5Aopen(file_id, attrblist[varIndex].name, H5P_DEFAULT);
        attributeH5type = H5Aget_type(attribute_id);
        attributeH5type = H5Tget_native_type(attributeH5type, H5T_DIR_DEFAULT);

        status = H5Aread(attribute_id, attributeH5type, attrblist[varIndex].data);
        status = H5Aclose(attribute_id);

        /// status = H5Sclose(dataspace_id);
        if (attrblist[varIndex].H5type == H5T_C_S1) {
            status = H5Tclose(attributeH5type);
        }
    }
    star_props->RotType = getDiffRotTypeFromString(BM_RotType_str);
    DM_props->RotType = getDiffRotTypeFromString(DM_RotType_str);

    /* ============================================= */
    /*  Read the data set for GRID points            */
    /* ============================================= */
    double s_toread[*sdiv], mu_toread[*mdiv];

    dataset_id = H5Dopen(file_id, "/s_qp", H5P_DEFAULT);

    datasetH5type = H5Dget_type(dataset_id);
    datasetH5type = H5Tget_native_type(datasetH5type, H5T_DIR_DEFAULT);
    status = H5Dread(dataset_id, datasetH5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, s_toread);
    status = H5Dclose(dataset_id);

    dataset_id = H5Dopen(file_id, "/mu", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mu_toread);
    status = H5Dclose(dataset_id);

    for (int s = 1; s <= *sdiv; s++) {
        s_gp[s] = s_toread[s - 1];
    }
    for (int m = 1; m <= *mdiv; m++) {
        mu[m] = mu_toread[m - 1];
    }


    /* ============================================= */
    /*  Read the data set for variables              */
    /* ============================================= */
    for (varIndex = 0; varIndex < varlistLEN; varIndex++) {
        dset_data = malloc(*sdiv * *mdiv * sizeof(double));
        dims[0] = *sdiv;
        dims[1] = *mdiv;
        dataset_id = H5Dopen(file_id, varlist[varIndex].name, H5P_DEFAULT);
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);

        var = varlist[varIndex].data;
        for (i = 0; i < *sdiv; i++) {
            for (j = 0; j < *mdiv; j++) {
                var[i + 1][j + 1] = dset_data[i * (*mdiv) + j];
            }
        }

        status = H5Dclose(dataset_id);
        free(dset_data);
    }

    /* Terminate access to the file. */
    status = H5Fclose(file_id);
}

void write_profile(double **variable, const char *filename) {
    if (access(filename, F_OK) == 0) {
        if (remove(filename) != 0) {
            fprintf(stderr, "Error deleting file: %s\n", filename);
            return;
        }
        printf("Deleted existing file: %s\n", filename);
    }
    FILE *fp;
    fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file\n");
        return;
    }
    fprintf(fp, "Radius\tVariable\n");
    for (int i = 1; i <= SDIV; i++) {
        fprintf(fp, "%.9f\t%.9f\n", SMAX * (i - 1.0) / (SDIV - 1.0), variable[i][1]);
    }
    fclose(fp);
}