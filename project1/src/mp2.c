// Add important libraries

#include <trexio.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Allocation of 4D arrays

double**** malloc_4d(size_t m, size_t n, size_t o, size_t p){
	double**** a = malloc(m*sizeof(double***));
	for (int i=0; i<m; i++){
		a[i] = malloc(n*sizeof(double**));
		for (int j=0; j<n; j++){
			a[i][j] = malloc(o*sizeof(double*));
			for (int k=0; k<o; k++){
				a[i][j][k] = malloc(p*sizeof(double));
			}
		}
	}

	return a;
}


// Free memory


void deallocate_4D(double ****array, size_t m, size_t n, size_t o) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < o; k++) {
                free(array[i][j][k]); // Free the 1D arrays
            }
            free(array[i][j]); // Free the 2D arrays
        }
        free(array[i]); // Free the 3D arrays
    }
    free(array); // Free the 4D array pointers
}


// Open the trexio file
trexio_t* open_file(trexio_exit_code rc, char* filename){
	trexio_t* trexio_file;
	

	trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);

	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error: %s \n", trexio_string_of_error(rc));
		exit(1);
	}

	return trexio_file;
}

// Close the trexio file

void close_file(trexio_exit_code rc, trexio_t* trexio_file){

	rc = trexio_close(trexio_file);

	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error: %s \n", trexio_string_of_error(rc));
		exit(1);
	}
	trexio_file = NULL;

}

// Calculate repulsion energy between nuclei (E_NN)

double get_repulsion_energy(trexio_exit_code rc, trexio_t* trexio_file){
	double repulsion_energy;

	rc = trexio_read_nucleus_repulsion(trexio_file, &repulsion_energy);
	
	// Check return code 
	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error reading nuclear repulsion energy: \n %s \n", trexio_string_of_error(rc));
		exit(1);
	}
	
	return repulsion_energy;
}


// Get number of occupied orbitals (RHF) == number of spin-up electrons

int get_Nocc(trexio_exit_code rc, trexio_t* trexio_file){
	int spin_up;

	rc = trexio_read_electron_up_num(trexio_file, &spin_up);

	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error reading spin-up electrons: \n %s \n", trexio_string_of_error(rc));
		exit(1);
	}

	return spin_up;
}

// Obtain number of canonical molecular orbitals for the given basis

int get_mo(trexio_exit_code rc, trexio_t* trexio_file){
	int mo;

	rc = trexio_read_mo_num(trexio_file, &mo);

	
	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error reading molecular orbitals: \n %s \n", trexio_string_of_error(rc));
		exit(1);
	}
	return mo;
}


// Obtain an array with all 1e integrals (core Hamiltonian) between mos


void get_1e_integrals(trexio_exit_code rc, trexio_t* trexio_file, double* data){
	
	rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);

	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error: %s \n", trexio_string_of_error(rc));
		exit(1);
	}
}

// Number of electron repulsion integrals 

int64_t get_n_integrals(trexio_exit_code rc, trexio_t* trexio_file){
	
	int64_t n_integrals;

	rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);

	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error: %s \n", trexio_string_of_error(rc));
		exit(1);
	}

	return n_integrals;
}

// Returns array of non-zero 4-center integrals and the corresponding orbital indexes (i,j,k,l)

void get_2e_integrals(trexio_exit_code rc, trexio_t* trexio_file,int64_t n_integrals, int32_t* index, double* value){
	int64_t offset_file = 0;
	int64_t buffer_size = n_integrals;

	rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, &buffer_size, index, value);

	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error: %s \n", trexio_string_of_error(rc));
		exit(1);
	}
}

void F_matrix(double**** F, int n_integrals, int32_t* index, double* value){
	for (int n=0; n<n_integrals; n++){		 
		int i = index[4*n+0];
		int j = index[4*n+1];
		int k = index[4*n+2];
		int l = index[4*n+3]; 
		
		// Need to account for all 8-fold permutations

		F[i][j][k][l] = value[n];
		F[i][l][k][j] = value[n];
		F[k][l][i][j] = value[n];
		F[k][j][i][l] = value[n];
		F[j][i][l][k] = value[n];
		F[l][i][j][k] = value[n];
		F[l][k][j][i] = value[n];
		F[j][k][l][i] = value[n];
		}
}



// Energy of molecular orbitals for MP2 calculation

void get_mo_energy(trexio_exit_code rc, trexio_t* trexio_file, double* mo_energy){
	
	rc = trexio_read_mo_energy(trexio_file, mo_energy);

	if (rc != TREXIO_SUCCESS){
		printf("TREXIO Error: %s \n", trexio_string_of_error(rc));
		exit(1);
	}

}



int main(int argc, char **argv){
	// Start by defining some variables
	
	trexio_exit_code rc;
	trexio_t* trexio_file; 
	char* filename; 
	double repulsion_energy;
	int spin_up;
	int Nocc;
	int mo;
	int64_t n_integrals;
	double core_energy = 0.0;
	double eri = 0.0;
	double hf_energy;
	double mp2_energy = 0.0;

	// First check that TREXIO file is specified in the input

	if (argc != 2){
		printf("Input error. Correct usage as: ./mp2 [file.h5] \n");
	}else{
		filename = argv[1];
	}



	// Open the file
	trexio_file = open_file(rc,filename);

	
	// Read nuclear repulsion energy

	repulsion_energy = get_repulsion_energy(rc,trexio_file);
	printf("Nuclear repulsion energy = %f a.u. \n", repulsion_energy);

	
	// Read number of up-spin electrons == number of occupied MOs

	Nocc = get_Nocc(rc,trexio_file);	
	printf("Number of spin-up electrons: %d \n", Nocc);


	// Read number of MOs
		
	mo = get_mo(rc,trexio_file);
	printf("Number of molecular orbitals: %d \n", mo);



	// Obtain core hamiltomian

	double* h = malloc(mo*mo*sizeof(double));

	get_1e_integrals(rc,trexio_file,h);


	// Obtain number of non-zero electron repulsion integrals
	
	n_integrals = get_n_integrals(rc, trexio_file);
	
	printf("Number of non-zero 2-electron integrals: %ld \n", n_integrals);
	

	// Allocate indices and values
	
	int32_t* const index = malloc(4*n_integrals*sizeof(int32_t));
	if (index == NULL){
		fprintf(stderr, "Malloc failed for index");
		exit(1);
	}
	
	double* const value = malloc(n_integrals*sizeof(double));
	if (value == NULL){
		fprintf(stderr, "Malloc failed for value");
		exit(1);
	}

	
	// Obtain electron repulsion integrals	

	
	get_2e_integrals(rc,trexio_file,n_integrals,index,value);

	
	// Generate a 4D array in which every dimension corresponds to an index and matrix elements are the integral values
	
	double**** F = malloc_4d(mo,mo,mo,mo);
	for (int i=0; i<mo;i++){
		for (int j=0; j<mo;j++){
			for (int k=0; k<mo;k++){
				for (int l=0; l<mo;l++){
					F[i][j][k][l] = 0.0;
				}
			}
		}
	}

	F_matrix(F,n_integrals,index,value);

	// Loop for core Hamiltonian energies (only diagonal elements <i|h|i>)

	for (int i=0; i<Nocc; i++){
		core_energy += h[i+i*mo];
	}

	// Loop over occupied orbitals for electron repulsion integrals (eri). IN this case, only two indexes are needed

	for (int i=0; i<Nocc; i++){
		for (int j=0; j<Nocc; j++){
			eri += 2*F[i][j][i][j];
			eri -= F[i][j][j][i];
		}
	}
	
	// Collect energies and use the formula

	hf_energy = repulsion_energy + 2*core_energy + eri;

	printf("Hartree-Fock energy = %lf a.u. \n", hf_energy);

	// Obtain orbital energies
	
	double* mo_energy = malloc(mo*sizeof(double));

	get_mo_energy(rc,trexio_file,mo_energy);

	// Calculate the MP2 energy correction 

	for (int i=0; i<Nocc; i++){
		for (int j=0; j<Nocc; j++){
			for (int a=Nocc; a<mo; a++){
				for (int b=Nocc; b<mo; b++){
					mp2_energy += (F[i][j][a][b]*(2*F[i][j][a][b] - F[i][j][b][a])) / (mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b]);
				}
			}
		}
	}


	printf("MP2 correction energy = %lf a.u. \n", mp2_energy);

	printf("MP2 total energy = %lf a.u. \n", hf_energy+mp2_energy);


	// Free allocated memory 
	
	free(h);
	free(index);
	free(value);
	deallocate_4D(F,mo,mo,mo);
	free(mo_energy);

	// Close the trexio file

	close_file(rc,trexio_file);
	
	return 0;
}


