#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpi.h>
#include <stdbool.h>

// Function to check if a number is prime using GMP
bool is_prime(mpz_t n) {
    return mpz_probab_prime_p(n, 25) > 1; // Mpz returns 2 when it is certain n is a prime
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double t1, t2;
    t1 = MPI_Wtime();

    const unsigned long long N = 1000000000; // Upper limit
    unsigned long long local_max_gap = 0;
    unsigned long long global_max_gap = 0;
    unsigned long long previous_prime = 2;
    unsigned long long first_prime_in_chunk = 0;
    unsigned long long last_prime_in_chunk = 0;
    unsigned long long gap_prime1 = 0, gap_prime2 = 0;

    // Divide the range among processes
    unsigned long long chunk_size = (N - 2) / size;
    unsigned long long start = 2 + rank * chunk_size;
    unsigned long long end = (rank == size - 1) ? N : start + chunk_size;

    // GMP variables
    mpz_t n;
    mpz_init(n);

    // Find primes in the local range and calculate the largest gap
    for (unsigned long long i = start; i < end; ++i) {
        mpz_set_ui(n, i); // Set GMP variable to the current number
        
        // Check if current number is prime
        if (is_prime(n)) {
            if (first_prime_in_chunk == 0) {
                first_prime_in_chunk = i; // Record the first prime in the chunk
                previous_prime = i; // First prime to start tracking gaps from
            }
            last_prime_in_chunk = i; // Update the last prime in the chunk

            // Find gap, compare to largest local gap and replace if larger
            unsigned long long gap = i - previous_prime;
            if (gap > local_max_gap) {
                local_max_gap = gap;
                gap_prime1 = previous_prime;
                gap_prime2 = i;
            }
            previous_prime = i; // Start looking for next gap
        }
    }

    // Worker processes send their results to the master
    if (rank != 0) {
        // Pack data into a buffer for easy reading
        unsigned long long data[5];
        data[0] = local_max_gap;
        data[1] = gap_prime1;
        data[2] = gap_prime2;
        data[3] = first_prime_in_chunk;
        data[4] = last_prime_in_chunk;

        // Send the data to the master process
        MPI_Send(data, 5, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);
    } else {
        // Master process collects data from all workers
        
        
        global_max_gap = local_max_gap; // Global max starts as the local max calculated by master
        unsigned long long master_gap_prime1 = gap_prime1;
        unsigned long long master_gap_prime2 = gap_prime2;

        // Check for gaps between chunks
        unsigned long long previous_last_prime = last_prime_in_chunk;

        for (int i = 1; i < size; ++i) {
            unsigned long long worker_data[5];
            MPI_Recv(worker_data, 5, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Check if the gap between chunks is larger than the current global max
            if (worker_data[3] - previous_last_prime > global_max_gap) {
                global_max_gap = worker_data[3] - previous_last_prime;
                master_gap_prime1 = previous_last_prime;
                master_gap_prime2 = worker_data[3];
            }

            // Check if the worker's local max gap is larger than the current global max
            if (worker_data[0] > global_max_gap) {
                global_max_gap = worker_data[0];
                master_gap_prime1 = worker_data[1];
                master_gap_prime2 = worker_data[2];
            }

            // Update the previous last prime for the next chunk
            previous_last_prime = worker_data[4];
        }

        // End timer
        t2 = MPI_Wtime();

        // Print the result
        if (rank == 0) {
            printf("The largest gap between two consecutive primes in the range 0 to %llu is %llu\n", N, global_max_gap);
            printf("The primes flanking this gap are %llu and %llu\n", master_gap_prime1, master_gap_prime2);
            printf("This was found in %f ticks using %i processes\n", t2-t1, size);
            printf("ticks: %f\n", t2-t1);
        }
    }

    // Clean up variables
    mpz_clear(n);
    MPI_Finalize();
    return 0;
}
