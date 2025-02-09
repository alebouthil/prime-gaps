#include <stdio.h>
// #include <stdlib.h>
#include <gmp.h>
#include <mpi.h>
#include <stdbool.h>

// Function to check if a number is prime using GMP
bool is_prime(mpz_t n) {
    return mpz_probab_prime_p(n, 25) > 1; // Mpz returns 2 when it is certain n is a prime
}

unsigned long long mpz2ull(mpz_t z) // From stack overflow - C. K. Young https://stackoverflow.com/questions/6248723
{
    unsigned long long result = 0;
    mpz_export(&result, 0, -1, sizeof result, 0, 0, z);
    return result;
}

void process_chunk(unsigned long long *start, unsigned long long *end, unsigned long long local_data[5]) {
    // Processes a chunk, returns relavant prime information
    mpz_t curr, next; // Init mpz variables for chunk
    mpz_init(curr);
    mpz_init(next);
    mpz_set_ui(curr, *start);
    
    // Find and set first prime, we will process chunk from this number
    if (!is_prime(curr)) {
        mpz_nextprime(next, curr);
        local_data[3] = mpz2ull(next);
        mpz_set_ui(curr, mpz_get_ui(next)); // Sets curr to the next prime in sequence
        mpz_nextprime(next, curr);
    } else {
        local_data[3] = mpz2ull(curr);
        mpz_nextprime(next, curr);
    }

    while (mpz_get_ui(next) <= *end) {
        if ((mpz2ull(next) - mpz2ull(curr)) > local_data[0]) {
            // Larger gap found, set variables and continue
            local_data[0] = mpz2ull(next) - mpz2ull(curr);
            local_data[1] = mpz2ull(curr);
            local_data[2] = mpz2ull(next);
            mpz_set_ui(curr, mpz_get_ui(next)); // Start looking for gap from next prime
            mpz_nextprime(next,curr); // Set now so we can check for oob using the while condition
        } else {
            // Gap not larger, check next
            mpz_set_ui(curr, mpz_get_ui(next));
            mpz_nextprime(next,curr);
        }
    }

    // Set last prime for chunk, clear mpz variables
    local_data[4]= mpz2ull(curr);
    mpz_clear(curr);
    mpz_clear(next);
}

void master_recieve(unsigned long long global_data[], unsigned long long local_data[], int size) {
    // Master process collects data from all workers
    // Local_data refers to gap info collected by proc 0


    // Check if gap found by proc 0 is higher than previous max gap
    if (global_data[0] > local_data[0]) {
        global_data[0] = local_data[0];
        global_data[1] = local_data[1];
        global_data[2] = local_data[2];
    }
    // Check for gaps between chunks, initialises to last prime from proc 0
    unsigned long long previous_last_prime = local_data[4];

    for (int i = 1; i < size; ++i) {
        unsigned long long worker_data[5];
        MPI_Recv(worker_data, 5, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Check if the gap between chunks is larger than the current global max
        if (worker_data[3] - previous_last_prime > global_data[0]) {
            global_data[0] = worker_data[3] - previous_last_prime;
            global_data[1] = previous_last_prime;
            global_data[2] = worker_data[3];
        }

        // Check if the worker's local max gap is larger than the current global max
        if (worker_data[0] > global_data[0]) {
            global_data[0] = worker_data[0];
            global_data[1] = worker_data[1];
            global_data[2] = worker_data[2];
        }

        // Update the previous last prime for the next chunk
        previous_last_prime = worker_data[4];
    }
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Start timer
    double t1, t2;
    t1 = MPI_Wtime();

    // We use unsigned long long as we do not care for negative numbers
    const unsigned long long N = 1000000000; // Upper limit of search
    unsigned long long up_chunk_size = N / 3; // Size of the three upper chunks to be searched
    unsigned long long global_max_gap = 0;
    unsigned long long chunk_size = up_chunk_size / size;

    // Tell process where to start and stop searching from, based on process rank
    unsigned long long start = rank * chunk_size;
    unsigned long long end = (rank == size - 1) ? up_chunk_size : start + chunk_size;

    // We use an array like this to make sending and receiving more readable
    unsigned long long local_data[5];
    local_data[0] = 0; // local_gap;
    local_data[1] = 0; // lower_prime;
    local_data[2] = 0; // upper_prime;
    local_data[3] = 0; // first_prime;
    local_data[4] = 0; // last_prime;

    unsigned long long global_data[5];
    global_data[0] = 0; // Elements represent same as above, but global
    global_data[1] = 0; // We will write output from this array at end
    global_data[2] = 0;
    global_data[3] = 0;
    global_data[4] = 0;

    // First third
    process_chunk(&start, &end, local_data);

    // MPI communication, check gaps and store largest gap found so far
    if (rank != 0) { MPI_Send(local_data, 5, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);}
    else { master_recieve(global_data, local_data, size); }

    // Synch all processes before moving to next chunk
    // Unsure if needed as Send and Recieve are blocking calls?
    // Probably makes no difference in run time but feels safe to have
    // Faster processes getting to start early before slower ones have a chance to report their result has no effect I think,
    // as we have to wait for the slowest process anyway to output result
    MPI_Barrier(MPI_COMM_WORLD);

    // Set up for second third
    start += up_chunk_size;
    end += up_chunk_size;

    // Second third
    process_chunk(&start, &end , local_data);
    if (rank != 0) { MPI_Send(local_data, 5, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);}
    else { master_recieve(global_data, local_data, size); }
    MPI_Barrier(MPI_COMM_WORLD);

    // Third third
    start += up_chunk_size; 
    end += up_chunk_size;
    process_chunk(&start, &end , local_data);
    if (rank != 0) { MPI_Send(local_data, 5, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);}
    else { master_recieve(global_data, local_data, size); }
    MPI_Barrier(MPI_COMM_WORLD);

    // End timer
    t2 = MPI_Wtime();

    // Clean up variables and print result
    MPI_Finalize();
    printf("The largest gap between two consecutive primes in the range 0 to %llu is %llu\n", N, global_data[0]);
    printf("The primes flanking this gap are %llu and %llu\n", global_data[1], global_data[2]);
    printf("This was found in %f ticks using %i processes\n", t2-t1, size);
    printf("ticks: %f\n", t2-t1);

    return 0;
}
