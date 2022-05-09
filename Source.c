#include<stdio.h>
#include <time.h>
#include <mpi.h>
int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL); //start mpi


    //get rank of thread and number of threads there are
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    int numproc;
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    

    //vars for backwards sub
    float ratio;
    float values[10];
    float backwardsSum = 0.0;

    //begin clock for performance analysis
    double time_spent = 0.0;
    clock_t begin = clock();



    // augmented matrix in form A|B
    float augmentedMatrix[3][4] = {
        {2,1,1,10},
        {3,2,3,18},
        {1,4,9,16}
    };

    int sizeOfArray = (sizeof(augmentedMatrix) / sizeof(augmentedMatrix[0])); //order of matrix (amount of rows)
    int partition = sizeOfArray / numproc; //splitting the workload for MPI

    int j;
    int i;
    int k;

    //forward elim
    for (k = myRank * partition; k < myRank * partition + partition; k++) { //splits into partitions
        for (i = myRank * partition; i < myRank * partition + partition; i++) {
            if (i > k) { //for upper triangle only
                printf("my rank is: %d\n", myRank);
                ratio = augmentedMatrix[i][k] / augmentedMatrix[k][k]; //ratio for multiplier
                for (j = myRank*partition; j <= myRank*partition+partition; j++) { //does forward elim calcs with the ratio
                    augmentedMatrix[i][j] = augmentedMatrix[i][j] - ratio * augmentedMatrix[k][j];
                }

            }
        }
    }


    //backwards sub
    values[sizeOfArray - 1] = augmentedMatrix[sizeOfArray - 1][sizeOfArray] / augmentedMatrix[sizeOfArray - 1][sizeOfArray - 1];

    for (int i = sizeOfArray - 2; i >= 0; i--) {

        backwardsSum = 0;

        for (int j = i; j <= sizeOfArray - 1; j++) {
            backwardsSum = backwardsSum + augmentedMatrix[i][j] * values[j];
        }

        values[i] = (augmentedMatrix[i][sizeOfArray] - backwardsSum) / augmentedMatrix[i][i];
    }


    //prints out the unknown values that were solved
    printf("\nThe solution for unknowns are: \n");
    for (int i = 0; i <= sizeOfArray - 1; i++) {
        printf("\n%f\t", values[i]); 
    }

    //calculates time taken and prints it
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\n\nThe elapsed time is %f seconds", time_spent);


    // Finalize the MPI environment.
    MPI_Finalize();
    return(0);
}