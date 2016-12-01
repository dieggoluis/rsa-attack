#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <ctime>

#include "batchGCD.hpp"
#include "mpi.h"

using namespace std;

//character that divide number in a string: ex: 123#23#455, numbers 123, 23, 455
//useful to pass saveral numbers as char* buffer for MPI functions 
const char divisor = '#';

//root will processor with taskid = 0
const int root = 0;

//convert vector of mpz_class to char* (ex: 123, 23, 455 would be {123#23#455})
char* parseToChar (vector<mpz_class>& v, int* sizeInput) {
    string aux;
    for (int i = 0; i < v.size(); i++) {
        aux += v[i].get_str() + divisor;
    }
    *sizeInput = aux.size();
    char * buff = new char[aux.size()+1];
    strcpy (buff, aux.c_str());
    return buff;
}

//convert char* to vector of mpz_class (ex: 123#23#455 would be {123, 23, 455})
vector<mpz_class> parseToMpz (char * keys, int size) {
    vector<mpz_class> localMpzKeys;
    for (int i = 0; i < size; i++) {
        int j = i;
        string s;
        while (j < size && keys[j] != divisor) {
            s += keys[j];
            j++;
        }
        mpz_class x(s);
        localMpzKeys.push_back(x);
        i = j;
    }
    return localMpzKeys;
}

//read input
vector<mpz_class> readInput (const char *filename) {
    vector<mpz_class> input;
    ifstream file(filename);
    int index;
    while (file >> index) {
        mpz_class x;
        file >> x;
        input.push_back(x);
    }
    return input;
}

//print vector for debug
void printVector(int * v, int n) {
    for (int i = 0; i<n; i++)
        cout << v[i] << " ";
    cout << endl;
}

//main function
//we distributed the algorithm in the main function
int main (int argc, char* argv[]) {
    // number of processors and the id
    int numtasks, taskid;

    //MPI settings
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    char * keys; //buffer with the keys
    int sizeInput = 0; //size of the buffer
    int numberKeys = 0; //number of the keys

    int * sendcounts = 0; //size sent to each processor
    int * displs = 0 ; //displacements (see the MPI_Scatter documentation)

    //if root read input
    vector<mpz_class> input;
    if (taskid == root) {
        input = readInput(argv[1]);
        numberKeys = input.size();
        keys = parseToChar(input, &sizeInput);

        //divide tasks among processor
        int q = numberKeys / numtasks;
        sendcounts = new int[numtasks];
        displs = new int[numtasks];

        int bef = 0;
        int t = 0;
        for (int i = 0; i < numtasks - 1; i++) {
            displs[i] = bef;
            int j = 0;
            while (j < q) {
                if (keys[t] == divisor)
                    j++;
                t++;
            }
            sendcounts[i] = t - bef;
            bef = t;
        }

        displs[numtasks-1] = bef;
        int r = numberKeys - (numtasks - 1) * q;
        int j = 0;
        while (j < r) {
            if (keys[t] == divisor)
                j++;
            t++;
        }
        sendcounts[numtasks-1] = t - bef;
    }

    clock_t begin;
    if (taskid == root) begin = clock();

    //send size of each buffer
    int localSize;
    MPI_Scatter(sendcounts, 1, MPI_INT, &localSize, 1, MPI_INT, root, MPI_COMM_WORLD);

    //local buffers
    char * localKeys;
    localKeys = new char[localSize+1];

    //root (master) send to each processor its local keys
    MPI_Scatterv(keys, sendcounts, displs, MPI_CHAR, localKeys, localSize, MPI_CHAR, root, MPI_COMM_WORLD);

    //parse to mpz vector
    vector<mpz_class> localMpzKeys = parseToMpz(localKeys, localSize);
    //local solution
    BatchGCD * localSolution = new BatchGCD(localMpzKeys);      

    //get root value (local tree multiplication. If the local entries are N1, N2, N3 the root would be N1 * N2 * N3)
    string s = (localSolution->getValueRoot()).get_str();
    s += divisor; //add divisor

    //array with the size of each number (number of digits + divisor)
    int * arraySizeSequential; 
    if (taskid == root) arraySizeSequential = new int[numtasks];

    //each processor send to root
    int sizeSequential = s.size();
    MPI_Gather(&sizeSequential, 1, MPI_INT, arraySizeSequential, 1, MPI_INT, root, MPI_COMM_WORLD);

    //buffer used by root to store the intermediate result
    char * buff = new char[sizeSequential + 1];
    strcpy (buff, s.c_str());

    char * arraySequential; //array of the intermediates products 
    int * displsSequential; //displacements (see MPI_Gatherv documentation)
    if (taskid == root) { //if master fill displacements and reserve memory fo arraySequential
        int sum = 0;
        for (int j = 0; j < numtasks; j++) {
            sum += arraySizeSequential[j];
        }
        arraySequential = new char[sum*sizeof(char) + 1];
        displsSequential = new int[numtasks];
        displsSequential[0] = 0;
        for (int j = 0; j < numtasks-1; j++) 
            displsSequential[j+1] = displsSequential[j] + arraySizeSequential[j];

    }

    //data sent to root (master)
    MPI_Gatherv(buff, sizeSequential, MPI_CHAR, arraySequential, arraySizeSequential, displsSequential, MPI_CHAR, root, MPI_COMM_WORLD);

    //root gets intermediates product and use the sequential algorithm to calculate the ramainders corresponding 
    //to each root of the intermediate trees 
    char * rootBuffer;
    int * rootCount;
    int * rootDispls;
    if (taskid == root) {
        int sum = arraySizeSequential[numtasks-1] + displsSequential[numtasks-1];
        vector<mpz_class> localRootKeys = parseToMpz(arraySequential, sum);
        BatchGCD * localRootSolution = new BatchGCD(localRootKeys);
        vector<mpz_class> remaindersRoot = localRootSolution->getRemainders();

        rootCount = new int[numtasks];
        rootDispls = new int[numtasks];
        rootBuffer = parseToChar(remaindersRoot, &sizeInput);

        int bef = 0;
        for (int i = 0; i < numtasks; i++) {
            rootDispls[i] = bef;
            int j = bef;
            while (rootBuffer[j++] != divisor);
            rootCount[i] = j - bef;
            bef = j;
        }

    }

    //send size of each buffer that will be sent to each processor
    int rootLocalSize;
    MPI_Scatter(rootCount, 1, MPI_INT, &rootLocalSize, 1, MPI_INT, root, MPI_COMM_WORLD);

    char * rootLocalKeys = new char[rootLocalSize+1];

    //root (master) sends data to each processor (slave)
    MPI_Scatterv(rootBuffer, rootCount, rootDispls, MPI_CHAR, rootLocalKeys, rootLocalSize, MPI_CHAR, root, MPI_COMM_WORLD);

    string m;
    for(int i = 0; i < rootLocalSize-1; i++)
        m = m +  rootLocalKeys[i];
    mpz_class x(m);

    //finaly each processor slave print the factorization of the local numbers 
    localSolution->changeValueRoot(x);
    vector<Factor> factors = localSolution->getFactorization();

    clock_t end;
    if (taskid == 0) end = clock();

    for (int i = 0; i < factors.size(); i++) {
        int index = factors[i].getIndex();
        mpz_class p = factors[i].getP();
        mpz_class q = factors[i].getQ();
        mpz_class N = p * q;
        //if (p!= 1 && q!=1)
        cout << "N = " << N << " P = " << p << " Q = " << q << "\n";
    }

    if (taskid == 0) {
        double exec_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "TIME EXECUTION " << exec_secs << "\n";
    }

    //finish MPI :-)
    MPI_Finalize();

    return 0;
}

