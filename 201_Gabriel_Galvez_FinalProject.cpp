// Gabriel Galvez

// Modules needed for this file 
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <omp.h>

// Namespace for consistency 
using namespace std;

// 1a)
//
// Base @class { BaseCadd }
// BASic Computer-Aided Drug Design
class BaseCadd {
public:
    /* Public members */
    // Constructor 
    BaseCadd(){}

    // Public member variables
    int const num_features = 185; 
    std::vector<float> S; // Used for storing similarity values between chemical compounds
 
    // Member functions
    // Sim computes the similarity vector 
    // TODO: Sij = exp(-Dij * gamma) | gamma = 1.00 \ num_features
    void sim(std::vector<float> distance_vector)
    {
        // Column
        float similarity_column;
        for(float col : distance_vector)
        {
            // Calculate Similarity value & append to vector
            similarity_column = exp(-col * (1.00 / this->num_features));
            S.push_back(similarity_column);
        }   
    }
};

// 1b)
//
// Descendant @class { Euclid } <- BaseCadd
class Euclid : public BaseCadd {
public:
    /* Public members */
    // Constructor
    Euclid(){}

    // Distance Vector public variable
    std::vector<float> D; 

    // Public member functions 
    // Euclidean distance 
    // Compute Euclid distance 
    // TODO: dij = sum((Xi - Xj)^(1/2))
    void dist(vector<vector<float>> info_matrix, vector<float> known_matrix)
    {
        float dij;
        // Iterate through each row
        for(vector<float> row : info_matrix)
        {
            // Counter for summation
            dij = 0.0;
            for(int i = 0; i < num_features; i++)
            {
                dij += pow(row.at(i) - known_matrix.at(i), 2.0);
            }
            // After summation, push back the Euclidean distance to the vector
            D.push_back(pow(dij, 0.5));
        }
    }
}; 


// 1c)
//
// Descendant @class { Mink } <- Euclidean
class Mink : public Euclid {
public:
    // Constructor 
    Mink()
    {
        euclid = new Euclid;
    }

    // Destructor
    ~Mink()
    {
        delete euclid;
    }

    /* Public members */
    Euclid* euclid;

    // Member Function
    // { dist } Computes the Minkowski distance 
    // Mink.dist overrides Euclid.dist via direct call 
    // TODO: dij = sum((Xi - Xj)^p)^(1 / p)) | p = 3
    void dist(vector<vector<float>> info_matrix, vector<float> known_matrix)
    {
        float dij;
        // Iterate through each row
        for(vector<float> row : info_matrix)
        {
            // Counter for the summation
            dij = 0.00;
            for(int i = 0; i < num_features; i++)
            {
                dij += pow((row.at(i) - known_matrix.at(i)), 3.0);
            }
            // After the summation, push back Minkowski distance to the vector
            D.push_back(pow(dij, (1.0 / 3.0)));
        }
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////

// 1d)
//
// TODO: Read biological responce feature dataset from {bioresponse_descriptors_matrix.txt }
// TODO: Store the data in 2D array of floats or matrix called A, noting the number of columns
// TODO: Read in the biological responce feature data for the known drug from { known_drug.txt }

// File paths:

// Data used to compute program 
string  bioresponse_descriptors_matrix_path = "bioresponse_descriptors_matrix.txt";
string  known_drug_path                     = "known_drug.txt";

// Output files
string const result1 = "result1.txt";
string const result2 = "result2.txt";

// Result sorted
string const sorted_result1 = "sorted_result1.txt";
string const sorted_result2 = "sorted_result2.txt";

// 2)
//
// TODO::Sort similarity vector S from largest to smallest, printing out the top 10 most similar entries
// TODO::Apply to both S's in MK1 && MK2, read into respective files {  DSI\\C++\\Final Project\\sorted_result }

// To accomplish this tasks we define method { sortSimilarityVectors }
// Processing both in parallel
void sortSimilarityVectors(Mink& mk, string filename)
{
    // Vector used to perform calculations 
    vector<pair<int, float>> indexed_similarity(mk.S.size());

    // Populate vector with value / index pairs
    for(size_t i = 0; i < mk.S.size(); i++)
    {
        indexed_similarity[i] = {i, mk.S[i]};
    }

    // Sort the vector in descending order of similarity
    std::sort(indexed_similarity.begin(), indexed_similarity.end(), [](const pair<int, float>& a, const pair<int, float>& b){
        return a.second > b.second;
    });

    // Open save file path
    ofstream sorted_save(filename);

    // Output top 10, & stream into file
    for(int i = 0; (i < 10) && (i < indexed_similarity.size()); i++)
    {
        // Save to file
        sorted_save << indexed_similarity[i].first << " ";

        // Output message
        cout << "Index:" << indexed_similarity[i].first << " Value:" << indexed_similarity[i].second << " ";
    }
    cout << endl;

    // Close save file
    sorted_save.close();
}

// 1e), 1f)
int main()
{
    try
    {
        // Const Variables 
        int const columns = 185; // 185 numerical features are computed per chemical
        int const rows    = 299; // description_matrix has 299 rows 
        
        // Declare a matrix for known drug & description_matrix
        std::vector<std::vector<float>> A;

        std::vector<float> known_drug_vector;

        // Parse data from file's into matrix data structure 
        // Open file & handle any issues
        ifstream description_matrix(bioresponse_descriptors_matrix_path);
        if(!description_matrix.is_open()) throw runtime_error("failed to open file.");
        
        // Get buffer & parse data into datastructure
        string line;
        while(getline(description_matrix, line))
        {
            // Stream buffer for given line
            std::stringstream ss(line);
            // Used to append rows into the matrix
            std::vector<float> row;

            float val;
            while(ss >> val)
            {
                row.push_back(val);
            }

            A.push_back(row);
        }
        // Close file after parsing 
        description_matrix.close();


        // Parse data from known drug matrix 
        // Open & handle file 
        ifstream known_drug_matrix(known_drug_path);
        if(!known_drug_matrix.is_open()) throw runtime_error("failed to open file.");

        // Get file line buffer
        string line2;
        getline(known_drug_matrix, line2);
        stringstream ss2(line2);

        float val2;
        while(ss2 >> val2)
        {
            known_drug_vector.push_back(val2);
        }

        // End loop && close file 
        known_drug_matrix.close();

        // Declare two @class { Mink } objects
        Mink MK1; Mink MK2;

        // MK1 will compute the Minkowski distances between known_drug_matrix & each row within { A }
        MK1.dist(A, known_drug_vector);
        // MK2 will compute the Euclidean distances alternatively
        MK2.euclid->dist(A, known_drug_vector);

        // Get Similarity vectors for MK1 && MK2
        MK1.sim(MK1.D); MK2.sim(MK2.euclid->D);

        // Open Output Files && handle errors
        ofstream output1(result1); ofstream output2(result2);
        if(!output1.is_open() || !output2.is_open()) throw runtime_error("failed to open output files.");

        // Parse data into output files for saving
        // Start with { MK1 } -> output1
        for(float v : MK1.S) output1 << v << " ";

        // End with { MK2 } -> output2
        for(float v : MK2.S) output2 << v << "";

        // Assure to close files to avoid memory leaks
        output1.close(); output2.close();

        // Sort && output top 10 similar values 
        sortSimilarityVectors(MK1, sorted_result1);  // MK1 
        sortSimilarityVectors(MK2, sorted_result2); // MK2
    }
    catch(const runtime_error& e)
    {
        // Error occured during program's life
        cout << "Error: " << e.what() << endl;
    }
    return 0;
}


// 3) Run in parallel?
//
// TODO: Running program in parallel using OpenMP
/*
class BaseCadd {
public:
    int const num_features = 185;
    vector<float> S;

    void sim(const vector<float>& distance_vector) {
        S.clear();
        S.resize(distance_vector.size());
        #pragma omp parallel for
        for (int i = 0; i < distance_vector.size(); ++i) {
            S[i] = exp(-distance_vector[i] * (1.00 / num_features));
        }
    }    
};
class Euclid : public BaseCadd {
public:
    vector<float> D;

    void dist(const vector<vector<float>>& info_matrix, const vector<float>& known_matrix) {
        D.resize(info_matrix.size());
        #pragma omp parallel for
        for (int i = 0; i < info_matrix.size(); ++i) {
            float dij = 0.0;
            for (size_t j = 0; j < num_features; ++j) {
                dij += pow(info_matrix[i][j] - known_matrix[j], 2.0);
            }
            D[i] = sqrt(dij);
        }
    }
};

class Mink : public Euclid {
public:
    void dist(const vector<vector<float>>& info_matrix, const vector<float>& known_matrix) {
        D.resize(info_matrix.size());
        #pragma omp parallel for
        for (int i = 0; i < info_matrix.size(); ++i) {
            float dij = 0.0;
            for (size_t j = 0; j < num_features; ++j) {
                dij += pow(info_matrix[i][j] - known_matrix[j], 3.0);
            }
            D[i] = pow(dij, 1.0/3.0);
        }
    }
};

*/