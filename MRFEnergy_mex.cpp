/*
 *  Filename: MRFEnergy_mex.cpp
 *  Project: TRW-S
 *  Created Date: Sunday January 20th 2019
 *  Author: Feng Panhe
 *  -----
 *  Last Modified:
 *  Modified By:
 *  -----
 *  Copyright (c) 2019 Feng Panhe
 */

#include "mex.hpp"
#include "mexAdapter.hpp"

using matlab::mex::ArgumentList;
using namespace matlab::data;
#include "MRFEnergy.h"
#include <stdio.h>

void DefaultErrorFn(char* msg)
{
    fprintf(stderr, "%s\n", msg);
    exit(1);
}
const int LABEL_ABC = 0;
const int LABEL_ACB = 1;
const int LABEL_BAC = 2;
const int LABEL_BCA = 3;
const int LABEL_CAB = 4;
const int LABEL_CBA = 5;

// ORDER_LABEL LABEL;
struct TJunction {
    int imid;
    int tid;
    double score[3];
    int eids[3];
    int rank[3];
};

double* D;
double* V;

// double** infoMatrix = NULL;
struct TJunction* tjs = NULL;
double penalties = 0;

double dCost(int pix, int label);
double fnCost(int pix1, int pix2, int label1, int label2);
bool existSameEdge(int pix1, int pix2);

class MexFunction : public matlab::mex::Function {
    void operator()(matlab::mex::ArgumentList outputs,
        matlab::mex::ArgumentList inputs)
    {
        checkArguments(outputs, inputs);

        penalties = inputs[0][0];
        matlab::data::TypedArray<double> inMatrix = std::move(inputs[1]);

        int row_num = inputs[2][0];
        int col_num = inputs[2][1];

        int t_junction_num = row_num / 3;
        tjs = new TJunction[t_junction_num];

        for (int i = 0; i < t_junction_num; i++) {
            int index0 = i * 3;
            int index1 = index0 + 1;
            int index2 = index0 + 2;
            if (inMatrix[index0][1] != inMatrix[index1][1] || inMatrix[index0][1] != inMatrix[index2][1]) {
                matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                    0, std::vector<matlab::data::Array>({ factory.createScalar("Tid is mismatch") }));
            }
            if (inMatrix[index0][0] != inMatrix[index1][0] || inMatrix[index0][0] != inMatrix[index2][0]) {
                matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                    0, std::vector<matlab::data::Array>({ factory.createScalar("Imageid is mismatch") }));
            }
            tjs[i].imid = inMatrix[index0][0];
            tjs[i].tid = inMatrix[index0][1];
            tjs[i].eids[0] = inMatrix[index0][2];
            tjs[i].eids[1] = inMatrix[index1][2];
            tjs[i].eids[2] = inMatrix[index2][2];
            tjs[i].score[0] = inMatrix[index0][3];
            tjs[i].score[1] = inMatrix[index1][3];
            tjs[i].score[2] = inMatrix[index2][3];
        }

        int imid = tjs[0].imid;
        int im_start_index = 0;
        for (int i = 1; i < t_junction_num; i++) {
            if (tjs[i].imid != imid) {
                this->use_TRWS(im_start_index, i - 1);
                im_start_index = i;
            }
        }
        for (int i = 0; i < t_junction_num; i++) {
            int index0 = i * 3;
            int index1 = index0 + 1;
            int index2 = index0 + 2;
            inMatrix[index0][3] = tjs[i].rank[0];
            inMatrix[index1][3] = tjs[i].rank[1];
            inMatrix[index2][3] = tjs[i].rank[2];
        }
        delete[] tjs;
        outputs[0] = std::move(inMatrix);

        // infoMatrix = new double*[row_num + 1];
        // for (int i = 0; i < row_num; i++) {
        //     infoMatrix[i] = new double[col_num + 1];
        // }

        // for (int i = 0; i < row_num; i++) {
        //     for (int j = 0; j < col_num; j++) {
        //         infoMatrix[i][j] = inMatrix[i][j];
        //     }
        // }
        // this->minimize_use_TRWS();
        // for (int i = 0; i < row_num; i++) {
        //     delete[] infoMatrix[i];
        // }
        // delete[] infoMatrix;
    }
    void use_TRWS(int im_start_index, int im_end_index)
    {
        MRFEnergy<TypePotts2>* mrf;
        MRFEnergy<TypePotts2>::NodeId* nodes;
        MRFEnergy<TypePotts2>::Options options;
        TypePotts2::REAL energy, lowerBound;

        const int nodeNum = im_end_index - im_start_index + 1; // number of nodes
        const int K = 6; // number of labels
        TypePotts2::REAL D[K];
        TypePotts2::REAL V[K * K];
        int x, y, z, k;

        mrf = new MRFEnergy<TypePotts2>(TypePotts2::GlobalSize(K), DefaultErrorFn);
        nodes = new MRFEnergy<TypePotts2>::NodeId[nodeNum];

        // construct energy
        for (int i = 0; i < nodeNum; i++) {
            for (int l = 0; l < k; l++) {
                D[l] = dCost(im_start_index + i, l);
            }
            nodes[i] = mrf->AddNode(TypePotts2::LocalSize(), TypePotts2::NodeData(D));
        }

        int pix1 = 0;
        int pix2 = 0;
        for (int i = 0; i < nodeNum - 1; i++) {
            for (int j = i + 1; j < nodeNum; j++) {
                pix1 = im_start_index + i;
                pix2 = im_start_index + j;
                if (existSameEdge(pix1, pix2)) {

                    for (int m = 0; m < K; m++) {
                        for (int n = 0; n < K; n++) {
                            V[m * k + n] = fnCost(pix1, pix2, m, n);
                        }
                    }
                    mrf->AddEdge(nodes[i], nodes[j], TypePotts2::EdgeData(V));
                }
            }
        }

        // TypePotts2::REAL aa[36] = { 0, 0.1, 0, 0, 0, 0, 0.1 };
        // mrf->AddEdge(nodes[0], nodes[2], TypePotts2::EdgeData(aa));
        // mrf->AddEdge(nodes[0], nodes[1], TypePotts2::EdgeData(aa));
        //     mrf->AddEdge(nodes[0], nodes[1], TypePotts2::EdgeData(2));
        //     mrf->AddEdge(nodes[3], nodes[0], TypePotts2::EdgeData(aa));
        // mrf->AddEdge(nodes[0], nodes[3], TypePotts2::EdgeData(aa));
        //     printf("%s\n", "hello");
        // Function below is optional - it may help if, for example, nodes are added in a random order
        // mrf->SetAutomaticOrdering();

        options.m_iterMax = 30; // maximum number of iterations
        mrf->Minimize_TRW_S(options, lowerBound, energy);

        // read solution
        x = mrf->GetSolution(nodes[0]);
        for (int i = 0; i < nodeNum; i++) {

            switch (mrf->GetSolution(nodes[i])) {
            case LABEL_ABC:
                tjs[im_start_index + i].rank[0] = 3;
                tjs[im_start_index + i].rank[2] = 2;
                tjs[im_start_index + i].rank[3] = 1;
                break;
            case LABEL_ACB:
                tjs[im_start_index + i].rank[0] = 3;
                tjs[im_start_index + i].rank[2] = 1;
                tjs[im_start_index + i].rank[3] = 2;
                break;
            case LABEL_BAC:
                tjs[im_start_index + i].rank[0] = 2;
                tjs[im_start_index + i].rank[2] = 3;
                tjs[im_start_index + i].rank[3] = 1;
                break;
            case LABEL_CAB:
                tjs[im_start_index + i].rank[0] = 2;
                tjs[im_start_index + i].rank[2] = 1;
                tjs[im_start_index + i].rank[3] = 3;
                break;
            case LABEL_CBA:
                tjs[im_start_index + i].rank[0] = 1;
                tjs[im_start_index + i].rank[2] = 2;
                tjs[im_start_index + i].rank[3] = 3;
                break;
            case LABEL_BCA:
                tjs[im_start_index + i].rank[0] = 1;
                tjs[im_start_index + i].rank[2] = 3;
                tjs[im_start_index + i].rank[3] = 2;
                break;
            default:
                break;
            }
        }

        // done
        delete nodes;
        delete mrf;
    }
    void checkArguments(matlab::mex::ArgumentList outputs,
        matlab::mex::ArgumentList inputs)
    {
        // std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        // matlab::data::ArrayFactory factory;

        if (inputs.size() != 3) {
            matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                0, std::vector<matlab::data::Array>({ factory.createScalar("Two inputs required") }));
        }

        if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE || inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE || inputs[0].getNumberOfElements() != 1) {
            matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input penalties must be a scalar") }));
        }

        if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE || inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input matrix must be type double") }));
        }

        if (inputs[1].getDimensions().size() != 2) {
            matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input must be m-by-n dimension") }));
        }
    }

    void displayOnMATLAB(std::ostringstream& stream)
    {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }

    // void minimize_use_TRWS(sizeX, sizeY, numLabels, energy)
    // {
    //     MRF* mrf;
    //     EnergyFunction* energy;
    //     MRF::EnergyVal E;
    //     double lowerBound;
    //     float t, tot_t;
    //     int iter;

    //     energy = generate_DataFUNCTION_SmoothGENERAL_FUNCTION();

    //     // printf("\n*******Started TRW-S *****\n");
    //     stream << "\n*******Started TRW-S *****\n";
    //     this.displayOnMATLAB(stream);

    //     mrf = new TRWS(sizeX, sizeY, numLabels, energy);

    //     // can disable caching of values of general smoothness function:
    //     // mrf->dontCacheSmoothnessCosts();

    //     mrf->initialize();
    //     mrf->clearAnswer();

    //     E = mrf->totalEnergy();
    //     printf("Energy at the Start= %g (%g,%g)\n", (float)E,
    //         (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    //     tot_t = 0;
    //     for (iter = 0; iter < 10; iter++) {
    //         mrf->optimize(10, t);

    //         E = mrf->totalEnergy();
    //         lowerBound = mrf->lowerBound();
    //         tot_t = tot_t + t;
    //         printf("energy = %g, lower bound = %f (%f secs)\n", (float)E,
    //             lowerBound, tot_t);
    //     }
    //     for (int pix = 0; pix < sizeX * sizeY; pix++)
    //         printf("Label of pixel %d is %d", pix, mrf->getLabel(pix));

    //     delete mrf;
    // }

private:
    // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    matlab::data::ArrayFactory factory;

    // Create an output stream
    std::ostringstream stream;
};

double dCost(int pix, int label)
{
    double cost = 0;
    double a = tjs[pix].score[0];
    double b = tjs[pix].score[1];
    double c = tjs[pix].score[2];
    switch (label) {
    case LABEL_ABC:
        cost = a * 2 - c * 2;
        break;
    case LABEL_ACB:
        cost = a * 2 - b * 2;
        break;
    case LABEL_BAC:
        cost = b * 2 - c * 2;
        break;
    case LABEL_BCA:
        cost = b * 2 - a * 2;
        break;
    case LABEL_CAB:
        cost = c * 2 - b * 2;
        break;
    case LABEL_CBA:
        cost = c * 2 - a * 2;
        break;

    default:
        break;
    }
    return -cost;
}

bool existSameEdge(int pix1, int pix2)
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (tjs[pix1].eids[i] == tjs[pix2].eids[j]) {
                return true;
            }
        }
    }
    return false;
}

double fnCost(int pix1, int pix2, int label1, int label2)
{

    double answer = 0;

    int equal_edgeid_index[3];

    for (int i = 0; i < 3; i++) {
        equal_edgeid_index[i] = -1;
        for (int j = 0; j < 3; j++) {
            if (tjs[pix1].eids[i] == tjs[pix2].eids[j]) {
                equal_edgeid_index[i] = j;
                break;
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        if (equal_edgeid_index[i] != -1) {
            bool flag1 = false;
            bool flag2 = false;

            switch (i) {
            case 0:
                if (label1 == LABEL_ABC || label1 == LABEL_ACB || label1 == LABEL_CAB) {
                    flag1 = true;
                }
                break;
            case 1:
                if (label1 == LABEL_ABC || label1 == LABEL_BAC || label1 == LABEL_BCA) {
                    flag1 = true;
                }
                break;
            case 2:
                if (label1 == LABEL_BCA || label1 == LABEL_CBA || label1 == LABEL_CAB) {
                    flag1 = true;
                }
                break;

            default:
                break;
            }

            switch (equal_edgeid_index[i]) {
            case 0:
                if (label2 == LABEL_ABC || label2 == LABEL_ACB || label2 == LABEL_CAB) {
                    flag2 = true;
                }
                break;
            case 1:
                if (label2 == LABEL_ABC || label2 == LABEL_BAC || label2 == LABEL_BCA) {
                    flag2 = true;
                }
                break;
            case 2:
                if (label2 == LABEL_BCA || label2 == LABEL_CBA || label2 == LABEL_CAB) {
                    flag2 = true;
                }
                break;

            default:
                break;
            }

            if (flag1 == flag2) {
                answer += penalties;
            }
        }
    }

    return answer;
}