/******************************************************************
typePotts.h

Energy function with Potts interactions:
   E(x)   =   \sum_i D_i(x_i)   +   \sum_ij V_ij(x_i,x_j)
   where x_i \in {0, 1, ..., K-1},
   V_ij(ki, kj) = edgePotential[k * ki + kj], edgePotential is a matrix k*k.

Example usage:

Minimize function E(x,y) = Dx(x) + Dy(y) + lambda*[x != y] where 
  x,y \in {0,1,2},
  Dx(0) = 0, Dx(1) = 1, Dx(2) = 2,
  Dy(0) = 3, Dy(1) = 4, Dy(2) = 5,
  lambda = 6,
  [.] is 1 if it's argument is true, and 0 otherwise.



#include "MRFEnergy.h"
#include <stdio.h>

void testPotts()
{
	MRFEnergy<TypePotts2>* mrf;
	MRFEnergy<TypePotts2>::NodeId* nodes;
	MRFEnergy<TypePotts2>::Options options;
	TypePotts2::REAL energy, lowerBound;

	const int nodeNum = 2; // number of nodes
	const int K = 3; // number of labels
	TypePotts2::REAL D[K];
	int x, y;

	mrf = new MRFEnergy<TypePotts2>(TypePotts2::GlobalSize(K));
	nodes = new MRFEnergy<TypePotts2>::NodeId[nodeNum];

	// construct energy
	D[0] = 0; D[1] = 1; D[2] = 2;
	nodes[0] = mrf->AddNode(TypePotts2::LocalSize(), TypePotts2::NodeData(D));
	D[0] = 3; D[1] = 4; D[2] = 5;
	nodes[1] = mrf->AddNode(TypePotts2::LocalSize(), TypePotts2::NodeData(D));
	mrf->AddEdge(nodes[0], nodes[1], TypePotts2::EdgeData(6));

	// Function below is optional - it may help if, for example, nodes are added in a random order
	// mrf->SetAutomaticOrdering();

	/////////////////////// TRW-S algorithm //////////////////////
	options.m_iterMax = 30; // maximum number of iterations
	mrf->Minimize_TRW_S(options, lowerBound, energy);

	// read solution
	x = mrf->GetSolution(nodes[0]);
	y = mrf->GetSolution(nodes[1]);

	printf("Solution: %d %d\n", x, y);

	//////////////////////// BP algorithm ////////////////////////
	mrf->ZeroMessages(); // in general not necessary - it may be faster to start 
	                     // with messages computed in previous iterations

	options.m_iterMax = 30; // maximum number of iterations
	mrf->Minimize_BP(options, energy);

	// read solution
	x = mrf->GetSolution(nodes[0]);
	y = mrf->GetSolution(nodes[1]);

	printf("Solution: %d %d\n", x, y);

	// done
	delete nodes;
	delete mrf;
}

*******************************************************************/

#ifndef __TYPEPOTTS2_H__
#define __TYPEPOTTS2_H__

#include <assert.h>
#include <string.h>

template <class T>
class MRFEnergy;

const int LABLE_NUM = 6;

class TypePotts2 {
private:
    struct Vector; // node parameters and messages
    struct Edge; // stores edge information and either forward or backward message

public:
    // types declarations
    typedef int Label;
    typedef double REAL;
    struct GlobalSize; // global information about number of labels
    struct LocalSize; // local information about number of labels (stored at each node)
    struct NodeData; // argument to MRFEnergy::AddNode()
    struct EdgeData; // argument to MRFEnergy::AddEdge()

    struct GlobalSize {
        GlobalSize(int K);

    private:
        friend struct Vector;
        friend struct Edge;
        int m_K; // number of labels
    };

    struct LocalSize // number of labels is stored at MRFEnergy::m_Kglobal
    {
    };

    struct NodeData {
        NodeData(REAL* data); // data = pointer to array of size MRFEnergy::m_Kglobal

    private:
        friend struct Vector;
        friend struct Edge;
        REAL* m_data;
    };

    struct EdgeData {
        EdgeData(REAL* edgePotential);

    private:
        friend struct Vector;
        friend struct Edge;
        REAL* m_edgePotential;
    };

    //////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Visible only to MRFEnergy /////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

private:
    friend class MRFEnergy<TypePotts2>;

    struct Vector {
        static int GetSizeInBytes(GlobalSize Kglobal, LocalSize K); // returns -1 if invalid K's
        void Initialize(GlobalSize Kglobal, LocalSize K, NodeData data); // called once when user adds a node
        void Add(GlobalSize Kglobal, LocalSize K, NodeData data); // called once when user calls MRFEnergy::AddNodeData()

        void SetZero(GlobalSize Kglobal, LocalSize K); // set this[k] = 0
        void Copy(GlobalSize Kglobal, LocalSize K, Vector* V); // set this[k] = V[k]
        void Add(GlobalSize Kglobal, LocalSize K, Vector* V); // set this[k] = this[k] + V[k]
        REAL GetValue(GlobalSize Kglobal, LocalSize K, Label k); // return this[k]
        REAL ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin); // return vMin = min_k { this[k] }, set kMin
        REAL ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K); // same as previous, but additionally set this[k] -= vMin (and kMin is not returned)

        static int GetArraySize(GlobalSize Kglobal, LocalSize K);
        REAL GetArrayValue(GlobalSize Kglobal, LocalSize K, int k); // note: k is an integer in [0..GetArraySize()-1].
            // For Potts functions GetArrayValue() and GetValue() are the same,
            // but they are different for, say, 2-dimensional labels.
        void SetArrayValue(GlobalSize Kglobal, LocalSize K, int k, REAL x);

    private:
        friend struct Edge;
        REAL m_data[1]; // actual size is MRFEnergy::m_Kglobal
    };

    struct Edge {
        static int GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data); // returns -1 if invalid data
        static int GetBufSizeInBytes(int vectorMaxSizeInBytes); // returns size of buffer need for UpdateMessage()
        void Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj); // called once when user adds an edge
        Vector* GetMessagePtr();
        void Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj); // if the client calls this function, then the meaning of 'dir'
            // in distance transform functions is swapped

        // When UpdateMessage() is called, edge contains message from dest to source.
        // The function must replace it with the message from source to dest.
        // The update rule is given below assuming that source corresponds to tail (i) and dest corresponds
        // to head (j) (which is the case if dir==0).
        //
        // 1. Compute Di[ki] = gamma*source[ki] - message[ki].  (Note: message = message from j to i).
        // 2. Compute distance transform: set
        //       message[kj] = min_{ki} (Di[ki] + V(ki,kj)). (Note: message = message from i to j).
        // 3. Compute vMin = min_{kj} m_message[kj].
        // 4. Set m_message[kj] -= vMin.
        // 5. Return vMin.
        //
        // If dir==1 then source corresponds to j, sink corresponds to i. Then the update rule is
        //
        // 1. Compute Dj[kj] = gamma*source[kj] - message[kj].  (Note: message = message from i to j).
        // 2. Compute distance transform: set
        //       message[ki] = min_{kj} (Dj[kj] + V(ki,kj)). (Note: message = message from j to i).
        // 3. Compute vMin = min_{ki} m_message[ki].
        // 4. Set m_message[ki] -= vMin.
        // 5. Return vMin.
        //
        // If Edge::Swap has been called odd number of times, then the meaning of dir is swapped.
        //
        // Vector 'source' must not be modified. Function may use 'buf' as a temporary storage.
        REAL UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* buf);

        // If dir==0, then sets dest[kj] += V(ksource,kj).
        // If dir==1, then sets dest[ki] += V(ki,ksource).
        // If Swap() has been called odd number of times, then the meaning of dir is swapped.
        void AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir);

    private:
        // edge information
        REAL m_edgePotential[LABLE_NUM][LABLE_NUM];

        // message
        Vector m_message;
    };
};

//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Implementation ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

inline TypePotts2::GlobalSize::GlobalSize(int K)
{
    m_K = K;
}

///////////////////// NodeData and EdgeData ///////////////////////

inline TypePotts2::NodeData::NodeData(REAL* data)
{
    m_data = data;
}

inline TypePotts2::EdgeData::EdgeData(REAL* edgePotential)
{
    m_edgePotential = edgePotential;
}

///////////////////// Vector ///////////////////////

inline int TypePotts2::Vector::GetSizeInBytes(GlobalSize Kglobal, LocalSize K)
{
    if (Kglobal.m_K < 1) {
        return -1;
    }
    return Kglobal.m_K * sizeof(REAL);
}
inline void TypePotts2::Vector::Initialize(GlobalSize Kglobal, LocalSize K, NodeData data)
{
    memcpy(m_data, data.m_data, Kglobal.m_K * sizeof(REAL));
}

inline void TypePotts2::Vector::Add(GlobalSize Kglobal, LocalSize K, NodeData data)
{
    for (int k = 0; k < Kglobal.m_K; k++) {
        m_data[k] += data.m_data[k];
    }
}

inline void TypePotts2::Vector::SetZero(GlobalSize Kglobal, LocalSize K)
{
    memset(m_data, 0, Kglobal.m_K * sizeof(REAL));
}

inline void TypePotts2::Vector::Copy(GlobalSize Kglobal, LocalSize K, Vector* V)
{
    memcpy(m_data, V->m_data, Kglobal.m_K * sizeof(REAL));
}

inline void TypePotts2::Vector::Add(GlobalSize Kglobal, LocalSize K, Vector* V)
{
    for (int k = 0; k < Kglobal.m_K; k++) {
        m_data[k] += V->m_data[k];
    }
}

inline TypePotts2::REAL TypePotts2::Vector::GetValue(GlobalSize Kglobal, LocalSize K, Label k)
{
    assert(k >= 0 && k < Kglobal.m_K);
    return m_data[k];
}

inline TypePotts2::REAL TypePotts2::Vector::ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin)
{
    REAL vMin = m_data[0];
    kMin = 0;
    for (int k = 1; k < Kglobal.m_K; k++) {
        if (vMin > m_data[k]) {
            vMin = m_data[k];
            kMin = k;
        }
    }

    return vMin;
}

inline TypePotts2::REAL TypePotts2::Vector::ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K)
{
    REAL vMin = m_data[0];
    for (int k = 1; k < Kglobal.m_K; k++) {
        if (vMin > m_data[k]) {
            vMin = m_data[k];
        }
    }
    for (int k = 0; k < Kglobal.m_K; k++) {
        m_data[k] -= vMin;
    }

    return vMin;
}

inline int TypePotts2::Vector::GetArraySize(GlobalSize Kglobal, LocalSize K)
{
    return Kglobal.m_K;
}

inline TypePotts2::REAL TypePotts2::Vector::GetArrayValue(GlobalSize Kglobal, LocalSize K, int k)
{
    assert(k >= 0 && k < Kglobal.m_K);
    return m_data[k];
}

inline void TypePotts2::Vector::SetArrayValue(GlobalSize Kglobal, LocalSize K, int k, REAL x)
{
    assert(k >= 0 && k < Kglobal.m_K);
    m_data[k] = x;
}

///////////////////// EdgeDataAndMessage implementation /////////////////////////

inline int TypePotts2::Edge::GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data)
{
    if (data.m_edgePotential == nullptr) {
        return -1;
    }
    return sizeof(Edge) - sizeof(Vector) + Kglobal.m_K * sizeof(REAL);
}

inline int TypePotts2::Edge::GetBufSizeInBytes(int vectorMaxSizeInBytes)
{
    return 0;
}

inline void TypePotts2::Edge::Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj)
{
    //     m_lambdaPotts = data.m_lambdaPotts;
    //     int len = Kglobal.m_K * Kglobal.m_K;
    //     memset(m_edgePotential, 0, len * sizeof(REAL));
    //     memcpy();
    for (int i = 0; i < LABLE_NUM; i++) {
        for (int j = 0; j < LABLE_NUM; j++)
            m_edgePotential[i][j] = data.m_edgePotential[i * LABLE_NUM + j];
    }

    memset(m_message.m_data, 0, Kglobal.m_K * sizeof(REAL));
}

inline TypePotts2::Vector* TypePotts2::Edge::GetMessagePtr()
{
    return &m_message;
}

inline void TypePotts2::Edge::Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj)
{
}

inline TypePotts2::REAL TypePotts2::Edge::UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* buf)
{
    int k;
    Vector* _buf = (Vector*)buf;
    REAL vMin;

    int ksource, kdest;

    for (ksource = 0; ksource < Kglobal.m_K; ksource++) {
        _buf->m_data[ksource] = gamma * source->m_data[ksource] - m_message.m_data[ksource];
    }

    for (kdest = 0; kdest < Kglobal.m_K; kdest++) {
        vMin = _buf->m_data[0] + this->m_edgePotential[0][kdest];
        for (ksource = 1; ksource < Kglobal.m_K; ksource++) {
            if (vMin > _buf->m_data[ksource] + this->m_edgePotential[ksource][kdest]) {
                vMin = _buf->m_data[ksource] + this->m_edgePotential[ksource][kdest];
            }
        }
        m_message.m_data[kdest] = vMin;
    }
    vMin = m_message.m_data[0];
    for (kdest = 1; kdest < Kglobal.m_K; kdest++) {
        if (vMin > m_message.m_data[kdest]) {
            vMin = m_message.m_data[kdest];
        }
    }

    for (kdest = 0; kdest < Kglobal.m_K; kdest++) {
        m_message.m_data[kdest] -= vMin;
    }

    return vMin;
}

inline void TypePotts2::Edge::AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir)
{
    for (int k = 0; k < Kglobal.m_K; k++) {
        dest->m_data[k] += this->m_edgePotential[ksource][k];
    }
}

//////////////////////////////////////////////////////////////////////////////////

#endif
