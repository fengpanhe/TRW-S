#include "MRFEnergy.h"
#include <stdio.h>

// Example: minimizing an energy function with Potts terms.
// See type*.h files for other types of terms.
// template MRFEnergy<TypePotts>;

// void testPotts()
// {
//     MRFEnergy<TypePotts>* mrf;
//     MRFEnergy<TypePotts>::NodeId* nodes;
//     MRFEnergy<TypePotts>::Options options;
//     TypePotts::REAL energy, lowerBound;

//     const int nodeNum = 4; // number of nodes
//     const int K = 3; // number of labels
//     TypePotts::REAL D[K];
//     int x, y, z, k;

//     mrf = new MRFEnergy<TypePotts>(TypePotts::GlobalSize(K), nullptr);
//     nodes = new MRFEnergy<TypePotts>::NodeId[nodeNum];

//     // construct energy
//     D[0] = 10;
//     D[1] = 3;
//     D[2] = 4;
//     nodes[0] = mrf->AddNode(TypePotts::LocalSize(), TypePotts::NodeData(D));
//     D[0] = 2;
//     D[1] = 10;
//     D[2] = 5;
//     nodes[1] = mrf->AddNode(TypePotts::LocalSize(), TypePotts::NodeData(D));
//     D[0] = 3;
//     D[1] = 10;
//     D[2] = 5;
//     nodes[2] = mrf->AddNode(TypePotts::LocalSize(), TypePotts::NodeData(D));
//     D[0] = 3;
//     D[1] = 10;
//     D[2] = 5;
//     nodes[3] = mrf->AddNode(TypePotts::LocalSize(), TypePotts::NodeData(D));

//     mrf->AddEdge(nodes[0], nodes[2], TypePotts::EdgeData(1));
//     //     mrf->AddEdge(nodes[0], nodes[1], TypePotts::EdgeData(2));
//     mrf->AddEdge(nodes[0], nodes[3], TypePotts::EdgeData(10));

//     // Function below is optional - it may help if, for example, nodes are added in a random order
//     // mrf->SetAutomaticOrdering();

//     /////////////////////// TRW-S algorithm //////////////////////
//     options.m_iterMax = 30; // maximum number of iterations
//     mrf->Minimize_TRW_S(options, lowerBound, energy);

//     // read solution
//     x = mrf->GetSolution(nodes[0]);
//     y = mrf->GetSolution(nodes[1]);
//     z = mrf->GetSolution(nodes[2]);
//     k = mrf->GetSolution(nodes[3]);

//     printf("Solution: %d %d %d %d\n", x, y, z, k);

//     //////////////////////// BP algorithm ////////////////////////
//     mrf->ZeroMessages(); // in general not necessary - it may be faster to start
//         // with messages computed in previous iterations.
//         // NOTE: in most cases, immediately after creating the energy
//         // all messages are zero. EXCEPTION: typeBinary and typeBinaryFast.
//         // So calling ZeroMessages for these types will NOT transform
//         // the energy to the original state.

//     options.m_iterMax = 30; // maximum number of iterations
//     mrf->Minimize_BP(options, energy);

//     // read solution
//     x = mrf->GetSolution(nodes[0]);
//     y = mrf->GetSolution(nodes[1]);
//     z = mrf->GetSolution(nodes[2]);
//     k = mrf->GetSolution(nodes[3]);

//     printf("Solution: %d %d %d %d\n", x, y, z, k);

//     // done
//     delete nodes;
//     delete mrf;
// }

void testPotts2()
{
    MRFEnergy<TypePotts2>* mrf;
    MRFEnergy<TypePotts2>::NodeId* nodes;
    MRFEnergy<TypePotts2>::Options options;
    TypePotts2::REAL energy, lowerBound;

    const int nodeNum = 4; // number of nodes
    const int K = 6; // number of labels
    TypePotts2::REAL D[K];
    int x, y, z, k;

    mrf = new MRFEnergy<TypePotts2>(TypePotts2::GlobalSize(K), nullptr);
    nodes = new MRFEnergy<TypePotts2>::NodeId[nodeNum];

    // construct energy
    D[0] = 1;
    D[1] = 0;
    D[2] = 3;
    D[3] = 4;
    D[4] = 5;
    D[5] = 6;
    nodes[0] = mrf->AddNode(TypePotts2::LocalSize(), TypePotts2::NodeData(D));
    D[0] = 1;
    D[1] = 2;
    D[2] = 3;
    D[3] = 4;
    D[4] = 5;
    D[5] = 6;
    nodes[1] = mrf->AddNode(TypePotts2::LocalSize(), TypePotts2::NodeData(D));
    D[0] = 1;
    D[1] = 2;
    D[2] = 3;
    D[3] = 4;
    D[4] = 5;
    D[5] = 6;
    nodes[2] = mrf->AddNode(TypePotts2::LocalSize(), TypePotts2::NodeData(D));
    D[0] = 1;
    D[1] = 2;
    D[2] = 3;
    D[3] = 4;
    D[4] = 5;
    D[5] = 6;
    nodes[3] = mrf->AddNode(TypePotts2::LocalSize(), TypePotts2::NodeData(D));

    TypePotts2::REAL aa[36] = { 0, 0.1, 0, 0, 0, 0, 0.1 };
    mrf->AddEdge(nodes[0], nodes[2], TypePotts2::EdgeData(aa));
    mrf->AddEdge(nodes[0], nodes[1], TypePotts2::EdgeData(aa));
    //     mrf->AddEdge(nodes[0], nodes[1], TypePotts2::EdgeData(2));
    //     mrf->AddEdge(nodes[3], nodes[0], TypePotts2::EdgeData(aa));
    mrf->AddEdge(nodes[0], nodes[3], TypePotts2::EdgeData(aa));
    //     printf("%s\n", "hello");
    // Function below is optional - it may help if, for example, nodes are added in a random order
    // mrf->SetAutomaticOrdering();

    /////////////////////// TRW-S algorithm //////////////////////
    options.m_iterMax = 30; // maximum number of iterations
    mrf->Minimize_TRW_S(options, lowerBound, energy);

    // read solution
    x = mrf->GetSolution(nodes[0]);
    y = mrf->GetSolution(nodes[1]);
    z = mrf->GetSolution(nodes[2]);
    k = mrf->GetSolution(nodes[3]);

    printf("Solution: %d %d %d %d\n", x, y, z, k);

    //////////////////////// BP algorithm ////////////////////////
    mrf->ZeroMessages(); // in general not necessary - it may be faster to start
        // with messages computed in previous iterations.
        // NOTE: in most cases, immediately after creating the energy
        // all messages are zero. EXCEPTION: typeBinary and typeBinaryFast.
        // So calling ZeroMessages for these types will NOT transform
        // the energy to the original state.

    options.m_iterMax = 30; // maximum number of iterations
    mrf->Minimize_BP(options, energy);

    // read solution
    x = mrf->GetSolution(nodes[0]);
    y = mrf->GetSolution(nodes[1]);
    z = mrf->GetSolution(nodes[2]);
    k = mrf->GetSolution(nodes[3]);

    printf("Solution: %d %d %d %d\n", x, y, z, k);

    // done
    delete nodes;
    delete mrf;
}
// typedef struct
// {
//         int a;
//         int b;
//         char buf[0];    // 或者char buf[];
// }Node;

// typedef struct T{
//     int aa;
//     int aaa;
//     int aaaa;
//     char m1[10];
//     Node *m2;
// }TT;

int main()
{
    // printf("%s\n", "hello");
    // char* a_str = "a";
    // char* b_str = "zxc";
    // TT* t = (TT*)new char(sizeof(TT) + sizeof(char) * strlen(a_str) + sizeof(char) * strlen(b_str) + 10);
    // t->m2 = (Node*)((char*)t + strlen(a_str) * sizeof(char) + 3 * sizeof(int) + 1);
    // //     memset(t->m2->a, 0, strlen(b_str));
    // strcpy(t->m1, a_str);
    // memcpy(t->m2->buf, b_str, strlen(b_str));
    // printf("t->a = %s\n", t->m1);
    // printf("t->b = %s\n", t->m2->buf);
    // printf("%d\n", (char*)t);
    // printf("%d\n", (char*)t->m1);
    // printf("%d\n", (char*)t->m2);
    // delete t;

    //      printf("%d\n", sizeof(Node));
    //
    //         Node *p = (Node *)malloc(sizeof(Node) + 16);
    //         p->a = 1;
    //         p->b = 2;
    //         strcpy(p->buf, "hello");
    //
    //         printf("node : %p\n", p);
    //         printf("node::a : %p, %d\n", &p->a, p->a);
    //         printf("node::b : %p, %d\n", &p->b, p->b);
    //         printf("node::buf : %p, %s\n", p->buf, p->buf);
    //
    //         free(p);

    // 	testPotts();
    testPotts2();
}
