# TRW-S
Tree-reweighted max-product message passing algorithm

[代码原有位置](https://www.microsoft.com/en-us/download/details.aspx?id=52499&from=http%3A%2F%2Fresearch.microsoft.com%2Fen-us%2Fdownloads%2Fdad6c31e-2c04-471f-b724-ded18bf70fe3%2F)



原有声明：

```
This software implements two algorithms for minimizing energy functions of the form 

   E(x)   =   \sum_i D_i(x_i)   +   \sum_ij V_ij(x_i,x_j)
   where x_i are discrete variables.

The two algorithms are max-product belief propagation (BP, Pearl'88) and
sequential tree-reweighted max-product message passing (TRW-S, Kolmogorov'05).

For example usage look at one of the type*.h files
(typeBinary.h, typeBinaryFast.h, typePotts.h, typeGeneral.h,
typeTruncatedLinear.h, typeTruncatedQuadratic.h, typeTruncatedLinear2D.h, typeTruncatedQuadratic2D.h).
If your energy function does not belong to the classes defined
in these files but terms V_ij allow fast distance transforms, then
it should be possible to extend the algorithms to your functions
using files type*.h as examples.

Written by Vladimir Kolmogorov (vnk@microsoft.com), 2005.
Tested under Microsoft Visual Studio .NET (Windows)
and GNU c++ compiler version 2.96 (Red Hat Linux 7.1).

(c) Microsoft Corporation. All rights reserved. 
```

