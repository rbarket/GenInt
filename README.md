# GenInt
Maple Implementation of the paper [Generating Elementary Integrable Expressions](https://link.springer.com/chapter/10.1007/978-3-031-41724-5_2).

## How to use

This repository contains code for the data generation method of generating (integrand, integral) pairs associated with the paper along with a curated dataset that is generated from the method.
`gen_risch_data.mpl` will generate the (integrand, integral) pairs For a given list of extensions and their arguments. If you wish to run this yourself with custom extensions and arguments, simply change the lists in the file.

The main function that actually generates the pair is `TR_lin_gen()` in `risch_funcs.mpl`. It takes the following arguments:
- num_terms: how many denominators to create for the partial fraction representation
- theta_sub: Which elementary extension we are working with (log(u) or exp(u), u in the field we are working in)
- m: the multiplicity of the denominators of each partial fraction. If m=1, then the denominator is square-free and only the TR-algorithm is used. If m>1, perform hermite reduction until the denominator is square-free and then use the TR-algorithm

The current implementation only generates linear denominators for each of the partial fractions generated (num_terms in `TR_lin_gen()`). Quadratic denominators are possible and the code in the function `lin_random_denom()` can be altered to include this case. The paper discusses the drawbacks of this and higher degree polynomials will not work due to the amount of time the `solve` command takes.

## Dataset

A curated dataset using the methodology of the paper is available in the `Dataset` folder. The name of the JSON for each file corresponds to the elementary extension used to generate the (integrand, integral) pair. Each JSON contains 10,000 examples. These can be directly read as Maple expressions using the `parse` command in Maple or to Python SymPy Expressions using the `sympify` command. 

## Reference
[Generating Elementary Integrable Expressions](https://link.springer.com/chapter/10.1007/978-3-031-41724-5_2) - Barket et al.

```
@InProceedings{Barket2023_generation,
author={Barket, Rashid
and England, Matthew
and Gerhard, J{\"u}rgen},
title={Generating Elementary Integrable Expressions},
series = {Lecture Notes in Computer Science},
volume={14139},
booktitle={Computer Algebra in Scientific Computing},
year={2023},
publisher={Springer Nature Switzerland},
pages={21--38},
doi={10.1007/978-3-031-41724-5_2}
}
```




