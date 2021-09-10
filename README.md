### Supplementary matlab code for the paper *Modeling asset allocations and a new portfolio performance score*, Digital Finance, Springer 2021.

- **Installation**

To install the dependencies, set the working space to the root of this repository and run 

```matlab
install
```



- **Execution of the code**

To run the code, the main function is `compute_scores`. The inputs of the function are:

1. `sigma`: The estimation of the covariance matrix of the distribution of the assets' returns.
2. `mu`: The estimation of the mean of the distribution of the assets' returns.
3. `asset_returns`: The vector of the assets' returns to evaluate the performance of the portfolio.
4. `Ptf`: The portfolio to evaluate its performance.
5. `num_risk_levels`: The number of levels of risk of the allocation strategies in the stock market.
6. `num_dispersion_levels`: The number of levels of dispersion for each level of risk of the allocation strategies in the stock market.
7. `risk_behavioral_function`: The behavioral function for the risk of the allocation strategies.
8. `dispersion_behavioral_function`: The behavioral function for the dispersion of each risk level.



The outputs of the function are:

1. `parametric_score`: The parametric score of the input portfolio.
2. `Ws`: The sequence of weight vectors corresponding to each score in the vector of parametric score.
3. `bias_vector`: The bias vector computed by the input behavioral functions.
4. `mean_score`: The mean score.
5. `a_sequence`: The sequence of the parameters that control the level of dispersion for each risk level.
6. `q_sequence`: The sequence of the parameters that control the level of risk of the allocation strategies.
7. `volatility_sequence`: The sequence of portfolio volatilities corresponding to the `q_sequence`.



You can also see the code in the script `example`.
