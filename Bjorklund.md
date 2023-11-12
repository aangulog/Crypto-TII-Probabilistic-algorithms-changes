# _T 

The new method _T was created analyzing the pseudo-code (and the one which was implemented in the Björklund code) of [[BBSV22](https://eprint.iacr.org/2021/913.pdf)] and contrasted with the original paper [[BKW19](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2019.26)] they coincide in big $O$ notation. Our new _T method is more precise as the hidden terms of the big $O$ notation are at most $O(1)$.

```Python
    @staticmethod
    def _T(n, m, λ):
        # uses brute force
        if n <= 1:            
            return 1
        T = 0  # Time complexity
        l = floor(λ * n)  # n1, l from Björklund.
        s = 48 * n + 1  # number of iterations
        sumbin_n_2 = sum_of_binomial_coefficients(n, 2)
        sumbin_B = sum_of_binomial_coefficients(n - l, l + 4)
        T += 2 ** (n-l)  # Line 5, initialize array of size 2^{n-l}
        # Lines 6-12: Enters for loop, so each summand is multiplied by s
        T += s * (l+2) * m * sumbin_n_2 # Line 7, generate alphas, multiply and add coeficients
        # Line 8: No need to fill with zeros V at this point
        # Line 9-10: enters a for loop, so each term inside is multiplied by sumbin_B
        # Line 10:        
        T += s * sumbin_B * sumbin_n_2 * (l+2) # Partially evaluate R_i's        
        T += s * sumbin_B * Bjorklund._T(l, l+2, λ) # Recursive call
        # Line 11: Interpolation: this function makes two calls to the Z-transform
        T += s * (n-l) * sumbin_B # The first Z-transform takes (n-l) * sumbin_B * O(1)        
        T += s * (2 ** (n-l) - sumbin_B)  # Fill with zeros the rest of the table
        T += s * (n-l) * 2 ** (n-l) # The second Z-transform takes over all space        
        T += s * 2 ** (n-l) # Line 12: update the score table of size 2 ** (n-l)
        T += 2 ** (n-l) # line 14-17 does a for of size 2 ** (n-l), inside the for takes O(1)
        return T
```

# Optimal λ
The next graph shows the values of logaritmic time complexity with different values of lambdas for 20, 40, 60, 80,..., 1000 variables. This bahaviour is also seen for the intermediate values of the number of variables.

From this graph we can see two important facts:
  1. We observe that for values of lambda near 1 the time complexity grows rapidly, this means that computing these values is much more slower than for values near   0.
  2. The time complexity function for any value of variables is unimodal with respect to lambda.

![20,40,60, ,1000 var](https://github.com/aangulog/Crypto-TII-Probabilistic-algorithms-changes/assets/101427877/1e3a961c-0879-4f0d-8865-2e4287c55a0d)

Checking the old code for computing the optimal lambda

```Python
    def λ(self):
        n, m = self.nvariables_reduced(), self.npolynomials_reduced()
        k = self._k # Does not do anything
        min_complexity = Infinity
        optimal_λ = None

        for l in range(1, min(m, n - 1)):
            λ_ = l / n
            complexity = self._time_complexity_(λ_)

            if complexity < min_complexity:
                min_complexity = complexity
                optimal_λ = λ_

        self._λ = optimal_λ
        return self._λ
```

We can see the problem with this naïve approach of cheacking all possible lambdas in the interval from the first observation we made. Therefore thanks to the unimodality of the function we can implement an early abort for the method. The moment the function starts increasing we know we already passed the optimal lambda and we return this optimal value.

```Python
    def λ_early_abort(self):

        n, m = self.nvariables_reduced(), self.npolynomials_reduced()
        min_complexity = Infinity
        optimal_λ = None

        for l in range(1, min(m, n - 1)):

            λ_ = l / n
            complexity = self._time_complexity_(λ_)

            if complexity < min_complexity:
                min_complexity = complexity
                optimal_λ = λ_

            else: #early abort if the function starts increasing
                self._λ = optimal_λ
                return self._λ

        self._λ = optimal_λ
        return self._λ
```
Let us see how this two methods compare in performance

|  Number of variables   | Early abort (seconds) | No early abort (seconds )| 
|------------------------| --------------------- | -------------------------|
|100|0.006|0.2| x33.3
|200|0.04|1.6| x40
|300|0.1|5.5| x55
|400|0.3|14.1| x47
|500|0.7|30.4| x43.4
|600|1.4|59.3| x42.4
|700|2.6|109.8| x42.2
|800|4.4|190.3| x43.3
|900|7.0|301.5| x43.1
|1000|11.1|452.7| x40,8

The early abort implementation makes an optimization of about x43. And they, in fact, provide the same lambda (tested up to 530 variables) as we can see in the next graph:

![early aborts vs no early abort](https://github.com/aangulog/Crypto-TII-Probabilistic-algorithms-changes/assets/101427877/b4e1fe9b-3812-48a6-b7a9-a65d11b57cd8)

Thanks to the unimodality of the function, early abort gives us an optimal upper bound of stop for searching the optimal lambda. Our problem now can be restated as: How can we find this upper limit as quickly as possible?

Remember that we are starting this search of lambda at $\dfrac{1}{n}$, where $n$ is the number of variables, we would like to find a lower bound of searching, letting us start further away from 0.

[[BKW19](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2019.26)] demonstrated that the function $T(n,m)$ (for us the method _T) is exponentially asymptotic to $2^{0.803225n}$, this value is achieved when $\lambda =  0.196774680497$, we can start from this value and work towards 0, insted of starting from 0 and works toward the optimimum. But why are we able to do this? We are in fact studying the optimization of the __time_complexity_ method.

```Python
def _time_complexity_(self, λ):

    n, m = self.nvariables_reduced(), self.npolynomials_reduced()
    k = self._k

    return 8*k * n * sum([Bjorklund._T(n - i, m + k + 2, λ) for i in range(1, n)])
```

This method relies in calculating the value of $T(n,m)$ for all $n\in[1,2,...,n-1]$, and the exponential nature of _T tells us that the higher values of $n$ will be the ones that impact the most this method, therefore is reasonable to assume that a $\lambda$ that optimizes the values of $T(n,m)$ for higher n's is also going to be the $\lambda$ that optimizes __time_complexity_.
Our strategy will be: insted of starting at $\dfrac{1}{n}$, we will start at a $\dfrac{l}{n}$ that is closest to the optimal value of lambda. Then we will check if we are in the increasing or decreasing part of the function. If we are in the increasing part we will start working towards 0, and thanks to the unimodality find the minimum. In the other hand if we are in the decreasing part our current code of early abort will work fine.

```Python

    def λearly_abort(self):
        
        n, m = self.nvariables_reduced(), self.npolynomials_reduced()
        min_complexity = Infinity
        optimal_λ = None
        asymptotic_λ = 0.196774680497

        #first we will find the l such that l/n is closest to the upper_bound
        l = floor(n*asymptotic_λ)

        #now we want to check if moving towards 0 decreases or increases the value of the function

        if  self._time_complexity_((l-1)/n) <= self._time_complexity_(l / n): #i.e. decreses towards 0

            while l>=0:
                λ_ = l / n
                complexity = self._time_complexity_(λ_)
                l -= 1

                if complexity < min_complexity:
                    min_complexity = complexity
                    optimal_λ = λ_

                else:
                    self._λ = optimal_λ
                    return self._λ
            
        else: #i.e. increases towards 0, this means decreases towards 1
            while l<min(m,n):
                λ_ = l / n
                complexity = self._time_complexity_(λ_)
                l += 1

                if complexity < min_complexity:
                    min_complexity = complexity
                    optimal_λ = λ_

                else:
                    self._λ = optimal_λ
                    return self._λ

        self._λ = optimal_λ
        return self._λ
```
![inverse_early_abort](https://github.com/aangulog/Crypto-TII-Probabilistic-algorithms-changes/assets/101427877/ec3d4c17-7be8-4365-829d-b43ddd27e001)

We can see how for large values of $n$ this implementation cuts out a lot of unnecesary calculations, since our value is much closer to the opitmal lambda than to 0. Just to be safe we still calculate if our optimal lambda is above or below this line. But the long term behaviour hints us that this value is always lower.

Finally, let us see how this new implementation changes the calculation time.

|  Number of variables   | Early abort starting from the asymptotic lambda (seconds ) | Early abort starting form 0 (seconds) | 
|------------------------| --------------------- | -------------------------|
|600|0.5|1.4| x2.8
|800|1.1|4.4| x4
|1000|2.6|11.1| x4.3
|1200|5.4|26.2| x4.9
|1500|12.4|71.4| x5.7
|2000|36.9|266.8| x7.2
|3000|178.1|1725.1| x9.6
