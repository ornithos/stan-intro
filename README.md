# stan-intro
Initial exploration of stan, including user defined functions

## Bayesian Logistic Regression
Playing around using sklearn wisconsin breast cancer dataset. Basic logistic regression using logit link function.

## User defined functions
Incremental attempts to incoporate more interesting operations and arguments. To date, the best I have is a `vector`,`vector` --> `double` function: a user defined dot product. 
  
I'm currently using Dan Foreman-Mackey's [hacked pystan](http://dfm.io/posts/stan-c++/) which allows `includes` as an argument to `pystan.StanModel` to permit custom user headers. Of course, this would be better pushed into a forked version of `stan::math` when I have a more useful function.

Now to code up a Bspline basis, function and derivative...
