# Presidential Plinko / Election Galton Boards

This repository contains code to produce Election Galton Boards: Galton boards
that roughly approximate the predictive distributions for the 2020 US
President electoral college vote, according to different poll aggregators.
You can see them in action at [presidential-plinko.com](http://presidential-plinko.com/).

The Galton boards are based on the models of [538](https://projects.fivethirtyeight.com/2020-election-forecast/)
and [the Economist](https://projects.economist.com/us-2020-forecast/president), both of
which are nice enough to open their data.

The Galton boards currently look like this:

![](galton_board-538.gif)
![](galton_board-the_economist.gif)

To download the latest versions of the Economist's and 538's models and re-render
everything, run [build.R](build.R).

The code for determining the binomial approximations to each modeler's predictions 
is in [binomial_approx_538.md](binomial_approx_538.md) ([source Rmd](binomial_approx_538.Rmd))
and [binomial_approx_economist.md](binomial_approx_538.md) ([source Rmd](binomial_approx_538.Rmd)).

The code for building the Galton boards is in [galton_board_quantile_ragg.Rmd](galton_board_quantile_ragg.Rmd).
This code is called from the above `binomial_approx*.Rmd` files.

The datasets are copyright their respective owners (see links above) and the
rest of the code in this repo is licensed under the MIT license.
