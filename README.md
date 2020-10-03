# Presidential Plinko / Election Galton Boards

This repository contains code to produce Election Galton Boards: Galton boards
that roughly approximate the predictive distributions for the 2020 US
President electoral college vote, according to different poll aggregators.
You can see them in action at [presidential-plinko.com](http://presidential-plinko.com/).

The Galton boards are based on the models of [538](https://projects.fivethirtyeight.com/2020-election-forecast/)
and [the Economist](https://projects.economist.com/us-2020-forecast/president), both of
which are nice enough to open their data.

The Galton boards currently look like this:

![](boards/galton_board-538.gif)
![](boards/galton_board-the_economist.gif)

To download the latest versions of the Economist's and 538's models and re-render
everything, run [build.R](build.R).

The code for determining the binomial approximations to each modeler's predictions 
and then building the Plinko boards is in [binomial_approx_both.md](binomial_approx_both.md)
([source Rmd](binomial_approx_both.Rmd)). It uses my experimental [plinko](https://mjskay.github.io/plinko/) 
R package for constructing animated Plinko boards.

An older version of the code (pre-[plinko](https://mjskay.github.io/plinko/)) 
for building the Galton boards is in [galton_board_quantile_ragg.Rmd](galton_board_quantile_ragg.Rmd).

The datasets are copyright their respective owners (see links above) and the
rest of the code in this repo is licensed under the MIT license.
