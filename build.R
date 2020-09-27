# Update data sources
# the economist
download.file(
  "https://cdn.economistdatateam.com/us-2020-forecast/data/president/economist_model_output.zip",
  "data/economist/model_output.zip"
)
unzip(
  "data/economist/model_output.zip",
  files = "output/site_data//electoral_college_simulations.csv",
  exdir = "data/economist",
  junkpaths = TRUE,
)
unlink("data/economist/model_output.zip")

# 538
download.file(
  "https://projects.fivethirtyeight.com/2020-general-data/presidential_ev_probabilities_2020.csv",
  "data/538/presidential_ev_probabilities_2020.csv"
)


# Build the animations
rmarkdown::render("binomial_approx_economist.Rmd")
rmarkdown::render("binomial_approx_538.Rmd")


# Rebuild index.html
rmarkdown::render("index.Rmd")
