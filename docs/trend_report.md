# Performance Trend Report Generator

## Goal

Provide a tool allowing to compare performance between two versions:

- Using a designated list of site with references to compare to
- Running automaticaly over these sites
- Extracting performance scores from logs
- Comparing to available previous runs (stored as jsons?)
- Copiling results in a graphical way into a report (stored as markdown with figures / color-coded tables)

All the results and reports should be added to the repo eventually. Best case scenario, the performance trend evaluation is automatically run by the CI.

## Approach

At the moment, the data is sored in test/data/... Generator code will be stored in test for the time being. Pending refactoring may alter these locations.

The generator is ran (by CI), collect everything from subdirectories (one by megasite since dem and wmap vrts have been defined by megasite). The .env file will export the input/extracts localisations as environment variables in order to run the whole generator with megasite specific information seamlessly. The goal is to heavily rely on the cachnig mecanism for dem/wmap and contourline since their generation is time consumming. A consolidated collection of the aforementioned files will have to be generated and maintained properly to insure consistency in the tren reports.

After completion (model estimation and validation on all sites), a "score scrapper will collect all the relevant metric from the logs and store them in a timestamped json (also containing the revision used by the run) that should be added to the repository for versionning and future reference.

The resulting json file will be interpreted and compared to available past jsons to generate a markdown report including multiple graph and table to illustrate the trends both globally and by site/mega-site. Timestamped markdown and associated figures should be added to the repository for versionning.
