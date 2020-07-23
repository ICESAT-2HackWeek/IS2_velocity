# Surface Velocity

Scripts for estimating ice velocity from ICESat-2 data.

Surface ice velocity from repeat passes using the across slope between beams and the slope between points on the same track. Cross overs could also be used. Which product?, maybe ATL06 or ATL03. It should be interesting to compare different areas (e.g. Antarctica, Greenland, Patagonia). Surface ice velocity from laser altimetry  has been done with ICESat 1 over the Ross Ice Shelf. The smaller footprint could improve the accuracy and the precise repeat passes could account for seasonality over certain regions.

## Files

* `.gitignore`
<br> Globally ignored files by `git` for the project.
* `environment.yml`
<br> `conda` environment description needed to run this project.
* `README.md`
<br> Description of the project. [Sample](https://geohackweek.github.io/wiki/github_project_management.html#project-guidelines)

## Folders

### `contributors`
Each team member has it's own folder under contributors, where he/she can
work on their contribution. Having a dedicated folder for one-self helps to 
prevent conflicts when merging with master.

### `notebooks`
Notebooks that are considered delivered results for the project should go in
here.

### `scripts`
Helper utilities that are shared with the team

