# demographic-scaling-model

`demographic-scaling-model` provides the R source code used in the paper: 

Bohk-Ewald, C. and M. Myrskylä (2020). A demographic scaling model to estimate COVID-19 infections. Preprint.

## How to run

The code needs to be executed in four steps, as defined by `step-*.R` files in the root directory.

### Prerequisites

- The `info-input-data.txt` file contains information about required input data.
  
- Due to copyrights you will need to download and save input data yourself. Original file names and URLs where to download them are given in the `info-input-data.txt` file. 

- Input data include confirmed cases and reported deaths attributable to COVID-19 of JHU CSSE (2020), population counts and abridged life tables of UNWPP 2019, infection fatality rates by 10-year age groups as, for example, published in Verity et al. (2020), and global age distribution of COVID-19 deaths as, for example, calculated based on data of Dudel et al. (2020).

- Make sure you have set the correct working directories before you start. 

### Execution

Run the `step-*.R` scripts in the root directory in the prescribed order.

## How to cite

If you use this code for academic research, please cite this GitHub repository as well as the paper noted above. 

## How to contribute

Please note that this source code is an academic project. We welcome any issues and pull requests.

## License

The source code of `demographic-scaling-model` is published under the [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). 
